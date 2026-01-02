#include <map>
#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include "types.h"
#include "debug.h"
#include "utils.h"
#include "modint.h"
#include "modint_avx.h"
#include "independent_set/square_grid_mod.h"
#include <random>
#include <mpi.h>
#include <boost/multiprecision/miller_rabin.hpp>



// Module to save partial results
/*
	Format (per line):
		n mod M result_0 result_1 ... result_M
	where result_i is the result for n x i grid(modulo `mod`).
*/
struct partial_results_writer {
	std::string fname;
	std::ofstream file;
	
	partial_results_writer () = default;
	partial_results_writer(const std::string &fname) : fname(fname), file(fname, std::ios_base::app) {
		if (!file) {
			std::cerr << "Failed to open results file in mode 'a': " << fname << std::endl;
			exit(1);
		}
	}
	// T is the type of calculation(u16 or u32)
	template<typename T> void write(size_t n, T mod, const std::vector<T> &results) {
		assert(results.size() > 0);
		size_t M = results.size() - 1;
		if (file) {
			file << n << " " << mod << " " << M;
			for (auto val : results) file << " " << val;
			file << std::endl;
		}
	}
};

partial_results_writer partial_results_file;

std::ofstream log_file;
std::mutex log_lock;



// Load results file, and return cache map: (n, mod) -> results vector
// T: calculation type
template<typename T>
std::map<std::pair<size_t, T>, std::vector<T> > load_saved_results(std::string txt_fname_prefix, bool verbose) {
	std::map<std::pair<size_t, T>, std::vector<T> > res_cache;
	
	size_t max_mod = 1ULL << (sizeof(T) * 8 - 1);
	
	// load from {txt_fname_prefix}0.txt, {txt_fname_prefix}1.txt, ... while they exist
	for (size_t i = 0; ; i++) {
		std::string fname = txt_fname_prefix + std::to_string(i) + ".txt";
		
		// Open file
		std::ifstream file(fname);
		if (!file) break;
		
		if (verbose) std::cerr << "Loading " << fname << "..." << std::endl;
		
		size_t n;
		size_t mod_;
		size_t M;
		std::vector<T> results;
		
		while (file >> n >> mod_ >> M) {
			// mod too large(happens when the recorded results are from u32 calculation and 
			// the current calculation is done in u16)
			if (mod_ > max_mod) continue;
			T mod = mod_;
			
			results.resize(M + 1);
			for (size_t i = 0; i <= M; i++) if (!(file >> results[i])) {
				std::cerr << "Failed to read " << i << "-th element of results vector (n=" << n <<
					", mod=" << mod << ")" << std::endl;
				break;
			}
			if (!file) break;
			
			if (res_cache.count({n, mod})) {
				auto &old_res = res_cache[{n, mod}];
				for (size_t i = 0; i < std::min(old_res.size(), results.size()); i++) if (old_res[i] != results[i]) {
					std::cerr << "Warning: inconsistent cache results(n=" << n << ", mod=" << mod << ")" << std::endl;
					break;
				}
				if (old_res.size() < results.size()) res_cache[{n, mod}] = results;
			} else res_cache[{n, mod}] = results;
		}
	}
	
	return res_cache;
}



/*
	calculate the number of independent set in P_n x P_m modulo `mod` for every 0 <= m <= M
	and return them in vector of length M+1.
*/
template<typename transitioner_type, typename T = typename transitioner_type::T>
std::vector<T> calc_for_one_mod(transitioner_type *transitioner, int n, int M, T mod, int n_threads, bool verbose) {
	int n_stages = M / 2;
	std::vector<T> res{ 1, (T) (transitioner->fib[n] % mod) };
	
	// prepare buffers
	auto setup_start = Timer::get();
	std::vector<T> dp = transitioner->get_initial_buffer(n);
	std::vector<T> dp_next(transitioner->get_buf_size(n));
	auto setup_end = Timer::get();
	if (verbose) std::cerr << "setup: " << Timer::diff_ms(setup_start, setup_end) << " ms" << std::endl;
	
	for (int i = 0; i < n_stages; i++) {
		transitioner->transition_one_row(n, dp, dp_next, mod, n_threads, verbose);
		
		auto inner_prod_start = Timer::get();
		std::pair<T, T> t = inner_product_mod_ab_bb(dp, dp_next, mod, -1, n_threads);
		res.push_back(t.first);
		res.push_back(t.second);
		auto inner_prod_end = Timer::get();
		if (verbose) std::cerr << "inner: " << Timer::diff_ms(inner_prod_start, inner_prod_end) << " ms" << std::endl;
		
		std::swap(dp, dp_next); // dp = dp_next
	}
	
	res.resize(M + 1);
	
	return res;
}


/*
	calculate the number of independent set in P_n x P_m for every 0 <= m <= M, taking cache into account,
	and return them in vector of length M+1.
	
	p: process id
	P: # of processes
	cache[{n, mod}][m]: cached results for those parameters
*/
template<typename mintNxM, typename T = typename mintNxM::T>
std::vector<cpp_int> calc(size_t n, size_t M, int n_threads,
	const std::map<std::pair<size_t, T>, std::vector<T> > &cache, size_t p, size_t P, bool verbose) {
	
	using T_double =  // if T is u16, T_double is u32. if T is u32, T_double is u64
		std::conditional_t<
			std::is_same_v<T, u16>, u32,
			std::conditional_t<
				std::is_same_v<T, u32>, u64,
				void
			>
		>;
	static_assert(!std::is_same_v<T_double, void>, "T must be u16 or u32");
	
	
	/*
		Generate enough number of mods(of type T) to recover the results.
		A mod must be less than 2^((# of bits in T) - 1) in order to avoid overflow in addition in type T
	*/
	size_t T_bits = sizeof(T) * 8;
	// Generate mods close to but less than 2^(T_bits-1)
	size_t n_mods = (n + 1) * (M + 1) * std::log2(1.5031) // hard square entropy constant
		/ (T_bits - 2) // one mod = ((T_bits - 1) - eps) bit
		+ 3; // to be safe
	
	using namespace boost::multiprecision;
	std::vector<T> mods;
	for (T x = (1U << (T_bits - 1)) - 1; x > 0; x--) {
		if (mods.size() >= n_mods) break;
		if (miller_rabin_test(x, 25)) mods.push_back(x);
	}
	assert(mods.size() == n_mods);
	
	
	// Show stastistics
	transitioner_IS_square_grid_mod<mintNxM> transitioner(n);
	if (p == 0) {
		std::cerr << "----------------------------" << std::endl;
		std::cerr << "# of mods: " << n_mods << std::endl;
		std::cerr << "Cost per mod: " << put_SI_prefix(transitioner.get_buf_size(n) * n * M) << "OP" << std::endl;
		std::cerr << "Total cost: " << put_SI_prefix(transitioner.get_buf_size(n) * n * M * n_mods) << "OP" << std::endl;
		std::cerr << "Memory Consumption: " << put_SI_prefix(transitioner.get_buf_size(n) * sizeof(T) * 2) << "B" << std::endl;
		std::cerr << "----------------------------" << std::endl;
	}
	
	// Actual Calculation
	bool incomplete = false;
	std::vector<std::vector<T> > results_mod;
	for (size_t i = 0; i < n_mods; i++) {
		T mod = mods[i];
		
		if (cache.count({n, mod}) && cache.at({n, mod}).size() >= M + 1) {
			// Cache exists && contains enough elements -> skip
			auto cached_res = cache.at({n, mod});
			cached_res.resize(M + 1);
			results_mod.push_back(cached_res);
			if (p == 0) std::cerr << "#" << i << " (mod=" << mod << "): cached  " << std::endl;
		} else {
			if (i % P == p) { // only calculate my task
				auto t0 = Timer::get();
				auto cur_res = calc_for_one_mod(&transitioner, n, M, mod, n_threads, verbose);
				auto t1 = Timer::get();
				assert(cur_res.size() == M + 1);
				partial_results_file.write(n, mod, cur_res);
				results_mod.push_back(cur_res);
				auto time_s = Timer::diff_s(t0, t1);
				auto ops = transitioner.get_buf_size(n) * n * M / time_s;
				std::cerr << "#" << i << " (mod=" << mod << "): done  ";
				std::cerr << time_s << " s  " << put_SI_prefix(ops) << "OP/s" << std::endl;
				
				log_lock.lock();
				log_file << "[" << timestamp_now() << "]: Finish (n,M)=(" << n << "," << M << ") mod_id=" << i
					<< " mod=" << mod << " time=" << time_s << " s(" << put_SI_prefix(ops) << " OP/s)" << std::endl;
				log_lock.unlock();
			} else incomplete = true; // there are mods for which the results the process doesn't know
		}
	}
	std::cerr << std::endl;
	
	// If `incomplete`, there are fresh new results calculated by other processes.
	// We will have the user to run the program again to get the final results
	//  (instead of writing complicated code of interprocess communication).
	if (incomplete) return {};
	
	// Recover actual values from modulo results by Garner algorithm
	std::vector<cpp_int> res;
	for (size_t i = 0; i <= M; i++) {
		std::vector<T> vals;
		for (size_t j = 0; j < results_mod.size(); j++) vals.push_back(results_mod[j][i]);
		res.push_back(garner<T, T_double>(vals, mods));
	}
	
	return res;
}

/*
	In this program, a **task** means calculating the target number(# of independent sets) in a single modulo.
	We could calculate the target number directly using cpp_int(multiprecision integer),
		but it would take too much RAM so we instead calculate it with a small(16 or 32 bit) modulo.
	By solving roughly (# of bits in target number) / (# of bits of a mod) tasks, we can restore the original number.
	
	If this MPI program is run in multiprocess(e.g. mpirun -n 2), they will independently solve tasks.
	Therefore, the required amount of RAM is multipled in the case.
	They only coordinate in that they share the task queue.
	
	On the other hand, you can specify the number of threads(T) in the in-program prompt.
	They will work on a single task, so T does not affect the required amount of RAM.
	(In fact, the threads are only spawn in transitioner::transition_one_row and inner_product_mod*_ab_bb)
*/


#define CACHE_FNAME_BASE "square_partial"
#define LOG_FNAME_BASE "square_log"


int main(int argc, char **argv) {
	using mintNxM = mint16x32; // change this to change avx type & calculation type(u32 or u16)
	std::string backend_name = "mint" + std::to_string(8 * sizeof(typename mintNxM::T)) +
		"x" + std::to_string(mintNxM::LEN);
	
	MPI_Init(&argc, &argv);
	
	int p; // my process id
	int P; // total # of processes
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Comm_rank(MPI_COMM_WORLD, &p);
	
	std::string my_cache_fname = CACHE_FNAME_BASE + std::to_string(p) + ".txt";
	std::string my_log_fname = LOG_FNAME_BASE + std::to_string(p) + ".txt";
	
	// Open log file and write startup log
	log_file = std::move(std::ofstream(my_log_fname, std::fstream::out | std::fstream::app));
	if (!log_file) {
		std::cerr << "Failed to open log file: " << my_log_fname << std::endl;
		return 1;
	}
	log_lock.lock();
	log_file << "[" << timestamp_now() << "] -------------------------------- Program startup --------------------------------" << std::endl;
	log_lock.unlock();
	
	
	// Input n, M, T
	int n;
	int M; // max m
	int T; // number of threads
	if (p == 0) {
		std::cout << "Independent set in P_n x P_m calculator(backend=" << backend_name << std::endl;
		std::cout << "Enter n M T (n, max m, # of threads, respectively)" << std::endl;
		std::cout << "Example: 24 30 2" << std::endl;
		std::cout << "Note that n <= 45 is recommended since the calculation time and the memory consumption is exponential to n." << std::endl;
		std::cout << "> " << std::flush;
		if (!(std::cin >> n >> M >> T)) {
			std::cerr << "Invalid input" << std::endl;
			n = M = T = -1;
		}
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	if (n >= 0) {
		log_lock.lock();
		log_file << "[" << timestamp_now() << "] Config: (n, M)=(" << n << "," << M << ") (p, P, T)=("
			<< p << ", " << P << ", " << T << ") backend=" << backend_name << std::endl;
		log_lock.unlock();
		
		// Load cache
		auto cache = load_saved_results<typename mintNxM::T>(CACHE_FNAME_BASE, p == 0); // load all cache files(incl. other than mine)
		MPI_Barrier(MPI_COMM_WORLD);
		partial_results_file = std::move(partial_results_writer(my_cache_fname)); // open my log file for writing
		
		// Actual calculation
		auto res = calc<mintNxM>(n, M, T, cache, p, P, false);
		if (p == 0) {
			if (!res.size()) std::cout << "Fresh results are not shared between processes. Rerun the program to get the global results." << std::endl;
			else {
				for (int i = 0; i <= M; i++) std::cout << n << " " << i << " " << res[i] << std::endl;
			}
		}
	}
	
	MPI_Finalize();
	
	return 0;
}

