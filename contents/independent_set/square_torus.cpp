#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <queue>
#include <mutex>
#include<thread>
#include<mpi.h>
#include "types.h"
#include "debug.h"
#include "utils.h"
#include "type_stepped_calculate.h"
#include "modint.h"
#include "fibonacci.h"
#include "inner_product.h"
#include "independent_set/circular.h"



struct task_t {
	size_t n; // the target graph is C_n x C_m (a single task calculates the values for m <= M for some M)
	size_t index; // compressed index: [0, compress_indices.size())
	u64 bitpattern; // (one of) bit pattern(s) of row 0: [0, 1 << n)
	int multiplier; // how many equivalent bitpatterns were merged into this single task: [1, n]
};

// Module to save partial results(per task)
/*
	Format (per line):
		n (task index) (task bitpattern) (task multiplier) M result_0 result_1 ... result_M
	where result_i is the result(without multiplier applied) for n x i grid.
	Task bitpattern(and task multiplier) can be calculated from task index, but they are recorded anyway for better human readability
*/
struct results_writer_t {
	std::string fname;
	std::ofstream file;
	
	results_writer_t () = default;
	results_writer_t(const std::string &fname) : fname(fname), file(fname, std::ios_base::app) {
		if (!file) {
			std::cerr << "Failed to open results file in mode 'a': " << fname << std::endl;
			exit(1);
		}
	}
	void write(const task_t &task, const std::vector<cpp_int> &results_unmultiplied) {
		assert(results_unmultiplied.size());
		size_t M = results_unmultiplied.size() - 1;
		if (file) {
			file << task.n << " " << task.index << " " << task.bitpattern << " " << task.multiplier << " " << M;
			for (auto &result : results_unmultiplied) file << " " << result;
			file << std::endl;
		}
	}
};

std::ofstream log_file;
std::mutex log_lock;



/*
	A struct for
	 - holding precalculated data necessary for the main calculation
	 - task enumeration/distribution mechanism
*/
struct Runner {
	size_t n;
	size_t M;
	size_t p, P;
	// transformer provides the core function transform_one_row which calculates a DP table for the next row from the DP table for the current row,
	// or equivalently multiplies the transfer matrix T_n to a given vector
	circular_no_adjacent_transformer transformer;
	std::vector<u64> &fib = transformer.fib;
	// compress_indices[i]: {minimum of (uncompressed) fibonacci indices, multiplicity}
	std::vector<std::pair<u64, int> > &compress_indices = transformer.compress_indices;
	
	std::vector<cpp_int> calc_task(task_t task) {
		int total_stages = (M + 1) >> 1;
		std::vector<cpp_int> res = {1};
		type_stepped_calculate(total_stages,
			// init function
			[&] (auto &dp) {
				dp.resize(transformer.get_buf_size(n));
				dp[compress_indices[task.index].first] = 1;
			},
			// upgrade-if-necessary function
			[&] (auto &cur_dp, auto &next_dp) {
				using cur_t = typename std::remove_reference<decltype(cur_dp[0])>::type;
				using next_t = typename std::remove_reference<decltype(next_dp[0])>::type;
				// std::cerr << " upgrade to: u" << sizeof(next_t) * 8 << std::endl;
				// (values in next vector) <= (sum of values in cur vector), so upgrade if (sum of current vector) exceeds type limit
				if (std::accumulate(cur_dp.begin(), cur_dp.end(), (next_t) 0) > std::numeric_limits<cur_t>::max()) {
					next_dp = std::vector<next_t>(cur_dp.begin(), cur_dp.end());
					return true;
				} else return false;
			},
			// transfer function
			[&] (int stage, auto &dp) {
				(void) stage;
				// std::cerr << "stage : " << stage << std::endl;
				auto prev = dp;
				transformer.transform_one_row(n, prev, dp);
				auto tmp = inner_product_ab_bb<u2048>(prev, dp);
				res.push_back(tmp.first);
				res.push_back(tmp.second);
			}
		);
		assert(res.size() >= M + 1);
		res.resize(M + 1);
		return res;
	}
	void thread_func(int thread_index, int n_threads) {
		(void) thread_index;
		while (1) {
			auto t0 = Timer::get();
			bool exit = false;
			// fetch a task
			task_t cur_task;
			{
				std::lock_guard<std::mutex> locking(thread_lock);
				if (!tasks.size()) exit = true;
				else cur_task = tasks.front(), tasks.pop();
			}
			if (exit) break;
			
			// (core) calculate
			auto cur_res = calc_task(cur_task);
			assert(cur_res.size() - 1 == M); // cur_res[0...M]
			
			auto t1 = Timer::get();
			// record the results and show the progress
			{
				std::lock_guard<std::mutex> locking(thread_lock);
				for (size_t i = 0; i <= M; i++) results[i] += cur_res[i] * cur_task.multiplier;
				result_writer.write(cur_task, cur_res);
				register_task_time(t1 - t0);
				double avg_gops = transform_cost / get_avg_time() / 1000000000.0;
				if (p == 0) fprintf(stderr, "\r%d / %d   avg time per task: %.2f s   avg GOP/s:  ST: %.4f  MT: %.4f   ",
					++completed_tasks_num, total_tasks_count, get_avg_time(), avg_gops, avg_gops * n_threads);
				// write log
				log_lock.lock();
				log_file << "[" << timestamp_now() << "] Thread #" << thread_index << ": Finish (n,M)=(" << n << "," << M << ") task_index=" << cur_task.index <<
					" bitpattern=" << to_binary_str(cur_task.bitpattern, n) << " time=" << Timer::s(t1 - t0) << " s" << std::endl;
				log_lock.unlock();
			}
		}
	}
	
	const double AVG_MULTIPLIER = 0.95;
	double thread_time_weighted_sum{};
	double thread_time_weight_sum = 0;
	Timer::duration_type thread_time_sum{};
	
	results_writer_t result_writer;
	
	std::queue<task_t> tasks;
	std::mutex thread_lock;
	std::vector<cpp_int> results;
	u32 completed_tasks_num = 0;
	u32 total_tasks_count;
	u64 transform_cost;
	
	void register_task_time(Timer::duration_type duration) {
		thread_time_weighted_sum *= AVG_MULTIPLIER;
		thread_time_weighted_sum += Timer::s(duration);
		thread_time_weight_sum *= AVG_MULTIPLIER;
		thread_time_weight_sum += 1.0;
		thread_time_sum += duration;
	}
	double get_avg_time() {
		return thread_time_weighted_sum / thread_time_weight_sum;
	}
	
	Runner (size_t n, size_t M, size_t p, size_t P) : n(n), M(M), p(p), P(P), transformer(n), results(M + 1) {
		for (size_t i = 0; i < compress_indices.size(); i++) if (i % P == p) tasks.push({n, i,
			get_bitpattern_from_circular_fibonacci_index<u64>(compress_indices[i].first, n), compress_indices[i].second});
		total_tasks_count = compress_indices.size();
		transform_cost = (u64) n * M * transformer.get_buf_size(n);
	}
	size_t get_workload() { return transform_cost * total_tasks_count; }
	size_t get_max_buffer_bytes_per_thread() { return transformer.get_buf_size(n) * (1024 / 8) * 2; } // u1024
};


// Load results file, and get rid of unnecessary tasks from runner
void load_saved_results_for_runner(std::string txt_fname_prefix, Runner *runner) {
	log_file << "[" << timestamp_now() << "] Loading saved results..." << std::endl;
	
	// statistics
	size_t insufficient_m = 0;
	
	std::set<u32> already_calculated_task_indices;
	
	// load from {txt_fname_prefix}0.txt, {txt_fname_prefix}1.txt, ... while they exist
	for (size_t i = 0; ; i++) {
		std::string fname = txt_fname_prefix + std::to_string(i) + ".txt";
		
		log_file << "[" << timestamp_now() << "]   Loading: " << fname << std::endl;
		
		// Open file
		std::ifstream file(fname);
		if (!file) break;
		
		task_t cur_task;
		size_t M;
		std::vector<cpp_int> results_unmultiplied;
		
		while (file >> cur_task.n >> cur_task.index >> cur_task.bitpattern >> cur_task.multiplier >> M) {
			results_unmultiplied.resize(M + 1);
			for (size_t i = 0; i <= M; i++) if (!(file >> results_unmultiplied[i])) {
				std::cerr << "Error: unexpected EOF or invalid input when reading results vector" << std::endl;
				break;
			}
			
			if (runner->n != cur_task.n) continue;
			if (M < runner->M) { // M not enough
				insufficient_m++;
				continue;
			}
			
			// validity check
			if (cur_task.index >= runner->compress_indices.size()) {
				std::cerr << "Warning: ignoring invalid results data: task index out of range(n=" << cur_task.n << ", task index=" << cur_task.index << ")" << std::endl;
				continue;
			}
			if (cur_task.bitpattern != get_bitpattern_from_circular_fibonacci_index<u64>(runner->compress_indices[cur_task.index].first, cur_task.n)) {
				std::cerr << "Warning: ignoring invalid results data: inconsistent bitpattern(n=" << cur_task.n << ", index="
					<< cur_task.index << ", bitpattern=" << cur_task.bitpattern << ")" << std::endl;
				continue;
			}
			if (already_calculated_task_indices.count(cur_task.index)) continue; // already loaded
			
			// add to runner->results
			for (size_t i = 0; i <= runner->M; i++) runner->results[i] += results_unmultiplied[i] * cur_task.multiplier;
			already_calculated_task_indices.insert(cur_task.index);
		}
	}
	// remove loaded tasks by filtering
	size_t prev_tasks_count = runner->total_tasks_count;
	std::string prev_workload_str = put_SI_prefix(runner->get_workload()) + "OP";
	
	std::queue<task_t> new_tasks;
	while (runner->tasks.size()) {
		auto task = runner->tasks.front();
		runner->tasks.pop();
		if (!already_calculated_task_indices.count(task.index)) new_tasks.push(task);
	}
	runner->tasks = new_tasks;
	runner->total_tasks_count = runner->tasks.size();
	
	
	log_file << "[" << timestamp_now() << "] " << already_calculated_task_indices.size() << " tasks cached out of " << prev_tasks_count <<
		"(" << prev_workload_str << ")" << std::endl;
	log_file << "[" << timestamp_now() << "] Remaining tasks: " << runner->total_tasks_count <<
		"(" << put_SI_prefix(runner->get_workload()) << "OP)" << std::endl;
	if (insufficient_m) log_file << "[" << timestamp_now() << "] " << insufficient_m << " items rejected because of insufficient M" << std::endl;
}

#define LOG_FNAME_BASE "square_torus_log"
#define CACHE_FNAME_BASE "square_torus_partial"

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	
	int p, P;
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
	
	
	int n;
	int M;
	int n_threads;
	if (p == 0) {
		std::cout << "Independent set in C_n x C_m calculator" << std::endl;
		std::cout << "Enter n M T  (n, max m, # of threads, respectively)" << std::endl;
		std::cout << "Example: 18 24 2" << std::endl;
		std::cout << "> ";
		if (!(std::cin >> n >> M >> n_threads)) {
			std::cerr << "Invalid input" << std::endl;
			n = M = n_threads = -1;
		}
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_threads, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (n >= 0) {
		log_lock.lock();
		log_file << "[" << timestamp_now() << "] Config: (n, M)=(" << n << "," << M << ") n_threads=" << n_threads << std::endl;
		log_lock.unlock();
		
		if (n == 0) { // needs special handling
			auto fib = fibonacci_sequence<cpp_int>(M);
			auto res = circular_fibonacci_sequence<cpp_int>(M, fib);
			for (int i = 0; i <= M; i++) std::cout << n << " " << i << " " << res[i] << std::endl;
		} else {
			Runner runner(n, M, p, P);
			load_saved_results_for_runner(CACHE_FNAME_BASE, &runner);
			runner.result_writer = std::move(results_writer_t(my_cache_fname));
			
			if (p == 0) std::cerr << "Memory consumption: " << put_SI_prefix(runner.get_max_buffer_bytes_per_thread() * n_threads) << "B" << std::endl;
			
			log_lock.lock();
			log_file << "[" << timestamp_now() << "] Finished loading cache file" << std::endl;
			log_lock.unlock();
			
			std::vector<std::thread> threads;
			for (int i = 0; i < n_threads; i++) threads.push_back(std::thread([&,i]() { runner.thread_func(i, n_threads); }));
			for (auto &thread : threads) thread.join();
			
			if (p == 0) {
				std::cerr << std::endl;
				for (size_t i = 0; i < runner.results.size(); i++) std::cout << n << " " << i << " " << runner.results[i] << std::endl;
			}
		}
	}
	
	MPI_Finalize();
	
	return 0;
}


