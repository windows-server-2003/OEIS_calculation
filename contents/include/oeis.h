#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <deque>
#include <thread>
#include <mutex>
#include <type_traits>
#include "utils.h"
#include "types.h"


struct OEISGenerator {
	int thread_num = 1;
	
	std::vector<cpp_int> result;
	std::mutex thread_lock;
	std::deque<int> tasks;
	int num_thread_running;
	std::vector<std::pair<int, int> > thread_limits;
	template<typename Sol> void thread_func(int thread_index, Sol *sol) {
		(void) thread_index;
		while (1) {
			bool exit = false;
			int next_n = -1;
			{
				std::lock_guard<std::mutex> locking(thread_lock);
				if (tasks.size()) {
					next_n = tasks.front();
					for (auto limit : thread_limits) if (next_n >= limit.first && num_thread_running > limit.second) exit = true;
					if (exit) num_thread_running--;
					else tasks.pop_front();
				} else exit = true;
			}
			if (exit) break;
			auto res = sol->template calc<cpp_int>(next_n);
			{
				std::lock_guard<std::mutex> locking(thread_lock);
				result[next_n] = res;
				std::cerr << next_n << " " << res << std::endl;
				// std::cerr << "Thread #" << thread_index << std::endl;
			}
		}
		// std::cerr << "Thread #" << thread_index << " exit" << std::endl;
	}
	
	struct general_tag {};
	struct special_tag : general_tag {};
	template<typename> struct int_ { typedef int type; };
	
	template<typename Sol, typename int_<decltype(Sol(0, 0))>::type = 0> Sol get_default_sol_(int n_min, int n_max, special_tag) { return Sol(n_min, n_max); }
	template<typename Sol> Sol get_default_sol_(int, int, general_tag) { return Sol(); }
	template<typename Sol> Sol get_default_sol(int n_min, int n_max) { return get_default_sol_<Sol>(n_min, n_max, special_tag()); }
	
	
	template<typename Sol> auto run_multithread(int n_min, int n_max, Sol &sol) -> decltype(std::declval<Sol>().template calc<cpp_int>(0), void()) {
		for (int n = n_min; n <= n_max; n++) tasks.push_back(n);
		result.resize(n_max + 1);
		
		num_thread_running = thread_num;
		std::vector<std::thread> threads;
		for (int i = 0; i < thread_num; i++) threads.push_back(std::thread(&OEISGenerator::thread_func<Sol>, this, i, &sol));
		for (auto &thread : threads) thread.join();
	}
	template<typename Sol> auto run_multithread(int n_min, int n_max, Sol &sol) -> decltype(std::declval<Sol>().template calc_all<cpp_int>(0), void()) {
		(void) n_min;
		(void) n_max;
		result = sol.template calc_all<cpp_int>(n_max);
	}
	template<typename Sol> auto is_solver_individual() -> decltype(std::declval<Sol>().template calc<cpp_int>(0), bool()) { return true; }
	template<typename Sol> auto is_solver_individual() -> decltype(std::declval<Sol>().template calc_all<cpp_int>(0), bool()) { return false; }
	
	// finalize function
	template<typename Sol, typename int_<decltype(std::declval<Sol>().finalize(0, 0, std::declval<std::vector<cpp_int> &>()))>::type = 0> void finalize_(int n_min, int n_max, Sol &sol, special_tag) {
		sol.finalize(n_min, n_max, result);
	}
	template<typename Sol> void finalize_(int, int, Sol &, general_tag) {}
	template<typename Sol> void finalize(int n_min, int n_max, Sol &sol) { finalize_<Sol>(n_min, n_max, sol, special_tag()); }

	struct GenOption {
		std::string outfile_name;
		bool outfile_append;
		
		int n_min;
		int n_max;
		int thread_num;
		std::vector<std::pair<int, int> > thread_limits;
	};
	template<typename Sol> void solve(const GenOption &opt) {
		std::ofstream file_stream(opt.outfile_name, opt.outfile_append ? (std::ios_base::out | std::ios_base::app) : std::ios_base::out);
		if (opt.outfile_name != "" && !file_stream) {
			std::cerr << "Failed to open " << opt.outfile_name << std::endl;
			return;
		}
		this->thread_num = opt.thread_num;
		this->thread_limits = opt.thread_limits;
		
		auto t0 = Timer::get();
		Sol sol = get_default_sol<Sol>(opt.n_min, opt.n_max);
		run_multithread(opt.n_min, opt.n_max, sol);
		finalize<Sol>(opt.n_min, opt.n_max, sol);
		auto t1 = Timer::get();
		
		std::string time_str = Timer::diff_str(t0, t1);
		std::string date_str = Timer::get_date_str();
		if (date_str.size()) date_str = date_str.substr(0, date_str.size() - 1); // ignore the last linebreak
		std::string thread_str = std::to_string(is_solver_individual<Sol>() ? thread_num : 1) + " threads";
		if (thread_limits.size()) {
			thread_str += " (limit";
			for (auto limit : thread_limits) thread_str += " (" + std::to_string(limit.first) + ", " + std::to_string(limit.second) + ")";
			thread_str += ")";
		}
		
		if (opt.outfile_name == "") {
			for (int i = opt.n_min; i <= opt.n_max; i++) std::cout << i << " " << result[i] << std::endl;
			std::cerr << time_str << std::endl;
		} else {
			file_stream << "# " << date_str << std::endl;
			file_stream << "# " << time_str << std::endl;
			file_stream << "# " << thread_str << std::endl;
			for (int i = opt.n_min; i <= opt.n_max; i++) file_stream << i << " " << result[i] << std::endl;
		}
	}
	
	int get_sequence_number(int argc, char **argv) {
		for (int i = 1; i < argc; i++) {
			if (argv[i][0] == 'A') {
				int res;
				try { res = std::stoi(argv[i] + 1); }
				catch (const std::exception &except) { res = -1; }
				
				if (res != -1) return res;
			}
		}
		return -1;
	}
	bool get_use_file(int argc, char **argv) {
		for (int i = 1; i < argc; i++) if (strcmp(argv[i], "-f") == 0) return true;
		return false;
	}
	bool get_append(int argc, char **argv) {
		for (int i = 1; i < argc; i++) if (strcmp(argv[i], "-a") == 0) return true;
		return false;
	}
	GenOption get_cmdline_options(int argc, char **argv, int seq_num) {
		GenOption res;
		
		bool use_file = get_use_file(argc, argv);
		char outfile_name_buf[64] = { 0 };
		if (use_file) snprintf(outfile_name_buf, 64, "b%06d.txt", seq_num);
		res.outfile_name = outfile_name_buf;
		res.outfile_append = get_append(argc, argv);
		return res;
	}
};

