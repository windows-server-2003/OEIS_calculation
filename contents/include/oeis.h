#pragma once
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

template<typename task_t> struct ThreadPool {
	std::deque<task_t> tasks;
	ThreadPool () {}
	ThreadPool (const std::vector<task_t> &tasks) : tasks(tasks.begin(), tasks.end()) {}
	
	void add_task(const task_t &task) { tasks.push_back(task); }
	
	std::mutex lock;
	int num_thread_running;
	template<typename process_func_t, typename check_limit_func_t> void thread_func(int thread_id,
		const process_func_t &process_func, const check_limit_func_t &check_limit_func) {
		
		while (1) {
			bool exit = false;
			task_t next_task{};
			{
				std::lock_guard<std::mutex> locking(lock);
				if (tasks.size()) {
					next_task = tasks.front();
					if (!check_limit_func(next_task, num_thread_running)) exit = true;
					else tasks.pop_front();
				} else exit = true;
			}
			if (exit) break;
			process_func(next_task, thread_id);
		}
		{
			std::lock_guard<std::mutex> locking(lock);
			num_thread_running--;
		}
	}
	
	template<typename process_func_t, typename check_limit_func_t> void run(int thread_num,
		const process_func_t &process_func, const check_limit_func_t &check_limit_func) {
		
		num_thread_running = thread_num;
		std::vector<std::thread> threads;
		for (int i = 0; i < thread_num; i++) threads.push_back(std::thread([&] (int thread_id) {
			thread_func(thread_id, process_func, check_limit_func);
		}, i));
		for (auto &thread : threads) thread.join();
	}
	
};

struct OEISGenerator {
	int thread_num = 1;
	
	std::vector<cpp_int> result;
	std::vector<std::pair<int, int> > thread_limits;
	
	struct general_tag {};
	struct special_tag : general_tag {};
	template<typename> struct int_ { typedef int type; };
	
	template<typename Sol, typename int_<decltype(Sol(0, 0))>::type = 0> Sol get_default_sol_(int n_min, int n_max, special_tag) { return Sol(n_min, n_max); }
	template<typename Sol> Sol get_default_sol_(int, int, general_tag) { return Sol(); }
	template<typename Sol> Sol get_default_sol(int n_min, int n_max) { return get_default_sol_<Sol>(n_min, n_max, special_tag()); }
	
	enum class RunType {
		TERM_THREADED,
		ALL,
		CUSTOM_THREADED
	};
	template<typename Sol> auto run_multithread(int n_min, int n_max, Sol &sol) -> decltype(std::declval<Sol>().template calc<cpp_int>(0), RunType()) {
		std::vector<int> tasks;
		for (int n = n_min; n <= n_max; n++) tasks.push_back(n);
		result.resize(n_max + 1);
		
		ThreadPool<int> pool(tasks);
		pool.run(thread_num, [&] (int task, int) {
			result[task] = sol.template calc<cpp_int>(task);
			{
				std::lock_guard<std::mutex> locking(pool.lock);
				std::cerr << task << " " << result[task] << std::endl;
			}
		}, [&] (int task, int thread_running) {
			for (auto limit : thread_limits) if (task >= limit.first && thread_running > limit.second) return false;
			return true;
		});
		
		return RunType::TERM_THREADED;
	}
	template<typename Sol> auto run_multithread(int n_min, int n_max, Sol &sol) -> decltype(std::declval<Sol>().template calc_all<cpp_int>(0), RunType()) {
		(void) n_min;
		result = sol.template calc_all<cpp_int>(n_max);
		return RunType::ALL;
	}
	template<typename Sol> auto run_multithread(int n_min, int n_max, Sol &sol) -> decltype(std::declval<Sol>().get_tasks(0, 0), RunType()) {
		using task_t = typename std::remove_reference<decltype(sol.get_tasks(0, 0)[0])>::type;
		ThreadPool<task_t> pool(sol.get_tasks(n_min, n_max));
		pool.run(thread_num, [&] (const task_t &task, int) {
			sol.calc_custom(task, pool.lock);
		}, [] (int, int) { return true; });
		auto tmp = sol.get_result(n_min, n_max);
		result = std::vector<cpp_int>(tmp.begin(), tmp.end());
		return RunType::CUSTOM_THREADED;
	}
	
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
		auto run_type = run_multithread(opt.n_min, opt.n_max, sol);
		finalize<Sol>(opt.n_min, opt.n_max, sol);
		auto t1 = Timer::get();
		
		std::string time_str = Timer::diff_str(t0, t1);
		std::string date_str = Timer::get_date_str();
		if (date_str.size()) date_str = date_str.substr(0, date_str.size() - 1); // ignore the last linebreak
		std::string thread_str = std::to_string(run_type != RunType::ALL ? thread_num : 1) + " threads";
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

