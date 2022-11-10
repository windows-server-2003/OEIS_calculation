#pragma once
#include <random>
#include <chrono>
#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <utility>

struct Random {
	std::mt19937 rnd; // intentionally default-constructed for consistent results across executions
	
	int rnd_int(int l, int r) { assert(l <= r); return std::uniform_int_distribution<int64_t>(l, r)(rnd); }
	int64_t rnd_int64(int64_t l, int64_t r) { assert(l <= r); return std::uniform_int_distribution<int64_t>(l, r)(rnd); }
	std::pair<int, int> rnd_pair(int l, int r) {
		int low = rnd_int(l, r), high = rnd_int(l, r);
		if (high < low) std::swap(low, high);
		return {low, high + 1};
	}
	template<typename T> void shuffle(T begin, T end) {
		int n = end - begin;
		for (int i = 1; i < n; i++) std::swap(begin[rnd_int(0, i)], begin[i]);
	}
};

struct Timer {
	using clock_type = decltype(std::chrono::high_resolution_clock::now());
	using duration_type = decltype(clock_type() - clock_type());
	static clock_type get() { return std::chrono::high_resolution_clock::now(); }
	static double ms(duration_type time) {
		return std::chrono::duration_cast<std::chrono::nanoseconds>(time).count() / 1000000.0;
	}
	static double s(duration_type time) {
		return std::chrono::duration_cast<std::chrono::nanoseconds>(time).count() / 1000000000.0;
	}
	static double diff_ms(clock_type start, clock_type end) { return ms(end - start); }
	static double diff_s(clock_type start, clock_type end) { return s(end - start); }
	static std::string get_date_str() {
		auto cur_time = std::chrono::system_clock::to_time_t(get());
		return std::ctime(&cur_time);
	}
	static std::string diff_str(duration_type time) {
		char res[64];
		if (s(time) < 10) snprintf(res, sizeof(res), "%d ms", (int) ms(time));
		else snprintf(res, sizeof(res), "%.3f s", s(time));
		return res;
	}
	static std::string diff_str(clock_type start, clock_type end) { return diff_str(end - start); }
	
	template<typename T> static void measure(const T &func) {
		auto t0 = get();
		func();
		auto t1 = get();
		std::cerr << Timer::diff_ms(t0, t1) << " ms" << std::endl;
	}
};
