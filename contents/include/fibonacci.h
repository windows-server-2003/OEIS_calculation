#pragma once
#include "types.h"
#include <type_traits>
#include <cassert>
#include <vector>

template<typename T> std::vector<T> fibonacci_sequence(int n) {
	assert(n >= 1);
	std::vector<T> res(n + 1);
	res[0] = 1;
	res[1] = 2;
	for (int i = 2; i <= n; i++) res[i] = res[i - 1] + res[i - 2];
	return res;
}
template<typename T> std::vector<T> circular_fibonacci_sequence(int n, const std::vector<T> &fib) {
	assert(n >= 1);
	if (n == 1) return {1, 1};
	std::vector<T> res{1, 1, 3};
	for (int i = 3; i <= n; i++) res.push_back(fib[i - 1] + fib[i - 3]);
	return res;
}

template<typename T, typename Func> void for_each_fibonacci(int n, const Func &func) {
	static_assert(std::is_unsigned<T>::value);
	assert(0 <= n && n < (int) sizeof(T) * 8);
	if (n <= 3) {
		func(0, 0);
		if (n >= 1) func(1, 1);
		if (n >= 2) func(2, 2);
		if (n >= 3) func(3, 4), func(4, 5);
		return;
	}
	T x = 0;
	for (size_t i = 0; x < T(1) << n; ) {
		auto call = [&] (int j) { func(i++, x | j); };
		call(0); call(1); call(2); call(4); call(5);
		if (x & 16) x |= 5;
		else call(8), call(9), call(10), x |= 10;
		T y = x | x >> 1;
		y = ~y & (y + 1);
		x = (x & ~(y - 1)) | y;
	}
}
template<typename T, typename Func0, typename Func1> void for_each_fibonacci_mod2(int n, const Func0 &func0, const Func1 &func1) {
	static_assert(std::is_unsigned<T>::value);
	assert(0 <= n && n < (int) sizeof(T) * 8);
	if (n <= 4) {
		func0(0, 0);
		if (n >= 1) func1(1, 1);
		if (n >= 2) func0(2, 1);
		if (n >= 3) func0(3, 2), func1(4, 3);
		if (n >= 4) func0(5, 3), func1(6, 4), func0(7, 4);
		return;
	}
	T x = 0;
	size_t sum = 0;
	for (size_t i = 0; x < T(1) << n; ) {
		auto call0 = [&] () { func0(i++, sum++); };
		auto call1 = [&] () { func1(i++, sum); };
		call0(); call1(); call0(); call0(); call1(); call0(); call1(); call0();
		if (x & 32) x |= 10;
		else call0(), call1(), call0(), call0(), call1(), x |= 21;
		T y = x | x >> 1;
		y = ~y & (y + 1);
		x = (x & ~(y - 1)) | y;
	}
}


template<typename T, typename Func> void for_each_circular_fibonacci(int n, const Func &func) {
	if (n <= 1) {
		func(0, 0);
		return;
	}
	for_each_fibonacci<T>(n - 1, func);
	size_t offset = fibonacci_sequence<size_t>(n - 1).back();
	T set_bit = T(1) << (n - 1);
	for_each_fibonacci<T>(n - 3, [&] (size_t i, T t) { func(offset + i, set_bit | t << 1); });
}

template<int K, int L> struct fibonacci_decoder {
	u64 offset[L][1 << K];
	
	fibonacci_decoder () {
		auto fib = fibonacci_sequence<u64>(K * L);
		for (int i = 0; i < L; i++) {
			u64 cur_offset = 0;
			for_each_fibonacci<u64>(K, [&] (size_t, u64 t) {
				offset[i][t] = cur_offset;
				cur_offset += fib[std::max(0, i * K - (int)(t & 1))];
			});
		}
	}
	
	u64 get_fibonacci_index(u64 x) {
		u64 res = 0;
		for (int i = 0; i < L; i++) res += offset[i][x >> (i * K) & ((1 << K) - 1)];
		return res;
	}
};

