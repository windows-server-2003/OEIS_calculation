#pragma once
#include <vector>
#include <cassert>
#include "types.h"
#include "fibonacci.h"

struct no_adjacent_transformer {
	std::vector<u64> fib;
	
	no_adjacent_transformer (int n) : fib(fibonacci_sequence<u64>(n)) {}
	
	template<typename val_t> void transform_one_square_first(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		u64 dst_head0 = 0;
		u64 dst_head1 = fib[n - 1];
		u64 src_head = 0;
		for_each_fibonacci<u64>(n - 2, [&] (size_t, u64 x) {
			if (!(x & 1)) {
				auto &d0 = next[dst_head0];
				auto &d1 = next[dst_head1];
				auto &d2 = next[dst_head0 + 1];
				
				d0 = dp[src_head + 0] + dp[src_head + 1];
				d1 = dp[src_head + 0] + dp[src_head + 2];
				d2 = dp[src_head + 2];
				dst_head0 += 2;
				dst_head1 += 1;
				src_head += 3;
			} else {
				auto &d0 = next[dst_head0];
				auto &d1 = next[dst_head1];
				d0 = dp[src_head + 0] + dp[src_head + 1];
				d1 = dp[src_head + 0];
				dst_head0 += 1;
				dst_head1 += 1;
				src_head += 2;
			}
		});
		assert(dst_head0 == fib[n - 1]);
		assert(dst_head1 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transform_one_square_last(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		u64 dst_head0 = 0;
		u64 dst_head1 = fib[n - 1];
		u64 src_head = 0;
		for (u64 i = 0; i < fib[n - 2]; i++) {
			next[dst_head0] = dp[src_head + 0] + dp[src_head + 1];
			next[dst_head1] = dp[src_head + 0];
			dst_head0 += 1;
			dst_head1 += 1;
			src_head += 2;
		}
		while (src_head < fib[n]) next[dst_head0++] = dp[src_head++];
		assert(dst_head0 == fib[n - 1]);
		assert(dst_head1 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transform_two_squares(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		u64 dst_head00 = 0;
		u64 dst_head10 = fib[column] * fib[n - 2 - column];
		u64 dst_head01 = fib[column + 1] * fib[n - 2 - column];
		u64 src_head = 0;
		{
			for (u64 i = 0; i < fib[column - 1]; i++) {
				for_each_fibonacci<u64>(n - column - 3, [&] (size_t, u64 right) {
					if (right & 1) {
						auto &d0 = next[dst_head00];
						auto &d1 = next[dst_head10];
						auto &d2 = next[dst_head01];
						d2 = dp[src_head + 0] + dp[src_head + 1] ;
						d1 = dp[src_head + 0] + dp[src_head + 2];
						d0 = d1 + dp[src_head + 1];
						dst_head00 += 1;
						dst_head10 += 1;
						dst_head01 += 1;
						src_head += 3;
					} else {
						auto &d0 = next[dst_head00];
						auto &d1 = next[dst_head10];
						auto &d2 = next[dst_head01];
						auto &d3 = next[dst_head00 + 1];
						auto &d4 = next[dst_head10 + 1];
						d4 = dp[src_head + 3];
						d3 = dp[src_head + 3] + dp[src_head + 4];
						d2 = dp[src_head + 0] + dp[src_head + 1] + d3;
						d1 = dp[src_head + 0] + dp[src_head + 2];
						d0 = d1 + dp[src_head + 1];
						dst_head00 += 2;
						dst_head10 += 2;
						dst_head01 += 1;
						src_head += 5;
					}
				});
			}
		}
		{
			for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
				for_each_fibonacci<u64>(n - column - 3, [&] (size_t, u64 right) {
					if (right & 1) {
						auto &d0 = next[dst_head00];
						auto &d1 = next[dst_head01];
						d1 = dp[src_head + 0];
						d0 = dp[src_head + 0] + dp[src_head + 1];
						dst_head00 += 1;
						dst_head01 += 1;
						src_head += 2;
					} else {
						auto &d0 = next[dst_head00];
						auto &d1 = next[dst_head01];
						auto &d2 = next[dst_head00 + 1];
						d2 = dp[src_head + 2];
						d1 = dp[src_head + 0] + dp[src_head + 2];
						d0 = dp[src_head + 0] + dp[src_head + 1];
						dst_head00 += 2;
						dst_head01 += 1;
						src_head += 3;
					}
				});
			}
			
		}
		assert(dst_head00 == fib[column] * fib[n - 2 - column]);
		assert(dst_head10 == fib[column + 1] * fib[n - 2 - column]);
		assert(dst_head01 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transform_three_squares(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		u64 dst_head000 = 0;
		u64 dst_head100 = fib[column] * fib[n - 3 - column];
		u64 dst_head010 = fib[column + 1] * fib[n - 3 - column];
		u64 dst_head001 = fib[column + 2] * fib[n - 3 - column];
		u64 dst_head101 = fib[column + 2] * fib[n - 3 - column] + fib[column] * fib[n - 4 - column];
		u64 src_head = 0;
		{
			for (u64 i = 0; i < fib[column - 1]; i++) {
				for_each_fibonacci<u64>(n - column - 4, [&] (size_t, u64 right) {
					if (right & 1) {
						auto &d0 = next[dst_head000];
						auto &d1 = next[dst_head100];
						auto &d2 = next[dst_head010];
						auto &d3 = next[dst_head001];
						auto &d4 = next[dst_head101];
						
						d4 = dp[src_head + 0] + dp[src_head + 2];
						d3 = d4 + dp[src_head + 1];
						val_t tmp = d3 + dp[src_head + 3];
						d1 = tmp - dp[src_head + 1];
						tmp += dp[src_head + 4];
						d0 = tmp;
						d2 = tmp - dp[src_head + 2];
						
						src_head += 5;
						dst_head000 += 1;
						dst_head100 += 1;
						dst_head010 += 1;
						dst_head001 += 1;
						dst_head101 += 1;
					} else {
						auto &d0 = next[dst_head000];
						auto &d1 = next[dst_head100];
						auto &d2 = next[dst_head010];
						auto &d3 = next[dst_head001];
						auto &d4 = next[dst_head101];
						auto &d5 = next[dst_head000 + 1];
						auto &d6 = next[dst_head100 + 1];
						auto &d7 = next[dst_head010 + 1];
						d7 = dp[src_head + 5] + dp[src_head + 6];
						d6 = dp[src_head + 5] + dp[src_head + 7];
						d5 = d7 + dp[src_head + 7];
						
						val_t tmp = dp[src_head + 0] + dp[src_head + 2];
						d4 = tmp + d6;
						tmp += dp[src_head + 1];
						d3 = tmp + d5;
						tmp += dp[src_head + 3];
						d1 = tmp - dp[src_head + 1];
						tmp += dp[src_head + 4];
						d0 = tmp;
						d2 = tmp - dp[src_head + 2];
						
						src_head += 8;
						dst_head000 += 2;
						dst_head100 += 2;
						dst_head010 += 2;
						dst_head001 += 1;
						dst_head101 += 1;
					}
				});
			}
		}
		{
			for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
				for_each_fibonacci<u64>(n - column - 4, [&] (size_t, u64 right) {
					if (right & 1) {
						auto &d0 = next[dst_head000];
						auto &d1 = next[dst_head010];
						auto &d2 = next[dst_head001];
						
						d2 = dp[src_head + 0] + dp[src_head + 1];
						d1 = dp[src_head + 0] + dp[src_head + 2];
						d0 = d2 + dp[src_head + 2];
						src_head += 3;
						dst_head000 += 1;
						dst_head010 += 1;
						dst_head001 += 1;
					} else {
						auto &d0 = next[dst_head000];
						auto &d1 = next[dst_head010];
						auto &d2 = next[dst_head001];
						auto &d3 = next[dst_head000 + 1];
						auto &d4 = next[dst_head010 + 1];
						
						d3 = dp[src_head + 3] + dp[src_head + 4];
						d4 = dp[src_head + 3];
						val_t tmp = dp[src_head + 0] + dp[src_head + 1];
						d2 = tmp + d3;
						d1 = dp[src_head + 0] + dp[src_head + 2];
						d0 = tmp + dp[src_head + 2];
						src_head += 5;
						dst_head000 += 2;
						dst_head010 += 2;
						dst_head001 += 1;
					}
				});
			}
			
		}
		assert(dst_head000 == fib[column] * fib[n - 3 - column]);
		assert(dst_head100 == fib[column + 1] * fib[n - 3 - column]);
		assert(dst_head010 == fib[column + 2] * fib[n - 3 - column]);
		assert(dst_head001 == fib[column + 2] * fib[n - 3 - column] + fib[column] * fib[n - 4 - column]);
		assert(dst_head101 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transform_four_squares(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 4;
		u64 dst_head0000 = 0;
		u64 dst_head1000 = fib[column] * fib[n - L - column];
		u64 dst_head0100 = fib[column + 1] * fib[n - L - column];
		u64 dst_head0010 = fib[column + 2] * fib[n - L - column];
		u64 dst_head1010 = (fib[column] + fib[column + 2]) * fib[n - L - column];
		u64 dst_head0001 = fib[column + 3] * fib[n - L - column];
		u64 dst_head1001 = fib[column + 3] * fib[n - L - column] + fib[column] * fib[n - L - 1 - column];
		u64 dst_head0101 = fib[column + 3] * fib[n - L - column] + fib[column + 1] * fib[n - L - 1 - column];
		u64 src_head = 0;
		{
			for (u64 i = 0; i < fib[column - 1]; i++) {
				for_each_fibonacci_mod2<u64>(n - column - L - 1, [&] (size_t j, size_t j2) {
						j2 += j;
						auto &d0 = next[dst_head0000 + j2];
						auto &d1 = next[dst_head1000 + j2];
						auto &d2 = next[dst_head0100 + j2];
						auto &d3 = next[dst_head0010 + j2];
						auto &d4 = next[dst_head1010 + j2];
						auto &d5 = next[dst_head0001 + j];
						auto &d6 = next[dst_head1001 + j];
						auto &d7 = next[dst_head0101 + j];
						auto &d8 = next[dst_head0000 + j2 + 1];
						auto &d9 = next[dst_head1000 + j2 + 1];
						auto &d10 = next[dst_head0100 + j2 + 1];
						auto &d11 = next[dst_head0010 + j2 + 1];
						auto &d12 = next[dst_head1010 + j2 + 1];
						d12 = dp[src_head + 8] + dp[src_head + 10];
						d9 = d12 + dp[src_head + 11];
						d11 = d12 + dp[src_head + 9];
						d8 = d9 + dp[src_head + 9] + dp[src_head + 12];
						d10 = d8 - dp[src_head + 10];
						
						auto tmp = dp[src_head + 0] + dp[src_head + 2];
						d4 = tmp + dp[src_head + 5] + dp[src_head + 7];
						tmp += dp[src_head + 3];
						d6 = tmp + d9;
						d3 = d4 + dp[src_head + 1] + dp[src_head + 6];
						d1 = d4 + dp[src_head + 3];
						tmp += dp[src_head + 1] + dp[src_head + 4];
						d5 = tmp + d8;
						tmp -= dp[src_head + 2];
						d7 = tmp + d10;
						d2 = tmp + dp[src_head + 5] + dp[src_head + 6];
						d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
						
						src_head += 13;
					}, [&] (size_t j, size_t j2) {
						j2 += j;
						auto &d0 = next[dst_head0000 + j2];
						auto &d1 = next[dst_head1000 + j2];
						auto &d2 = next[dst_head0100 + j2];
						auto &d3 = next[dst_head0010 + j2];
						auto &d4 = next[dst_head1010 + j2];
						auto &d5 = next[dst_head0001 + j];
						auto &d6 = next[dst_head1001 + j];
						auto &d7 = next[dst_head0101 + j];
						auto tmp = dp[src_head + 0] + dp[src_head + 2];
						d6 = tmp + dp[src_head + 3];
						d4 = tmp + dp[src_head + 5] + dp[src_head + 7];
						d3 = d4 + dp[src_head + 1] + dp[src_head + 6];
						d1 = d4 + dp[src_head + 3];
						d5 = d6 + dp[src_head + 1] + dp[src_head + 4];
						d7 = d5 - dp[src_head + 2];
						d2 = d7 + dp[src_head + 5] + dp[src_head + 6];
						d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
						
						src_head += 8;
					}
				);
				dst_head0001 += fib[n - column - L - 1];
				dst_head1001 += fib[n - column - L - 1];
				dst_head0101 += fib[n - column - L - 1];
				dst_head0000 += fib[n - column - L];
				dst_head1000 += fib[n - column - L];
				dst_head0100 += fib[n - column - L];
				dst_head0010 += fib[n - column - L];
				dst_head1010 += fib[n - column - L];
			}
		}
		{
			for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
				for_each_fibonacci_mod2<u64>(n - column - L - 1, [&] (size_t j, size_t j2) {
					j2 += j;
					auto &d0 = next[dst_head0000 + j2];
					auto &d1 = next[dst_head0100 + j2];
					auto &d2 = next[dst_head0010 + j2];
					auto &d3 = next[dst_head0001 + j];
					auto &d4 = next[dst_head0101 + j];
					auto &d5 = next[dst_head0000 + j2 + 1];
					auto &d6 = next[dst_head0100 + j2 + 1];
					auto &d7 = next[dst_head0010 + j2 + 1];
					d7 = dp[src_head + 5] + dp[src_head + 6];
					d6 = dp[src_head + 5] + dp[src_head + 7];
					d5 = d7 + dp[src_head + 7];
					
					auto tmp = dp[src_head + 0] + dp[src_head + 2];
					d4 = tmp + d6;
					d1 = tmp + dp[src_head + 3];
					tmp += dp[src_head + 1];
					d3 = tmp + d5;
					d0 = tmp + dp[src_head + 3] + dp[src_head + 4];
					d2 = d0 - dp[src_head + 2];
					
					src_head += 8;
				}, [&] (size_t j, size_t j2) {
					j2 += j;
					auto &d0 = next[dst_head0000 + j2];
					auto &d1 = next[dst_head0100 + j2];
					auto &d2 = next[dst_head0010 + j2];
					auto &d3 = next[dst_head0001 + j];
					auto &d4 = next[dst_head0101 + j];
					
					d4 = dp[src_head + 0] + dp[src_head + 2];
					d1 = d4 + dp[src_head + 3];
					d3 = d4 + dp[src_head + 1];
					d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
					d2 = d0 - dp[src_head + 2];
					
					src_head += 5;
				});
				dst_head0000 += fib[n - column - L];
				dst_head0100 += fib[n - column - L];
				dst_head0010 += fib[n - column - L];
				dst_head0001 += fib[n - column - L - 1];
				dst_head0101 += fib[n - column - L - 1];
			}
			
		}
		assert(dst_head0000 == fib[column] * fib[n - L - column]);
		assert(dst_head1000 == fib[column + 1] * fib[n - L - column]);
		assert(dst_head0100 == fib[column + 2] * fib[n - L - column]);
		assert(dst_head0010 == (fib[column] + fib[column + 2]) * fib[n - L - column]);
		assert(dst_head1010 == fib[column + 3] * fib[n - L - column]);
		assert(dst_head0001 == fib[column + 3] * fib[n - L - column] + fib[column] * fib[n - L - 1 - column]);
		assert(dst_head1001 == fib[column + 3] * fib[n - L - column] + fib[column + 1] * fib[n - L - 1 - column]);
		assert(dst_head0101 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transpose4(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int k = 4;
		constexpr u64 fib[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073, 4807526976, 7778742049, 12586269025,};
		u64 head[fib[k]], head_last[fib[k]];
		head[0] = 0;
		for_each_fibonacci<u64>(k, [&] (size_t i, u64 t) {
			if (i + 1 < fib[k]) {
				head[i + 1] = head[i] + fib[column - (t & 1)] * fib[n - k - column - (t >> (k - 1) & 1)];
				head_last[i] = head[i + 1];
			}
		});
		head_last[fib[k] - 1] = fib[n];
		
		u64 src_head = 0;
		
		u64 r0 = 0;
		u64 r1 = 0;
		for (u64 i = 0; i < fib[column - 1]; i++) {
			for_each_fibonacci<u64>(n - column - k, [&] (size_t, u64 right) {
				if (right & 1) {
					for (u64 j = 0; j < fib[k - 1]; j++) next[head[j] + r0] = dp[src_head + j];
					src_head += fib[k - 1];
					r0++;
				} else {
					for (u64 j = 0; j < fib[k - 1]; j++) next[head[j] + r0] = dp[src_head + j];
					for (u64 j = fib[k - 1]; j < fib[k]; j++) next[head[j] + r1] = dp[src_head + j];
					src_head += fib[k];
					r0++;
					r1++;
				}
			});
		}
		for (u64 j = 0; j < fib[k - 1]; j++) head[j] += r0;
		for (u64 j = fib[k - 1]; j < fib[k]; j++) head[j] += r1;
		for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
			for_each_fibonacci<u64>(n - column - k, [&] (size_t, u64 right) {
				if (right & 1) {
					for_each_fibonacci<u64>(k - 1, [&] (size_t j, u64 t) {
						if (!(t & 1)) next[head[j]++] = dp[src_head++];
					});
				} else {
					for_each_fibonacci<u64>(k, [&] (size_t j, u64 t) {
						if (!(t & 1)) next[head[j]++] = dp[src_head++];
					});
				}
			});
		}
		assert(src_head == fib[n]);
		for (u64 i = 0; i < fib[k]; i++) assert(head[i] == head_last[i]);
	}
	template<typename val_t, int k> void transpose(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr u64 fib[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073, 4807526976, 7778742049, 12586269025,};
		u64 head[fib[k]], head_last[fib[k]];
		head[0] = 0;
		for_each_fibonacci<u64>(k, [&] (size_t i, u64 t) {
			if (i + 1 < fib[k]) {
				head[i + 1] = head[i] + fib[column - (t & 1)] * fib[n - k - column - (t >> (k - 1) & 1)];
				head_last[i] = head[i + 1];
			}
		});
		head_last[fib[k] - 1] = fib[n];
		
		u64 src_head = 0;
		
		u64 r0 = 0;
		u64 r1 = 0;
		for (u64 i = 0; i < fib[column - 1]; i++) {
			for_each_fibonacci<u64>(n - column - k, [&] (size_t, u64 right) {
				if (right & 1) {
					for (u64 j = 0; j < fib[k - 1]; j++) next[head[j] + r0] = dp[src_head + j];
					src_head += fib[k - 1];
					r0++;
				} else {
					for (u64 j = 0; j < fib[k - 1]; j++) next[head[j] + r0] = dp[src_head + j];
					for (u64 j = fib[k - 1]; j < fib[k]; j++) next[head[j] + r1] = dp[src_head + j];
					src_head += fib[k];
					r0++;
					r1++;
				}
			});
		}
		for (u64 j = 0; j < fib[k - 1]; j++) head[j] += r0;
		for (u64 j = fib[k - 1]; j < fib[k]; j++) head[j] += r1;
		for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
			for_each_fibonacci<u64>(n - column - k, [&] (size_t, u64 right) {
				if (right & 1) {
					for_each_fibonacci<u64>(k - 1, [&] (size_t j, u64 t) {
						if (!(t & 1)) next[head[j]++] = dp[src_head++];
					});
				} else {
					for_each_fibonacci<u64>(k, [&] (size_t j, u64 t) {
						if (!(t & 1)) next[head[j]++] = dp[src_head++];
					});
				}
			});
		}
		assert(src_head == fib[n]);
		for (u64 i = 0; i < fib[k]; i++) assert(head[i] == head_last[i]);
	}
	template<typename val_t> void transpose_one(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int k = 1;
		constexpr u64 fib[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073, 4807526976, 7778742049, 12586269025,};
		u64 head[fib[k]], head_last[fib[k]];
		head[0] = 0;
		for_each_fibonacci<u64>(k, [&] (size_t i, u64 t) {
			if (i + 1 < fib[k]) {
				head[i + 1] = head[i] + fib[column - (t & 1)] * fib[n - k - column - (t >> (k - 1) & 1)];
				head_last[i] = head[i + 1];
			}
		});
		head_last[fib[k] - 1] = fib[n];
		
		u64 src_head = 0;
		
		for (u64 i = 0; i < fib[column - 1]; i++) {
			for_each_fibonacci<u64>(n - column, [&] (size_t i, u64 right) {
				next[head[right & 1]++] = dp[src_head + i];
			});
			src_head += fib[n - column];
		}
		memcpy(&next[head[0]], &dp[src_head], fib[std::max(0, column - 2)] * fib[n - column - k] * sizeof(val_t));
		head[0] += fib[std::max(0, column - 2)] * fib[n - column - k];
		src_head += fib[std::max(0, column - 2)] * fib[n - column - k];
		
		assert(src_head == fib[n]);
		for (u64 i = 0; i < fib[k]; i++) assert(head[i] == head_last[i]);
	}
	template<typename val_t> void transform_five_squares(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 5;
		u64 dst_head00000 = 0;
		u64 dst_head10000 = fib[column] * fib[n - L - column];
		u64 dst_head01000 = fib[column + 1] * fib[n - L - column];
		u64 dst_head00100 = fib[column + 2] * fib[n - L - column];
		u64 dst_head10100 = (fib[column] + fib[column + 2]) * fib[n - L - column];
		u64 dst_head00010 = fib[column + 3] * fib[n - L - column];
		u64 dst_head10010 = (fib[column] + fib[column + 3]) * fib[n - L - column];
		u64 dst_head01010 = (fib[column + 1] + fib[column + 3]) * fib[n - L - column];
		u64 dst_head00001 = fib[column + 4] * fib[n - L - column];
		u64 dst_head10001 = fib[column + 4] * fib[n - L - column] + fib[column] * fib[n - L - 1 - column];
		u64 dst_head01001 = fib[column + 4] * fib[n - L - column] + fib[column + 1] * fib[n - L - 1 - column];
		u64 dst_head00101 = fib[column + 4] * fib[n - L - column] + fib[column + 2] * fib[n - L - 1 - column];
		u64 dst_head10101 = fib[column + 4] * fib[n - L - column] + (fib[column] + fib[column + 2]) * fib[n - L - 1 - column];
		
		u64 src_head = 0;
		{
			for (u64 i = 0; i < fib[std::max(0, column - 1)]; i++) {
				for_each_fibonacci_mod2<u64>(n - column - L - 1, [&] (size_t j, size_t j2) {
					j2 += j;
					
					auto &d0 = next[dst_head00000 + j2];
					auto &d1 = next[dst_head10000 + j2];
					auto &d2 = next[dst_head01000 + j2];
					auto &d3 = next[dst_head00100 + j2];
					auto &d4 = next[dst_head10100 + j2];
					auto &d5 = next[dst_head00010 + j2];
					auto &d6 = next[dst_head10010 + j2];
					auto &d7 = next[dst_head01010 + j2];
					auto &d8 = next[dst_head00001 + j];
					auto &d9 = next[dst_head10001 + j];
					auto &d10 = next[dst_head01001 + j];
					auto &d11 = next[dst_head00101 + j];
					auto &d12 = next[dst_head10101 + j];
					auto &d13 = next[dst_head00000 + j2 + 1];
					auto &d14 = next[dst_head10000 + j2 + 1];
					auto &d15 = next[dst_head01000 + j2 + 1];
					auto &d16 = next[dst_head00100 + j2 + 1];
					auto &d17 = next[dst_head10100 + j2 + 1];
					auto &d18 = next[dst_head00010 + j2 + 1];
					auto &d19 = next[dst_head10010 + j2 + 1];
					auto &d20 = next[dst_head01010 + j2 + 1];
					
					auto tmp = dp[src_head + 13] + dp[src_head + 15];
					d19 = tmp + dp[src_head + 16];
					d17 = tmp + dp[src_head + 18] + dp[src_head + 20];
					d14 = d17 + dp[src_head + 16];
					d18 = d19 + dp[src_head + 14] + dp[src_head + 17];
					d20 = d18 - dp[src_head + 15];
					d15 = d20 + dp[src_head + 18] + dp[src_head + 19];
					d16 = d17 + dp[src_head + 14] + dp[src_head + 19];
					d13 = d16 + dp[src_head + 16] + dp[src_head + 17];
					
					
					auto l02 = dp[src_head + 0] + dp[src_head + 2];
					auto r02 = dp[src_head + 8] + dp[src_head + 10];
					auto r023 = r02 + dp[src_head + 11];
					auto l0257 = l02 + dp[src_head + 5] + dp[src_head + 7];
					d12 = l0257 + d17;
					d4 = l0257 + r02;
					auto l02357 = l0257 + dp[src_head + 3];
					d1 = l02357 + r023;
					d9 = l02357 + d14;
					auto l33 = l0257 + dp[src_head + 1] + dp[src_head + 6];
					auto r3 = r02 + dp[src_head + 9];
					d3 = l33 + r3;
					d11 = l33 + d16;
					auto l023 = l02 + dp[src_head + 3];
					d6 = l023 + r023;
					auto r5 = r023 + dp[src_head + 9] + dp[src_head + 12];
					auto l5 = l023 + dp[src_head + 1] + dp[src_head + 4];
					d5 = l5 + r5;
					auto l5_1 = l5 - dp[src_head + 2];
					auto l8 = l33 + dp[src_head + 3] + dp[src_head + 4];
					d8 = l8 + d13;
					d0 = l8 + r5;
					auto r22 = r5 - dp[src_head + 10];
					d7 = l5_1 + r22;
					auto l24 = l5_1 + dp[src_head + 5] + dp[src_head + 6];
					d10 = l24 + d15;
					d2 = l24 + r22;
					
					src_head += 21;
				}, [&] (size_t j, size_t j2) {
					j2 += j;
					auto &d0 = next[dst_head00000 + j2];
					auto &d1 = next[dst_head10000 + j2];
					auto &d2 = next[dst_head01000 + j2];
					auto &d3 = next[dst_head00100 + j2];
					auto &d4 = next[dst_head10100 + j2];
					auto &d5 = next[dst_head00010 + j2];
					auto &d6 = next[dst_head10010 + j2];
					auto &d7 = next[dst_head01010 + j2];
					auto &d8 = next[dst_head00001 + j];
					auto &d9 = next[dst_head10001 + j];
					auto &d10 = next[dst_head01001 + j];
					auto &d11 = next[dst_head00101 + j];
					auto &d12 = next[dst_head10101 + j];
					
					auto l02 = dp[src_head + 0] + dp[src_head + 2];
					auto r02 = dp[src_head + 8] + dp[src_head + 10];
					auto r023 = r02 + dp[src_head + 11];
					auto l0257 = l02 + dp[src_head + 5] + dp[src_head + 7];
					d12 = l0257;
					d4 = l0257 + r02;
					auto l02357 = l0257 + dp[src_head + 3];
					d1 = l02357 + r023;
					d9 = l02357;
					auto l33 = l0257 + dp[src_head + 1] + dp[src_head + 6];
					auto r3 = r02 + dp[src_head + 9];
					d3 = l33 + r3;
					d11 = l33;
					auto l023 = l02 + dp[src_head + 3];
					d6 = l023 + r023;
					auto r5 = r023 + dp[src_head + 9] + dp[src_head + 12];
					auto l5 = l023 + dp[src_head + 1] + dp[src_head + 4];
					d5 = l5 + r5;
					auto l5_1 = l5 - dp[src_head + 2];
					auto l8 = l33 + dp[src_head + 3] + dp[src_head + 4];
					d8 = l8;
					d0 = l8 + r5;
					auto r22 = r5 - dp[src_head + 10];
					d7 = l5_1 + r22;
					auto l24 = l5_1 + dp[src_head + 5] + dp[src_head + 6];
					d10 = l24;
					d2 = l24 + r22;
					
					src_head += 13;
				});
				dst_head00000 += fib[n - L - column];
				dst_head10000 += fib[n - L - column];
				dst_head01000 += fib[n - L - column];
				dst_head00100 += fib[n - L - column];
				dst_head10100 += fib[n - L - column];
				dst_head00010 += fib[n - L - column];
				dst_head10010 += fib[n - L - column];
				dst_head01010 += fib[n - L - column];
				dst_head00001 += fib[n - L - column - 1];
				dst_head10001 += fib[n - L - column - 1];
				dst_head01001 += fib[n - L - column - 1];
				dst_head00101 += fib[n - L - column - 1];
				dst_head10101 += fib[n - L - column - 1];
			}
		}
		if (column) {
			for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
				for_each_fibonacci_mod2<u64>(n - column - L - 1, [&] (size_t j, size_t j2) {
					j2 += j;
					auto &d0 = next[dst_head00000 + j2];
					auto &d1 = next[dst_head01000 + j2];
					auto &d2 = next[dst_head00100 + j2];
					auto &d3 = next[dst_head00010 + j2];
					auto &d4 = next[dst_head01010 + j2];
					auto &d5 = next[dst_head00001 + j];
					auto &d6 = next[dst_head01001 + j];
					auto &d7 = next[dst_head00101 + j];
					auto &d8 = next[dst_head00000 + j2 + 1];
					auto &d9 = next[dst_head01000 + j2 + 1];
					auto &d10 = next[dst_head00100 + j2 + 1];
					auto &d11 = next[dst_head00010 + j2 + 1];
					auto &d12 = next[dst_head01010 + j2 + 1];
					
					d12 = dp[src_head + 8] + dp[src_head + 10];
					d11 = d12 + dp[src_head + 9];
					d9 = d12 + dp[src_head + 11];
					d8 = d11 + dp[src_head + 11] + dp[src_head + 12];
					d10 = d8 - dp[src_head + 10];
					
					auto ltmp = dp[src_head + 0] + dp[src_head + 2];
					auto rtmp = dp[src_head + 5] + dp[src_head + 7];
					auto rall = rtmp + dp[src_head + 6];
					d4 = ltmp + rtmp;
					d3 = ltmp + dp[src_head + 1] + rall;
					ltmp += dp[src_head + 3];
					d6 = ltmp + d9;
					d1 = ltmp + rtmp;
					ltmp += dp[src_head + 1] + dp[src_head + 4];
					d5 = ltmp + d8;
					d0 = ltmp + rall;
					ltmp -= dp[src_head + 2];
					d2 = ltmp + dp[src_head + 5] + dp[src_head + 6];
					d7 = ltmp + d10;
					
					src_head += 13;
				}, [&] (size_t j, size_t j2) {
					j2 += j;
					
					auto &d0 = next[dst_head00000 + j2];
					auto &d1 = next[dst_head01000 + j2];
					auto &d2 = next[dst_head00100 + j2];
					auto &d3 = next[dst_head00010 + j2];
					auto &d4 = next[dst_head01010 + j2];
					auto &d5 = next[dst_head00001 + j];
					auto &d6 = next[dst_head01001 + j];
					auto &d7 = next[dst_head00101 + j];
					
					auto ltmp = dp[src_head + 0] + dp[src_head + 2];
					auto rtmp = dp[src_head + 5] + dp[src_head + 7];
					auto rall = rtmp + dp[src_head + 6];
					d4 = ltmp + rtmp;
					d3 = ltmp + dp[src_head + 1] + rall;
					ltmp += dp[src_head + 3];
					d6 = ltmp;
					d1 = ltmp + rtmp;
					ltmp += dp[src_head + 1] + dp[src_head + 4];
					d5 = ltmp;
					d0 = ltmp + rall;
					ltmp -= dp[src_head + 2];
					d2 = ltmp + dp[src_head + 5] + dp[src_head + 6];
					d7 = ltmp;
					
					src_head += 8;
				});
				dst_head00000 += fib[n - column - L];
				dst_head01000 += fib[n - column - L];
				dst_head00100 += fib[n - column - L];
				dst_head00010 += fib[n - column - L];
				dst_head01010 += fib[n - column - L];
				dst_head00001 += fib[n - column - L - 1];
				dst_head01001 += fib[n - column - L - 1];
				dst_head00101 += fib[n - column - L - 1];
			}
		}
		
		assert(dst_head00000 == fib[column] * fib[n - L - column]);
		assert(dst_head10000 == fib[column + 1] * fib[n - L - column]);
		assert(dst_head01000 == fib[column + 2] * fib[n - L - column]);
		assert(dst_head00100 == (fib[column] + fib[column + 2]) * fib[n - L - column]);
		assert(dst_head10100 == fib[column + 3] * fib[n - L - column]);
		assert(dst_head00010 == (fib[column] + fib[column + 3]) * fib[n - L - column]);
		assert(dst_head10010 == (fib[column + 1] + fib[column + 3]) * fib[n - L - column]);
		assert(dst_head01010 == fib[column + 4] * fib[n - L - column]);
		assert(dst_head00001 == fib[column + 4] * fib[n - L - column] + fib[column] * fib[n - L - 1 - column]);
		assert(dst_head10001 == fib[column + 4] * fib[n - L - column] + fib[column + 1] * fib[n - L - 1 - column]);
		assert(dst_head01001 == fib[column + 4] * fib[n - L - column] + fib[column + 2] * fib[n - L - 1 - column]);
		assert(dst_head00101 == fib[column + 4] * fib[n - L - column] + (fib[column] + fib[column + 2]) * fib[n - L - 1 - column]);
		assert(dst_head10101 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transform_last_five_squares(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 5;
		int column = n - L;
		u64 dst_head00000 = 0;
		u64 dst_head10000 = fib[column] * fib[n - L - column];
		u64 dst_head01000 = fib[column + 1] * fib[n - L - column];
		u64 dst_head00100 = fib[column + 2] * fib[n - L - column];
		u64 dst_head10100 = (fib[column] + fib[column + 2]) * fib[n - L - column];
		u64 dst_head00010 = fib[column + 3] * fib[n - L - column];
		u64 dst_head10010 = (fib[column] + fib[column + 3]) * fib[n - L - column];
		u64 dst_head01010 = (fib[column + 1] + fib[column + 3]) * fib[n - L - column];
		u64 dst_head00001 = fib[column + 4] * fib[n - L - column];
		u64 dst_head10001 = fib[column + 4] * fib[n - L - column] + fib[column] * fib[std::max(n - L - 1 - column, 0)];
		u64 dst_head01001 = fib[column + 4] * fib[n - L - column] + fib[column + 1] * fib[std::max(n - L - 1 - column, 0)];
		u64 dst_head00101 = fib[column + 4] * fib[n - L - column] + fib[column + 2] * fib[std::max(n - L - 1 - column, 0)];
		u64 dst_head10101 = fib[column + 4] * fib[n - L - column] + (fib[column] + fib[column + 2]) * fib[std::max(n - L - 1 - column, 0)];
		
		u64 src_head = 0;
		{
			for (u64 i = 0; i < fib[std::max(0, column - 1)]; i++) {
				auto &d0 = next[dst_head00000 + i];
				auto &d1 = next[dst_head10000 + i];
				auto &d2 = next[dst_head01000 + i];
				auto &d3 = next[dst_head00100 + i];
				auto &d4 = next[dst_head10100 + i];
				auto &d5 = next[dst_head00010 + i];
				auto &d6 = next[dst_head10010 + i];
				auto &d7 = next[dst_head01010 + i];
				auto &d8 = next[dst_head00001 + i];
				auto &d9 = next[dst_head10001 + i];
				auto &d10 = next[dst_head01001 + i];
				auto &d11 = next[dst_head00101 + i];
				auto &d12 = next[dst_head10101 + i];
				
				auto l02 = dp[src_head + 0] + dp[src_head + 2];
				auto r02 = dp[src_head + 8] + dp[src_head + 10];
				auto r023 = r02 + dp[src_head + 11];
				auto l0257 = l02 + dp[src_head + 5] + dp[src_head + 7];
				d12 = l0257;
				d4 = l0257 + r02;
				auto l02357 = l0257 + dp[src_head + 3];
				d1 = l02357 + r023;
				d9 = l02357;
				auto l33 = l0257 + dp[src_head + 1] + dp[src_head + 6];
				auto r3 = r02 + dp[src_head + 9];
				d3 = l33 + r3;
				d11 = l33;
				auto l023 = l02 + dp[src_head + 3];
				d6 = l023 + r023;
				auto r5 = r023 + dp[src_head + 9] + dp[src_head + 12];
				auto l5 = l023 + dp[src_head + 1] + dp[src_head + 4];
				d5 = l5 + r5;
				auto l5_1 = l5 - dp[src_head + 2];
				auto l8 = l33 + dp[src_head + 3] + dp[src_head + 4];
				d8 = l8;
				d0 = l8 + r5;
				auto r22 = r5 - dp[src_head + 10];
				d7 = l5_1 + r22;
				auto l24 = l5_1 + dp[src_head + 5] + dp[src_head + 6];
				d10 = l24;
				d2 = l24 + r22;
				
				src_head += 13;
			}
			dst_head00000 += fib[std::max(0, column - 1)];
			dst_head10000 += fib[std::max(0, column - 1)];
			dst_head01000 += fib[std::max(0, column - 1)];
			dst_head00100 += fib[std::max(0, column - 1)];
			dst_head10100 += fib[std::max(0, column - 1)];
			dst_head00010 += fib[std::max(0, column - 1)];
			dst_head10010 += fib[std::max(0, column - 1)];
			dst_head01010 += fib[std::max(0, column - 1)];
			dst_head00001 += fib[std::max(0, column - 1)];
			dst_head10001 += fib[std::max(0, column - 1)];
			dst_head01001 += fib[std::max(0, column - 1)];
			dst_head00101 += fib[std::max(0, column - 1)];
			dst_head10101 += fib[std::max(0, column - 1)];
		}
		{
			for (u64 i = 0; i < fib[std::max(0, column - 2)]; i++) {
				auto &d0 = next[dst_head00000 + i];
				auto &d1 = next[dst_head01000 + i];
				auto &d2 = next[dst_head00100 + i];
				auto &d3 = next[dst_head00010 + i];
				auto &d4 = next[dst_head01010 + i];
				auto &d5 = next[dst_head00001 + i];
				auto &d6 = next[dst_head01001 + i];
				auto &d7 = next[dst_head00101 + i];
				
				auto ltmp = dp[src_head + 0] + dp[src_head + 2];
				auto rtmp = dp[src_head + 5] + dp[src_head + 7];
				auto rall = rtmp + dp[src_head + 6];
				d4 = ltmp + rtmp;
				d3 = ltmp + dp[src_head + 1] + rall;
				ltmp += dp[src_head + 3];
				d6 = ltmp;
				d1 = ltmp + rtmp;
				ltmp += dp[src_head + 1] + dp[src_head + 4];
				d5 = ltmp;
				d0 = ltmp + rall;
				ltmp -= dp[src_head + 2];
				d2 = ltmp + dp[src_head + 5] + dp[src_head + 6];
				d7 = ltmp;
				
				src_head += 8;
			}
			dst_head00000 += fib[std::max(0, column - 2)];
			dst_head01000 += fib[std::max(0, column - 2)];
			dst_head00100 += fib[std::max(0, column - 2)];
			dst_head00010 += fib[std::max(0, column - 2)];
			dst_head01010 += fib[std::max(0, column - 2)];
			dst_head00001 += fib[std::max(0, column - 2)];
			dst_head01001 += fib[std::max(0, column - 2)];
			dst_head00101 += fib[std::max(0, column - 2)];
		}
		
		assert(dst_head00000 == fib[column] * fib[n - L - column]);
		assert(dst_head10000 == fib[column + 1] * fib[n - L - column]);
		assert(dst_head01000 == fib[column + 2] * fib[n - L - column]);
		assert(dst_head00100 == (fib[column] + fib[column + 2]) * fib[n - L - column]);
		assert(dst_head10100 == fib[column + 3] * fib[n - L - column]);
		assert(dst_head00010 == (fib[column] + fib[column + 3]) * fib[n - L - column]);
		assert(dst_head10010 == (fib[column + 1] + fib[column + 3]) * fib[n - L - column]);
		assert(dst_head01010 == fib[column + 4] * fib[n - L - column]);
		assert(dst_head00001 == fib[column + 4] * fib[n - L - column] + fib[column] * fib[std::max(n - L - 1 - column, 0)]);
		assert(dst_head10001 == fib[column + 4] * fib[n - L - column] + fib[column + 1] * fib[std::max(n - L - 1 - column, 0)]);
		assert(dst_head01001 == fib[column + 4] * fib[n - L - column] + fib[column + 2] * fib[std::max(n - L - 1 - column, 0)]);
		assert(dst_head00101 == fib[column + 4] * fib[n - L - column] + (fib[column] + fib[column + 2]) * fib[std::max(n - L - 1 - column, 0)]);
		assert(dst_head10101 == fib[n]);
		assert(src_head == fib[n]);
	}
	template<typename val_t> void transform_one_square(int n, int column, std::vector<val_t> &dp, std::vector<val_t> &next) {
		u64 dst_head0 = 0;
		u64 dst_head1 = fib[column] * fib[n - 1 - column];
		u64 src_head = 0;
		for (u64 i = 0; i < fib[column - 1]; i++) {
			for_each_fibonacci<u64>(n - column - 2, [&] (size_t, u64 x) {
				if (!(x & 1)) {
					auto &d0 = next[dst_head0];
					auto &d1 = next[dst_head1];
					auto &d2 = next[dst_head0 + 1];
					
					d0 = dp[src_head + 0] + dp[src_head + 1];
					d1 = dp[src_head + 0] + dp[src_head + 2];
					d2 = dp[src_head + 2];
					dst_head0 += 2;
					dst_head1 += 1;
					src_head += 3;
				} else {
					auto &d0 = next[dst_head0];
					auto &d1 = next[dst_head1];
					d0 = dp[src_head + 0] + dp[src_head + 1];
					d1 = dp[src_head + 0];
					dst_head0 += 1;
					dst_head1 += 1;
					src_head += 2;
				}
			});
			assert(dst_head0 <= fib[column] * fib[n - 1 - column]);
			assert(dst_head1 <= fib[n]);
		}
		{
			u64 size = fib[std::max(0, column - 2)] * fib[n - 1 - column];
			for (u64 i = 0; i < size; i++) next[dst_head0 + i] = dp[src_head + i];
			dst_head0 += size;
			src_head += size;
		}
		assert(dst_head0 == fib[column] * fib[n - 1 - column]);
		assert(dst_head1 == fib[n]);
	}
	template<typename val_t> void transform_one_row(int n, std::vector<val_t> &dp) {
		assert(dp.size() == fib[n]);
		
		if (n == 0) return;
		if (n == 1) {
			dp = {dp[0] + dp[1], dp[0]};
			return;
		}
		
		std::vector<val_t> dp_swap(fib[n]);
		for (int i = 0; i < n; ) {
			int advance = 1;
			if (i == 0) {
				if (n > 5) transform_five_squares(n, i, dp, dp_swap), advance = 5;
				else transform_one_square_first(n, dp, dp_swap);
			} else if (i == n - 1) transform_one_square_last(n, dp, dp_swap);
			else if (i > n - 5 || i == n - 6) transform_one_square(n, i, dp, dp_swap);
			else if (i == n - 5) transform_last_five_squares(n, dp, dp_swap), advance = 5;
			else if (i == n - 7) transform_two_squares(n, i, dp, dp_swap), advance = 2;
			else if (i == n - 8) transform_three_squares(n, i, dp, dp_swap), advance = 3;
			else if (i == n - 9) transform_four_squares(n, i, dp, dp_swap), advance = 4;
			else                 transform_five_squares(n, i, dp, dp_swap), advance = 5;
			std::swap(dp, dp_swap);
			i += advance;
		}
	}
};