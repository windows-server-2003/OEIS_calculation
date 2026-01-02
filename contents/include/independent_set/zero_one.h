#pragma once
#include <vector>
#include <cassert>
#include "types.h"
#include "fibonacci.h"

struct no_adjacent_zero_one_transformer {
	std::vector<u64> fib;
	
	no_adjacent_zero_one_transformer (int n) : fib(fibonacci_sequence<u64>(n)) {}
	
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
				d1 = dp[src_head + 0];
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
						d1 = dp[src_head + 0];
						d2 = dp[src_head + 0] + dp[src_head + 1];
						d0 = d2 + dp[src_head + 2];
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
						d1 = dp[src_head + 0];
						d2 = dp[src_head + 0] + dp[src_head + 1];
						d0 = d2 + dp[src_head + 2];
						
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
						d1 = dp[src_head + 0];
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
						
						d4 = dp[src_head + 0];
						d2 = dp[src_head + 0] + dp[src_head + 1];
						d3 = d2 + dp[src_head + 2];
						d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
						d1 = dp[src_head + 0] + dp[src_head + 3];
						
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
						d6 = dp[src_head + 5];
						d7 = dp[src_head + 5] + dp[src_head + 6];
						d5 = d7 + dp[src_head + 7];
						
						d4 = dp[src_head + 0];
						d2 = dp[src_head + 0] + dp[src_head + 1];
						d3 = d2 + dp[src_head + 2];
						d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
						d1 = dp[src_head + 0] + dp[src_head + 3];
						
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
						
						d1 = dp[src_head + 0];
						d2 = dp[src_head + 0] + dp[src_head + 1];
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
						
						d1 = dp[src_head + 0];
						d2 = dp[src_head + 0] + dp[src_head + 1];
						d0 = d2 + dp[src_head + 2];
						
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
						
						d12 = dp[src_head + 8];
						d10 = dp[src_head + 8] + dp[src_head + 9];
						d11 = d10 + dp[src_head + 10];
						d8 = d11 + dp[src_head + 11] + dp[src_head + 12];
						d9 = dp[src_head + 8] + dp[src_head + 11];
						
						d4 = dp[src_head + 0];
						d7 = dp[src_head + 0] + dp[src_head + 1];
						auto tmp = dp[src_head + 5] + dp[src_head + 6];
						d2 = d7 + tmp;
						d6 = dp[src_head + 0] + dp[src_head + 3];
						d1 = d6 + dp[src_head + 5];
						d3 = d7 + dp[src_head + 2];
						d5 = d3 + dp[src_head + 3] + dp[src_head + 4];
						d0 = d5 + tmp + dp[src_head + 7];

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
						
						d4 = dp[src_head + 0];
						d7 = dp[src_head + 0] + dp[src_head + 1];
						auto tmp = dp[src_head + 5] + dp[src_head + 6];
						d2 = d7 + tmp;
						d6 = dp[src_head + 0] + dp[src_head + 3];
						d1 = d6 + dp[src_head + 5];
						d3 = d7 + dp[src_head + 2];
						d5 = d3 + dp[src_head + 3] + dp[src_head + 4];
						d0 = d5 + tmp + dp[src_head + 7];
						
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
					
					d6 = dp[src_head + 5];
					d7 = dp[src_head + 5] + dp[src_head + 6];
					d5 = d7 + dp[src_head + 7];
					
					d4 = dp[src_head + 0];
					d1 = dp[src_head + 0] + dp[src_head + 3];
					d2 = dp[src_head + 0] + dp[src_head + 1];
					d3 = d2 + dp[src_head + 2];
					d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
					
					src_head += 8;
				}, [&] (size_t j, size_t j2) {
					j2 += j;
					auto &d0 = next[dst_head0000 + j2];
					auto &d1 = next[dst_head0100 + j2];
					auto &d2 = next[dst_head0010 + j2];
					auto &d3 = next[dst_head0001 + j];
					auto &d4 = next[dst_head0101 + j];
					
					d4 = dp[src_head + 0];
					d1 = dp[src_head + 0] + dp[src_head + 3];
					d2 = dp[src_head + 0] + dp[src_head + 1];
					d3 = d2 + dp[src_head + 2];
					d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
					
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
					
					d17 = dp[src_head + 13];
					d20 = dp[src_head + 13] + dp[src_head + 14];
					auto tmp = dp[src_head + 18] + dp[src_head + 19];
					d15 = d20 + tmp;
					d19 = dp[src_head + 13] + dp[src_head + 16];
					d14 = d19 + dp[src_head + 18];
					d16 = d20 + dp[src_head + 15];
					d18 = d16 + dp[src_head + 16] + dp[src_head + 17];
					d13 = d18 + tmp + dp[src_head + 20];
					
					d12 = dp[src_head + 0];
					d7 = dp[src_head + 0] + dp[src_head + 1];
					d4 = dp[src_head + 0] + dp[src_head + 8];
					auto mtmp = dp[src_head + 5] + dp[src_head + 6];
					d11 = d7 + dp[src_head + 2];
					d10 = d7 + mtmp;
					d6 = dp[src_head + 0] + dp[src_head + 3];
					d9 = d6 + dp[src_head + 5];
					d5 = d11 + dp[src_head + 3] + dp[src_head + 4];
					d8 = d5 + mtmp + dp[src_head + 7];
					auto rtmp = dp[src_head + 8] + dp[src_head + 9];
					d2 = d10 + rtmp;
					d1 = d9 + dp[src_head + 8] + dp[src_head + 11];
					rtmp += dp[src_head + 10];
					d3 = d11 + rtmp;
					d0 = d8 + rtmp + dp[src_head + 11] + dp[src_head + 12];
					
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
					
					d12 = dp[src_head + 0];
					d7 = dp[src_head + 0] + dp[src_head + 1];
					d4 = dp[src_head + 0] + dp[src_head + 8];
					auto mtmp = dp[src_head + 5] + dp[src_head + 6];
					d11 = d7 + dp[src_head + 2];
					d10 = d7 + mtmp;
					d6 = dp[src_head + 0] + dp[src_head + 3];
					d9 = d6 + dp[src_head + 5];
					d5 = d11 + dp[src_head + 3] + dp[src_head + 4];
					d8 = d5 + mtmp + dp[src_head + 7];
					auto rtmp = dp[src_head + 8] + dp[src_head + 9];
					d2 = d10 + rtmp;
					d1 = d9 + dp[src_head + 8] + dp[src_head + 11];
					rtmp += dp[src_head + 10];
					d3 = d11 + rtmp;
					d0 = d8 + rtmp + dp[src_head + 11] + dp[src_head + 12];
					
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
					
					d12 = dp[src_head + 8];
					d10 = dp[src_head + 8] + dp[src_head + 9];
					d11 = d10 + dp[src_head + 10];
					d8 = d11 + dp[src_head + 11] + dp[src_head + 12];
					d9 = dp[src_head + 8] + dp[src_head + 11];
					
					d4 = dp[src_head + 0];
					d7 = dp[src_head + 0] + dp[src_head + 1];
					auto tmp = dp[src_head + 5] + dp[src_head + 6];
					d2 = d7 + tmp;
					d6 = dp[src_head + 0] + dp[src_head + 3];
					d1 = d6 + dp[src_head + 5];
					d3 = d7 + dp[src_head + 2];
					d5 = d3 + dp[src_head + 3] + dp[src_head + 4];
					d0 = d5 + tmp + dp[src_head + 7];

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
					
					d4 = dp[src_head + 0];
					d7 = dp[src_head + 0] + dp[src_head + 1];
					auto tmp = dp[src_head + 5] + dp[src_head + 6];
					d2 = d7 + tmp;
					d6 = dp[src_head + 0] + dp[src_head + 3];
					d1 = d6 + dp[src_head + 5];
					d3 = d7 + dp[src_head + 2];
					d5 = d3 + dp[src_head + 3] + dp[src_head + 4];
					d0 = d5 + tmp + dp[src_head + 7];
					
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
				
				d12 = dp[src_head + 0];
				d7 = dp[src_head + 0] + dp[src_head + 1];
				d4 = dp[src_head + 0] + dp[src_head + 8];
				auto mtmp = dp[src_head + 5] + dp[src_head + 6];
				d11 = d7 + dp[src_head + 2];
				d10 = d7 + mtmp;
				d6 = dp[src_head + 0] + dp[src_head + 3];
				d9 = d6 + dp[src_head + 5];
				d5 = d11 + dp[src_head + 3] + dp[src_head + 4];
				d8 = d5 + mtmp + dp[src_head + 7];
				auto rtmp = dp[src_head + 8] + dp[src_head + 9];
				d2 = d10 + rtmp;
				d1 = d9 + dp[src_head + 8] + dp[src_head + 11];
				rtmp += dp[src_head + 10];
				d3 = d11 + rtmp;
				d0 = d8 + rtmp + dp[src_head + 11] + dp[src_head + 12];
				
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
				
				d4 = dp[src_head + 0];
				d7 = dp[src_head + 0] + dp[src_head + 1];
				auto tmp = dp[src_head + 5] + dp[src_head + 6];
				d2 = d7 + tmp;
				d6 = dp[src_head + 0] + dp[src_head + 3];
				d1 = d6 + dp[src_head + 5];
				d3 = d7 + dp[src_head + 2];
				d5 = d3 + dp[src_head + 3] + dp[src_head + 4];
				d0 = d5 + tmp + dp[src_head + 7];

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
					d1 = dp[src_head + 0];
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
		}
		{
			u64 size = fib[std::max(0, column - 2)] * fib[n - 1 - column];
			for (u64 i = 0; i < size; i++) next[dst_head0 + i] = dp[src_head + i];
			dst_head0 += size;
			src_head += size;
		}
		assert(src_head == fib[n]);
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
			if (i == 0) transform_one_square_first(n, dp, dp_swap);
			else if (i == n - 1) transform_one_square_last(n, dp, dp_swap);
			else if (i > n - 5 || i == n - 6) transform_one_square(n, i, dp, dp_swap);
			else if (i == n - 5) transform_last_five_squares(n, dp, dp_swap), advance = 5;
			else if (i == n - 7) transform_two_squares(n, i, dp, dp_swap), advance = 2;
			else if (i == n - 8) transform_three_squares(n, i, dp, dp_swap), advance = 3;
			else if (i == n - 9) transform_four_squares(n, i, dp, dp_swap), advance = 4;
			else transform_five_squares(n, i, dp, dp_swap), advance = 5;
			std::swap(dp, dp_swap);
			i += advance;
		}
	}
	template<typename res_t, typename val_t> res_t reverse_match_product(int n, const std::vector<val_t> &a, const std::vector<val_t> &b) {
		assert(a.size() == fib[n]);
		assert(b.size() == fib[n]);
		
		fibonacci_decoder<15, 3> decoder;
		
		int low_bits = (n + 1) / 2;
		int high_bits = n / 2;
		
		auto bit_reverse = [] (u64 x, int n) {
			u64 res = 0;
			for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
			return res;
		};
		
		u64 high_reverse[fib[high_bits]];
		for_each_fibonacci<u64>(high_bits, [&] (size_t i, u64 t) {
			t = bit_reverse(t, high_bits);
			high_reverse[i] = decoder.get_fibonacci_index(t);
		});
		
		u64 high_reverse1[fib[high_bits - 1]];
		for_each_fibonacci<u64>(high_bits - 1, [&] (size_t i, u64 t) {
			t = bit_reverse(t, high_bits - 1);
			high_reverse1[i] = decoder.get_fibonacci_index(t);
		});
		
		res_t res = 0;
		for_each_fibonacci<u64>(low_bits, [&] (size_t i, u64 low) {
			u64 b_head = decoder.get_fibonacci_index(bit_reverse(low, low_bits) << high_bits);
			u64 a_head = i;
			if (low >> (low_bits - 1) & 1) {
				for_each_fibonacci<u64>(high_bits - 1, [&] (size_t j, u64 high) {
					res += (res_t) a[a_head] * b[b_head + high_reverse1[j]];
					a_head += fib[low_bits + 1 - (high & 1)];
				});
			} else {
				for_each_fibonacci<u64>(high_bits, [&] (size_t j, u64 high) {
					res += (res_t) a[a_head] * b[b_head + high_reverse[j]];
					a_head += fib[low_bits - (high & 1)];
				});
			}
		});
		return res;
	}
	template<typename res_t, typename val_t> res_t reverse_match_product2(int n, const std::vector<val_t> &a, const std::vector<val_t> &b) {
		assert(a.size() == fib[n]);
		assert(b.size() == fib[n]);
		
		constexpr int K = 15;
		constexpr int L = 3;
		u64 offset[L][1 << K];
		{
			auto bit_reverse = [] (u64 x, int n) {
				u64 res = 0;
				for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
				return res;
			};
			auto fib = fibonacci_sequence<u64>(K * L);
			for (int i = 0; i < L; i++) {
				if (i * K >= n) continue;
				u64 cur_offset = 0;
				for_each_fibonacci<u64>(std::min(K, n - i * K), [&] (size_t, u64 t) {
					offset[i][bit_reverse(t, std::min(K, n - i * K))] = cur_offset;
					cur_offset += fib[std::max(0, n - (i + 1) * K - (int)(t & 1))];
				});
			}
		}
		
		int stage = (n + K - 1) / K;
		
		res_t res = 0;
		for_each_fibonacci<u64>(n, [&] (size_t i, u64 t) {
			u64 rev_i = 0;
			for (int j = 0; j < stage; j++) rev_i += offset[j][t >> (j * K) & ((1 << K) - 1)];
			res += (res_t) a[i] * b[rev_i];
		});
		return res;
	}
};