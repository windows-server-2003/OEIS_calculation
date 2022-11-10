#pragma once
#include <vector>
#include <cassert>
#include "types.h"
#include "debug.h"
#include "fibonacci.h"
#include "plain.h"
#include "utils.h"

template<typename T> size_t get_circular_compressed_size(int n) {
	size_t res = 1 + !(n & 1);
	for_each_fibonacci<T>(n - 4, [&] (size_t, T i) {
		auto bit_reverse = [&] (T x) {
			T res = 0;
			for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
			return res;
		};
		auto bit_rotate = [&] (T x, int k) { return (x >> k) | (x & (((T) 1 << k) - 1)) << (n - k); };
		
		i = i << 2 | 1;
		
		T tmp = i;
		for (int j = 1; j < n; j++) if (bit_rotate(tmp, j) < i) return;
		T rev = bit_reverse(i);
		if (rev < i) return;
		for (int j = 1; j < n; j++) if (bit_rotate(rev, j) < i) return;
		
		res++;
	});
	return res;
}
template<typename T> std::vector<std::pair<T, int> > get_circular_compress_indices(int n) {
	std::vector<std::pair<T, int> > res;
	
	for_each_circular_fibonacci<T>(n, [&] (size_t index, T i) {
		auto bit_reverse = [&] (T x) {
			T res = 0;
			for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
			return res;
		};
		auto bit_rotate = [&] (T x, int k) { return (x >> k) | (x & (((T) 1 << k) - 1)) << (n - k); };
		
		T tmp = i;
		for (int j = 1; j < n; j++) if (bit_rotate(tmp, j) < i) return;
		T rev = bit_reverse(i);
		if (rev < i) return;
		for (int j = 1; j < n; j++) if (bit_rotate(rev, j) < i) return;
		
		T buf[128];
		int buf_head = 0;
		buf[buf_head++] = tmp;
		for (int j = 1; j < n; j++) {
			auto cur = bit_rotate(tmp, j);
			if (cur != tmp) buf[buf_head++] = cur;
			else break;
		}
		if (std::find(buf, buf + buf_head, rev) == buf + buf_head) buf_head *= 2;
		
		res.push_back({index, buf_head});
	});
	return res;
}
template<typename T> std::vector<std::tuple<T, T, int> > get_circular_compress_indices_and_values(int n) {
	std::vector<std::tuple<T, T, int> > res;
	
	for_each_circular_fibonacci<T>(n, [&] (size_t index, T i) {
		auto bit_reverse = [&] (T x) {
			T res = 0;
			for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
			return res;
		};
		auto bit_rotate = [&] (T x, int k) { return (x >> k) | (x & (((T) 1 << k) - 1)) << (n - k); };
		
		T tmp = i;
		for (int j = 1; j < n; j++) if (bit_rotate(tmp, j) < i) return;
		T rev = bit_reverse(i);
		if (rev < i) return;
		for (int j = 1; j < n; j++) if (bit_rotate(rev, j) < i) return;
		
		T buf[128];
		int buf_head = 0;
		buf[buf_head++] = tmp;
		for (int j = 1; j < n; j++) {
			auto cur = bit_rotate(tmp, j);
			if (cur != tmp) buf[buf_head++] = cur;
			else break;
		}
		if (std::find(buf, buf + buf_head, rev) == buf + buf_head) buf_head *= 2;
		
		res.push_back(std::make_tuple(index, i, buf_head));
	});
	return res;
}
/*
template<typename T> struct no_adjacent_circular_transformer_fast_dp_t {
	std::vector<u64> fib, circular_fib;
	std::vector<T> data;
	T zero;
	T zero_one;
	no_adjacent_circular_transformer_fast_dp_t (int n) : fib(fibonacci_sequence<u64>(n)), circular_fib(circular_fibonacci_sequence<u64>(n, fib)),
		data(fib[n - 4]), zero(0), zero_one(0) { assert(n >= 5); }
	
	void fill_all(T val) {
		for_each_fibonacci<u64>(n - 4, [&] (size_t i, u64 t) {
			data[i] = ((t & 1) ? 1 : 2) * ((t >> (n - 5) & 1) ? 2 : 3);
		});
	}
};
struct no_adjacent_circular_transformer_fast {
	no_adjacent_transformer transformer;
	
	no_adjacent_circular_transformer (int n) : transformer(n) {}
	
	void transform_one_row(no_adjacent_circular_transformer_fast_dp_t &dp) {
		
	}
};*/

struct no_adjacent_circular_transformer_fast {
	int n;
	no_adjacent_transformer transformer;
	std::vector<u64> &fib = transformer.fib;
	std::vector<u64> circular_fib;
	std::vector<std::tuple<u64, u64, int> > compress_indices_and_values;
	
	std::vector<u32> compressed_index;
	
	static constexpr int K = 17;
	static constexpr int L = 3;
	u64 offset[L][1 << K];
	
	u64 get_fibonacci_index(u64 x) {
		u64 res = 0;
		for (int i = 0; i < L; i++) res += offset[i][x >> (i * K) & ((1 << K) - 1)];
		return res;
	}
	u64 get_circular_fib_index(u64 x) {
		if (x >> (n - 1) & 1) return fib[n - 1] + get_fibonacci_index((x ^ (1ULL << (n - 1))) >> 1);
		return get_fibonacci_index(x);
	}

	no_adjacent_circular_transformer_fast (int n, bool build_compressed_index_table = true) : n(n), transformer(n), circular_fib(circular_fibonacci_sequence<u64>(n, fib)),
		compress_indices_and_values(get_circular_compress_indices_and_values<u64>(n)) {
		
		assert(n >= 5);
		assert(n <= K * L);
		{
			auto fib = fibonacci_sequence<u64>(L * K);
			for (int i = 0; i < L; i++) {
				u64 cur_offset = 0;
				for_each_fibonacci<u64>(K, [&] (size_t, u64 t) {
					offset[i][t] = cur_offset;
					cur_offset += fib[std::max(0, i * K - (int)(t & 1))];
				});
			}
		}
		
		auto bit_reverse = [&] (u64 x) {
			u64 res = 0;
			for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
			return res;
		};
		auto bit_rotate = [&] (u64 x, int k) { return (x >> k) | (x & (((u64) 1 << k) - 1)) << (n - k); };
		
		if (build_compressed_index_table) {
			compressed_index.resize(fib[n - 3]); // compressed_index[0] = 0
			for (size_t i = 0; i < compress_indices_and_values.size(); i++) {
				auto val = std::get<1>(compress_indices_and_values[i]);
				auto process = [&] (u64 x) {
					if (x & 1) compressed_index[get_fibonacci_index(x >> 2)] = i;
				};
				for (int j = 0; j < n; j++) process(bit_rotate(val, j));
				val = bit_reverse(val);
				for (int j = 0; j < n; j++) process(bit_rotate(val, j));
			}
		}
	}
	
	template<typename val_t> void prepare_internal_dp(const std::vector<val_t> &dp, std::vector<val_t> &internal_dp) {
		if (compressed_index.size()) { // faster, but the table eats more memory
			u64 bit0 = 1ULL << (n - 2);
			u64 bit1 = 1ULL << (n - 1);
			for_each_fibonacci<u64>(n - 4, [&] (size_t i, u64 t) {
				auto process = [&] (u64 x) {
					u64 index = 0;
					if (x) index = compressed_index[get_fibonacci_index(x >> (__builtin_ctzll(x) + 2))];
					internal_dp[i] += dp[index];
				};
				process(t << 2);
				process(t << 2 | bit1);
				if (!(t >> (n - 5) & 1)) process(t << 2 | bit0);
				if (!(t & 1)) {
					process(2 | t << 2);
					process(2 | t << 2 | bit1);
					if (!(t >> (n - 5) & 1)) process(2 | t << 2 | bit0);
				};
			});
		} else {
			for (size_t i = 0; i < dp.size(); i++) {
				auto tmp = std::get<1>(compress_indices_and_values[i]);
				int num = std::get<2>(compress_indices_and_values[i]);
				auto process = [&] (u64 x) {
					if (x & 1) return;
					internal_dp[get_fibonacci_index(x >> 2 & ((1ULL << (n - 4)) - 1))] += dp[i];
				};
				auto bit_reverse = [&] (u64 x) {
					u64 res = 0;
					for (int i = 0; i < n; i++) res |= (x >> i & 1) << (n - 1 - i);
					return res;
				};
				auto bit_rotate = [&] (u64 x, int k) { return (x >> k) | (x & (((u64) 1 << k) - 1)) << (n - k); };
				
				if (num == 2 * n) {
					for (int j = 0; j < n; j++) process(bit_rotate(tmp, j));
					tmp = bit_reverse(tmp);
					for (int j = 0; j < n; j++) process(bit_rotate(tmp, j));
				} else {
					std::vector<u64> indices;
					for (int j = 0; j < n; j++) indices.push_back(bit_rotate(tmp, j));
					tmp = bit_reverse(tmp);
					for (int j = 0; j < n; j++) indices.push_back(bit_rotate(tmp, j));
					std::sort(indices.begin(), indices.end());
					indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
					for (auto j : indices) process(j);
				}
			}
		}
	}
	template<typename val_t> void write_back_dp(const std::vector<val_t> &internal_dp, std::vector<val_t> &dp) {
		val_t new_dp_front_val = 0;
		
		for (size_t i = 0; i < dp.size(); i++) new_dp_front_val += dp[i] * std::get<2>(compress_indices_and_values[i]);
		if ((n & 1) == 0) {
			auto zero_one_val = std::get<1>(compress_indices_and_values.back());
			val_t new_dp_back_val = 0;
			for (size_t i = 0; i < dp.size(); i++) {
				auto val = std::get<1>(compress_indices_and_values[i]);
				int num = std::get<2>(compress_indices_and_values[i]);
				int ok_num = 0;
				if (val == 0) ok_num = num;
				else if (!(val & zero_one_val) || !(val & zero_one_val << 1)) ok_num = num / 2, assert(!(num & 1));
				new_dp_back_val += dp[i] * ok_num;
			}
			dp.back() = new_dp_back_val;
		}
		dp[0] = new_dp_front_val;
		for (size_t i = 1; i < dp.size(); i++) {
			if ((n & 1) == 0 && i + 1 == dp.size()) continue;
			auto val = std::get<1>(compress_indices_and_values[i]);
			assert((val & 3) == 1);
			assert((val >> (n - 2)) == 0);
			dp[i] = internal_dp[get_fibonacci_index(val >> 2 & ((1ULL << (n - 4)) - 1))];
		}
	}
	template<typename val_t> void transform_one_row(std::vector<val_t> &dp) {
		assert(dp.size() == compress_indices_and_values.size());
		
		// auto r0 = Timer::get();
		std::vector<val_t> internal_dp(fib[n - 4]);
		prepare_internal_dp(dp, internal_dp);
		// auto r1 = Timer::get();
		transformer.transform_one_row(n - 4, internal_dp);
		// auto r2 = Timer::get();
		write_back_dp(internal_dp, dp);
		// auto r3 = Timer::get();
		/*
		std::cerr << Timer::diff_ms(r0, r1) << " ms" << std::endl;
		std::cerr << Timer::diff_ms(r1, r2) << " ms" << std::endl;
		std::cerr << Timer::diff_ms(r2, r3) << " ms" << std::endl;
		std::cerr << std::endl;*/
	}
};

struct no_adjacent_circular_transformer {
	std::vector<u64> fib, circular_fib;
	std::vector<std::pair<u64, int> > compress_indices;
	
	no_adjacent_circular_transformer (int n) : fib(fibonacci_sequence<u64>(n)), circular_fib(circular_fibonacci_sequence<u64>(n, fib)),
		compress_indices(get_circular_compress_indices<u64>(n)) {}
	
	template<typename A, typename B> void compress(const std::vector<A> &a, std::vector<B> &b) {
		b.resize(compress_indices.size());
		u64 head = 0;
		for (auto i : compress_indices) b[head++] = (B) a[i.first] * i.second;
	}
	// a : compressed
	// b : uncompressed
	template<typename T, typename S> T uncompress_inner_product(const std::vector<S> &a, const std::vector<S> &b) {
		assert(a.size() == compress_indices.size());
		T res = 0;
		for (u64 i = 0; i < compress_indices.size(); i++) res += (T) a[i] * b[compress_indices[i].first];
		return res;
	}
	
	#define BIT(x, i) ((x) >> (i) & 1)
	template<typename val_t> void transform_one_square(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 2;
		u64 dst_head0 = 0;
		u64 dst_head1 = fib[n - 1];
		u64 src_head = 0;
		for_each_fibonacci<u64>(std::max(0, n - L - 1), [&] (size_t, u64 x) {
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
		if (n >= 3) {
			u64 size = fib[std::max(0, n - 3)];
			for (u64 i = 0; i < size; i++) next[dst_head0 + i] = dp[src_head + i];
			dst_head0 += size;
			src_head += size;
		}
		assert(dst_head0 == fib[n - 1]);
		assert(dst_head1 == circular_fib[n]);
		assert(src_head == circular_fib[n]);
	}
	template<typename val_t> void transform_one_square_first(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 2;
		u64 head0 = 0;
		u64 head1 = fib[n - 1];
		u64 dst_head00 = 0;
		u64 dst_head01 = fib[n - 2];
		u64 dst_head10 = fib[n - 1];
		for_each_fibonacci<u64>(n - L - 1, [&] (size_t, u64 t) {
			bool left = BIT(t, 0);
			bool right = BIT(t, n - L - 2);
			if (!left && !right) {
				auto &d0 = next[dst_head00];
				auto &d1 = next[dst_head01];
				auto &d2 = next[dst_head10];
				auto &d3 = next[dst_head00 + 1];
				auto &d4 = next[dst_head01 + 1];
				auto &s0 = dp[head0];
				auto &s1 = dp[head1];
				auto &s2 = dp[head0 + 1];
				auto &s3 = dp[head0 + 2];
				auto &s4 = dp[head1 + 1];
				
				d0 = s0 + s2;
				d1 = s1;
				d2 = s0 + s1 + s3 + s4;
				d3 = s3;
				d4 = s4;
				head0 += 3;
				head1 += 2;
				dst_head00 += 2;
				dst_head01 += 2;
				dst_head10 += 1;
			} else if (!right) {
				auto &d0 = next[dst_head00];
				auto &d1 = next[dst_head01];
				auto &d2 = next[dst_head10];
				auto &s0 = dp[head0];
				auto &s1 = dp[head1];
				auto &s2 = dp[head0 + 1];
				d0 = s0 + s2;
				d1 = s1;
				d2 = s0 + s1;
				head0 += 2;
				head1 += 1;
				dst_head00 += 1;
				dst_head01 += 1;
				dst_head10 += 1;
			} else if (!left) {
				auto &d0 = next[dst_head00];
				auto &d1 = next[dst_head10];
				auto &d2 = next[dst_head00 + 1];
				auto &s0 = dp[head0];
				auto &s1 = dp[head0 + 1];
				auto &s2 = dp[head0 + 2];
				d0 = s0 + s1;
				d1 = s0 + s2;
				d2 = s2;
				head0 += 3;
				dst_head00 += 2;
				dst_head10 += 1;
			} else {
				auto &d0 = next[dst_head00];
				auto &d1 = next[dst_head10];
				auto &s0 = dp[head0];
				auto &s1 = dp[head0 + 1];
				d0 = s0 + s1;
				d1 = s0;
				head0 += 2;
				dst_head00 += 1;
				dst_head10 += 1;
			}
		});
	}
	template<typename val_t> void transform_one_square_last(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 1;
		u64 head = 0;
		u64 dst_head0 = 0;
		u64 dst_head1 = fib[n - 1];
		for_each_fibonacci<u64>(n - L - 1, [&] (size_t, u64 t) {
			bool left = BIT(t, 0);
			bool right = BIT(t, n - L - 2);
			if (!left && !right) {
				auto &d0 = next[dst_head0];
				auto &d1 = next[dst_head1];
				auto &d2 = next[dst_head0 + 1];
				d0 = dp[head + 0] + dp[head + 1];
				d1 = dp[head + 0];
				d2 = dp[head + 2];
				head += 3;
				dst_head0 += 2;
				dst_head1 += 1;
			} else if (!right) {
				auto &d0 = next[dst_head0];
				auto &d1 = next[dst_head1];
				d0 = dp[head + 0] + dp[head + 1];
				d1 = dp[head + 0];
				head += 2;
				dst_head0++;
				dst_head1++;
			} else if (!left) {
				auto &d0 = next[dst_head0];
				auto &d1 = next[dst_head0 + 1];
				d0 = dp[head + 0];
				d1 = dp[head + 1];
				head += 2;
				dst_head0 += 2;
			} else {
				auto &d0 = next[dst_head0];
				d0 = dp[head + 0];
				head += 1;
				dst_head0++;
			}
		});
		assert(dst_head0 == fib[n - 1]);
		assert(dst_head1 == circular_fib[n]);
		assert(head == circular_fib[n]);
	}
	template<typename val_t> void transform_two_squares(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 2;
		assert(n >= L + 3);
		u64 src_head = 0;
		u64 dst_head00 = 0;
		u64 dst_head10 = fib[n - 2];
		u64 dst_head01 = fib[n - 1];
		
		for_each_fibonacci<u64>(n - L - 2, [&] (size_t, u64 t) {
			if (t & 1) {
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
		for_each_fibonacci<u64>(n - L - 3, [&] (size_t, u64 t) {
			if (t & 1) {
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
		assert(src_head == circular_fib[n]);
		assert(dst_head00 == fib[n - 2]);
		assert(dst_head10 == fib[n - 1]);
		assert(dst_head01 == circular_fib[n]);
	}
	template<typename val_t> void transform_three_squares(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 3;
		assert(n >= L + 3);
		u64 src_head = 0;
		u64 dst_head000 = 0;
		u64 dst_head100 = fib[n - 3];
		u64 dst_head010 = fib[n - 2];
		u64 dst_head001 = fib[n - 1];
		u64 dst_head101 = fib[n - 1] + fib[n - 4];
		
		for_each_fibonacci<u64>(n - L - 2, [&] (size_t, u64 right) {
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
		for_each_fibonacci<u64>(n - L - 3, [&] (size_t, u64 right) {
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
		assert(src_head == circular_fib[n]);
		assert(dst_head000 == fib[n - 3]);
		assert(dst_head100 == fib[n - 2]);
		assert(dst_head010 == fib[n - 1]);
		assert(dst_head001 == fib[n - 1] + fib[n - 4]);
		assert(dst_head101 == circular_fib[n]);
	}
	template<typename val_t> void transform_four_squares(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 4;
		assert(n >= L + 3);
		u64 src_head = 0;
		u64 dst_head0000 = 0;
		u64 dst_head1000 = fib[n - 4];
		u64 dst_head0100 = fib[n - 3];
		u64 dst_head0010 = fib[n - 2];
		u64 dst_head1010 = fib[n - 2] + fib[n - 4];
		u64 dst_head0001 = fib[n - 1];
		u64 dst_head1001 = fib[n - 1] + fib[n - 5];
		u64 dst_head0101 = fib[n - 1] + fib[n - 4];
		
		for_each_fibonacci_mod2<u64>(n - L - 2, [&] (size_t j, size_t j2) {
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
		dst_head0001 += fib[n - L - 2];
		dst_head1001 += fib[n - L - 2];
		dst_head0101 += fib[n - L - 2];
		dst_head0000 += fib[n - L - 1];
		dst_head1000 += fib[n - L - 1];
		dst_head0100 += fib[n - L - 1];
		dst_head0010 += fib[n - L - 1];
		dst_head1010 += fib[n - L - 1];
		
		for_each_fibonacci_mod2<u64>(n - L - 3, [&] (size_t j, size_t j2) {
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
		dst_head0000 += fib[n - L - 2];
		dst_head0100 += fib[n - L - 2];
		dst_head0010 += fib[n - L - 2];
		dst_head0001 += fib[n - L - 3];
		dst_head0101 += fib[n - L - 3];
		
		assert(src_head == circular_fib[n]);
		assert(dst_head0000 == fib[n - 4]);
		assert(dst_head1000 == fib[n - 3]);
		assert(dst_head0100 == fib[n - 2]);
		assert(dst_head0010 == fib[n - 2] + fib[n - 4]);
		assert(dst_head1010 == fib[n - 1]);
		assert(dst_head0001 == fib[n - 1] + fib[n - 5]);
		assert(dst_head1001 == fib[n - 1] + fib[n - 4]);
		assert(dst_head0101 == circular_fib[n]);
	}
	template<typename val_t> void transform_five_squares(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 5;
		assert(n >= L + 3);
		u64 src_head = 0;
		
		u64 dst_head00000 = 0;
		u64 dst_head10000 = fib[n - 5];
		u64 dst_head01000 = fib[n - 4];
		u64 dst_head00100 = fib[n - 3];
		u64 dst_head10100 = fib[n - 3] + fib[n - 5];
		u64 dst_head00010 = fib[n - 2];
		u64 dst_head10010 = fib[n - 2] + fib[n - 5];
		u64 dst_head01010 = fib[n - 2] + fib[n - 4];
		u64 dst_head00001 = fib[n - 1];
		u64 dst_head10001 = fib[n - 1] + fib[n - 6];
		u64 dst_head01001 = fib[n - 1] + fib[n - 5];
		u64 dst_head00101 = fib[n - 1] + fib[n - 4];
		u64 dst_head10101 = fib[n - 1] + fib[n - 4] + fib[n - 6];
		
		for_each_fibonacci_mod2<u64>(n - L - 2, [&] (size_t j, size_t j2) {
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
		dst_head00000 += fib[n - L - 1];
		dst_head10000 += fib[n - L - 1];
		dst_head01000 += fib[n - L - 1];
		dst_head00100 += fib[n - L - 1];
		dst_head10100 += fib[n - L - 1];
		dst_head00010 += fib[n - L - 1];
		dst_head10010 += fib[n - L - 1];
		dst_head01010 += fib[n - L - 1];
		dst_head00001 += fib[n - L - 2];
		dst_head10001 += fib[n - L - 2];
		dst_head01001 += fib[n - L - 2];
		dst_head00101 += fib[n - L - 2];
		dst_head10101 += fib[n - L - 2];
		
		for_each_fibonacci_mod2<u64>(n - L - 3, [&] (size_t j, size_t j2) {
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
		dst_head00000 += fib[n  - L - 2];
		dst_head01000 += fib[n  - L - 2];
		dst_head00100 += fib[n  - L - 2];
		dst_head00010 += fib[n  - L - 2];
		dst_head01010 += fib[n  - L - 2];
		dst_head00001 += fib[n  - L - 3];
		dst_head01001 += fib[n  - L - 3];
		dst_head00101 += fib[n  - L - 3];
		
		assert(src_head == circular_fib[n]);
		assert(dst_head00000 == fib[n - 5]);
		assert(dst_head10000 == fib[n - 4]);
		assert(dst_head01000 == fib[n - 3]);
		assert(dst_head00100 == fib[n - 3] + fib[n - 5]);
		assert(dst_head10100 == fib[n - 2]);
		assert(dst_head00010 == fib[n - 2] + fib[n - 5]);
		assert(dst_head10010 == fib[n - 2] + fib[n - 4]);
		assert(dst_head01010 == fib[n - 1]);
		assert(dst_head00001 == fib[n - 1] + fib[n - 6]);
		assert(dst_head10001 == fib[n - 1] + fib[n - 5]);
		assert(dst_head01001 == fib[n - 1] + fib[n - 4]);
		assert(dst_head00101 == fib[n - 1] + fib[n - 4] + fib[n - 6]);
		assert(dst_head10101 == circular_fib[n]);
	}
	template<typename val_t> void transform_last_five_squares(int n, std::vector<val_t> &dp, std::vector<val_t> &next) {
		constexpr int L = 5;
		assert(n >= L + 2);
		u64 src_head = 0;
		
		u64 dst_head00000 = 0;
		u64 dst_head10000 = fib[n - 5];
		u64 dst_head01000 = fib[n - 4];
		u64 dst_head00100 = fib[n - 3];
		u64 dst_head10100 = fib[n - 3] + fib[n - 5];
		u64 dst_head00010 = fib[n - 2];
		u64 dst_head10010 = fib[n - 2] + fib[n - 5];
		u64 dst_head01010 = fib[n - 2] + fib[n - 4];
		u64 dst_head00001 = fib[n - 1];
		u64 dst_head10001 = fib[n - 1] + fib[n - 6];
		u64 dst_head01001 = fib[n - 1] + fib[n - 5];
		u64 dst_head00101 = fib[n - 1] + fib[n - 4];
		u64 dst_head10101 = fib[n - 1] + fib[n - 4] + fib[n - 6];
		
		for_each_fibonacci_mod2<u64>(n - L - 1, [&] (size_t j, size_t j2) {
			auto &d0 = next[dst_head00000 + j];
			auto &d1 = next[dst_head10000 + j];
			auto &d2 = next[dst_head01000 + j];
			auto &d3 = next[dst_head00100 + j];
			auto &d4 = next[dst_head10100 + j];
			auto &d5 = next[dst_head00010 + j];
			auto &d6 = next[dst_head10010 + j];
			auto &d7 = next[dst_head01010 + j];
			auto &d8 = next[dst_head00001 + j2];
			auto &d9 = next[dst_head10001 + j2];
			auto &d10 = next[dst_head01001 + j2];
			auto &d11 = next[dst_head00101 + j2];
			auto &d12 = next[dst_head10101 + j2];
			
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
		}, [&] (size_t j, size_t) {
			auto &d0 = next[dst_head00000 + j];
			auto &d1 = next[dst_head10000 + j];
			auto &d2 = next[dst_head01000 + j];
			auto &d3 = next[dst_head00100 + j];
			auto &d4 = next[dst_head10100 + j];
			auto &d5 = next[dst_head00010 + j];
			auto &d6 = next[dst_head10010 + j];
			auto &d7 = next[dst_head01010 + j];
			
			auto l02 = dp[src_head + 0] + dp[src_head + 2];
			d4 = l02 + dp[src_head + 5] + dp[src_head + 7];
			d1 = d4 + dp[src_head + 3];
			d3 = d4 + dp[src_head + 1] + dp[src_head + 6];
			d6 = l02 + dp[src_head + 3];
			d5 = d6 + dp[src_head + 1] + dp[src_head + 4];
			d0 = d3 + dp[src_head + 3] + dp[src_head + 4];
			d7 = d5 - dp[src_head + 2];
			d2 = d7 + dp[src_head + 5] + dp[src_head + 6];
			
			src_head += 8;
		});
		dst_head00000 += fib[n - L - 1];
		dst_head10000 += fib[n - L - 1];
		dst_head01000 += fib[n - L - 1];
		dst_head00100 += fib[n - L - 1];
		dst_head10100 += fib[n - L - 1];
		dst_head00010 += fib[n - L - 1];
		dst_head10010 += fib[n - L - 1];
		dst_head01010 += fib[n - L - 1];
		dst_head00001 += fib[n - L - 2];
		dst_head10001 += fib[n - L - 2];
		dst_head01001 += fib[n - L - 2];
		dst_head00101 += fib[n - L - 2];
		dst_head10101 += fib[n - L - 2];
		
		for_each_fibonacci_mod2<u64>(n - L - 2, [&] (size_t j, size_t j2) {
			auto &d0 = next[dst_head00000 + j];
			auto &d1 = next[dst_head01000 + j];
			auto &d2 = next[dst_head00100 + j];
			auto &d3 = next[dst_head00010 + j];
			auto &d4 = next[dst_head01010 + j];
			auto &d5 = next[dst_head00001 + j2];
			auto &d6 = next[dst_head01001 + j2];
			auto &d7 = next[dst_head00101 + j2];
			
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
		}, [&] (size_t j, size_t) {
			auto &d0 = next[dst_head00000 + j];
			auto &d1 = next[dst_head01000 + j];
			auto &d2 = next[dst_head00100 + j];
			auto &d3 = next[dst_head00010 + j];
			auto &d4 = next[dst_head01010 + j];
			
			d4 = dp[src_head + 0] + dp[src_head + 2];
			d3 = d4 + dp[src_head + 1];
			d1 = d4 + dp[src_head + 3];
			d0 = d1 + dp[src_head + 1] + dp[src_head + 4];
			d2 = d0 - dp[src_head + 2];
			
			src_head += 5;
		});
		dst_head00000 += fib[n  - L - 2];
		dst_head01000 += fib[n  - L - 2];
		dst_head00100 += fib[n  - L - 2];
		dst_head00010 += fib[n  - L - 2];
		dst_head01010 += fib[n  - L - 2];
		dst_head00001 += fib[std::max(0, n  - L - 3)];
		dst_head01001 += fib[std::max(0, n  - L - 3)];
		dst_head00101 += fib[std::max(0, n  - L - 3)];
		
		assert(src_head == circular_fib[n]);
		assert(dst_head00000 == fib[n - 5]);
		assert(dst_head10000 == fib[n - 4]);
		assert(dst_head01000 == fib[n - 3]);
		assert(dst_head00100 == fib[n - 3] + fib[n - 5]);
		assert(dst_head10100 == fib[n - 2]);
		assert(dst_head00010 == fib[n - 2] + fib[n - 5]);
		assert(dst_head10010 == fib[n - 2] + fib[n - 4]);
		assert(dst_head01010 == fib[n - 1]);
		assert(dst_head00001 == fib[n - 1] + fib[n - 6]);
		assert(dst_head10001 == fib[n - 1] + fib[n - 5]);
		assert(dst_head01001 == fib[n - 1] + fib[n - 4]);
		assert(dst_head00101 == fib[n - 1] + fib[n - 4] + fib[n - 6]);
		assert(dst_head10101 == circular_fib[n]);
	}
	template<typename val_t> void transform_one_row(int n, std::vector<val_t> &dp) {
		assert(dp.size() == circular_fib[n]);
		
		if (n <= 1) return;
		if (n == 2) {
			dp = {dp[0] + dp[1] + dp[2], dp[0] + dp[2], dp[0] + dp[1]};
			return;
		}
		if (n == 3) {
			dp = {
				dp[0] + dp[1] + dp[2] + dp[3],
				dp[0] + dp[2] + dp[3],
				dp[0] + dp[1] + dp[3],
				dp[0] + dp[1] + dp[2],
			};
			return;
		}
		
		std::vector<val_t> dp_swap(circular_fib[n]);
		for (int i = 0; i < n; ) {
			int advance = 1;
			if (i == 0) transform_one_square_first(n, dp, dp_swap);
			else if (i == n - 1) transform_one_square_last(n, dp, dp_swap);
			else if (n <= 6 || i >= n - 4 || i == n - 6) transform_one_square(n, dp, dp_swap);
			else if (i == n - 5) transform_last_five_squares(n, dp, dp_swap), advance = 5;
			else if (i == n - 7) transform_two_squares(n, dp, dp_swap), advance = 2;
			else if (i == n - 8) transform_three_squares(n, dp, dp_swap), advance = 3;
			else if (i == n - 9) transform_four_squares(n, dp, dp_swap), advance = 4;
			else transform_five_squares(n, dp, dp_swap), advance = 5;
			
			std::swap(dp, dp_swap);
			i += advance;
		}
	}
};