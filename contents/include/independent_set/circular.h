#pragma once
#include <vector>
#include <map>
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
	
	size_t get_buf_size(int n) {
		assert(n == this->n);
		return compress_indices_and_values.size();
	}
	
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
		
		auto r0 = Timer::get();
		std::vector<val_t> internal_dp(fib[n - 4]);
		prepare_internal_dp(dp, internal_dp);
		auto r1 = Timer::get();
		transformer.transform_one_row(n - 4, internal_dp);
		auto r2 = Timer::get();
		write_back_dp(internal_dp, dp);
		auto r3 = Timer::get();
		std::cerr << Timer::diff_ms(r0, r1) << " ms" << std::endl;
		std::cerr << Timer::diff_ms(r1, r2) << " ms" << std::endl;
		std::cerr << Timer::diff_ms(r2, r3) << " ms" << std::endl;
		std::cerr << std::endl;
		/*
			true:
			871.774 ms
			356.555 ms
			17.1513 ms
			
			false:
			2351.17 ms
			382.738 ms
			19.2428 ms
		*/
	}
};

struct circular_no_adjacent_transformer {
	std::vector<u64> fib;
	std::vector<std::pair<u64, int> > compress_indices;
	
	circular_no_adjacent_transformer (int n) : fib(fibonacci_sequence<u64>(n)),
		compress_indices(get_circular_compress_indices<u64>(n)) {}
	
	u64 get_buf_size(u32 n) {
		if (n == 1) return 1;
		if (n == 2) return fib[2];
		assert(n < fib.size());
		return fib[n - 1] + fib[n - 3];
	}
	
	
	/*
		Executes k-bit transform:
			buf: fib[k]-long
			output buf[i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[j]
				where addition is done via function `add`
		
		Algorithm:
			Let
				A = buf[fibindex(0000...0)...fibindex(0011...1)]
				B = buf[fibindex(0100...0)...fibindex(0111...1)]
				C = buf[fibindex(1000...0)...fibindex(1000...1)].
			Note that
				len(A) = fib[k - 2], A = buf[0...fib[k-2]-1]
				len(B) = fib[k - 3], B = buf[fib[k-2]..fib[k-1]-1]
				len(C) = fib[k - 2], B = buf[fib[k-1]..fib[k]-1]
			and 
				buf = concat(A, B, C).
			Then deciding the uppermost one or two bits either 0 or 10 corresponds:
				buf -> [A + C, B, A + B] (additions are vector addition, treating out-of-bound element as zero)
				 : can be done by swapping A and C, adding A to C and adding B to A
			We can recursively do this on [A + C, B] and [A + B] (of size fib[n - 1] and fib[n - 2], resp.)
	*/
	template<typename T, typename AddOpType> void execute_abstract_transform(int k, T *buf, const AddOpType &add) {
		static_assert(std::is_convertible<AddOpType, std::function<T (T, T)> >::value ||
					  std::is_convertible<AddOpType, std::function<T (const T&, const T&)> >::value);
		if (k <= 0) return;
		if (k == 1) {
			std::swap(buf[0], buf[1]);
			buf[0] = add(buf[0], buf[1]);
			return;
		}
		if (k == 2) { // due to fib[-1] = 1 undefined
			T ab = add(buf[0], buf[1]);
			T ac = add(buf[0], buf[2]);
			buf[0] = add(ab, buf[2]);
			buf[1] = ac;
			buf[2] = ab;
			return;
		}
		if (k == 3) { // 7 additions
			T a0 = buf[0];
			T a1 = buf[1];
			T a2 = buf[2];
			T a3 = buf[3];
			T a4 = buf[4];
			T a03 = add(a0, a3);
			T a14 = add(a1, a4);
			T a02 = add(a0, a2);
			T a032 = add(a03, a2);
			
			buf[0] = add(a032, a14); // 01234
			buf[1] = a032;
			buf[2] = add(a03, a14); // 0134
			buf[3] = add(a02, a1); // 012
			buf[4] = a02;
			return;
		}
		if (k == 4) {
			T a0 = buf[0];
			T a1 = buf[1];
			T a2 = buf[2];
			T a3 = buf[3];
			T a4 = buf[4];
			T a5 = buf[5];
			T a6 = buf[6];
			T a7 = buf[7];
			T a05 = add(a0, a5);
			T a16 = add(a1, a6);
			T a27 = add(a2, a7);
			T a03 = add(a0, a3);
			T a14 = add(a1, a4);
			T a053 = add(a05, a3);
			T a164 = add(a16, a4);
			T a0527 = add(a05, a27);
			T a05327 = add(a053, a27);
			T a032 = add(a03, a2);
			
			buf[0] = add(a05327, a164); // 01234567
			buf[1] = a05327; // 02357
			buf[2] = add(a053, a164); // 013456
			buf[3] = add(a0527, a16); // 012567
			buf[4] = a0527;
			buf[5] = add(a032, a14);
			buf[6] = a032;
			buf[7] = add(a03, a14);
			return;
		}
		for (u64 i = 0; i < fib[k - 3]; i++) {
			T a = buf[i];
			T b = buf[i + fib[k - 2]];
			T c = buf[i + fib[k - 1]];
			buf[i] = add(a, c);
			buf[i + fib[k - 1]] = add(a, b);
		}
		for (u64 i = fib[k - 3]; i < fib[k - 2]; i++) {
			T a = buf[i];
			T c = buf[i + fib[k - 1]];
			buf[i] = add(a, c);
			buf[i + fib[k - 1]] = a;
		}
		
		execute_abstract_transform(k - 1, buf, add);
		execute_abstract_transform(k - 2, buf + fib[k - 1], add);
	}
	/*
		Executes k-bit transform but with additional uppermost bit that turns off if (k-1)-th bit turns on by transformation:
			buf: fib[k+1]-long
			output buf[0i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[0j] +
							 sum_{j: (fibstr(i) & fibstr(j)) == 0, fibstr(j)[k-1] == 1} buf[1j]
			output buf[1i] = sum_{j: (fibstr(i) & fibstr(j)) == 0, fibstr(j)[k-1] == 0} buf[1j]
				where addition is done via function `add`
		
	*/
	template<typename T, typename AddOpType> void execute_abstract_transform_plus_one(int k, T *buf, const AddOpType &add) {
		execute_abstract_transform(k, buf, add);
		execute_abstract_transform(k - 1, buf + fib[k], add);
		
		for (u64 i = 0; i < fib[k - 2]; i++) buf[i + fib[k - 1]] = add(buf[i + fib[k - 1]], buf[i + fib[k]]);
	}
	
	
	
	
	template<typename T> void execute_abstract_transform_vecs(int k, T *buf, u64 n_vecs) {
		if (k <= 0) return;
		using vec_t = typename std::remove_reference<decltype(buf[0][0])>::type;
		
		// std::cerr << "!!!!!!!! " << k << " !!!!!!!!!!!!!!"  << std::endl;
		if (k == 1) { // 1 addition
			for (u64 i = 0; i < n_vecs; i++) {
				vec_t a0 = buf[0][i];
				vec_t a1 = buf[1][i];
				buf[0][i] = a0 + a1;
				buf[1][i] = a0;
			}
			return;
		}
		if (k == 2) { // 3 additions
			for (u64 i = 0; i < n_vecs; i++) {
				vec_t a0 = buf[0][i];
				vec_t a1 = buf[1][i];
				vec_t a2 = buf[2][i];
				vec_t a01 = a0 + a1;
				buf[0][i] = a01 + a2;
				buf[1][i] = a0 + a2;
				buf[2][i] = a01;
			}
			return;
		}
		if (k == 3) { // 7 additions
			for (u64 i = 0; i < n_vecs; i++) {
				vec_t a0 = buf[0][i];
				vec_t a1 = buf[1][i];
				vec_t a2 = buf[2][i];
				vec_t a3 = buf[3][i];
				vec_t a4 = buf[4][i];
				vec_t a03 = a0 + a3;
				vec_t a14 = a1 + a4;
				vec_t a02 = a0 + a2;
				vec_t a032 = a03 + a2;
				
				buf[0][i] = a032 + a14; // 01234
				buf[1][i] = a032;
				buf[2][i] = a03 + a14; // 0134
				buf[3][i] = a02 + a1; // 012
				buf[4][i] = a02;
			}
			return;
		}
		
		u64 unit5 = fib[std::max(0, k - 5)];
		u64 unit4 = fib[k - 4];
		u64 unit3 = fib[k - 3];
		u64 unit2 = fib[k - 2];
		u64 unit1 = fib[k - 1];
		for (u64 i = 0; i < unit5; i++) {
			T &t0 = buf[i];
			T &t1 = buf[unit4 + i];
			T &t2 = buf[unit3 + i];
			T &t3 = buf[unit2 + i];
			T &t4 = buf[unit2 + unit4 + i];
			T &t5 = buf[unit1 + i];
			T &t6 = buf[unit1 + unit4 + i];
			T &t7 = buf[unit1 + unit3 + i];
			for (u64 j = 0; j < n_vecs; j++) {
				vec_t a05 = t0[j] + t5[j];
				vec_t a16 = t1[j] + t6[j];
				vec_t a27 = t2[j] + t7[j];
				vec_t a03 = t0[j] + t3[j];
				vec_t a14 = t1[j] + t4[j];
				vec_t a053 = a05 + t3[j];
				vec_t a164 = a16 + t4[j];
				vec_t a0527 = a05 + a27;
				vec_t a05327 = a053 + a27;
				vec_t a032 = a03 + t2[j];
				
				t0[j] = a05327 + a164; // 01234567
				t1[j] = a05327; // 02357
				t2[j] = a053 + a164; // 013456
				t3[j] = a0527 + a16; // 012567
				t4[j] = a0527;
				t5[j] = a032 + a14;
				t6[j] = a032;
				t7[j] = a03 + a14;
			}
		}
		for (u64 i = unit5; i < unit4; i++) {
			T &t0 = buf[i];
			T &t1 = buf[unit4 + (i-unit5)]; // destination only
			T &t2 = buf[unit3 + i];
			T &t3 = buf[unit2 + i];
			T &t4 = buf[unit2 + unit4 + (i-unit5)]; // destination only
			T &t5 = buf[unit1 + i];
			T &t6 = buf[unit1 + unit4 + (i-unit5)]; // destination only
			T &t7 = buf[unit1 + unit3 + i];
			
			for (u64 j = 0; j < n_vecs; j++) {
				vec_t a05 = t0[j] + t5[j];
				vec_t a27 = t2[j] + t7[j];
				vec_t a03 = t0[j] + t3[j];
				vec_t a053 = a05 + t3[j];
				vec_t a032 = a03 + t2[j];
				vec_t a05327 = a053 + a27;
				vec_t a0527 = a05 + a27;
				
				t0[j] = a05327;
				t1[j] += a05327;
				t2[j] = a053;
				t3[j] = a0527;
				t4[j] += a0527;
				t5[j] = a032;
				t6[j] += a032;
				t7[j] = a03;
			}
		}
		
		
		execute_abstract_transform_vecs(k - 4, buf, n_vecs);
		execute_abstract_transform_vecs(k - 5, buf + unit4, n_vecs);
		execute_abstract_transform_vecs(k - 4, buf + unit3, n_vecs);
		execute_abstract_transform_vecs(k - 4, buf + unit2, n_vecs);
		execute_abstract_transform_vecs(k - 5, buf + unit2 + unit4, n_vecs);
		execute_abstract_transform_vecs(k - 4, buf + unit1, n_vecs);
		execute_abstract_transform_vecs(k - 5, buf + unit1 + unit4, n_vecs);
		execute_abstract_transform_vecs(k - 4, buf + unit1 + unit3, n_vecs);
	}
	/*
		Executes k-bit transform but with additional uppermost bit that turns off if (k-1)-th bit turns on by transformation:
			buf: fib[k+1]-long
			output buf[0i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[0j] +
							 sum_{j: (fibstr(i) & fibstr(j)) == 0, fibstr(j)[k-1] == 1} buf[1j]
			output buf[1i] = sum_{j: (fibstr(i) & fibstr(j)) == 0, fibstr(j)[k-1] == 0} buf[1j]
				where addition is done via function `add`
		
	*/
	template<typename T> void execute_abstract_transform_vecs_plus_one(int k, T *buf, u64 n_vecs) {
		execute_abstract_transform_vecs(k, buf, n_vecs);
		execute_abstract_transform_vecs(k - 1, buf + fib[k], n_vecs);
		
		for (u64 i = 0; i < fib[k - 2]; i++) {
			for (u64 j = 0; j < n_vecs; j++) {
				buf[i + fib[k - 1]][j] += buf[i + fib[k]][j];
			}
		}
	}
	
	
	template<typename T> void transform_first_C_squares_naive(int n, int C, const std::vector<T> &dp, std::vector<T> &next) {
		std::map<u64, size_t> rev;
		for_each_circular_fibonacci<u64>(n, [&] (size_t i, u64 str) { rev[str] = i; });
		
		next.assign(fib[n - 1] + fib[n - 3], 0);
		for_each_circular_fibonacci<u64>(n, [&] (size_t i, u64 cur_str) {
			for_each_fibonacci<u64>(C, [&] (size_t j, u64 new_low) {
				(void) j;
				if (cur_str & new_low) return;
				
				u64 new_val = (cur_str & ~((1 << C) - 1)) | new_low;
				if ((new_low >> (C - 1) & 1) && (cur_str >> C & 1)) new_val ^= 1ULL << C;
				if ((new_low & 1) && (cur_str >> (n-1) & 1)) new_val ^= 1ULL << (n - 1);
				
				assert(rev.count(new_val));
				next[rev[new_val]] += dp[i];
			});
		});
	}
	template<typename T>  void transform_first_C_squares(int n, int C, const std::vector<T> &dp, std::vector<T> &next) {
		assert(n - C - 2 >= 0);
		assert(C >= 2);
		// enumerate upper n - C - 1 digits
		u64 head0 = 0;          // head of *******0
		u64 head1 = fib[n - 1]; // head of *******1
		
		for_each_fibonacci<u64>(std::max(0, n - C - 2), [&] (size_t i, u64 upper) {
			(void) i;
			
			if (upper >> (n - C - 3) & 1) {
				for (u64 j = 0; j < fib[C + 1 - (upper & 1)]; j++) next[head0 + j] = dp[head0 + j];
				
				if (upper & 1) execute_abstract_transform(C, next.data() + head0, [] (const T &x, const T &y) { return x + y; });
				else execute_abstract_transform_plus_one(C, next.data() + head0, [] (const T &x, const T &y) { return x + y; });
				
				head0 += fib[C + 1 - (upper & 1)];
			} else {
				for (u64 j = 0; j < fib[C + 1 - (upper & 1)]; j++) next[head0 + j] = dp[head0 + j];
				for (u64 j = 0; j < fib[C - (upper & 1)]; j++) next[head1 + j] = dp[head1 + j];
				
				if (upper & 1) {
					execute_abstract_transform(C,   next.data() + head0, [] (const T &x, const T &y) { return x + y; });
					execute_abstract_transform(C-1, next.data() + head1, [] (const T &x, const T &y) { return x + y; });
				} else {
					execute_abstract_transform_plus_one(C  , next.data() + head0, [] (const T &x, const T &y) { return x + y; });
					execute_abstract_transform_plus_one(C-1, next.data() + head1, [] (const T &x, const T &y) { return x + y; });
				}
				
				for_each_fibonacci<u64>(C - 1 - (upper & 1), [&] (size_t, u64 middle) {
					// 10******0 += 00*******1
					next[head0 + 1] += next[head1 + 0];
					
					if (!(middle & 1)) {
						head0 += 3;
						head1 += 2;
					} else {
						head0 += 2;
						head1 += 1;
					}
				});
			}
		});
		assert(head0 == fib[n - 1]);
		assert(head1 == fib[n - 1] + fib[n - 3]);
	}
	
	template<typename T> void transform_C_squares_naive(int n, int column, int C, std::vector<T> &dp) {
		std::map<u64, size_t> rev;
		for_each_circular_fibonacci<u64>(n, [&] (size_t i, u64 str) { rev[str] = i; });
		
		std::vector<T> next(dp.size());
		
		for_each_circular_fibonacci<u64>(n, [&] (size_t i, u64 cur_str) {
			for_each_fibonacci<u64>(C, [&] (size_t j, u64 new_part) {
				(void) j;
				if (cur_str >> column & new_part) return;
				
				u64 replaced = (cur_str & ~(((1ULL << C) - 1) << column)) | (new_part << column);
				if (column + C == n) {
					if ((new_part >> (C - 1) & 1) && (replaced & 1)) return;
				} else {
					if ((new_part >> (C - 1) & 1) && (cur_str >> (C + column) & 1)) replaced ^= 1ULL << (C + column);
				}
				if ((new_part & 1) && (cur_str >> (column - 1) & 1)) return;
				
				assert(rev.count(replaced));
				next[rev[replaced]] += dp[i];
			});
		});
		dp = next;
	}
	// column -> column + C
	// in-place
	template<typename T>  void transform_C_squares_(int n, int column, int C, T *dp) {
		constexpr int D = 5;
		assert(column >= D + 1);
		assert(column + C <= n);
		assert(C >= 2);
		
		
		std::vector<T *> ptrs_buf(fib[C + 1]);
		/*
			  low <- [D] [column - D] [C][1] [n - column - C - 1] -> high
			or 
			  low <- [D] [column - D] [C][1](virtual, 0)  [1](virtual, 1) -> high
			if column + C == n
		*/
		
		u64 head = 0;
		for_each_fibonacci<u64>(std::max(0, n - column - C - 1), [&] (size_t, u64 up_upper) {
			// virtually set up_upper ([n-column-C-1]) to one to fix [1] to zero
			if (column + C == n) up_upper = 1;
			
			u64 j_head_offset = 0;
			for_each_fibonacci<u64>(column - D, [&] (size_t, u64 low_upper) {
				int c = C;
				bool do_plus_one_transform = true;
				if (up_upper & 1) do_plus_one_transform = false;
				if (low_upper >> (column - D - 1) & 1) c--; // the lowest bit of [C] is fixed to zero in this case
				int d = D;
				if (low_upper & 1) d--;
				u64 n_vecs = fib[d]; // # of vecs that actually need addition (=ceil(fib[d] / 8))
				
				// load
				u64 head_local = head + j_head_offset;
				for_each_fibonacci<u64>(c + do_plus_one_transform, [&] (size_t k, u64 up_lower) {
					// load dp[in_head..in_head+fib[d]-1] by vector
					ptrs_buf[k] = dp + head_local;
					
					head_local += fib[column +
						(low_upper >> (column - D - 1) & 1) - // entire low_upper is offset by one in this case
						(up_lower & 1)];
				});
				// transform
				if (do_plus_one_transform) execute_abstract_transform_vecs_plus_one(c, ptrs_buf.data(), n_vecs);
				else execute_abstract_transform_vecs(c, ptrs_buf.data(), n_vecs);
				
				assert(head_local == head + j_head_offset + fib[column + C + 1 - (up_upper & 1)]);
				j_head_offset += fib[d];
			});
			head += fib[column + C + 1 - (up_upper & 1)];
		});
	}
	// equivalent to transform_C_squares_naive(n, column, C, dp)
	template<typename T>  void transform_C_squares(int n, int column, int C, std::vector<T> &dp) {
		assert(column + C <= n - 3);
		transform_C_squares_(n - 1, column, C, dp.data());
		transform_C_squares_(n - 3, column - 1, C, dp.data() + fib[n - 1]);
	}
	// equivalent to transform_C_squares_naive(n, n - C, C, dp)
	template<typename T>  void transform_last_C_squares(int n, int C, std::vector<T> &dp) {
		transform_C_squares_(n - 1, n - C, C - 1, dp.data());
		transform_C_squares_(n - 3, n - C - 1, C - 2, dp.data() + fib[n - 1]);
		
		u64 head00 = 0;          // *******00: the uppermost digit can be altered
		u64 head10 = fib[n - 2]; // ******010: the uppermost digit can be altered
		u64 head01 = fib[n - 1]; // 0******01: the uppermost two digits can be altered
		// enumerate n-3 digits except the lowermost and the uppermost two
		for_each_fibonacci<u64>(n - 3, [&] (size_t, u64 middle) {
			if (!(middle & 1) && !(middle >> (n - 4) & 1)) {
				// a0_00: 0[middle]00
				// a1_10: 1[middle]10, and so on...
				T &a0_00 = dp[head00 + 0]; // -> a0_00 + a0_01
				// T &a1_00 = dp[head00 + 1]; // -> a1_00
				T &a0_10 = dp[head10 + 0]; // -> a0_10 + a0_01
				// T &a1_10 = dp[head10 + 1]; // -> a1_10
				T &a0_01 = dp[head01 + 0]; // -> a0_00
				
				T new_a0_00 = a0_00 + a0_01;
				T new_a0_10 = a0_10 + a0_01;
				T new_a0_01 = a0_00;
				a0_00 = new_a0_00;
				a0_10 = new_a0_10;
				a0_01 = new_a0_01;
				head00 += 2;
				head10 += 2;
				head01 += 1;
			} else if (!(middle & 1)) {
				T &a0_00 = dp[head00 + 0]; // -> a0_00 + a0_01
				// T &a1_00 = dp[head00 + 1]; // -> a1_00
				T &a0_01 = dp[head01 + 0]; // -> a0_00
				
				T new_a0_00 = a0_00 + a0_01;
				T new_a0_01 = a0_00;
				a0_00 = new_a0_00;
				a0_01 = new_a0_01;
				head00 += 2;
				head01 += 1;
			} else if (!(middle >> (n - 4) & 1)) {
				T &a0_00 = dp[head00 + 0]; // -> a0_00 + a0_01
				T &a0_10 = dp[head10 + 0]; // -> a0_10 + a0_01
				T &a0_01 = dp[head01 + 0]; // -> a0_00
				
				T new_a0_00 = a0_00 + a0_01;
				T new_a0_10 = a0_10 + a0_01;
				T new_a0_01 = a0_00;
				a0_00 = new_a0_00;
				a0_10 = new_a0_10;
				a0_01 = new_a0_01;
				head00 += 1;
				head10 += 1;
				head01 += 1;
			} else { // both end of `middle` are one
				T &a0_00 = dp[head00 + 0]; // -> a0_00 + a0_01
				T &a0_01 = dp[head01 + 0]; // -> a0_00
				
				T new_a0_00 = a0_00 + a0_01;
				T new_a0_01 = a0_00;
				a0_00 = new_a0_00;
				a0_01 = new_a0_01;
				head00 += 1;
				head01 += 1;
			}
		});
		assert(head00 == fib[n - 2]);
		assert(head10 == fib[n - 1]);
		assert(head01 == fib[n - 1] + fib[n - 3]);
	}
	 
	template<typename T> void transform_one_row(int n, const std::vector<T> &dp, std::vector<T> &next) {
		assert(dp.size() == get_buf_size(n));
		assert(next.size() == get_buf_size(n));
		
		if (n == 0) return;
		
		if (n <= 16) {
			next = dp;
			transform_C_squares_naive(n, 0, n, next);
			return;
		}
		constexpr int C = 12;
		int n_steps = (n + C - 1) / C;
		
		for (int i = 0; i < n_steps; i++) {
			int cur_column = n * i / n_steps;
			int next_column = n * (i + 1) / n_steps;
			// std::cerr << "step " << i << ": " << next_column - cur_column << std::endl;
			// Timer::measure([&] () {
				if (i == 0) transform_first_C_squares(n, next_column - cur_column, dp, next);
				else if (i != n_steps - 1) transform_C_squares(n, cur_column, next_column - cur_column, next);
				else transform_last_C_squares(n, next_column - cur_column, next);
			// });
		}
	}
};