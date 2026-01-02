#pragma once
#include "types.h"
#include "avx_mask_generator.h"
#include <immintrin.h>
#include <mutex>
#include <thread>

struct mint32x16 {
	using T = u32;
	static constexpr size_t LEN = 16; // # of T in a vector
	
	class no_mask_tag {};
	static constexpr no_mask_tag NO_MASK{};
	
	using vec_t = __m512i;
	using mask_t = __mmask16;
	static vec_t mod_vec;
	
	using prefix_mask_generator_t = avx512_prefix_mask_gen_t<mint32x16>;
	
	vec_t val;
	mint32x16 () : val(_mm512_setzero_si512()) {}
	mint32x16 (u32 x) : val(_mm512_set1_epi32(x)) {} // x should be less than any of the mods
	mint32x16 (const vec_t &val) : val(val) {}
	mint32x16 (const mint32x16 &rhs) = default;
	mint32x16 & operator += (const mint32x16 &rhs) {
		val = _mm512_add_epi32(val, rhs.val);
		val = _mm512_min_epu32(val, _mm512_sub_epi32(val, mod_vec));
		return *this;
	}
	mint32x16 & operator -= (const mint32x16 &rhs) {
		val = _mm512_sub_epi32(val, rhs.val);
		val = _mm512_min_epu32(val, _mm512_add_epi32(val, mod_vec));
		return *this;
	}
	mint32x16 operator + (const mint32x16 &rhs) const { return mint32x16(*this) += rhs; }
	mint32x16 operator - (const mint32x16 &rhs) const { return mint32x16(*this) -= rhs; }
	
	static mint32x16 loadu(const u32 *ptr) { return { _mm512_loadu_si512((vec_t *) ptr) }; }
	static mint32x16 loadu(const mint32x16 *ptr) { return { _mm512_loadu_si512(&ptr->val) }; }
	void storeu(u32 *ptr) const { _mm512_storeu_si512((vec_t *) ptr, val); }
	void storeu_masked(u32 *ptr, mask_t mask) const { _mm512_mask_storeu_epi32((int *) ptr, mask, val); }
	void storeu_masked(u32 *ptr, no_mask_tag) const { storeu(ptr); } // for common interface with actual masked store
	
	static mint32x16 zero() { return { _mm512_set1_epi32(0) }; }
	static void set_mod(u32 mod) { mod_vec = _mm512_set1_epi32(mod); }
	u32 operator [] (int i) const {
		u32 tmp[LEN];
		storeu(tmp);
		return tmp[i];
	}
	friend std::ostream & operator << (std::ostream &stream, const mint32x16 &rhs) {
		stream << "[";
		for (size_t i = 0; i < LEN; i++) {
			if (i) stream << ", ";
			stream << rhs[i];
		}
		stream << "]";
		return stream;
	}
};
mint32x16::vec_t mint32x16::mod_vec;



struct mint32x8 {
	using T = u32;
	static constexpr size_t LEN = 8; // # of T in a vector
	
	class no_mask_tag {};
	static constexpr no_mask_tag NO_MASK{};
	
	using vec_t = __m256i;
	using mask_t = mint32x8;
	static vec_t mod_vec;
	
	using prefix_mask_generator_t = avx2_prefix_mask_generator_t<mint32x8>;
	
	vec_t val;
	mint32x8 () : val(_mm256_setzero_si256()) {}
	mint32x8 (u32 x) : val(_mm256_set1_epi32(x)) {} // x should be less than any of the mods
	mint32x8 (const vec_t &val) : val(val) {}
	mint32x8 (const mint32x8 &rhs) = default;
	mint32x8 & operator += (const mint32x8 &rhs) {
		val = _mm256_add_epi32(val, rhs.val);
		val = _mm256_min_epu32(val, _mm256_sub_epi32(val, mod_vec));
		return *this;
	}
	mint32x8 & operator -= (const mint32x8 &rhs) {
		val = _mm256_sub_epi32(val, rhs.val);
		val = _mm256_min_epu32(val, _mm256_add_epi32(val, mod_vec));
		return *this;
	}
	mint32x8 operator + (const mint32x8 &rhs) const { return mint32x8(*this) += rhs; }
	mint32x8 operator - (const mint32x8 &rhs) const { return mint32x8(*this) -= rhs; }
	
	static mint32x8 loadu(const u32 *ptr) { return { _mm256_loadu_si256((vec_t *) ptr) }; }
	static mint32x8 loadu(const mint32x8 *ptr) { return { _mm256_loadu_si256(&ptr->val) }; }
	void storeu(u32 *ptr) const { _mm256_storeu_si256((vec_t *) ptr, val); }
	void storeu_masked(u32 *ptr, mask_t mask) const { _mm256_maskstore_epi32((int *) ptr, mask.val, val); }
	void storeu_masked(u32 *ptr, no_mask_tag) const { storeu(ptr); } // for common interface with actual masked store
	
	static mint32x8 zero() { return { _mm256_set1_epi32(0) }; }
	static void set_mod(u32 mod) { mod_vec = _mm256_set1_epi32(mod); }
	u32 operator [] (int i) const {
		u32 tmp[LEN];
		storeu(tmp);
		return tmp[i];
	}
	friend std::ostream & operator << (std::ostream &stream, const mint32x8 &rhs) {
		stream << "[";
		for (size_t i = 0; i < LEN; i++) {
			if (i) stream << ", ";
			stream << rhs[i];
		}
		stream << "]";
		return stream;
	}
};
mint32x8::vec_t mint32x8::mod_vec;





// returns {dot(a, b), dot(b, b)} mod `mod`, multithreaded by T threads
std::pair<u32, u32> inner_product_mod_ab_bb(const std::vector<u32> &a, const std::vector<u32> &b, u32 mod, size_t len = -1, int T = 1) {
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size());
	
	u32 res_ab = 0;
	u32 res_bb = 0;
	std::mutex lock;
	auto thread_func = [&] (size_t thread_id) {
		size_t start = thread_id * len / T;
		size_t end = (thread_id + 1) * len / T;
		size_t i = start;
		
		using vec_t = __m256i;
		vec_t zero = _mm256_set1_epi64x(0);
		vec_t sum_ab = _mm256_set1_epi64x(0);
		vec_t sum_bb = _mm256_set1_epi64x(0);
		vec_t mod_large_vec = _mm256_set1_epi64x((u64) mod * mod * 2);
		for (; i + 7 < end; i += 8) {
			vec_t ax = _mm256_loadu_si256((vec_t *) (a.data() + i));
			vec_t bx = _mm256_loadu_si256((vec_t *) (b.data() + i));
			vec_t prod_ab0 = _mm256_mul_epu32(ax, bx);
			vec_t prod_ab1 = _mm256_mul_epu32(_mm256_srli_epi64(ax, 32), _mm256_srli_epi64(bx, 32));
			vec_t prod_bb0 = _mm256_mul_epu32(bx, bx);
			vec_t prod_bb1 = _mm256_mul_epu32(_mm256_srli_epi64(bx, 32), _mm256_srli_epi64(bx, 32));
			
			auto reduce_large = [&] (vec_t a) {
				vec_t subed = _mm256_sub_epi64(a, mod_large_vec);
				return _mm256_blendv_epi8(a, subed, _mm256_cmpgt_epi64(zero, a));
			};
			
			sum_ab = reduce_large(_mm256_add_epi64(sum_ab, _mm256_add_epi64(prod_ab0, prod_ab1)));
			sum_bb = reduce_large(_mm256_add_epi64(sum_bb, _mm256_add_epi64(prod_bb0, prod_bb1)));
		}
		auto vec_sum = [&] (vec_t sum_vec) {
			u64 buf[4];
			_mm256_storeu_si256((vec_t *) buf, sum_vec);
			u32 res = 0;
			for (int i = 0; i < 4; i++) res = (res + buf[i]) % mod;
			return res;
		};
		// add up to global results
		lock.lock();
		res_ab = (res_ab + vec_sum(sum_ab)) % mod;
		res_bb = (res_bb + vec_sum(sum_bb)) % mod;
		for (; i < end; i++) {
			res_ab = (res_ab + (u64) a[i] * b[i]) % mod;
			res_bb = (res_bb + (u64) b[i] * b[i]) % mod;
		}
		lock.unlock();
	};
	
	std::vector<std::thread> threads;
	for (int i = 0; i < T; i++) threads.push_back(std::thread(thread_func, i));
	for (auto &thread : threads) thread.join();
	return {res_ab, res_bb};
}
u32 inner_product_mod32(const std::vector<u32> &a, const std::vector<u32> &b, u32 mod, size_t len = -1) {
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size());
	size_t i;
	using vec_t = __m256i;
	vec_t zero = _mm256_set1_epi64x(0);
	vec_t sum = _mm256_set1_epi64x(0);
	vec_t mod_large_vec = _mm256_set1_epi64x((u64) mod * mod * 2);
	for (i = 0; i + 7 < len; i += 8) {
		vec_t ax = _mm256_loadu_si256((vec_t *) (a.data() + i));
		vec_t bx = _mm256_loadu_si256((vec_t *) (b.data() + i));
		vec_t prod0 = _mm256_mul_epu32(ax, bx);
		vec_t prod1 = _mm256_mul_epu32(_mm256_srli_epi64(ax, 32), _mm256_srli_epi64(bx, 32));
		sum = _mm256_add_epi64(sum, _mm256_add_epi64(prod0, prod1));
		vec_t sum_subed = _mm256_sub_epi64(sum, mod_large_vec);
		sum = _mm256_blendv_epi8(sum, sum_subed, _mm256_cmpgt_epi64(zero, sum));
	}
	u64 buf[4];
	_mm256_storeu_si256((vec_t *) buf, sum);
	u32 res = 0;
	for (int i = 0; i < 4; i++) res = (res + buf[i]) % mod;
	for (; i < len; i++) res = (res + (u64) a[i] * b[i]) % mod;
	return res;
}

