#pragma once
#include "types.h"
#include "avx_mask_generator.h"
#include <immintrin.h>
#include <mutex>
#include <thread>

struct mint16x16 {
	using T = u16;
	static constexpr size_t LEN = 16; // # of T in a vector
	
	class no_mask_tag {};
	static constexpr no_mask_tag NO_MASK{};
	
	using vec_t = __m256i;
	using mask_t = mint16x16;
	static vec_t mod_vec;
	
	using prefix_mask_generator_t = avx2_prefix_mask_generator_t<mint16x16>;
	
	vec_t val;
	mint16x16 () : val(_mm256_setzero_si256()) {}
	mint16x16 (u16 x) : val(_mm256_set1_epi16(x)) {} // x should be less than any of the mods
	mint16x16 (const vec_t &val) : val(val) {}
	mint16x16 (const mint16x16 &rhs) = default;
	mint16x16 & operator += (const mint16x16 &rhs) {
		val = _mm256_add_epi16(val, rhs.val);
		val = _mm256_min_epu16(val, _mm256_sub_epi16(val, mod_vec));
		return *this;
	}
	mint16x16 & operator -= (const mint16x16 &rhs) {
		val = _mm256_sub_epi16(val, rhs.val);
		val = _mm256_min_epu16(val, _mm256_add_epi16(val, mod_vec));
		return *this;
	}
	mint16x16 operator + (const mint16x16 &rhs) const { return mint16x16(*this) += rhs; }
	mint16x16 operator - (const mint16x16 &rhs) const { return mint16x16(*this) -= rhs; }
	
	static mint16x16 loadu(const u16 *ptr) { return { _mm256_loadu_si256((vec_t *) ptr) }; }
	static mint16x16 loadu(const mint16x16 *ptr) { return { _mm256_loadu_si256(&ptr->val) }; }
	void storeu(u16 *ptr) const { _mm256_storeu_si256((vec_t *) ptr, val); }
	void storeu_masked(u16 *ptr, mask_t mask) const { // avx2 doesn't have mask-store for u16
		u16 buf[LEN], mask_buf[LEN];
		storeu(buf);
		mask.storeu(mask_buf);
		for (size_t i = 0; i < LEN; i++) if (mask_buf[i]) ptr[i] = buf[i];
	}
	void storeu_masked(u16 *ptr, no_mask_tag) const { storeu(ptr); } // for common interface with actual masked store
	
	static mint16x16 zero() { return { _mm256_set1_epi16(0) }; }
	static void set_mod(u16 mod) { mod_vec = _mm256_set1_epi16(mod); }
	u16 operator [] (int i) const {
		u16 tmp[LEN];
		storeu(tmp);
		return tmp[i];
	}
	friend std::ostream & operator << (std::ostream &stream, const mint16x16 &rhs) {
		stream << "[";
		for (size_t i = 0; i < LEN; i++) {
			if (i) stream << ", ";
			stream << rhs[i];
		}
		stream << "]";
		return stream;
	}
};

mint16x16::vec_t mint16x16::mod_vec;



struct mint16x32 {
	using T = u16;
	static constexpr size_t LEN = 32; // # of T in a vector
	
	class no_mask_tag {};
	static constexpr no_mask_tag NO_MASK{};
	
	using vec_t = __m512i;
	using mask_t = __mmask32;
	static vec_t mod_vec;
	
	using prefix_mask_generator_t = avx512_prefix_mask_gen_t<mint16x32>;
	
	vec_t val;
	mint16x32 () : val(_mm512_setzero_si512()) {}
	mint16x32 (u16 x) : val(_mm512_set1_epi16(x)) {} // x should be less than any of the mods
	mint16x32 (const vec_t &val) : val(val) {}
	mint16x32 (const mint16x32 &rhs) = default;
	mint16x32 & operator += (const mint16x32 &rhs) {
		val = _mm512_add_epi16(val, rhs.val);
		val = _mm512_min_epu16(val, _mm512_sub_epi16(val, mod_vec));
		return *this;
	}
	mint16x32 & operator -= (const mint16x32 &rhs) {
		val = _mm512_sub_epi16(val, rhs.val);
		val = _mm512_min_epu16(val, _mm512_add_epi16(val, mod_vec));
		return *this;
	}
	mint16x32 operator + (const mint16x32 &rhs) const { return mint16x32(*this) += rhs; }
	mint16x32 operator - (const mint16x32 &rhs) const { return mint16x32(*this) -= rhs; }
	
	static mint16x32 loadu(const u16 *ptr) { return { _mm512_loadu_si512((vec_t *) ptr) }; }
	static mint16x32 loadu(const mint16x32 *ptr) { return { _mm512_loadu_si512(&ptr->val) }; }
	void storeu(u16 *ptr) const { _mm512_storeu_si512((vec_t *) ptr, val); }
	void storeu_masked(u16 *ptr, mask_t mask) const { _mm512_mask_storeu_epi16(ptr, mask, val); }
	void storeu_masked(u16 *ptr, no_mask_tag) const { storeu(ptr); } // for common interface with actual masked store
	
	static mint16x32 zero() { return { _mm512_set1_epi16(0) }; }
	static void set_mod(u16 mod) { mod_vec = _mm512_set1_epi16(mod); }
	u16 operator [] (int i) const {
		u16 tmp[LEN];
		storeu(tmp);
		return tmp[i];
	}
	friend std::ostream & operator << (std::ostream &stream, const mint16x32 &rhs) {
		stream << "[";
		for (size_t i = 0; i < LEN; i++) {
			if (i) stream << ", ";
			stream << rhs[i];
		}
		stream << "]";
		return stream;
	}
};
mint16x32::vec_t mint16x32::mod_vec;





// returns {dot(a, b), dot(b, b)} mod `mod`, multithreaded by T threads
// assuming mod <= 2^15
std::pair<u16, u16> inner_product_mod_ab_bb(const std::vector<u16> &a, const std::vector<u16> &b, u16 mod, size_t len = -1, int T = 1) {
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size());
	
	
	u16 res_ab = 0;
	u16 res_bb = 0;
	std::mutex lock;
	auto thread_func = [&] (size_t thread_id) {
		size_t start = thread_id * len / T;
		size_t end = (thread_id + 1) * len / T;
		size_t i = start;
		
		using vec_t = __m256i;
		vec_t zero = _mm256_set1_epi32(0);
		vec_t sum_ab = _mm256_set1_epi32(0);
		vec_t sum_bb = _mm256_set1_epi32(0);
		vec_t mod_large_vec = _mm256_set1_epi32((u32) mod * mod * 2);
		for (; i + 15 < end; i += 16) {
			vec_t ax = _mm256_loadu_si256((vec_t *) (a.data() + i));
			vec_t bx = _mm256_loadu_si256((vec_t *) (b.data() + i));
			vec_t prod_ab = _mm256_madd_epi16(ax, bx); // madd: multiply and horizontally-add
			vec_t prod_bb = _mm256_madd_epi16(bx, bx);
			
			auto reduce_large = [&] (vec_t a) {
				// if a is negative(which means it overflew), subtract mod_large
				vec_t subed = _mm256_sub_epi32(a, mod_large_vec);
				return _mm256_blendv_epi8(a, subed, _mm256_cmpgt_epi32(zero, a));
			};
			
			sum_ab = reduce_large(_mm256_add_epi32(sum_ab, prod_ab));
			sum_bb = reduce_large(_mm256_add_epi32(sum_bb, prod_bb));
		}
		auto vec_sum = [&] (vec_t sum_vec) {
			u32 buf[8];
			_mm256_storeu_si256((vec_t *) buf, sum_vec);
			u16 res = 0;
			for (int i = 0; i < 8; i++) res = (res + buf[i]) % mod;
			return res;
		};
		// add up to global results
		lock.lock();
		res_ab = (res_ab + vec_sum(sum_ab)) % mod;
		res_bb = (res_bb + vec_sum(sum_bb)) % mod;
		for (; i < end; i++) {
			res_ab = (res_ab + (u32) a[i] * b[i]) % mod;
			res_bb = (res_bb + (u32) b[i] * b[i]) % mod;
		}
		lock.unlock();
	};
	
	std::vector<std::thread> threads;
	for (int i = 0; i < T; i++) threads.push_back(std::thread(thread_func, i));
	for (auto &thread : threads) thread.join();
	return {res_ab, res_bb};
}
u16 inner_product_mod16(const std::vector<u16> &a, const std::vector<u16> &b, u16 mod, size_t len = -1) {
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size());
	size_t i;
	using vec_t = __m256i;
	vec_t zero = _mm256_set1_epi32(0);
	vec_t sum = _mm256_set1_epi32(0);
	vec_t mod_large_vec = _mm256_set1_epi32((u64) mod * mod * 2);
	for (i = 0; i + 15 < len; i += 16) {
		vec_t ax = _mm256_loadu_si256((vec_t *) (a.data() + i));
		vec_t bx = _mm256_loadu_si256((vec_t *) (b.data() + i));
		vec_t prod_ab = _mm256_madd_epi16(ax, bx);
		auto reduce_large = [&] (vec_t a) {
			vec_t subed = _mm256_sub_epi32(a, mod_large_vec);
			return _mm256_blendv_epi8(a, subed, _mm256_cmpgt_epi32(zero, a));
		};
		sum = reduce_large(_mm256_add_epi32(sum, prod_ab));
	}
	u32 buf[8];
	_mm256_storeu_si256((vec_t *) buf, sum);
	u16 res = 0;
	for (int i = 0; i < 8; i++) res = (res + buf[i]) % mod;
	for (; i < len; i++) res = (res + (u64) a[i] * b[i]) % mod;
	return res;
}

