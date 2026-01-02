#pragma once
#include "types.h"
#include <immintrin.h>


// AVX2 uses __m256i as the mask: uppermost bit of each integer in the vector represents the mask
// LEN: vector length
// T: integer type (scalar type)
template<typename mintNxM> struct avx2_prefix_mask_generator_t {
	static constexpr size_t LEN = mintNxM::LEN;
	using T = typename mintNxM::T;
	
	T buf_for_mask[LEN * 2]; // the first half are all -1, the latter half 0
	avx2_prefix_mask_generator_t () {
		for (size_t i = 0; i < LEN; i++) buf_for_mask[i] = (T) -1; // 0b1111..
		for (size_t i = 0; i < LEN; i++) buf_for_mask[LEN + i] = 0; // 0b0000..
	}
	mintNxM mask_for_first_k_elements(int k) { // 0 <= k <= VEC_LEN
		return mintNxM::loadu(buf_for_mask + LEN - k);
	}
};

// AVX512 uses dedicated type (e.g. __mmask16) for masking
template<typename mintNxM> struct avx512_prefix_mask_gen_t {
	avx512_prefix_mask_gen_t () {}
	typename mintNxM::mask_t mask_for_first_k_elements(int k) { // 0 <= k <= VEC_LEN
		return (typename mintNxM::mask_t) ((1ULL << k) - 1); // supports up to VEC_LEN <= 32
	}
};
