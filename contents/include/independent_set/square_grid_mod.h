#include <vector>
#include <map>
#include <cassert>
#include <type_traits>
#include "types.h"
#include "fibonacci.h"
#include "thread_info.h"

/*
	Class for efficiently multiplying transition matrix for independent set counting on square grid graph.
	
	 - IS: independent set
	 - square: square grid; two squares are adjacent iff their Manhattan/L1 distance is 1
	 - mod: all values are Z/mZ where the size(# of bits) of mod `m` depends on the N of mintNxM
	
	mintNxM: class providing implementation of vectorized modint(Z/mZ)
		e.g. mintNxM contains eight Z/mZ elements where the mod is 32-bit (acutally up to 2^31).
			It is probably implemented using AVX2 because the vector size is 256 bit.
	
	The function `transition_one_row` is thread-unsafe, though it is okay to 
		call it from different *process*.
*/
template<typename mintNxM>
struct transitioner_IS_square_grid_mod {
	using T = typename mintNxM::T; // underlying scalar type (e.g. u32 if mintNxM is mintNxM)
	using mask_type = typename mintNxM::mask_t;
	static constexpr size_t VEC_LEN = mintNxM::LEN;
	
	T mod;
	std::vector<size_t> fib;
	T mint_add_scalar(T x, T y) {
		x += y;
		if (x >= mod) x -= mod;
		return x;
	}
	
	
	/*
		Masks are used when dealing with incomplete vectors(i.e. a long vector whose length not divisible by VEC_LEN)
		We often load the entire NxM-bit vector even if its partially out-of-bound instead of loading with a mask. (*1)
		Instead, we use masks when storing them(=storing only some number of elements from the beginning of the vector).
		
		(*1) This is why we allocate VEC_LEN more (useless) elements for vectors(see: get_buf_size())
	*/
	typename mintNxM::prefix_mask_generator_t prefix_gen;
	// returns the mask enabling exactly k elements from the beginning of the vector (like 11111000 if k=5)
	mask_type mask_for_first_k_elements(int k) { return prefix_gen.mask_for_first_k_elements(k); }
	
	
	// Constructor
	transitioner_IS_square_grid_mod (int n) : fib(fibonacci_sequence<size_t>(n)) {}
	
	void set_mod(T mod) {
		this->mod = mod;
		mintNxM::set_mod(mod);
	}
	size_t get_buf_size(size_t n) {
		assert(n < fib.size());
		return fib[n] + VEC_LEN; // add padding to avoid illegal access when loading a vector
	}
	std::vector<T> get_initial_buffer(int n) {
		std::vector<T> res(get_buf_size(n), 1);
		// the padding should be zero..?
		for (size_t i = fib[n]; i < get_buf_size(n); i++) res[i] = 0;
		return res;
	}
	
	
	/*
		Executes k-bit transition:
			input buf: fib[k]-long vector of T
			output buf[i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[j]
				where addition is done mod `mod`
		
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
				buf = concat([A, B, C]).
			Then deciding the new uppermost bit(s) to either 0 or 10 corresponds:
				buf <- concat([A + C, B, A + B]) (additions are vector addition, treating out-of-bound element as zero)
				 : can be done by swapping A and C, adding A to C and adding B to A
			We can recursively do this on [A + C, B] and [A + B] (of size fib[n - 1] and fib[n - 2], resp.)
		
		transform_one_row could just call this function with k=n and it should work, but it is
		 - cache-inefficient (access the entire buffer every 1 or 2 depths)
		 - not very easy to parallelize(make it multi-thread)
	*/
	void transition(int k, T *buf) {
		auto add = [this] (T x, T y) { return mint_add_scalar(x, y); };
		
		if (k <= 0) return;
		if (k == 1) {
			std::swap(buf[0], buf[1]);
			buf[0] = add(buf[0], buf[1]);
			return;
		}
		if (k == 2) {
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
		
		/* // before vectorization
		for (size_t i = 0; i < fib[k - 3]; i++) {
			T a = buf[i];
			T b = buf[i + fib[k - 2]];
			T c = buf[i + fib[k - 1]];
			buf[i] = add(a, c);
			buf[i + fib[k - 1]] = add(a, b);
		}
		for (size_t i = fib[k - 3]; i < fib[k - 2]; i++) {
			T a = buf[i];
			T c = buf[i + fib[k - 1]];
			buf[i] = add(a, c);
			buf[i + fib[k - 1]] = a;
		}
		*/
		
		/*
		{
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a = mintNxM::loadu(&buf[i]);
				mintNxM b = mintNxM::loadu(&buf[i + fib[k - 2]]);
				mintNxM c = mintNxM::loadu(&buf[i + fib[k - 1]]);
				(a + c).storeu_masked(&buf[i], mask);
				(a + b).storeu_masked(&buf[i + fib[k - 1]], mask);
			};
			
			size_t i = 0;
			for (; i + VEC_LEN < fib[k - 3]; i += VEC_LEN) kernel(i, mintNxM::NO_MASK);
			kernel(i, mask_for_first_k_elements(fib[k - 3] - i));
		}
		{
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a = mintNxM::loadu(&buf[i]); \
				mintNxM c = mintNxM::loadu(&buf[i + fib[k - 1]]); \
				(a + c).storeu_masked(&buf[i], mask); \
				(a    ).storeu_masked(&buf[i + fib[k - 1]], mask); \
			};
			
			size_t i = fib[k - 3];
			for (; i + VEC_LEN < fib[k - 2]; i += VEC_LEN) kernel(i, mintNxM::NO_MASK);
			kernel(i, mask_for_first_k_elements(fib[k - 2] - i));
		}
		*/
		
		size_t unit5 = fib[std::max(0, k - 5)];
		size_t unit4 = fib[k - 4];
		size_t unit3 = fib[k - 3];
		size_t unit2 = fib[k - 2];
		size_t unit1 = fib[k - 1];
		{
			T *t0 = buf;
			T *t1 = buf + unit4;
			T *t2 = buf + unit3;
			T *t3 = buf + unit2;
			T *t4 = buf + unit2 + unit4;
			T *t5 = buf + unit1;
			T *t6 = buf + unit1 + unit4;
			T *t7 = buf + unit1 + unit3;
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a0 = mintNxM::loadu(&t0[i]);
				mintNxM a1 = mintNxM::loadu(&t1[i]);
				mintNxM a2 = mintNxM::loadu(&t2[i]);
				mintNxM a3 = mintNxM::loadu(&t3[i]);
				mintNxM a4 = mintNxM::loadu(&t4[i]);
				mintNxM a5 = mintNxM::loadu(&t5[i]);
				mintNxM a6 = mintNxM::loadu(&t6[i]);
				mintNxM a7 = mintNxM::loadu(&t7[i]);
				mintNxM a05 = a0 + a5;
				mintNxM a16 = a1 + a6;
				mintNxM a27 = a2 + a7;
				mintNxM a03 = a0 + a3;
				mintNxM a14 = a1 + a4;
				mintNxM a053 = a05 + a3;
				mintNxM a164 = a16 + a4;
				mintNxM a0527 = a05 + a27;
				mintNxM a05327 = a053 + a27;
				mintNxM a032 = a03 + a2;
				(a05327 + a164).storeu_masked(&t0[i], mask);
				(a05327       ).storeu_masked(&t1[i], mask);
				(a053   + a164).storeu_masked(&t2[i], mask);
				(a0527  + a16 ).storeu_masked(&t3[i], mask);
				(a0527        ).storeu_masked(&t4[i], mask);
				(a032   + a14 ).storeu_masked(&t5[i], mask);
				(a032         ).storeu_masked(&t6[i], mask);
				(a03    + a14 ).storeu_masked(&t7[i], mask);
			};
			size_t i = 0;
			for (; i + VEC_LEN < unit5; i += VEC_LEN) kernel(i, mintNxM::NO_MASK);
			kernel(i, mask_for_first_k_elements(unit5 - i));
		}
		{
			T *t0 = buf;
			T *t1 = buf + unit4 - unit5;
			T *t2 = buf + unit3;
			T *t3 = buf + unit2;
			T *t4 = buf + unit2 + unit4 - unit5;
			T *t5 = buf + unit1;
			T *t6 = buf + unit1 + unit4 - unit5;
			T *t7 = buf + unit1 + unit3;
			
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a0 = mintNxM::loadu(&t0[i]);
				mintNxM a1 = mintNxM::loadu(&t1[i]);
				mintNxM a2 = mintNxM::loadu(&t2[i]);
				mintNxM a3 = mintNxM::loadu(&t3[i]);
				mintNxM a4 = mintNxM::loadu(&t4[i]);
				mintNxM a5 = mintNxM::loadu(&t5[i]);
				mintNxM a6 = mintNxM::loadu(&t6[i]);
				mintNxM a7 = mintNxM::loadu(&t7[i]);
				mintNxM a05 = a0 + a5;
				mintNxM a27 = a2 + a7;
				mintNxM a03 = a0 + a3;
				mintNxM a053 = a05 + a3;
				mintNxM a032 = a03 + a2;
				mintNxM a05327 = a053 + a27;
				mintNxM a0527 = a05 + a27;
				(a05327       ).storeu_masked(&t0[i], mask);
				(a05327 + a1  ).storeu_masked(&t1[i], mask);
				(a053         ).storeu_masked(&t2[i], mask);
				(a0527        ).storeu_masked(&t3[i], mask);
				(a0527  + a4  ).storeu_masked(&t4[i], mask);
				(a032         ).storeu_masked(&t5[i], mask);
				(a032   + a6  ).storeu_masked(&t6[i], mask);
				(a03          ).storeu_masked(&t7[i], mask);
			};
			// for i in [unit5, unit4)
			size_t i = unit5;
			for (; i + VEC_LEN < unit4; i += VEC_LEN) kernel(i, mintNxM::NO_MASK);
			kernel(i, mask_for_first_k_elements(unit4 - i));
		}
		
		
		transition(k - 4, buf                );
		transition(k - 5, buf + unit4        );
		transition(k - 4, buf + unit3        );
		transition(k - 4, buf + unit2        );
		transition(k - 5, buf + unit2 + unit4);
		transition(k - 4, buf + unit1        );
		transition(k - 5, buf + unit1 + unit4);
		transition(k - 4, buf + unit1 + unit3);
	}
	/*
		Executes k-bit transition but with additional uppermost bit that turns off if (k-1)-th bit turns on by transition:
			buf: fib[k+1]-long vector of T
			output buf[0i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[0j] +
							 sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[1j] (only if i[k - 1] == 1)
			output buf[1i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[1j] (i[k - 1] == 0)
				where addition is done mod `mod`
		
	*/
	void transition_plus_one(int k, T *buf) {
		transition(k, buf);
		transition(k - 1, buf + fib[k]);
		
		// buf[fib[k-1] : fib[k]] += fib[fib[k] : fib[k]+fib[k-2]]
		auto kernel = [&] (size_t i, auto mask) {
			mintNxM x = mintNxM::loadu(&buf[i + fib[k - 1]]);
			mintNxM y = mintNxM::loadu(&buf[i + fib[k]]);
			(x + y).storeu_masked(&buf[i + fib[k - 1]], mask);
		};
		size_t i = 0;
		for (; i + VEC_LEN < fib[k - 2]; i += VEC_LEN) kernel(i, mintNxM::NO_MASK);
		kernel(i, mask_for_first_k_elements(fib[k - 2] - i));
	}
	
	
	
	/*
		Same as `transition` but each buf[i] is a long vector of T instead of a single T
		All long vectors (buf[0], buf[1], ...) must have the same length `n_elements`.
		
		Algorithm:
			Same as `transition` but we further split the buffer.
			Let
				T0 = buf[fibindex(0000*)] (meaining buf[fibindex(00000...0)...fibindex(00010...0)-1])
				T1 = buf[fibindex(00010*)]
				T2 = buf[fibindex(0010*)]
				T3 = buf[fibindex(0100*)]
				T4 = buf[fibindex(01010*)]
				T5 = buf[fibindex(1000*)]
				T6 = buf[fibindex(10010*)]
				T7 = buf[fibindex(1010*)]
			, so that
				buf = concat([t0, t1, t2, ..., t7])
			
			Then, when we decide the upper 4 or 5 bits to one of the above(000, 00010, 0010, 0100, ...)
				new buf = concat([
					sum_{j in 01234567} Tj,
					sum_{j in 02357} Tj,
					sum_{j in 013456} Tj, // !
					sum_{j in 012567} Tj, // !
					sum_{j in 0257} Tj,
					sum_{j in 01234} Tj, // !
					sum_{j in 023} Tj,
					sum_{j in 0134} Tj // !
				])
			Note that T0, T2, T3, T5, T7 are fib[k - 4]-long whereas others are fib[k - 5]-long, and
				the out-of-bound elements are treated as zero in the above sums.
			Also, sumed array that should fit into a short fib[k - 5]-length section(lines with // !) are round wrapped:
				let tmp = sum_{j in 013456} Tj, and tmp[i] (i >= fib[k-5]) are round-wrapped and summed into tmp[i - fib[k-5]].
			This corresponds to the fact that when we decide the upper 5 bits to be 00010, for example, the fifth uppermost bit of i(old index)
				turns off in j whether it was 0 or 1.
			Each Tj is `n_elements`-long vector of T.
	*/
	void transition_vecs(int k, T **buf, u64 n_elements) {
		size_t n_full_vecs = (n_elements - 1) / VEC_LEN;
		size_t last_vec_n_elements = n_elements - n_full_vecs * VEC_LEN; // # of elements in the last vector
		transition_vecs_internal(k, buf, n_full_vecs, mask_for_first_k_elements(last_vec_n_elements));
	}
	void transition_vecs_internal(int k, T **buf, size_t n_full_vecs, mask_type last_vec_mask) {
		if (k <= 0) return;
		
		// std::cerr << "!!!!!!!! " << k << " !!!!!!!!!!!!!!"  << std::endl;
		if (k == 1) { // 1 addition
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a0 = mintNxM::loadu(&buf[0][i]);
				mintNxM a1 = mintNxM::loadu(&buf[1][i]);
				(a0 + a1).storeu_masked(&buf[0][i], mask);
				(a0     ).storeu_masked(&buf[1][i], mask);
			};
			for (size_t i = 0; i < n_full_vecs; i++) kernel(i * VEC_LEN, mintNxM::NO_MASK);
			kernel(n_full_vecs * VEC_LEN, last_vec_mask);
			return;
		}
		if (k == 2) { // 3 additions
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a0 = mintNxM::loadu(&buf[0][i]);
				mintNxM a1 = mintNxM::loadu(&buf[1][i]);
				mintNxM a2 = mintNxM::loadu(&buf[2][i]);
				mintNxM a01 = a0 + a1;
				(a01 + a2).storeu_masked(&buf[0][i], mask);
				(a0  + a2).storeu_masked(&buf[1][i], mask);
				(a01     ).storeu_masked(&buf[2][i], mask);
			};
			for (size_t i = 0; i < n_full_vecs; i++) kernel(i * VEC_LEN, mintNxM::NO_MASK);
			kernel(n_full_vecs * VEC_LEN, last_vec_mask);
			return;
		}
		if (k == 3) { // 7 additions
			auto kernel = [&] (size_t i, auto mask) {
				mintNxM a0 = mintNxM::loadu(&buf[0][i]);
				mintNxM a1 = mintNxM::loadu(&buf[1][i]);
				mintNxM a2 = mintNxM::loadu(&buf[2][i]);
				mintNxM a3 = mintNxM::loadu(&buf[3][i]);
				mintNxM a4 = mintNxM::loadu(&buf[4][i]);
				mintNxM a03 = a0 + a3;
				mintNxM a14 = a1 + a4;
				mintNxM a02 = a0 + a2;
				mintNxM a032 = a03 + a2;
				(a032 + a14).storeu_masked(&buf[0][i], mask);
				(a032      ).storeu_masked(&buf[1][i], mask);
				(a03  + a14).storeu_masked(&buf[2][i], mask);
				(a02  + a1 ).storeu_masked(&buf[3][i], mask);
				(a02       ).storeu_masked(&buf[4][i], mask);
			};
			for (size_t i = 0; i < n_full_vecs; i++) kernel(i * VEC_LEN, mintNxM::NO_MASK);
			kernel(n_full_vecs * VEC_LEN, last_vec_mask);
			return;
		}
		
		size_t unit5 = fib[std::max(0, k - 5)];
		size_t unit4 = fib[k - 4];
		size_t unit3 = fib[k - 3];
		size_t unit2 = fib[k - 2];
		size_t unit1 = fib[k - 1];
		for (size_t i = 0; i < unit5; i++) {
			T *t0 = buf[i];
			T *t1 = buf[unit4 + i];
			T *t2 = buf[unit3 + i];
			T *t3 = buf[unit2 + i];
			T *t4 = buf[unit2 + unit4 + i];
			T *t5 = buf[unit1 + i];
			T *t6 = buf[unit1 + unit4 + i];
			T *t7 = buf[unit1 + unit3 + i];
			auto kernel = [&] (size_t j, auto mask) {
				mintNxM a0 = mintNxM::loadu(&t0[j]);
				mintNxM a1 = mintNxM::loadu(&t1[j]);
				mintNxM a2 = mintNxM::loadu(&t2[j]);
				mintNxM a3 = mintNxM::loadu(&t3[j]);
				mintNxM a4 = mintNxM::loadu(&t4[j]);
				mintNxM a5 = mintNxM::loadu(&t5[j]);
				mintNxM a6 = mintNxM::loadu(&t6[j]);
				mintNxM a7 = mintNxM::loadu(&t7[j]);
				mintNxM a05 = a0 + a5;
				mintNxM a16 = a1 + a6;
				mintNxM a27 = a2 + a7;
				mintNxM a03 = a0 + a3;
				mintNxM a14 = a1 + a4;
				mintNxM a053 = a05 + a3;
				mintNxM a164 = a16 + a4;
				mintNxM a0527 = a05 + a27;
				mintNxM a05327 = a053 + a27;
				mintNxM a032 = a03 + a2;
				(a05327 + a164).storeu_masked(&t0[j], mask);
				(a05327       ).storeu_masked(&t1[j], mask);
				(a053   + a164).storeu_masked(&t2[j], mask);
				(a0527  + a16 ).storeu_masked(&t3[j], mask);
				(a0527        ).storeu_masked(&t4[j], mask);
				(a032   + a14 ).storeu_masked(&t5[j], mask);
				(a032         ).storeu_masked(&t6[j], mask);
				(a03    + a14 ).storeu_masked(&t7[j], mask);
			};
			for (size_t j = 0; j < n_full_vecs; j++) kernel(j * VEC_LEN, mintNxM::NO_MASK);
			kernel(n_full_vecs * VEC_LEN, last_vec_mask);
		}
		for (size_t i = unit5; i < unit4; i++) {
			T *t0 = buf[i];
			T *t1 = buf[unit4 + (i-unit5)]; // destination only
			T *t2 = buf[unit3 + i];
			T *t3 = buf[unit2 + i];
			T *t4 = buf[unit2 + unit4 + (i-unit5)]; // destination only
			T *t5 = buf[unit1 + i];
			T *t6 = buf[unit1 + unit4 + (i-unit5)]; // destination only
			T *t7 = buf[unit1 + unit3 + i];
			
			auto kernel = [&] (size_t j, auto mask) {
				mintNxM a0 = mintNxM::loadu(&t0[j]);
				mintNxM a1 = mintNxM::loadu(&t1[j]);
				mintNxM a2 = mintNxM::loadu(&t2[j]);
				mintNxM a3 = mintNxM::loadu(&t3[j]);
				mintNxM a4 = mintNxM::loadu(&t4[j]);
				mintNxM a5 = mintNxM::loadu(&t5[j]);
				mintNxM a6 = mintNxM::loadu(&t6[j]);
				mintNxM a7 = mintNxM::loadu(&t7[j]);
				mintNxM a05 = a0 + a5;
				mintNxM a27 = a2 + a7;
				mintNxM a03 = a0 + a3;
				mintNxM a053 = a05 + a3;
				mintNxM a032 = a03 + a2;
				mintNxM a05327 = a053 + a27;
				mintNxM a0527 = a05 + a27;
				(a05327       ).storeu_masked(&t0[j], mask);
				(a05327 + a1  ).storeu_masked(&t1[j], mask);
				(a053         ).storeu_masked(&t2[j], mask);
				(a0527        ).storeu_masked(&t3[j], mask);
				(a0527  + a4  ).storeu_masked(&t4[j], mask);
				(a032         ).storeu_masked(&t5[j], mask);
				(a032   + a6  ).storeu_masked(&t6[j], mask);
				(a03          ).storeu_masked(&t7[j], mask);
			};
			for (size_t j = 0; j < n_full_vecs; j++) kernel(j * VEC_LEN, mintNxM::NO_MASK);
			kernel(n_full_vecs * VEC_LEN, last_vec_mask);
		}
		
		
		transition_vecs_internal(k - 4, buf                , n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 5, buf + unit4        , n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 4, buf + unit3        , n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 4, buf + unit2        , n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 5, buf + unit2 + unit4, n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 4, buf + unit1        , n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 5, buf + unit1 + unit4, n_full_vecs, last_vec_mask);
		transition_vecs_internal(k - 4, buf + unit1 + unit3, n_full_vecs, last_vec_mask);
	}
	/*
		Executes k-bit transform but with additional uppermost bit that turns off if (k-1)-th bit turns on by transformation:
			buf: fib[k+1]-long
			output buf[0i] = sum_{j: (fibstr(i) & fibstr(j)) == 0} buf[0j] +
							 sum_{j: (fibstr(i) & fibstr(j)) == 0, fibstr(j)[k-1] == 1} buf[1j]
			output buf[1i] = sum_{j: (fibstr(i) & fibstr(j)) == 0, fibstr(j)[k-1] == 0} buf[1j]
				where addition is done via function `add`
		
	*/
	void transition_vecs_plus_one(int k, T **buf, u64 n_elements) {
		transition_vecs(k,     buf,          n_elements);
		transition_vecs(k - 1, buf + fib[k], n_elements);
		
		for (size_t i = 0; i < fib[k - 2]; i++) {
			T *t0 = buf[i + fib[k - 1]];
			T *t1 = buf[i + fib[k]];
			auto kernel = [&] (size_t j, auto mask) {
				mintNxM x = mintNxM::loadu(&t0[j]);
				mintNxM y = mintNxM::loadu(&t1[j]);
				(x + y).storeu_masked(&t0[j], mask);
			};
			size_t j = 0;
			for (; j + VEC_LEN < n_elements; j += VEC_LEN) kernel(j, mintNxM::NO_MASK);
			kernel(j, mask_for_first_k_elements(n_elements - j));
		}
	}
	
	
	
	void transform_first_C_squares_naive(int n, int C, const std::vector<T> &dp, std::vector<T> &next, thread_info_t thread_info) {
		if (thread_info.my_id != 0) return; // only work on thread 0(if there are multiple)
		assert(dp.size() == get_buf_size(n));
		std::map<u64, size_t> rev;
		for_each_fibonacci<u64>(n, [&] (size_t i, u64 str) { rev[str] = i; });
		
		next.assign(get_buf_size(n), 0);
		for_each_fibonacci<u64>(n, [&] (size_t i, u64 cur_str) {
			for_each_fibonacci<u64>(C, [&] (size_t j, u64 new_low) {
				(void) j;
				if (cur_str & new_low) return;
				
				u64 new_val = (cur_str & ~((1 << C) - 1)) | new_low;
				if ((new_low >> (C - 1) & 1) && (cur_str >> C & 1)) new_val ^= 1ULL << C;
				
				auto itr = rev.find(new_val);
				assert(itr != rev.end());
				size_t index = itr->second;
				next[index] += dp[i];
				if (next[index] >= mod) next[index] -= mod;
			});
		});
	}
	void transform_first_C_squares(int n, int C, const std::vector<T> &dp, std::vector<T> &next, thread_info_t thread_info) {
		// enumerate upper n - C - 1 digits
		u64 head = 0;
		
		for_each_fibonacci<u64>(std::max(0, n - C - 1), [&] (size_t i, u64 upper) {
			if (n == C) upper = 1; // to fix the middle 1 bit to zero
			
			size_t len = fib[(upper & 1) ? C : C + 1];
			if (i % thread_info.total_num == thread_info.my_id) {
				size_t i = 0;
				for (; i + VEC_LEN < len; i += VEC_LEN) mintNxM::loadu(&dp[head + i]).storeu(&next[head + i]);
				mintNxM::loadu(&dp[head + i]).storeu_masked(&next[head + i], mask_for_first_k_elements(len - i));
				
				if (upper & 1) transition(C, &next[head]);
				else transition_plus_one(C, &next[head]);
				
			}
			head += len;
			
		});
		assert(head == fib[n]);
	}
	
	
	void transform_C_squares_naive(int n, int column, int C, std::vector<T> &dp, thread_info_t thread_info) {
		if (thread_info.my_id != 0) return; // only work on thread 0(if any)
		std::map<u64, size_t> rev;
		for_each_fibonacci<u64>(n, [&] (size_t i, u64 str) { rev[str] = i; });
		
		std::vector<T> next(dp.size());
		
		for_each_fibonacci<u64>(n, [&] (size_t i, u64 cur_str) {
			for_each_fibonacci<u64>(C, [&] (size_t j, u64 new_low) {
				(void) j;
				if (cur_str >> column & new_low) return;
				
				u64 new_val = (cur_str & ~(((1ULL << C) - 1) << column)) | (new_low << column);
				if ((new_low >> (C - 1) & 1) && (cur_str >> (C + column) & 1)) new_val ^= 1ULL << (C + column);
				if ((new_low & 1) && (cur_str >> (column - 1) & 1)) return;
				
				auto itr = rev.find(new_val);
				assert(itr != rev.end());
				size_t index = itr->second;
				next[index] += dp[i];
				if (next[index] >= mod) next[index] -= mod;
			});
		});
		dp = next;
	}
	// column -> column + C
	// in-place
	void transform_C_squares(int n, int column, int C, std::vector<T> &dp, thread_info_t thread_info) {
		constexpr int D = 9;
		assert(column > D);
		assert(column + C <= n);
		
		
		std::vector<T *> ptrs_to_vecs(fib[C + 1]);
		
		
		/*
			  low <- [D] [column - D] [C][1] [n - column - C - 1] -> high
			or 
			  low <- [D] [column - D] [C][1](virtual, 0)  [1](virtual, 1) -> high
			if column + C == n
		*/
		
		u64 head_ = 0;
		for_each_fibonacci<u64>(std::max(0, n - column - C - 1), [&] (size_t i, u64 up_upper) {
			// virtually set up_upper ([n-column-C-1]) to one to fix [1] to zero
			if (column + C == n) up_upper = 1;
			
			u64 j_head_offset = 0;
			for_each_fibonacci<u64>(column - D, [&] (size_t j, u64 low_upper) {
				int d = D;
				if (low_upper & 1) d--;
				if ((i + j) % thread_info.total_num == thread_info.my_id) {
					int c = C;
					bool do_plus_one_transform = true;
					if (up_upper & 1) do_plus_one_transform = false;
					if (low_upper >> (column - D - 1) & 1) c--; // the lowest bit of [C] is fixed to zero in this case
					
					
					// prepare pointers to the long vectors
					u64 in_head = head_ + j_head_offset;
					for_each_fibonacci<u64>(c + do_plus_one_transform, [&] (size_t k, u64 up_lower) {
						// load dp[in_head..in_head+fib[d]-1] by vector
						ptrs_to_vecs[k] = &dp[in_head];
						
						in_head += fib[column +
							(low_upper >> (column - D - 1) & 1) - // entire low_upper is offset by one in this case
							(up_lower & 1)];
					});
					// transform
					if (do_plus_one_transform) transition_vecs_plus_one(c, ptrs_to_vecs.data(), fib[d]);
					else transition_vecs(c, ptrs_to_vecs.data(), fib[d]);
					
					assert(in_head == head_ + j_head_offset + fib[column + C + 1 - (up_upper & 1)]);
				}
				j_head_offset += fib[d];
			});
			head_ += fib[column + C + 1 - (up_upper & 1)];
		});
	}
	// apply the transition matrix to `dp` and write the resulting vector in `next`.
	void transition_one_row(int n, const std::vector<T> &dp, std::vector<T> &next, T mod, size_t n_threads, bool verbose = false) {
		this->set_mod(mod);
		
		assert(dp.size() == get_buf_size(n));
		assert(next.size() == get_buf_size(n));
		
		// for smaller n where main algorithm does not work, use less optimized way
		if (n <= 16) {
			if (verbose) {
				fprintf(stderr, "[INFO] transition_one_row: Using small-n fallback.\n");
			}
			next = dp;
			transition(n, next.data());
			return;
		}
		
		int s = std::max(12, (int) (n * 0.4));
		std::vector<int> t = {0, s, s + (n - s) / 2, n};
		int n_steps = t.size() - 1;
		for (int i = 0; i < n_steps; i++) {
			int cur_column = t[i];
			int next_column = t[i + 1];
			
			auto t0 = Timer::get();
			// Spawn threads
			std::vector<std::thread> threads;
			for (size_t my_tid = 0; my_tid < n_threads; my_tid++) threads.push_back(std::thread([&,my_tid] () {
				if (i == 0) transform_first_C_squares(n, next_column - cur_column, dp, next, {my_tid, n_threads});
				else transform_C_squares(n, cur_column, next_column - cur_column, next, {my_tid, n_threads});
			}));
			for (auto &thread : threads) thread.join();
			
			auto t1 = Timer::get();
			if (verbose) {
				fprintf(stderr, "#%d: %2d columns: %.1f ms\n", (int) i, next_column - cur_column, Timer::diff_ms(t0, t1));
			}
		}
	}
};
