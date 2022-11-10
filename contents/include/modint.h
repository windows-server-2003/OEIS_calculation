#pragma once
#include "types.h"

template<typename base_type, typename mult_type, base_type mod> struct mint {
	using base_t = base_type;
	static base_type get_mod() { return mod; }
	base_type val = 0;
	mint () = default;
	mint (base_type val) : val(val) {}
	mint (const mint &rhs) = default;
	mint & operator += (const mint &rhs) { val += rhs.val; if (val >= mod) val -= mod; return *this; }
	mint & operator -= (const mint &rhs) { val = val < rhs.val ? val - rhs.val + mod : val - rhs.val; return *this; }
	mint & operator *= (const mint &rhs) { val = (mult_type) val * rhs.val % mod; return *this; }
	mint & operator <<= (int x) {
		base_type cur = 2;
		for (; x; x >>= 1) {
			if (x & 1) val = (mult_type) val * cur % mod;
			cur = (mult_type) cur * cur % mod;
		}
		return *this;
	}
	mint  operator + (const mint &rhs) const { return mint(*this) += rhs; }
	mint  operator - (const mint &rhs) const { return mint(*this) -= rhs; }
	mint  operator * (const mint &rhs) const { return mint(*this) *= rhs; }
	mint  operator << (int x) const { return mint(*this) <<= x; }
	friend std::ostream & operator << (std::ostream &stream, const mint &rhs) { stream << rhs.val; return stream; }
	explicit operator bool () { return val; }
};

u64 mod_inv_(u64 x, u64 mod) {
	s64 a = x, b = mod, u = 1, v = 0, t;
	while (b > 0) {
		t = a / b;
		a -= t * b;
		std::swap(a, b);
		u -= t * v;
		std::swap(u, v);
	}
	if (u < 0) u += mod;
	assert(0 <= u && (u64) u < mod);
	assert((__uint128_t) u * x % mod == 1);
	return u;
}
template<typename base_type, typename mult_type> cpp_int garner(std::vector<base_type> a, std::vector<base_type> mod) {
	int m = a.size();
	cpp_int fact = 1;
	std::vector<base_type> fact_mod(m, 1);
	cpp_int res = 0;
	std::vector<base_type> cur_remainder(m);
	for (int i = 0; i < m; i++) {
		base_type add_times = a[i] < cur_remainder[i] ? a[i] - cur_remainder[i] + mod[i] : a[i] - cur_remainder[i];
		add_times = (mult_type) add_times * mod_inv_(fact_mod[i], mod[i]) % mod[i];
		res += add_times * fact;
		fact *= mod[i];
		for (int j = i + 1; j < m; j++) {
			cur_remainder[j] = (cur_remainder[j] + (mult_type) fact_mod[j] * add_times) % mod[j];
			fact_mod[j] = (mult_type) fact_mod[j] * mod[i] % mod[j];
		}
	}
	for (int i = 0; i < m; i++) assert(res % mod[i] == a[i]);
	return res;
}