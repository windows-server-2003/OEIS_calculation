#pragma once
#include <vector>

template<typename val_t> struct combination_precalc_t {
	int n;
	std::vector<std::vector<val_t> > ncr;
	std::vector<val_t> fact;
	combination_precalc_t (int n) : n(n), ncr(n + 1, std::vector<val_t>(n + 1)), fact(n + 1, 1) {
		ncr[0][0] = 1;
		for (int i = 1; i <= n; i++) {
			ncr[i][0] = ncr[i][i] = 1;
			for (int j = 1; j < i; j++) ncr[i][j] = ncr[i - 1][j] + ncr[i - 1][j - 1];
		}
		for (int i = 2; i <= n; i++) fact[i] = fact[i - 1] * i;
	}
	val_t nCr(int n, int r) {
		assert(0 <= r && r <= n && n <= this->n);
		return ncr[n][r];
	}
	val_t nHr(int n, int r) {
		if (!n && !r) return 1;
		return nCr(n + r - 1, r);
	}
};
template<typename func_t> bool for_each_permutation(int n, const func_t &func) {
	std::vector<int> a(n);
	std::iota(a.begin(), a.end(), 0);
	size_t index = 0;
	do {
		if (!func(index, a)) return false;
	} while (std::next_permutation(a.begin(), a.end()));
	return true;
}
