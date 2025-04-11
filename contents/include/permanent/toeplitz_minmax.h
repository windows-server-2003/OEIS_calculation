#include <vector>
#include <numeric>
#include <map>
#include <cmath>
#include "debug.h"
#include "comb.h"

template<typename val_t> val_t A358327(int n) {
	std::vector<int> a(n);
	std::iota(a.begin(), a.end(), 0);
	val_t res = 0;
	do {
		std::vector<std::vector<val_t> > m(n, std::vector<val_t>(n));
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) m[i][j] = a[n - 1 - abs(i - j)];
		val_t cur_perm = 0;
		for (int i = 0; i < 1 << n; i++) {
			val_t cur_prod = 1;
			for (int j = 0; j < n; j++) {
				val_t cur_sum = 0;
				for (int k = 0; k < n; k++) if (i >> k & 1) cur_sum += m[j][k];
				cur_prod *= cur_sum;
			}
			if ((n - __builtin_popcount(i)) & 1) cur_perm -= cur_prod;
			else cur_perm += cur_prod;
		}
		res = std::max(res, cur_perm);
	} while (std::next_permutation(a.begin(), a.end()));
	
	return res;
}
template<typename val_t> val_t A358327_2(int n) {
	int h = (n + 1) / 2;
	std::vector<int> a(n);
	std::iota(a.begin(), a.begin() + h, 0);
	std::iota(a.begin() + h, a.end() - 1, h + 1);
	a.back() = h;
	// std::cerr << a << std::endl;
	
	auto get_perm = [&] () {
		std::vector<std::vector<val_t> > m(n, std::vector<val_t>(n));
		for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) m[i][j] = a[n - 1 - abs(i - j)];
		val_t cur_perm = 0;
		for (int i = 0; i < 1 << n; i++) {
			val_t cur_prod = 1;
			for (int j = 0; j < n; j++) {
				val_t cur_sum = 0;
				for (int k = 0; k < n; k++) if (i >> k & 1) cur_sum += m[j][k];
				cur_prod *= cur_sum;
			}
			if ((n - __builtin_popcount(i)) & 1) cur_perm -= cur_prod;
			else cur_perm += cur_prod;
		}
		return cur_perm;
	};
	
	return get_perm();
}

template<typename val_t> struct A358327_3 {
	int n;
	val_t global_lower_bound;
	std::vector<std::vector<val_t> > power;
	std::vector<int> cnt;
	val_t dfs(int i, const std::map<std::vector<int>, val_t> &coefs) {
		val_t res = 0;
		
		
		// if (i < 4) std::cerr << i << " " << coefs.size() << std::endl;
		if (i == n - 1) {
			for (auto &j : coefs) res += power[i][j.first[0]] * j.second;
			return res;
		}
		if (0) {
			val_t upper_bound = 0;
			for (auto &j : coefs) upper_bound += power[n - 1][std::accumulate(j.first.begin(), j.first.end(), 0)] * j.second;
			if (upper_bound <= global_lower_bound) {
				// std::cerr << "ng" << std::endl;
				return 0;
			}//  else std::cerr << "ok" << std::endl;
		}
		if (1) {
			val_t upper_bound = 0;
			for (auto j : coefs) {
				auto powers = j.first;
				std::sort(powers.begin(), powers.end(), std::greater<>());
				val_t cur_prod = 1;
				for (int k = 0; k < (int) powers.size(); k++) cur_prod *= power[n - 1 - k][powers[k]];
				upper_bound += cur_prod * j.second;
			}
			if (upper_bound <= global_lower_bound) return 0;
		}
		cnt[i]++;
		
		int max = n - i;
		for (int j = 0; j < max; j++) {
			std::map<std::vector<int>, val_t> next_coefs;
			for (auto &k : coefs) {
				auto next_key = k.first;
				auto times = power[i][next_key[j]];
				next_key.erase(next_key.begin() + j);
				next_coefs[next_key] += k.second * times;
			}
			res = std::max(res, dfs(i + 1, next_coefs));
		}
		return res;
	}
	val_t solve(int n) {
		this->n = n;
		this->global_lower_bound = A358327_2<val_t>(n) * 9 / 10;
		this->power.assign(n, {});
		for (int i = 0; i < n; i++) {
			this->power[i].assign(n + 1, 1);
			for (int j = 1; j <= n; j++) this->power[i][j] = this->power[i][j - 1] * i;
		}
		this->cnt.resize(n);
		
		std::map<std::vector<int>, val_t> all;
		for_each_permutation(n, [&] (size_t, const std::vector<int> &a) {
			std::vector<int> cnt(n);
			for (int i = 0; i < n; i++) cnt[n - 1 - std::abs(i - a[i])]++;
			all[cnt]++;
			return true;
		});
		std::cerr << "size : " << all.size() << std::endl;
		
		auto res = dfs(0, all);
		
		std::cerr << cnt << std::endl;
		
		return res;
	}
};
