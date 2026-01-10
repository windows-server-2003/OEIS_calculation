#define _GLIBCXX_DEBUG
#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"
#include "debug.h"
#include "types.h"
#include "comb.h"
#include "oeis.h"
#include "avoid_001_110_row_col_nwse.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

template<typename T> void apply_mul(std::vector<T> &lhs, std::vector<T> rhs, int limit) {
	int m = rhs.size();
	assert(m);
	int n = lhs.size();
	lhs.resize(std::min(limit, n + m - 1));
	for (int i = std::min(limit, n); i--; ) {
		for (int j = 1; j < m && i + j < limit; j++) lhs[i + j] += lhs[i] * rhs[j];
		lhs[i] *= rhs[0];
	}
}

struct A244288 {
	template<typename T> std::vector<T> calc_all(size_t n) {
		combination_precalc_t<T> comb(n + 1);
		
		std::vector<T> poly{1};
		std::vector<T> res{ 0 };
		for (size_t i = 1; i <= n; i++) {
			std::vector<T> cur_poly((i + 1) / 2 + 1);
			for (size_t j = 0; j <= (i + 1) / 2; j++) cur_poly[j] = comb.nCr(i + 1 - j, j);
			for (int j = 0; j < 2; j++) {
				apply_mul(poly, cur_poly, n + 1);
				if (j == 0) res.push_back(poly[i]);
			}
		}
		return res;
	}
};
struct A182563 {
	template<typename T> T calc(size_t n) {
		if (n == 1) return 1;
		combination_precalc_t<T> comb(n + 1);
		
		std::map<int, int> orders;
		for (size_t i = 0; i < n - 2; i++) orders[1 + i / 2]++;
		for (size_t i = 0; i < n; i++) {
			orders[std::min(1 + (n - 2) / 2, n - i)]++;
			orders[std::min(1 + (n - 1) / 2, n - i)]++;
		}
		
		std::vector<T> res{1};
		for (auto &i : orders) {
			size_t m = i.first;
			size_t k = i.second;
			std::vector<T> poly(m + 1);
			for (size_t i = 0; i <= (m + 1) / 2; i++) poly[i] = comb.nCr(m + 1 - i, i);
			
			for (; k; k >>= 1) {
				if (k & 1) apply_mul(res, poly, n + 1);
				if (k >= 2) apply_mul(poly, poly, n + 1);
			}
		}
		return res[n];
	}
};
struct A197989 {
	template<typename T> T calc(size_t n) {
		combination_precalc_t<T> comb(n + 1);
		
		std::vector<T> poly(n + 1);
		for (size_t i = 0; i <= (n + 1) / 2; i++) poly[i] = comb.nCr(n + 1 - i, i);
		
		std::vector<T> res{1};
		for (size_t t = n; t; t >>= 1) {
			if (t & 1) apply_mul(res, poly, n + 1);
			apply_mul(poly, poly, n + 1);
		}
		return res[n];
	}
};

int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = gen.get_sequence_number(argc, argv);
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	if (seq_num == -1) {
		std::cerr << "No valid sequence number found" << std::endl;
		return 1;
	}
	
	opt.n_min = 1;
	opt.n_max = ri();
	opt.thread_num = 4;
	// opt.thread_limits = {{30, 2}, {35, 1}};
	
	switch (seq_num) {
		case 244288 : gen.solve<A244288>(opt); break;
		case 182563 : gen.solve<A182563>(opt); break;
		case 197989 : gen.solve<A197989>(opt); break;
		case 202447 : gen.solve<A202447>(opt); break;
		case 202439 : gen.solve<A202439>(opt); break;
		default :
			std::cerr << "Invalid seq number : " << seq_num << std::endl;
	}
	
	return 0;
}

