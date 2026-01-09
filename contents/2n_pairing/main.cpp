#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include "utils.h"
#include "debug.h"
#include "types.h"
#include "comb.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}
struct A320129 {
	template<typename val_t> val_t calc(int n) {
		std::vector<int> cnt(2 * (2 * n - 2) + 1, 0);
		std::vector<u64> dp;
		u64 res = 0;
		u64 cur_bits = 0;
		std::map<std::vector<int>, int> all;
		for (u64 i = 0; i < 1ULL << (2 * n - 1); i++) {
			if (i) {
				int t = __builtin_ctzll(i);
				if (cur_bits >> t & 1) {
					cur_bits ^= 1ULL << t;
					for (int j = 0; j < 2 * n; j++) if (cur_bits >> j & 1) cnt[t + j - 1]--;
					cnt[t + (2 * n - 1) - 1]--;
				} else {
					for (int j = 0; j < 2 * n; j++) if (cur_bits >> j & 1) cnt[t + j - 1]++;
					cnt[t + (2 * n - 1) - 1]++;
					cur_bits ^= 1ULL << t;
				}
			}
			auto cnt1 = cnt;
			std::sort(cnt1.begin(), cnt1.end());
			all[cnt1]++;
			dp.assign(2 * (2 * n - 2) + 2, 0);
			dp[0] = 1;
			
			int k_max = 0;
			for (auto j : cnt) if (j) {
				for (int k = k_max; k >= 0; k--) dp[k + 1] += dp[k] * j;
				k_max++;
			}
			dp[n] *= (cur_bits ? __builtin_ctzll(cur_bits) : 2 * n - 1) + 1;
			if (i & 1) res += dp[n];
			else res -= dp[n];
		}
		std::cerr << "total : " << (1ULL << (2 * n - 1)) << std::endl;
		std::cerr << "all   : " << all.size() << std::endl;
		
		
		return res;
	}
};
struct A060963 {
	template<typename val_t> val_t calc(int n) {
		std::vector<int> cnt(2 * n, 0);
		std::vector<u64> dp;
		u64 res = 0;
		u64 cur_bits = 0;
		for (u64 i = 0; i < 1ULL << (2 * n - 1); i++) {
			if (i) {
				int t = __builtin_ctzll(i);
				if (cur_bits >> t & 1) {
					cur_bits ^= 1ULL << t;
					for (int j = 0; j < t; j++) if (cur_bits >> j & 1) cnt[t - j]--;
					for (int j = t + 1; j < 2 * n - 1; j++) if (cur_bits >> j & 1) cnt[j - t]--;
					cnt[(2 * n - 1) - t]--;
				} else {
					for (int j = 0; j < t; j++) if (cur_bits >> j & 1) cnt[t - j]++;
					for (int j = t + 1; j < 2 * n - 1; j++) if (cur_bits >> j & 1) cnt[j - t]++;
					cnt[(2 * n - 1) - t]++;
					cur_bits ^= 1ULL << t;
				}
			}
			dp.assign(2 * (2 * n - 2) + 2, 0);
			dp[0] = 1;
			
			int k_max = 0;
			for (auto j : cnt) if (j) {
				for (int k = k_max; k >= 0; k--) dp[k + 1] += dp[k] * j;
				k_max++;
			}
			dp[n] *= (cur_bits ? __builtin_ctzll(cur_bits) : 2 * n - 1) + 1;
			if (i & 1) res += dp[n];
			else res -= dp[n];
		}
		
		return res;
	}
};

int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = gen.get_sequence_number(argc, argv);
	if (seq_num == -1) {
		std::cerr << "No valid sequence number found" << std::endl;
		return 1;
	}
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	opt.n_min = opt.n_max = ri();
	opt.thread_num = 1;
	
	switch (seq_num) {
		case 320129 : gen.solve<A320129>(opt); break;
		case  60963 : gen.solve<A060963>(opt); break;
	}
	
	
	return 0;
}
