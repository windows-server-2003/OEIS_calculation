#include <cstdio>
#include <iostream>
#include <vector>
#include "debug.h"
#include "types.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

// Number of partitions of n such that every orderless pair of distinct parts has a different sum
struct A325857 {
	using val_t = u64;
	
	int n;
	u256 sum_mask; // usage of sums (i...i+255)
	u256 used; // bit representation of usage of integers
	int sum; // sum of the integers that the dfs has so far decided to use; some of them may be used more than once but it counts towards the index of `dp`
	std::vector<val_t> res;
	// dp[x]: number of ways to choose additional integers up to current dfs depth such that sum of the additional numbers are x
	// note that the dfs determines the set of integers used, but their multiplicities are free so if the same number is used more than once, it counts toward x
	std::vector<val_t> dp;
	void dfs(int i) {
		for (int j = 0; j <= n - sum; j++) res[j + sum] += dp[j];
		
		auto sum_mask_bak = sum_mask;
		for (int j = i + 1; j <= n - sum; j++) {
			sum_mask >>= 1;
			if (!(sum_mask & used)) {
				sum += j; // base sum: j must be used at least once
				if (n - sum < j + 1) { // (to speed up) there is no room for yet another number to be added (except if n - sum == j, j can be used once more)
					for (int k = 0; k <= n - sum; k++) res[k + sum] += dp[k];
					if (n - sum == j) res[n] += dp[0];
				} else {
					sum_mask ^= used;
					used ^= (u256) 1 << j;
					for (int k = 0; k <= n - sum - j; k++) dp[k + j] += dp[k]; // j can be additionally used for any number of times
					dfs(j);
					for (int k = n - sum - j; k >= 0; k--) dp[k + j] -= dp[k]; // revert the above process
					used ^= (u256) 1 << j;
					sum_mask ^= used;
				}
				sum -= j;
			}
		}
		sum_mask = sum_mask_bak;
	}
	// this val_t_ is not actually used
	template<typename val_t_> std::vector<val_t_> calc_all(int n) {
		assert(n < 512);
		
		this->n = n;
		this->sum_mask = this->used = 0;
		this->sum = 0;
		this->res.assign(n + 1, 0);
		this->dp.assign(n + 1, 0);
		dp[0] = 1;
		
		dfs(0);
		return std::vector<val_t_>(res.begin(), res.end());
	}
};
// Number of strict partitions of n such that every orderless pair of distinct parts has a different sum
struct A325877 {
	using val_t = u64;
	
	int n;
	u256 sum_mask;
	u256 used;
	int remaining;
	std::vector<val_t> res;
	void dfs(int i) {
		res[n - remaining]++;
		
		auto sum_mask_bak = sum_mask;
		for (int j = i + 1; j <= remaining; j++) {
			sum_mask >>= 1;
			if (!(sum_mask & used)) {
				remaining -= j;
				if (remaining < j + 1) res[n - remaining]++;
				else {
					sum_mask ^= used;
					used ^= (u256) 1 << j;
					dfs(j);
					used ^= (u256) 1 << j;
					sum_mask ^= used;
				}
				remaining += j;
			}
		}
		sum_mask = sum_mask_bak;
	}
	
	// this val_t_ is not actually used
	template<typename val_t_> std::vector<val_t_> calc_all(int n) {
		assert(n < 512);
		
		this->n = n;
		this->sum_mask = this->used = 0;
		this->remaining = n;
		this->res.assign(n + 1, 0);
		
		dfs(0);
		return std::vector<val_t_> (res.begin(), res.end());
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
	opt.thread_num = 1;
	
	switch (seq_num) {
		case 325857 : gen.solve<A325857>(opt); break;
		case 325877 : gen.solve<A325877>(opt); break;
	}
	
	return 0;
}

