#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <bitset>
#include "utils.h"
#include "debug.h"
#include "types.h"
#include "prime.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

/*
	
*/

template<typename val_t> struct A339556 {
	std::vector<bool> is_prime;
	std::vector<int> primes;
	
	void init_primes(size_t num_primes) {
		constexpr int n = 20000;
		is_prime.assign(n + 1, true);
		primes.clear();
		for (int i = 2; i <= n; i++) if (is_prime[i]) {
			for (int j = i + i; j <= n; j += i) is_prime[j] = false;
			primes.push_back(i);
			if (primes.size() >= num_primes) {
				is_prime.resize(i + 1);
				return;
			}
		}
		assert(0);
	}
	int n_max;
	std::vector<val_t> res;
	
	A339556 (int n_min, int n_max) : n_max(n_max), res(n_max + 1, 2) {
		assert(n_min == 1);
		assert(n_max >= 2);
		init_primes(n_max);
		primes.erase(primes.begin()); // subsets contatining 2 and satisfying the condition is only { 2 }
		primes.erase(primes.begin()); // subsets contatining 3 and satisfying the condition is only { 3 }
	}
	
	std::vector<int> get_tasks(int, int) { return primes; }
	void calc_custom(int p, std::mutex &lock) { // more than one thread didn't make it faster, probably because of memory bandwidths
		size_t plus_sum = 0;
		for (auto q : primes) if (q > p) plus_sum += (q * q - p * p) / 24;
		
		std::vector<val_t> dp = { 1 };
		std::vector<val_t> cur_res(n_max + 1);
		for (size_t i = 0; i < primes.size(); i++) {
			int q = primes[i];
			int diff = (q * q - p * p) / 24;
			if (diff < 0) {
				diff *= -1;
				dp.resize(std::min(dp.size() + diff, plus_sum + 1));
				if ((int) dp.size() > diff) for (size_t j = dp.size() - diff; j--; ) dp[j + diff] += dp[j];
			} else {
				if (diff > 0) for (size_t j = 0; j + diff < dp.size(); j++) dp[j] += dp[j + diff];
				cur_res[3 + i] = dp[0] * 2 - 1;
			}
		}
		{
			std::lock_guard<std::mutex> locking(lock);
			for (int i = 0; i <= n_max; i++) res[i] += cur_res[i];
		}
	}
	std::vector<val_t> get_result(int, int) {
		res[0] = 0;
		res[1] = 1;
		res[2] = 2;
		return res;
	}
};

// A108018

template<typename val_t> struct A339485 {
	std::vector<bool> is_prime;
	std::vector<int> primes;
	
	void init_primes(size_t num_primes) {
		constexpr int n = 40000;
		is_prime.assign(n + 1, true);
		primes.clear();
		for (int i = 2; i <= n; i++) if (is_prime[i]) {
			for (int j = i + i; j <= n; j += i) is_prime[j] = false;
			primes.push_back(i);
			if (primes.size() >= num_primes) {
				is_prime.resize(i + 1);
				return;
			}
		}
		assert(0);
	}
	int n_max;
	std::vector<val_t> res;
	
	A339485 (int n_min, int n_max) : n_max(n_max), res(n_max + 1, 1) {
		assert(n_min == 1);
		init_primes(n_max);
		primes.erase(primes.begin()); // subsets contatining 2 and satisfying the condition is only { 2 }
	}
	
	std::vector<int> get_tasks(int, int) { return primes; }
	
	void calc_custom(int p, std::mutex &lock) {
		size_t plus_sum = 0;
		for (auto q : primes) if (q > p) plus_sum += (q - p) / 2;
		
		std::vector<val_t> dp = { 1 };
		std::vector<val_t> cur_res(n_max + 1);
		for (size_t i = 0; i < primes.size(); i++) {
			int q = primes[i];
			int diff = (q - p) / 2;
			if (diff < 0) {
				diff *= -1;
				dp.resize(std::min(dp.size() + diff, plus_sum + 1));
				if ((int) dp.size() > diff) for (size_t j = dp.size() - diff; j--; ) dp[j + diff] += dp[j];
			} else {
				if (diff > 0) for (size_t j = 0; j + diff < dp.size(); j++) dp[j] += dp[j + diff];
				cur_res[2 + i] += dp[0] * 2 - 1;
			}
		}
		{
			std::cerr << p << std::endl;
			std::lock_guard<std::mutex> locking(lock);
			for (int i = 0; i <= n_max; i++) res[i] += cur_res[i];
		}
	}
	std::vector<val_t> get_result(int, int) {
		res[0] = 0;
		res[1] = 1;
		return res;
	}
};

struct A115030 {
	
	template<typename val_t> std::vector<val_t> calc_all(int n) {
		constexpr size_t MAX = 10000000;
		auto primes = get_first_n_primes<u64>(n);
		assert(std::accumulate(primes.begin(), primes.end(), (u64) 0) < MAX);
		std::bitset<MAX> bitset;
		bitset[0] = 1;
		std::vector<val_t> res = { 0 };
		size_t sum = 0;
		for (auto p : primes) {
			bitset |= bitset << p;
			sum += p;
			
			for (s64 i = p/2; i >= 0; i--) if (!bitset[i]) {
				std::cerr << p << " : sum=" << sum << " -> " << p/2 - i << " " << i << std::endl;
				break;
			}
			
			res.push_back(bitset.count() - 1);
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
	
	opt.n_min = 1;
	opt.n_max = ri();
	opt.thread_num = 1;
	
	switch (seq_num) {
		case 339556 : gen.solve<A339556<u512> >(opt); break;
		case 339485 : gen.solve<A339485<cpp_int> >(opt); break;
		case 115030 : gen.solve<A115030>(opt); break;
		default : std::cerr << "Unknown sequence : A" << seq_num << " found" << std::endl;
	}
	
	
	return 0;
}
