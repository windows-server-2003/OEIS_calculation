#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include "utils.h"
#include "debug.h"
#include "types.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

/*
	A051424(n) = number of ways to partition n into n = a_1 + a_2 + ... + a_k such that a_i are coprime.
	
	We call primes smaller than sqrt(n) *small*, and ones larger than sqrt(n) *large*.
	We use dynamic programming. Looking from larger primes, we decide
	 - whether we use the prime p in any of a_i
	 - if so, what is a_i
	The DP state is (sum of currently decided a_i's, whether we used each of *small* primes).
	We only have to remember the usage of small primes, because
	 - when choosing whether to use a large prime p, it is not possible that another large prime q less than p has been used.
		because it would mean another prime r > q is used together with q but rq exceeds n.
	 - when choosing whether to use a small prime p, usage of any primes larger than p is no longer relevant.
*/
struct A051424 {
	int n;
	// divisor[i]: a prime factor of i; 0 if i is prime.
	// Prime factoring is easily done with this array.
	std::vector<int> divisor;
	std::vector<int> primes; // list of primes less than n
	int sq; // number of small primes
	
	void init_primes(int max) {
		divisor.assign(max + 1, 0);
		primes.clear();
		sq = 0;
		
		for (int i = 2; i <= max; i++) if (!divisor[i]) {
			primes.push_back(i);
			for (int j = i + i; j <= max; j += i) divisor[j] = i;
		}
		for (; sq + 1 < (int) primes.size() && primes[sq] * primes[sq + 1] <= n; sq++);
	}
	int prime_index(int p) { return std::lower_bound(primes.begin(), primes.end(), p) - primes.begin(); }
	
	// We manage the usage of small primes by u64 (where 2^i bit refers to the usage of primes[i])
	// get_factors(x): if x consists of small primes only, return its usage of small primes
	// if x has large prime factor, return (u64) -1
	u64 get_factors(int x) {
		if (x == 1) return 0;
		u64 res = 0;
		while (divisor[x]) {
			int index = prime_index(divisor[x]);
			if (index <= sq) res |= 1ULL << index; // primes[sq] * primes[sq] may be <= n
			else return (u64) -1;
			x /= divisor[x];
		}
		int index = prime_index(x);
		if (index <= sq) res |= 1ULL << index;
		else return (u64) -1;
		return res;
	}
	/*
		Replacement of map<size_t, T> but managers smaller keys with vector<T> to save map size and memory
	*/
	template<typename T, size_t ARRAY_SIZE_LIMIT> struct PartialMapArray {
		size_t size_;
		std::vector<T> array;
		std::map<size_t, T> map;
		PartialMapArray (size_t size) : size_(size), array(std::min(ARRAY_SIZE_LIMIT, size)) {}
		T & operator [] (size_t i) { return i < ARRAY_SIZE_LIMIT ? array[i] : map[i]; }
		template<typename F> void for_each_rev(size_t r, const F &f) {
			if (map.lower_bound(r) != map.begin()) {
				for (auto i = std::prev(map.lower_bound(r)); ; i--) {
					f(i->first, i->second);
					if (i == map.begin()) break;
				}
			}
			for (size_t i = std::min(r, array.size()); i--; ) f(i, array[i]);
		}
		template<typename F> void for_each_rev(const F &f) { for_each_rev((size_t) -1, f); }
		template<typename F> void for_each(size_t l, const F &f) {
			for (size_t i = l; i < array.size(); i++) f(i, array[i]);
			for (auto i = map.lower_bound(l); i != map.end(); i++) f(i->first, i->second);
		}
		void shrink(size_t size) {
			this->size_ = size;
			map.erase(map.lower_bound(size), map.end());
			if (array.size() > size) array.resize(size), array.shrink_to_fit();
		}
		void expand(size_t size) {
			this->size_ = size;
			array.resize(std::min(ARRAY_SIZE_LIMIT, size));
		}
		size_t size() { return size_; }
		size_t mem() { return sizeof(array) + array.size() * sizeof(T) + sizeof(map) + map.size() * (16 + sizeof(size_t) + sizeof(T)); }
	};
	template<typename val_t> std::vector<val_t> calc_all_(int n) {
		this->n = n;
		init_primes(n);
		
		assert(sq <= 60);
		
		// dp[i][j]: the number of ways to get small primes usage i and decided sum j.
		PartialMapArray<std::vector<val_t>, 10000000> dp(1);
		dp[0] = std::vector<val_t>(n + 1);
		dp[0][0] = 1;
		
		std::vector<u64> factors(n + 1);
		for (int i = 1; i <= n; i++) factors[i] = get_factors(i);
		
		std::cerr << sq << " " << (double) (1ULL << sq) * n * 16 / 1000000000 << " GB" << std::endl;
		
		for (int i = primes.size(); i--; ) {
			// after deciding on primes[i], DP only has to remember the usage of first `next_bits` primes
			//   primes[next_bits] * primes[i] > n   -> primes[next_bits] can never have been used (after deciding on primes[i])
			int next_bits = 0;
			while (next_bits < (int) primes.size() && next_bits < i && primes[next_bits] * primes[i] <= n) next_bits++;
			size_t next_size = 1ULL << next_bits;
			
			u64 forbidden_mask = 0;
			if (i <= 62) forbidden_mask = ~((1ULL << (i + 1)) - 1); // cannot use primes[>= i + 1]
			auto go = [&] (int i, int j, val_t &val) {
				if (!dp[i].size()) dp[i].assign(n + 1, 0);
				dp[i][j] += val;
			};
			if (i >= sq) { // dp table size(next_size) will grow or stay same
				assert(next_size >= dp.size());
				dp.expand(next_size);
				dp.for_each_rev([&] (size_t j, std::vector<val_t> &dp_j) {
					// j: old small primes usage;  k: old sum;  l*primes[i]: a_i that contains primes[i]
					if (dp_j.size()) for (int k = n; k >= 0; k--) if (dp_j[k]) for (int l = 1; k + l * primes[i] <= n; l++)
						if (!(factors[l] & forbidden_mask) && !(j & factors[l])) go((j | factors[l]) & (next_size - 1), k + l * primes[i], dp_j[k]);
				});
			} else { // dp table shrinks to half
				assert(next_size == dp.size() / 2);
				assert(next_size == (1ULL << i));
				dp.for_each_rev(next_size, [&] (size_t j, std::vector<val_t> &dp_j) {
					if (dp_j.size()) for (int k = n; k >= 0; k--) if (dp_j[k]) for (int l = 1; k + l * primes[i] <= n; l++)
						if (!(factors[l] & forbidden_mask) && !(j & factors[l])) go((j | factors[l]) & (next_size - 1), k + l * primes[i], dp_j[k]);
				});
				// convolve second half into first half(=forget the usage of prime[i])
				dp.for_each(next_size, [&] (size_t j, std::vector<val_t> &dp_j) {
					if (dp_j.size()) for (int k = 0; k <= n; k++) if (dp_j[k]) go(j - next_size, k, dp_j[k]);
				});
				dp.shrink(next_size);
			}
			size_t dp_size = 0;
			dp.for_each(0, [&] (size_t, std::vector<val_t> &i) { dp_size += i.size(); });
			if (next_bits >= 11 || i <= sq) {
				fprintf(stderr, "[MEMINFO] Processing prime #%d: std::vector mem = %10sB,  content mem = %10sB\n",
					i, put_SI_prefix(dp.mem()).c_str(), put_SI_prefix(dp_size * sizeof(val_t)).c_str());
			}
		}
		assert(dp.size() == 1);
		// any number of 1's
		for (int i = 0; i < n; i++) dp[0][i + 1] += dp[0][i];
		return dp[0];
	}
	// slightly faster version to calculate the maximum RAM usage in advance; not used unless manually coded to use it
	template<typename val_t> void calc_mem(int n) {
		this->n = n;
		init_primes(n);
		
		assert(sq <= 60);
		
		std::map<size_t, int> dp = { {0, 0} };
		
		std::vector<u64> factors(n + 1);
		for (int i = 1; i <= n; i++) factors[i] = get_factors(i);
		
		std::cerr << sq << " " << (double) (1ULL << sq) * n * 16 / 1000000000 << " GB" << std::endl;
		
		size_t content_max = 0;
		size_t dp_max = 0;
		size_t max = 0;
		
		for (int i = primes.size(); i--; ) {
			int next_bits = 0;
			while (next_bits < (int) primes.size() && next_bits < i && primes[next_bits] * primes[i] <= n) next_bits++;
			size_t next_size = 1ULL << next_bits;
			
			u64 forbidden_mask = 0;
			if (i <= 62) forbidden_mask = ~((1ULL << (i + 1)) - 1); // cannot use primes[>= i + 1]
			auto go = [&] (size_t i, int j) {
				if (!dp.count(i)) dp[i] = j;
				else dp[i] = std::min(dp[i], j);
			};
			if (i >= sq) {
				assert(dp.size());
				assert(next_size >= dp.size());
				for (auto j = std::prev(dp.end()); ; j--) {
					for (int l = 1; j->second + l * primes[i] <= n; l++) if (!(factors[l] & forbidden_mask) && !(j->first & factors[l])) go((j->first | factors[l]) & (next_size - 1), j->second + l * primes[i]);
					if (j == dp.begin()) break;
				}
			} else {
				assert(dp.size());
				assert(next_size == (1ULL << i));
				for (auto j = std::prev(dp.lower_bound(next_size)); ; j--) {
					for (int l = 1; j->second + l * primes[i] <= n; l++) if (!(factors[l] & forbidden_mask) && !(j->first & factors[l])) go((j->first | factors[l]) & (next_size - 1), j->second + l * primes[i]);
					if (j == dp.begin()) break;
				}
				for (auto j = dp.lower_bound(next_size); j != dp.end(); j++) go(j->first - next_size, j->second);
				dp.erase(dp.lower_bound(next_size), dp.end());
			}
			size_t content_size = dp.size() * (n + 1);
			content_size *= sizeof(val_t);
			size_t vec_size = (std::min<size_t>(10000000, next_size) + dp.size()) * 24;
			
			if (next_bits >= 12) 
			std::cerr << "#" << i << " " << (long double) content_size / 1000000000 << " GB    "
				<< (long double) vec_size / 1000000000 << " GB" << std::endl;
				
			content_max = std::max(content_max, content_size);
			dp_max = std::max(dp_max, vec_size);
			max = std::max(max, content_size + vec_size);
		}
		std::cerr << (long double) content_max / 1000000000 << " GB" << std::endl;
		std::cerr << (long double) dp_max / 1000000000 << " GB" << std::endl;
		std::cerr << (long double) max / 1000000000 << " GB" << std::endl;
	}
	template<typename val_t> std::vector<val_t> calc_all(int n) {
		// val_t is not used internally
		auto res = calc_all_<u256>(n);
		return std::vector<val_t>(res.begin(), res.end());
	}
};


int main(int argc, char **argv) {
	std::cerr << "Example: Max n = 5000 typically takes 10 seconds and around 1 GB of RAM" << std::endl;
	std::cerr << "Example: Max n = 10000 typically takes 1 minute and around 10 GB of RAM" << std::endl;
	std::cerr << "Example: Max n = 20000 typically takes 20 minutes and around 200 GB of RAM" << std::endl;
	std::cerr << "Max n > ";
	int n_max = ri();
	
	OEISGenerator gen;
	int seq_num = 51424;
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	opt.n_min = 0;
	opt.n_max = n_max;
	opt.thread_num = 1;
	
	gen.solve<A051424>(opt);
	
	return 0;
}
