#include <cstdio>
#include <iostream>
#include <vector>
#include <map>
#include <array>
#include "debug.h"
#include "types.h"
#include "comb.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}


/*
	Number of permutations p of [0, ..., n-1] such that
	 - |p(i+k) - p(i)| != k
	
	O(n^{2m})
*/
template<typename val_t> val_t solve(int n, int m) {
	// dp[key][flag]
	// key: list of {number count, violation count} in a segment; descending except that key[0] must be the segment
	//  that has the most recently inserted number.
	// flag: whether i-1 & i-2 are adjacent
	std::map<std::vector<std::pair<int, int> >, std::array<val_t, 2> > dp;
	std::vector<std::pair<int, int> > init_key(m);
	init_key[0].first = 1;
	dp[init_key][0] = m;
	
	std::vector<int> disable_violation_indices;
	{
		int cur = 0;
		for (int i = 0; i < m; i++) disable_violation_indices.push_back(cur += (n + m - (i + 1)) / m);
	}
	int max = (n + m - 1) / m;
	for (int i = 1; i < n; i++) {
		bool disable_violation = std::find(disable_violation_indices.begin(), disable_violation_indices.end(), i) != disable_violation_indices.end();
		
		std::map<std::vector<std::pair<int, int> >, std::array<val_t, 2> > next;
		for (auto &j : dp) {
			for (int k = 0; k < 2; k++) {
				auto &cur_val = j.second[k];
				if (cur_val == 0) continue;
				
				for (int l = 0; l < m; l++) {
					if (j.first[l].first == max) continue;
					auto go = [&] (int violation_change, int next_flag, int coef) {
						auto next_key = j.first;
						std::swap(next_key[0], next_key[l]); // next_key[l] now has the most recently inserted integer
						next_key[0].first++;
						next_key[0].second += violation_change;
						if (max - next_key[0].first < next_key[0].second) return; // no room for resolving all violations 
						if (l) { // maintain sorted
							for (int t = l; t + 1 < m && next_key[t] > next_key[t + 1]; t++) std::swap(next_key[t], next_key[t + 1]);
							for (int t = l; t > 1 && next_key[t - 1] > next_key[t]; t--) std::swap(next_key[t - 1], next_key[t]);
						}
						next[next_key][next_flag] += cur_val * coef;
					};
					int cnt = j.first[l].first;
					int violation = j.first[l].second;
					go(0, 0, cnt + 1 - violation - (l == 0 ? k ? 1 : 2 : 0)); // insert into non-violating gap
					if (violation) go(-1, 0, violation - (l == 0 && k)); // insert into violating gap
					if (l == 0) { // insert next to i-1
						go(!disable_violation, !disable_violation, 2 - k); // insert into non-violating gap
						if (k) go(!disable_violation - 1, !disable_violation, k); // insert into violating gap (between i-1 & i-2)
					}
				}
			}
		}
		std::swap(dp, next);
	}
	
	std::vector<std::pair<int, int> > res_key(m);
	for (int i = 0; i < m; i++) res_key[i].first = (n + i) / m;
	auto res = dp[res_key][0];
	if (n % m) {
		std::swap(res_key[0], res_key[m - n % m]);
		res += dp[res_key][0];
	}
	combination_precalc_t<val_t> comb(n);
	res /= comb.nCr(m, n % m);
	return res;
}

/*
int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = 189281;
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	opt.n_min = 0;
	opt.n_max = ri();
	opt.thread_num = 2;
	
	gen.solve<A189281>(opt);
	
	if (1) {
		std::cerr << "Checking for recurrence :";
		for (int i = 1; i + 8 <= opt.n_max; i++) {
			std::cerr << " " << i;
			check_recurrence(gen.result, i);
		}
		std::cerr << std::endl;
		std::cerr << "OK" << std::endl;
	}
	
	return 0;
}
*/

int main() {
	int n = ri();
	int k = ri();
	Timer::measure([&] () {
		std::cerr << solve<cpp_int>(n, k) << std::endl;
	});
	return 0;
}
