#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "utils.h"
#include "types.h"
#include "comb.h"
#include "debug.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

struct A137432 {
	template<typename val_t> std::vector<val_t> calc_all(int n) {
		combination_precalc_t<val_t> comb(2 * n);
		
		int max_lines = n / 2;
		
		// indices : whether we have already passed the start point, b_i, # of lines crossing per direction, # of vertices used
		std::vector<std::vector<std::vector<val_t> > > dp(2, std::vector<std::vector<val_t> >(max_lines + 1, std::vector<val_t>(n)));
		auto dp_initialized = dp;
		dp[0][0][0] = 1;
		
		std::vector<val_t> res = { 1 };
		for (int i = 0; i <= n; i++) {
			// allow closing
			for (int after_start = 0; after_start < 2; after_start++) {
				for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[after_start][k][l]) {
					int closable = after_start ? k : k - 1;
					for (int m = 1; m <= closable && l + m < n; m++) {
						dp[after_start][k - m][l + m] += dp[after_start][k][l] * comb.nCr(closable, m);
					}
				}
			}
			
			decltype(dp) dp_next = dp_initialized;
			{
				// allow knotting
				{ // when continuing normally
					for (int after_start = 0; after_start < 2; after_start++) {
						for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[after_start][k][l]) {
							for (int knot = 0; knot <= std::min({n - 1 - l, 2 * k, (n - 1 - l) - k + 1}); knot++) 
								dp_next[after_start][k][l + knot] += dp[after_start][k][l] * comb.nCr(2 * k, knot);
						}
					}
				}
				{ // when passing the starting point 
					for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[0][k][l]) {
						// (top|bottom)_dir : 0:left, 1:right
						for (int top_dir = 0; top_dir < 2; top_dir++) for (int bottom_dir = 0; bottom_dir < 2; bottom_dir++) {
							if (!k && (!top_dir || !bottom_dir)) continue;
							int left_knottable = k - !top_dir;
							int right_knottable = k - !bottom_dir;
							int k_diff = (top_dir + bottom_dir) - 1;
							if (k + k_diff > max_lines) continue;
							assert(k + k_diff >= 0);
							for (int knot = 0; knot <= std::min({left_knottable + right_knottable, n - 1 - l, (n - 1 - l) - (k + k_diff) + 1}); knot++) 
								dp_next[1][k + k_diff][l + knot] += dp[0][k][l] * comb.nCr(left_knottable + right_knottable, knot);
						}
					}
				}
			}
			if (i) {
				val_t cur_res = 0;
				int cur_n = i;
				for (int j = 0; j < cur_n; j++) cur_res += dp_next[1][0][cur_n - 1 - j] * comb.nHr(cur_n + 1 - j, j);
				cur_res += ((val_t) 1 << cur_n) * (cur_n + 1); // [i, i, ..., i] (0 <= i <= n)
				std::cerr << cur_n << " " << cur_res << std::endl;
				res.push_back(cur_res);
			}
			{
				// allow knotting
				{ // when continuing normally
					for (int after_start = 0; after_start < 2; after_start++) {
						for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[after_start][k][l]) {
							for (int knot = k; knot <= std::min({n - 1 - l, 2 * k, (n - 1 - l) - k + 1}); knot++) 
								dp_next[after_start][k][l + knot] += dp[after_start][k][l] * comb.nCr(k, knot - k);
						}
					}
				}
				{ // when passing the starting point 
					for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[0][k][l]) {
						// (top|bottom)_dir : 0:left, 1:right
						for (int top_dir = 0; top_dir < 2; top_dir++) for (int bottom_dir = 0; bottom_dir < 2; bottom_dir++) {
							if (!k && (!top_dir || !bottom_dir)) continue;
							int left_knottable = k - !top_dir;
							int right_knottable = k - !bottom_dir;
							int k_diff = (top_dir + bottom_dir) - 1;
							if (k + k_diff > max_lines) continue;
							assert(k + k_diff >= 0);
							for (int knot = left_knottable; knot <= std::min({left_knottable + right_knottable, n - 1 - l, (n - 1 - l) - (k + k_diff) + 1}); knot++) 
								dp_next[1][k + k_diff][l + knot] += dp[0][k][l] * comb.nCr(right_knottable, knot - left_knottable);
						}
					}
				}
			}
			
			// allow opening
			for (int after_start = 0; after_start < 2; after_start++) {
				for (int k = max_lines; k >= 0; k--) for (int l = 0; l < n; l++) if (dp_next[after_start][k][l]) {
					int openable_place_num = after_start ? k : k + 1;
					if (openable_place_num) {
						for (int m = 1; m <= std::min({n - 1 - l, max_lines - k, ((n - 1 - l) - k + 1) / 2}); m++) {
							dp_next[after_start][k + m][l + m] += dp_next[after_start][k][l] * comb.nHr(openable_place_num, m);
						}
					}
				}
			}
			dp = dp_next;
		}
		return res;
	}
};


int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = 137432;
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	opt.n_min = 1;
	opt.n_max = ri();
	opt.thread_num = 1;
	
	gen.solve<A137432>(opt);
	
	return 0;
}

