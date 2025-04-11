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
		
		/*
			dp[after_start][k][l] :
				after_start: whether we have already passed the start point
				k: number of pairs of (left-going and right-going) horizontal lines crossing the current vertical line
				l: # of vertices already created
			It can be proved that the dp values are same regardless of b_i = 0 or b_i = 1, so the same, single value is stored in dp[after_start][k][l]
		*/
		std::vector<std::vector<std::vector<val_t> > > dp(2, std::vector<std::vector<val_t> >(max_lines + 1, std::vector<val_t>(n)));
		auto dp_initialized = dp;
		dp[0][0][0] = 1;
		
		std::vector<val_t> res = { 1 };
		for (int i = 0; i <= n; i++) {
			// Step 1: allow closing pairs of horizontal lines on vertical line #i
			for (int after_start = 0; after_start < 2; after_start++) {
				// careful direction of loop for in-place update avoiding double-transition(k decreases by at least one in traisition)
				for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[after_start][k][l]) {
					// closable: number of pairs of horizontal lines that can be closed
					// 	if after_start: all the k pairs of horizontal lines may be closed
					// 	if !after_start: closing all pairs would leave no space for creating start point
					int closable = after_start ? k : k - 1;
					for (int m = 1; m <= closable && l + m < n; m++) { // actually closing `m` pairs out of `closable`
						// k: decreases by m because m pairs were closed
						// l: increases by m because closing a pair creates a vertex
						dp[after_start][k - m][l + m] += dp[after_start][k][l] * comb.nCr(closable, m);
					}
				}
			}
			
			decltype(dp) dp_next = dp_initialized;
			// Step 2: allow knotting on vertical line #i
			// NOTE that transitions in Step 2 and Step 3 are done simultaneously from `dp` to `dp_next`:
			//    values added in Step 2 does not affect Step 3
			// Step 2 only takes in account cases where b_i = b_{i - 1} (or, as if we could only choose b_0, b_1, ..., b_{i-1})
			// This helps calculate answers for lower n(Step 2.5), where b_i does not exist
			// Then Step 3 adds the case when b_i != b_{i-1}
			{ // case A: when continuing normally
				for (int after_start = 0; after_start < 2; after_start++) {
					for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[after_start][k][l]) {
						// create `knot` vertices on some of the 2*k horizontal lines
						// knot <= n - 1 - l: we already created l vertices, so only n - 1 - l remaining
						// knot <= (n - 1 - l) - k + 1: we need to create at least (k-1) vertex to close k pairs of horizontal lines
						//   (one pair might be closed when placing the starting point (case B below) and setting top_dir and bottom_dir left)
						//   this constraint is only for pruning (avoid creating unnecessary non-zero element in `dp`) so not mandatory 
						for (int knot = 0; knot <= std::min({n - 1 - l, 2 * k, (n - 1 - l) - k + 1}); knot++) 
							dp_next[after_start][k][l + knot] += dp[after_start][k][l] * comb.nCr(2 * k, knot);
					}
				}
			}
			{ // case B: when placing the starting point on vertical line #i 
				for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[0][k][l]) {
					// We can free choose the direction of horizontal line that extends from the start vertex (top_dir)
					// also we can choose the direction of horizontal line that comes into the start vertex at the bottom (bottom_dir)
					// These will affect the new k because either being right increases it by 0.5 and being left decreases it by 0.5
					// (top|bottom)_dir : 0:left, 1:right
					for (int top_dir = 0; top_dir < 2; top_dir++) for (int bottom_dir = 0; bottom_dir < 2; bottom_dir++) {
						// top_dir and bottom_dir cannot be left if there is currently no pairs of horizontal lines
						if (!k && (!top_dir || !bottom_dir)) continue;
						
						// if top_dir is left, that left-going horizontal line cannot be knotted because start vertex is already there
						// same thing for when bottom_dir is left
						int left_knottable = k - !top_dir;
						int right_knottable = k - !bottom_dir;
						
						// k changes as stated above
						int k_diff = (top_dir + bottom_dir) - 1;
						if (k + k_diff > max_lines) continue;
						assert(k + k_diff >= 0);
						
						// create `knot` vertices (apart from the start vertex which is not counted in l in the first place)
						// on some of (left_knottable + right_knottable) horizontal lines
						for (int knot = 0; knot <= std::min({left_knottable + right_knottable, n - 1 - l, (n - 1 - l) - (k + k_diff) + 1}); knot++) 
							dp_next[1][k + k_diff][l + knot] += dp[0][k][l] * comb.nCr(left_knottable + right_knottable, knot);
					}
				}
			}
			
			// Step 2.5: calculate answer for n = i
			if (i) {
				val_t cur_res = 0;
				int cur_n = i;
				// j remaining vertices: create j additional vertices at existing knots, creating multiple knots(corresponding to same values adjacent in c)
				// here we can also create additional vertices at the starting vertex(both one at the top and one at the bottom which are distinguished)
				for (int j = 0; j < cur_n; j++) cur_res += dp_next[1][0][cur_n - 1 - j] * comb.nHr(cur_n + 1 - j, j);
				// we did not count the case when c = [i, i, ..., i] (0 <= i <= n), so add here (b can be any of 2^n 01 sequences)
				cur_res += ((val_t) 1 << cur_n) * (cur_n + 1);
				std::cerr << cur_n << " " << cur_res << std::endl;
				res.push_back(cur_res);
			}
			
			
			// Step 3: knot on vertical line #i (b_i != b_{i-1} case)
			// either left-going or right-going horizontal lines should all be knotted(here 'either' does not mean we can choose which)
			{ // when continuing normally
				for (int after_start = 0; after_start < 2; after_start++) {
					for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[after_start][k][l]) {
						// Same as Step 2 but:
						// If all k left-going lines are knotted, we can only choose which (knot - k) lines out of k right-going lines will be knotted
						// Above case calculates when (b_{i-1}, b_i) = (0, 1); same for when (1, 0) (which results in the same number so not calculated twice)
						for (int knot = k; knot <= std::min({n - 1 - l, 2 * k, (n - 1 - l) - k + 1}); knot++) 
							dp_next[after_start][k][l + knot] += dp[after_start][k][l] * comb.nCr(k, knot - k);
					}
				}
			}
			{ // when passing the starting point 
				for (int k = 0; k <= max_lines; k++) for (int l = 0; l < n; l++) if (dp[0][k][l]) {
					// same as Step 2 but we can only choose which (knot - left_knottable) out of right_knottable right-going lines will be knotted
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
			
			// Step 4: allow opening pairs of horizontal lines on vertical line #i
			for (int after_start = 0; after_start < 2; after_start++) {
				// careful direction of loop for in-place update avoiding double-transition(k increases by at least one in traisition)
				for (int k = max_lines; k >= 0; k--) for (int l = 0; l < n; l++) if (dp_next[after_start][k][l]) {
					// after_start: new pair of horizontal lines can be opened anywhere in the middle of two horizontal lines currently paired
					// !after_start: new pair of horizontal lines can be opened anywhere between two existing pairs
					int openable_place_num = after_start ? k : k + 1;
					if (openable_place_num) {
						// actually open m pairs
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

