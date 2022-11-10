// #define _GLIBCXX_DEBUG
#include <cstdio>
#include <iostream>
#include <vector>
#include <array>
#include "debug.h"
#include "types.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

struct A189281 {
	/*
		calc<val_t>(n) :
			returns the n-th term of A189281 using val_t as the type used in calculation
	*/
	template<typename val_t> val_t calc(int n) {
		if (n == 0) return 1;
		
		int half = (n - 1) >> 1;
		
		/*
			dp[g_pos][left_violation_num][right_violation_num][prev_pos]
			g_pos : the number of element in the left of g
			left_violation_num : the number of violating left gaps
			right_violation_num : the number of violating right gaps
			prev_pos : the position of the last-inserted element(1 : the left element of g, 0 : anything left, 2 : anything right)
		*/
		std::vector<std::vector<std::vector<std::array<val_t, 3> > > > dp(n + 1,
			std::vector<std::vector<std::array<val_t, 3> > >(n,
				std::vector<std::array<val_t, 3> >(n)));
		std::vector<std::vector<std::vector<std::array<val_t, 3> > > > next;
		
		// start with [0]
		dp[0][0][0][2] = 1; // choose the left of 0 as g
		dp[1][0][0][1] = 1; // choose the right of 0 as g
		
		for (int i = 1; i < n; i++) {
			// insert i
			
			// clear next with zero
			next.assign(n + 1, std::vector<std::vector<std::array<val_t, 3> > >(n,
				std::vector<std::array<val_t, 3> >(n)));
			
			bool disable_violation = (i - 1 == half);
			for (int g_pos = 0; g_pos <= i; g_pos++) {
				for (int left_violation = 0; left_violation <= std::max(0, g_pos - 1); left_violation++) {
					for (int right_violation = 0; right_violation <= std::max(0, i - g_pos - 1); right_violation++) {
						for (int prev_pos = 0; prev_pos < 3; prev_pos++) {
							if (dp[g_pos][left_violation][right_violation][prev_pos] == 0) continue;
							auto &cur_val = dp[g_pos][left_violation][right_violation][prev_pos];
							
							int left_nonviolating = g_pos - left_violation;
							int right_nonviolating = (i - g_pos) - right_violation;
							// insert into left, non-violating gap
							if (left_nonviolating - (prev_pos == 0) > 0) next[g_pos + 1][left_violation][right_violation][0] += cur_val * (left_nonviolating - (prev_pos == 0));
							if (prev_pos == 0) // creating a new violation due to inserting i immediately after i - 1
								next[g_pos + 1][left_violation + !disable_violation][right_violation][0] += cur_val;
							
							// insert into left, violating gap
							if (left_violation) next[g_pos + 1][left_violation - 1][right_violation][0] += cur_val * left_violation;
							
							// insert into g; g will be split into two
							next[g_pos][left_violation][right_violation][2] += cur_val; // choose the left one as new g
							next[g_pos + 1][left_violation + (prev_pos == 1 && !disable_violation)][right_violation][1] += cur_val; // choose the right one as new g
							
							// insert into right, non-violating gap
							if (right_nonviolating - (prev_pos == 2) > 0) next[g_pos][left_violation][right_violation][2] += cur_val * (right_nonviolating - (prev_pos == 2));
							if (prev_pos == 2) // creating a new violation due to inserting i immediately after i - 1
								next[g_pos][left_violation][right_violation + !disable_violation][2] += cur_val;
							
							// insert into right, violating gap
							if (right_violation) next[g_pos][left_violation][right_violation - 1][2] += cur_val * right_violation;
						}
					}
				}
			}
			std::swap(dp, next);
		}
		val_t res = 0;
		for (int prev_pos = 0; prev_pos < 3; prev_pos++) res += dp[half + 1][0][0][prev_pos];
		
		return res;
	}
};

template<typename val_t> void check_recurrence(const std::vector<val_t> &a, int n_) {
	val_t n[12];
	n[1] = n_;
	for (int i = 2; i < 12; i++) n[i] = n[i - 1] * n_;
	val_t test_val = 
		(262711*n[1] + 1387742*n[2] - 824875*n[3] - 1855253*n[4] - 111530*n[5] + 680983*n[6] + 364242*n[7] + 84992*n[8] + 10332*n[9] + 640*n[10] + 16*n[11])*a[n_] +
		(-1050844*n[1] - 9705192*n[2] - 7414683*n[3] + 3536494*n[4] + 6459004*n[5] + 3326393*n[6] + 903534*n[7] + 144684*n[8] + 13756*n[9] + 720*n[10] + 16*n[11])*a[n_+1] +
		(3492344 - 2212342*n[1] - 8507169*n[2] - 11544227*n[3] - 12034116*n[4] - 8216995*n[5] - 3442049*n[6] - 890050*n[7] - 142300*n[8] - 13660*n[9] - 720*n[10] - 16*n[11])*a[n_+2] +
		(19817984 + 45323852*n[1] + 825228*n[2] - 57004661*n[3] - 57059306*n[4] - 28077270*n[5] - 8398637*n[6] - 1631510*n[7] - 207980*n[8] - 16828*n[9] - 784*n[10] - 16*n[11])*a[n_+3] +
		(9586160 + 6680237*n[1] - 13772613*n[2] - 27689586*n[3] - 22162455*n[4] - 9855085*n[5] - 2629562*n[6] - 427656*n[7] - 41332*n[8] - 2176*n[9] - 48*n[10])*a[n_+4] +
		(22192864 + 44710768*n[1] - 2924668*n[2] - 52385912*n[3] - 45161616*n[4] - 18784740*n[5] - 4549208*n[6] - 674256*n[7] - 60400*n[8] - 3008*n[9] - 64*n[10])*a[n_+5] +
		(557152 - 2032472*n[1] - 2937392*n[2] - 1594200*n[3] - 517688*n[4] - 122032*n[5] - 19856*n[6] - 1792*n[7] - 64*n[8])*a[n_+6] +
		(3786960 + 7105324*n[1] - 1191064*n[2] - 8059160*n[3] - 5938996*n[4] - 2073752*n[5] - 402736*n[6] - 44528*n[7] - 2624*n[8] - 64*n[9])*a[n_+7] +
		(-598208 - 943004*n[1] + 414196*n[2] + 1213772*n[3] + 728648*n[4] + 203584*n[5] + 29616*n[6] + 2176*n[7] + 64*n[8])*a[n_+8];
	assert(test_val == 0);
}

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

