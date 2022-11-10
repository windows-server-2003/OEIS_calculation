#include "partition.h"

/*
	Integer partition related sequence
	
	s = (s_1, s_2, ..., s_k) is a partition of n iff :
		 - s_i is an integer for every i
		 - s_1 >= s_2 >= ... >= s_k >= 1
		 - (sum s_i) = n
	A partition s of n is strict iff :
		 - s_i > s_{i + 1} for every i
	
	The conjugate partition of s is the sequence t = (t_1, t_2, ..., t_{max s_i}) such that :
		 - t_i = #{j | s_j >= i}
	Here, t is also a partition of n
*/

/*
	a(n) : number of partition s such that :
		 - s_i != 2 * s_{i + 1} for every i
*/
struct A350837 {
	template<typename T> std::vector<T> calc_all(int n) {
		assert(n >= 1);
		int half = n / 2;
		std::vector<std::vector<T> > dp_offset(half + 1, std::vector<T>(n + 1));
		std::vector<T> dp(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i++) {
			for (int j = 0; j <= n - i; j++) { // use i here and transit
				auto cur = dp[j];
				if (!(i & 1)) cur += dp_offset[i >> 1][j];
				
				if (i <= half) dp_offset[i][j + i] -= cur;
				dp[j + i] += cur;
			}
		}
		return dp;
	}
};


/*
	a(n) : number of partition s such that :
		 - s_i - 1 > s_{i + 1} holds for some i
		 - t_i - 1 > t_{i + 1} holds for some i (where t is the conjugate partition of s)
*/
struct A350839 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<std::vector<T> > dp1(6, std::vector<T>(n + 1));
		std::vector<std::vector<T> > dp2(6, std::vector<T>(n + 1));
		for (int i = 1; i <= n; i++) {
			std::vector<std::vector<T> > dp0(6, std::vector<T>(n + 1));
			
			dp0[0][i] += 1;
			for (int j = 0; j < 6; j++) for (int k = 0; k + i <= n; k++) {
				int next_j = (j >= 2 && j < 4) ? j + 2 : j;
				dp0[next_j][k + i] += dp1[j][k];
				dp0[next_j | 1][k + i] += dp2[j][k];
			}
			for (int j = 0; j < 6; j++) for (int k = 0; k + i <= n; k++) 
				dp0[j < 2 ? j + 2 : j][k + i] += dp0[j][k];
			
			for (int j = 0; j < 6; j++) for (int k = 0; k <= n; k++) dp2[j][k] += dp1[j][k];
			dp1 = dp0;
		}
		
		std::vector<T> res(n + 1);
		for (int i = 0; i <= n; i++) res[i] = dp1[5][i] + dp2[5][i];
		
		return res;
	}
};


/*
	a(n) : number of partition s such that :
		 - s is strict
		 - s_i != 2 * s_{i + 1} for every i
*/
struct A350840 {
	template<typename T> std::vector<T> calc_all(int n) {
		assert(n >= 1);
		int half = n / 2;
		std::vector<std::vector<T> > dp_offset(half + 1, std::vector<T>(n + 1));
		std::vector<T> dp(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i++) {
			for (int j = n - i; j >= 0; j--) { // use i here and transit
				auto cur = dp[j];
				if (!(i & 1)) cur += dp_offset[i >> 1][j];
				
				if (i <= half) dp_offset[i][j + i] -= cur;
				dp[j + i] += cur;
			}
		}
		
		return dp;
	}
};


/*
	a(n) : number of partition s such that :
		 - s_i - 2 != s_{i + 1} for every i
*/
struct A350842 {
	template<typename T> std::vector<T> calc_all(int n) {
		assert(n >= 1);
		std::vector<T> dp1(n + 1);
		std::vector<T> dp2(n + 1);
		std::vector<T> dp3(n + 1);
		dp3[0] = 1;
		
		for (int i = 1; i <= n; i++) {
			std::vector<T> dp0(n + 1);
			for (int j = 0; j + i <= n; j++) dp0[j + i] += dp1[j] + dp3[j] + dp0[j];
			
			for (int j = 0; j <= n; j++) dp3[j] += dp2[j];
			dp2 = std::move(dp1);
			dp1 = std::move(dp0);
		}
		
		std::vector<T> res(n + 1);
		for (int i = 0; i <= n; i++) res[i] = dp1[i] + dp2[i] + dp3[i];
		return res;
	}
};


/*
	a(n) : number of partition s such that :
		 - s is strict
		 - s_i - 2 != s_{i + 1} for every i
*/
struct A350844 {
	template<typename T> std::vector<T> calc_all(int n) {
		assert(n >= 1);
		std::vector<T> dp1(n + 1);
		std::vector<T> dp2(n + 1);
		std::vector<T> dp3(n + 1);
		dp3[0] = 1;
		
		for (int i = 1; i <= n; i++) {
			std::vector<T> dp0(n + 1);
			for (int j = 0; j + i <= n; j++) dp0[j + i] += dp1[j] + dp3[j];
			
			for (int j = 0; j <= n; j++) dp3[j] += dp2[j];
			dp2 = std::move(dp1);
			dp1 = std::move(dp0);
		}
		
		std::vector<T> res(n + 1);
		for (int i = 0; i <= n; i++) res[i] = dp1[i] + dp2[i] + dp3[i];
		return res;
	}
};


/*
	a(n) : number of partition s such that :
		 - s_i - 2 == s_{i + 1} for some i
*/
struct A350846 {
	template<typename T> std::vector<T> calc_all(int n) {
		auto res = get_partition_sequence<T>(n);
		auto complement = A350837().calc_all<T>(n);
		for (int i = 0; i <= n; i++) res[i] -= complement[i];
		return res;
	}
};


/*
	a(n) : number of partition s such that :
		 - s_i != 2 * s_j for every i, j
*/
struct A323092 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<T> dp(n + 1);
		std::vector<T> dp_lastused(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i += 2) {
			for (int t = i; t <= n; t <<= 1) {
				std::vector<T> dp_next_lastused(n + 1);
				for (int j = 0; j + t <= n; j++) dp_next_lastused[j + t] = dp[j] + dp_next_lastused[j];
				for (int j = 0; j <= n; j++) dp[j] += dp_lastused[j];
				std::swap(dp_lastused, dp_next_lastused);
			}
			for (int j = 0; j <= n; j++) dp[j] += dp_lastused[j];
			dp_lastused.assign(n + 1, 0);
		}
		
		return dp;
	}
};
/*
	a(n) : number of partition s such that :
		 - s_i != 2^k * s_j for every i, j, k>0
*/
struct A323093 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<T> dp(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i += 2) {
			std::vector<T> used(n + 1);
			for (int t = i; t <= n; t <<= 1) {
				std::vector<T> tmp(n + 1);
				for (int j = 0; j + t <= n; j++) tmp[j + t] = dp[j] + tmp[j];
				for (int j = t; j <= n; j++) used[j] += tmp[j];
			}
			for (int j = 0; j <= n; j++) dp[j] += used[j];
		}
		
		return dp;
	}
};
/*
	a(n) : number of partition s such that :
		 - s is strict
		 - s_i != 2^k * s_j for every i, j, k>0
*/
struct A323094 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<T> dp(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i += 2) {
			std::vector<T> used(n + 1);
			for (int t = i; t <= n; t <<= 1) 
				for (int j = t; j <= n; j++) used[j] += dp[j - t];
			for (int j = 0; j <= n; j++) dp[j] += used[j];
		}
		
		return dp;
	}
};


/*
	a(n) : number of partition s such that :
		 - s is strict
		 - s_i != 2 * s_j for every i, j
*/
struct A120641 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<T> dp(n + 1);
		std::vector<T> dp_lastused(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i += 2) {
			for (int t = i; t <= n; t <<= 1) {
				std::vector<T> dp_next_lastused(n + 1);
				for (int j = 0; j + t <= n; j++) dp_next_lastused[j + t] = dp[j];
				for (int j = 0; j <= n; j++) dp[j] += dp_lastused[j];
				std::swap(dp_lastused, dp_next_lastused);
			}
			for (int j = 0; j <= n; j++) dp[j] += dp_lastused[j];
			dp_lastused.assign(n + 1, 0);
		}
		
		return dp;
	}
};


/*
	a(n) : number of partition s such that :
		 - s_i is not a perfect power(integers that can be expressed m^a (a >= 2)) including 1
*/
struct A323089 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<bool> is_usable(n + 1, true);
		for (s64 i = 2; i * i <= n; i++) for (s64 j = i * i; j <= n; j *= i) is_usable[j] = false;
		
		std::vector<T> dp(n + 1);
		dp[0] = 1;
		
		for (int i = 1; i <= n; i++) if (is_usable[i]) for (int j = n - i; j >= 0; j--) dp[j + i] += dp[j];
		
		return dp;
	}
};

/*
	a(n) : number of partition s such that :
		 - s_i is not a perfect power(integers that can be expressed m^a (a >= 2)) excluding 1
*/
struct A323088 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<bool> is_usable(n + 1, true);
		for (s64 i = 2; i * i <= n; i++) for (s64 j = i * i; j <= n; j *= i) is_usable[j] = false;
		
		std::vector<T> dp(n + 1);
		dp[0] = 1;
		
		for (int i = 2; i <= n; i++) if (is_usable[i]) for (int j = n - i; j >= 0; j--) dp[j + i] += dp[j];
		
		return dp;
	}
};

/*
	a(n) : number of partition {s_k} such that :
		 - k >= 3
*/
struct A004250 {
	template<typename T> std::vector<T> calc_all(int n) {
		auto res = get_partition_sequence<T>(n);
		for (int i = 0; i <= n; i++) res[i] -= 1 + i / 2;
		return res;
	}
};



/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that
		 - s can be arranged so that it is a alternately up-down(starting with either) sequence
*/
struct A351012 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<std::vector<T> > dp(n + 1);
		for (int i = 0; i <= n; i++) dp[i].resize(i + 1);
		dp[0][0] = 1;
		for (int i = 1; i <= n; i++) for (int j = 1; j <= i; j++) 
			dp[i][j] = dp[i - 1][j - 1] + (i - j >= j ? dp[i - j][j] : 0);
		
		T res = 0;
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j * i <= n; j++) {
				int remaining = n - i * j;
				for (int k = 0; k <= j - 2 && k <= remaining; k++) res += dp[remaining][k];
				if (remaining >= i) {
					for (int k = 0; k <= j - 3 && k <= remaining - i; k++) res -= dp[remaining - i][k];
				}
				std::cerr << i << " x " << j << " : " << res << std::endl;
			}
		}
		
		return std::vector<T>(n + 1);
	}
};



