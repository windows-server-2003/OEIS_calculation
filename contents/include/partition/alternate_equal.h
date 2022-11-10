
/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that s_i == s_{i + 1} for every even i
*/
struct A351003 {
	template<typename T> std::vector<T> calc_all(int n) {
		int half = n / 2;
		std::vector<T> dp(n + 1);
		for (int i = half + 1; i <= n; i++) dp[i] = 1;
		
		std::vector<T> res_acc(n + 1);
		
		for (int i = half; i >= 1; i--) {
			dp[i] += 1;
			for (int j = 2 * i; j <= n; j++) {
				auto increase = dp[j - 2 * i];
				res_acc[j] += increase;
				if (j + i + 1 <= n) res_acc[j + i + 1] -= increase;
				
				dp[j] += increase;
			}
		}
		for (int i = 1; i <= n; i++) res_acc[i] += res_acc[i - 1];
		for (int i = 1; i <= n; i++) res_acc[i] += 1; // length 1
		for (int i = 2; i <= n; i++) res_acc[i] += i / 2; // length 2
		
		res_acc[0] = 1;
		return res_acc;
	}
};

/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that s_i == s_{i + 1} for every odd i
*/
struct A351004 {
	template<typename T> std::vector<T> calc_all(int n) {
		int half = n / 2;
		std::vector<T> dp(half + 1);
		dp[0] = 1;
		
		std::vector<T> res_acc(n + 1);
		res_acc[0] += 1;
		
		for (int i = half; i >= 1; i--) {
			for (int j = i; j <= half; j++) {
				auto increase = dp[j - i];
				res_acc[2 * j] += increase;
				if (2 * j + i + 1 <= n) res_acc[2 * j + i + 1] -= increase;
				
				dp[j] += increase;
			}
		}
		for (int i = 1; i <= n; i++) res_acc[i] += res_acc[i - 1];
		return res_acc;
	}
};


/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that
		 - s_i == s_{i + 1} for every odd i
		 - s_i != s_{i + 1} for every even i
*/
struct A351005 {
	template<typename T> std::vector<T> calc_all(int n) {
		int half = n / 2;
		std::vector<T> dp(half + 1);
		dp[0] = 1;
		
		std::vector<T> res_acc(n + 1);
		res_acc[0] += 1;
		
		for (int i = half; i >= 1; i--) {
			for (int j = i; j <= half; j++) {
				auto increase = dp[j - i];
				res_acc[2 * j] += increase;
				if (2 * j + i <= n) res_acc[2 * j + i] -= increase;
			}
			for (int j = half - i; j >= 0; j--) dp[j + i] += dp[j];
		}
		for (int i = 1; i <= n; i++) res_acc[i] += res_acc[i - 1];
		return res_acc;
	}
};


/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that
		 - s_i != s_{i + 1} for every odd i
		 - s_i == s_{i + 1} for every even i
*/
struct A351006 {
	template<typename T> std::vector<T> calc_all(int n) {
		int half = n / 2;
		std::vector<T> dp(n + 1);
		for (int i = half + 1; i <= n; i++) dp[i] = 1;
		
		std::vector<T> res_acc(n + 1);
		
		for (int i = half; i >= 1; i--) {
			for (int j = n; j >= 2 * i; j--) {
				dp[j] += dp[j - 2 * i];
				auto increase = dp[j - 2 * i];
				res_acc[j] += increase;
				if (j + i <= n) res_acc[j + i] -= increase;
			}
			dp[i] += 1;
		}
		for (int i = 1; i <= n; i++) res_acc[i] += res_acc[i - 1];
		for (int i = 1; i <= n; i++) res_acc[i] += 1; // length 1
		for (int i = 3; i <= n; i++) res_acc[i] += (i - 1) / 2; // length 2
		
		res_acc[0] = 1;
		return res_acc;
	}
};

/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that
		 - k is even
		 - s_i != s_{i + 1} for every odd i
		 - s_i == s_{i + 1} for every even i
*/
struct A351007 {
	template<typename T> std::vector<T> calc_all(int n) {
		int half = n / 2;
		std::vector<T> dp(n + 1);
		for (int i = half + 1; i <= n; i++) dp[i] = 1;
		
		std::vector<T> res_acc(n + 1);
		
		for (int i = half; i >= 1; i--) {
			for (int j = n; j >= 2 * i; j--) {
				dp[j] += dp[j - 2 * i];
				auto increase = dp[j - 2 * i];
				if (j + 1 <= n) res_acc[j + 1] += increase;
				if (j + i <= n) res_acc[j + i] -= increase;
			}
			dp[i] += 1;
		}
		for (int i = 1; i <= n; i++) res_acc[i] += res_acc[i - 1];
		for (int i = 3; i <= n; i++) res_acc[i] += (i - 1) / 2; // length 2
		
		res_acc[0] = 1;
		return res_acc;
	}
};


/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that
		 - k is even
		 - s_i != s_{i + 1} for every odd i
*/
struct A351008 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<T> dp0(n + 1);
		std::vector<T> dp1(n + 1);
		dp0[0] = 1;
		
		for (int i = n; i >= 1; i--) {
			for (int j = n - i; j >= 0; j--) {
				dp1[j + i] += dp0[j];
				dp0[j + i] += dp1[j];
				if (j + i + i <= n) dp1[j + i + i] += dp1[j];
			}
		}
		return dp0;
	}
};

/*
	a(n) : number of partition of n (s_1 >= s_2 >= ... >= s_k ((sum s_i) = n))
		such that
		 - k is even
		 - s_i == s_{i + 1} for every even i
*/
struct A351012 {
	template<typename T> std::vector<T> calc_all(int n) {
		int half = n / 2;
		std::vector<T> dp(n + 1);
		for (int i = half + 1; i <= n; i++) dp[i] = 1;
		
		std::vector<T> res_acc(n + 1);
		
		for (int i = half; i >= 1; i--) {
			dp[i] += 1;
			for (int j = 2 * i; j <= n; j++) {
				auto increase = dp[j - 2 * i];
				if (j + 1 <= n) res_acc[j + 1] += increase;
				if (j + i + 1 <= n) res_acc[j + i + 1] -= increase;
				
				dp[j] += increase;
			}
		}
		for (int i = 1; i <= n; i++) res_acc[i] += res_acc[i - 1];
		for (int i = 2; i <= n; i++) res_acc[i] += i / 2; // length 2
		
		res_acc[0] = 1;
		return res_acc;
	}
};
