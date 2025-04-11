#include <vector>
#include <queue>
#include "comb.h"

/*
	Calculates permanents of all principal submatrices of the following matrix m[i][j] in O(N^2) :
		m[i][j] = a(i)f(j) + b(i)
			where
				 - a(i) and b(i) are specified by row_func(i).first, row_func(i).second, respectively
				 - f(j) is specified by col_func(j)
	
	res_vec[i] contains the permanent of the top-left ixi submatrix
*/
template<typename T, typename row_func_t, typename col_func_t> std::vector<T> permanents_of_rowwise_linear_n_independent_matrix(int N, const row_func_t &row_func, const col_func_t &col_func) {
	std::vector<T> fact(N + 1, 1);
	for (int i = 1; i <= N; i++) fact[i] = fact[i - 1] * i;
	
	std::vector<T> poly0 = { 1 };
	std::vector<T> poly1 = { 1 };
	std::vector<T> res = { 1 };
	for (int i = 0; i < N; i++) {
		{
			poly0.push_back(0);
			auto coef = row_func(i);
			for (int j = i; j >= 0; j--) {
				poly0[j + 1] += poly0[j] * coef.first * (j + 1);
				poly0[j] *= coef.second;
			}
		}
		{
			poly1.push_back(0);
			auto val = col_func(i);
			for (int j = i; j >= 0; j--) poly1[j + 1] += poly1[j] * val;
		}
		int n = i + 1;
		T cur_res = 0;
		for (int j = 0; j <= n; j++) cur_res += poly0[j] * poly1[j] * fact[n - j];
		res.push_back(cur_res);
	}
	return res;
}
/*
	Same as permanents_of_rowwise_linear_n_independent_matrix except that this one does not calculate permanents of the submatrices and is slightly faster
*/
template<typename T, typename row_func_t, typename col_func_t> T permanent_of_rowwise_linear_matrix(int n, const row_func_t &row_func, const col_func_t &col_func) {
	std::vector<T> fact(n + 1, 1);
	for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i;
	
	std::vector<T> poly0 = { 1 };
	std::vector<T> poly1 = { 1 };
	for (int i = 0; i < n; i++) {
		{
			poly0.push_back(0);
			auto coef = row_func(n, i);
			for (int j = i; j >= 0; j--) {
				poly0[j + 1] += poly0[j] * coef.first;
				poly0[j] *= coef.second;
			}
		}
		{
			poly1.push_back(0);
			auto val = col_func(n, i);
			for (int j = i; j >= 0; j--) poly1[j + 1] += poly1[j] * val;
		}
	}
	T res = 0;
	for (int j = 0; j <= n; j++) res += poly0[j] * poly1[j] * fact[j] * fact[n - j];
	return res;
}


template<typename T, typename row_func_t, typename col_func_t> std::vector<T> permanents_of_rowwise_quadratic_n_independent_matrix(int N, const row_func_t &row_func, const col_func_t &col_func) {
	std::vector<T> fact(N + 1, 1);
	for (int i = 1; i <= N; i++) fact[i] = fact[i - 1] * i;
	
	std::vector<std::vector<T> > poly0 = { { 1 } };
	std::vector<std::vector<T> > poly1 = { { 1 } };
	std::vector<T> res = { 1 };
	for (int i = 0; i < N; i++) {
		std::cerr << "stage " << i << std::endl;
		{
			for (auto &j : poly0) j.push_back(0);
			poly0.push_back({ 0 });
			auto coefs = row_func(i);
			for (int j = i; j >= 0; j--) for (int k = i - j; k >= 0; k--) {
				poly0[j + 1][k] += poly0[j][k] * std::get<1>(coefs) * (j + 1);
				poly0[j][k + 1] += poly0[j][k] * std::get<2>(coefs) * (k + 1);
				poly0[j][k] *= std::get<0>(coefs);
			}
		}
		{
			for (auto &j : poly1) j.push_back(0);
			poly1.push_back({ 0 });
			auto val = col_func(i);
			auto val2 = val * val;
			for (int j = i; j >= 0; j--) for (int k = i - j; k >= 0; k--) {
				poly1[j + 1][k] += poly1[j][k] * val;
				poly1[j][k + 1] += poly1[j][k] * val2;
			}
		}
		
		int n = i + 1;
		T cur_res = 0;
		for (int sum = 0; sum <= n; sum++) {
			T cur_sum = 0;
			for (int j = 0; j <= sum; j++) cur_sum += poly0[j][sum - j] * poly1[j][sum - j];
			cur_res += cur_sum * fact[n - sum];
		}
		res.push_back(cur_res);
	}
	return res;
}

/*
	Calculates permanents of the principal submatrices of the following NxN matrices m[i][j] in O(N^2) :
		m[i][i] = f(i)
		m[i][j] = 1   (i != j)
			where f(i) is specified by val_func(i)
	
	res_vec[i] contains the permanent of the top-left ixi submatrix
*/
template<typename T, typename val_func_t> std::vector<T> permanents_of_diagonal_plus1_matrix(int N, const val_func_t &val_func) {
	std::vector<T> monmort(N + 1);
	monmort[0] = 1;
	for (int i = 2; i <= N; i++) monmort[i] = (monmort[i - 1] + monmort[i - 2]) * (i - 1);
	
	auto dp = monmort;
	std::vector<T> res = { 1 };
	for (int i = 0; i < N; i++) {
		auto val = val_func(i);
		for (int j = N; j--; ) dp[j + 1] += dp[j] * val;
		
		res.push_back(dp[i + 1]);
	}
	return res;
}

// left_func / right_func : (n, row, column cnt, column sum) -> element sum
// old function
template<typename T, typename left_func_t, typename right_func_t> T permanent_of_half_linear_matrix2(int n, const left_func_t &left_func, const right_func_t &right_func) {
	std::vector<std::vector<std::vector<std::vector<T> > > > dp(1, std::vector<std::vector<std::vector<T> > >(1, 
		std::vector<std::vector<T> >(n + 1)));
	for (int i = 0; i <= n; i++) dp[0][0][i].assign(i * (n - i) + 1, 1);
	
	for (int i = 0; i < n; i++) {
		dp.resize(i + 2);
		for (int j = 0; j <= i + 1; j++) dp[j].resize(j * (i + 1 - j) + 1);
		dp[i + 1][0].resize(n - i);
		for (int j = 0; j < n - i; j++) dp[i + 1][0][j].resize(j * (n - i - 1 - j) + 1);
		
		for (size_t j = dp.size(); j--; ) {
			for (size_t k = 0; k < dp[j].size(); k++) {
				int left_index_sum = j * (j - 1) / 2 + k;
				int left_val_sum = left_func(n, i, j, left_index_sum);
				for (size_t l = dp[j][k].size(); l--; ) {
					for (size_t m = 0; m < dp[j][k][l].size(); m++) {
						if (dp[j][k][l][m] == 0) continue;
						int right_sum = i * l + l * (l - 1) / 2 + m;
						// std::cerr << i << " " << j << " " << k << " " << l << " " << m << " : " << dp[j][k][l][m] << "x" << ((k + (n - i) * j) + (n * l - right_sum + i * l)) << std::endl;
						dp[j][k][l][m] *= left_val_sum + right_func(n, i, l, right_sum);
						// use
						if (l && m <= (l - 1) * (n - i - 1 - (l - 1))) dp[j + 1][k + (i - j)][l - 1][m] -= dp[j][k][l][m];
						// don't use
						if ((int) l < n - i && m >= l && m - l <= l * (n - i - 1 - l)) {
							if (l) dp[j][k][l][m - l] += dp[j][k][l][m], dp[j][k][l][m] = 0;
						} else dp[j][k][l][m] = 0;
					}
				}
				dp[j][k].resize(n - i);
				for (size_t l = 0; l < dp[j][k].size(); l++) dp[j][k][l].resize(l * (n - i - 1 - l) + 1);
			}
		}
	}
	T res = 0;
	for (size_t i = 0; i < dp.size(); i++) for (size_t j = 0; j < dp[i].size(); j++) 
		res += dp[i][j][0][0];
	return (n & 1) ? -res : res;
}

// left_func / right_func : (n, row, column cnt, column sum) -> element sum
template<typename T, typename left_func_t, typename right_func_t> T permanent_of_half_linear_matrix(int n, const left_func_t &left_func, const right_func_t &right_func) {
	std::vector<std::vector<std::vector<std::vector<T> > > > dp = { { std::vector<std::vector<T> >(n + 1, { 1 }) } };
	
	for (int i = 0; i < n; i++) {
		dp.resize(i + 2);
		for (int j = 0; j <= i + 1; j++) {
			dp[j].resize(j * (i + 1 - j) + 1, std::vector<std::vector<T> >(n - i, std::vector<T>(i + 2)));
			for (auto &k : dp[j]) for (auto &l : k) l.resize(i + 2);
		}
		
		for (size_t j = i + 1; j--; ) {
			for (size_t k = 0; k <= j * (i - j); k++) {
				int left_index_sum = j * (j - 1) / 2 + k;
				int left_val_sum = left_func(n, i, j, left_index_sum);
				for (size_t l = n - i + 1; l--; ) {
					auto r_coef0 = right_func(n, i, l, 0);
					auto r_coef1 = right_func(n, i, l, 1) - right_func(n, i, l, 0);
					auto coef0 = left_val_sum + r_coef0 - left_index_sum * r_coef1;
					auto coef1 = r_coef1;
					
					for (size_t m = i + 1; m--; ) {
						dp[j][k][l][m + 1] += dp[j][k][l][m] * coef1;
						dp[j][k][l][m] *= coef0;
					}
					auto add = [] (std::vector<T> &dst, const std::vector<T> &src) {
						for (size_t i = 0; i < src.size(); i++) dst[i] += src[i];
					};
					// use
					if (l) add(dp[j + 1][k + (i - j)][l - 1], dp[j][k][l]);
					// don't use : keep
				}
				dp[j][k].resize(n - i);
				for (auto &l : dp[j][k]) l.resize(i + 2);
			}
		}
	}
	
	T res = 0;
	for (size_t i = 0; i < dp.size(); i++) for (size_t j = 0; j < dp[i].size(); j++) {
		int sum = i * (i - 1) / 2 + j;
		T cur_res = 0;
		T cur_power = 1;
		for (auto &k : dp[i][j][0]) {
			cur_res += cur_power * k;
			cur_power *= sum;
		}
		if ((n - i) & 1) res -= cur_res;
		else res += cur_res;
	}
	return res;
}

template<typename T, typename val_func_t> std::vector<T> permanent_of_mod3_dependent_general(int n, const val_func_t &val_func) {
	constexpr int MOD = 3;
	int n0 = (n + 2) / MOD;
	int n1 = (n + 1) / MOD;
	std::vector<std::vector<T> > dp(n0 + 1, std::vector<T>(n1 + 1));
	dp[0][0] = 1;
	
	T cur_coef = 1;
	std::vector<T> res = { 1 };
	for (int i = 0; i < n; i++) {
		std::cerr << "stage " << i << std::endl;
		int mod = i % MOD;
		for (int j = n0; j >= 0; j--) for (int k = std::min(n1, i - j); k >= 0; k--) {
			if (j < n0) dp[j + 1][k] += dp[j][k] * val_func(0, mod);
			if (k < n1) dp[j][k + 1] += dp[j][k] * val_func(1, mod);
			dp[j][k] *= val_func(2, mod);
		}
		int cur_n = i + 1;
		auto cur_res = dp[(cur_n + 2) / MOD][(cur_n + 1) / MOD];
		cur_coef *= (i / MOD) + 1;
		res.push_back(cur_res * cur_coef);
	}
	return res;
}
template<typename T, typename val_func_t> std::vector<T> permanent_of_mod3_dependent(int n, const val_func_t &val_func) {
	constexpr int MOD = 3;
	std::vector<bool> has_zero_row(MOD), has_zero_column(MOD);
	for (int i = 0; i < MOD; i++) for (int j = 0; j < MOD; j++) if (val_func(i, j) == 0) has_zero_row[i] = has_zero_column[j] = true;
	
	auto all_one = [] (auto &a) { return std::accumulate(a.begin(), a.end(), 0) == (int) a.size(); };
	if (!all_one(has_zero_row) || !all_one(has_zero_column)) return permanent_of_mod3_dependent_general<T>(n, val_func);
	
	std::vector<std::vector<std::pair<int, T> > > value_coef(MOD, std::vector<std::pair<int, T> >(MOD));
	std::queue<std::pair<int, int> > que;
	for (int i = 0; i < MOD; i++) {
		for (int j = 0; j < MOD; j++) if (val_func(i, j)) {
			value_coef[i][j] = {1, 0};
			que.push({i, j});
			break;
		}
		if (que.size()) break;
	}
	int ns[MOD] = {(n + 2) / MOD, (n + 1) / MOD, n / MOD};
	while (que.size()) {
		auto i = que.front();
		auto cur_val = value_coef[i.first][i.second];
		que.pop();
		
		auto get_only_index = [] (int n, const auto &f) {
			int only = -1; // -1 : not yet assigned, -2 : multiple found
			for (int i = 0; i < n; i++) if (f(i)) {
				if (only == -1) only = i;
				else only = -2;
			}
			return only < 0 ? -1 : only;
		};
		{
			int only_nonzero = get_only_index(MOD, [&] (int j) { return j != i.first && val_func(j, i.second); });
			if (only_nonzero >= 0 && value_coef[only_nonzero][i.second].first == 0) {
				value_coef[only_nonzero][i.second] = {-cur_val.first, ns[i.second] - cur_val.second};
				que.push({only_nonzero, i.second});
			}
		}
		{
			int only_nonzero = get_only_index(MOD, [&] (int j) { return j != i.second && val_func(i.first, j); });
			if (only_nonzero >= 0 && value_coef[i.first][only_nonzero].first == 0) {
				value_coef[i.first][only_nonzero] = {-cur_val.first, ns[i.first] - cur_val.second};
				que.push({i.first, only_nonzero});
			}
		}
	}
	for (int j = 0; j < MOD; j++) for (int k = 0; k < MOD; k++) if (val_func(j, k) && !value_coef[j][k])
		return permanent_of_mod3_dependent_general<T>(n, val_func);
	
	int min = 0;
	int max = 1000000000;
	for (int j = 0; j < MOD; j++) for (int k = 0; k < MOD; k++) if (val_func(j, k) && !value_coef[j][k])
	
	for (int j = 0; j < MOD; j++) {
		for (int k = 0; k < MOD; k++) std::cerr << value_coef[j][k] << " ";
		std::cerr << std::endl;
	}
	return std::vector<T>(n + 1);
}
template<typename T, typename val_func_t> std::vector<T> permanent_of_mod4_dependent(int n, const val_func_t &val_func) {
	constexpr int MOD = 4;
	int n0 = (n + 3) / MOD;
	int n1 = (n + 2) / MOD;
	int n2 = (n + 1) / MOD;
	std::vector<std::vector<std::vector<T> > > dp(n0 + 1, std::vector<std::vector<T> >(n1 + 1, std::vector<T>(n2 + 1)));
	dp[0][0][0] = 1;
	
	T cur_coef = 1;
	std::vector<T> res = { 1 };
	for (int i = 0; i < n; i++) {
		std::cerr << "stage " << i << std::endl;
		int mod = i % MOD;
		for (int j = std::min(n0, i); j >= 0; j--) for (int k = std::min(n1, i - j); k >= 0; k--) for (int l = std::min(n2, i - j - k); l >= 0; l--) {
			if (j < n0) dp[j + 1][k][l] += dp[j][k][l] * val_func(0, mod);
			if (k < n1) dp[j][k + 1][l] += dp[j][k][l] * val_func(1, mod);
			if (l < n2) dp[j][k][l + 1] += dp[j][k][l] * val_func(2, mod);
			dp[j][k][l] *= val_func(3, mod);
		}
		int cur_n = i + 1;
		auto cur_res = dp[(cur_n + 3) / MOD][(cur_n + 2) / MOD][(cur_n + 1) / MOD];
		cur_coef *= (i / MOD) + 1;
		res.push_back(cur_res * cur_coef);
	}
	return res;
}


