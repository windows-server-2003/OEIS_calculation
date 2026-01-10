#include <vector>
#include <cassert>

struct A202447 {
	template<typename val_t> std::vector<val_t> calc_2d_all(int n, int m_max) {
		assert(n >= 3);
		std::vector<val_t> res{0, 0, 0};
		
		if (n == 3) {
			auto ok = [] (int x) { return x != 0b001 && x != 0b110; };
			std::vector<std::vector<val_t> > dp(1 << n, std::vector<val_t>(1 << n));
			for (int i = 0; i < 1 << n; i++) if (ok(i)) for (int j = 0; j < 1 << n; j++) if (ok(j)) dp[i][j] = 1;
			for (int i = 2; i < m_max; i++) {
				std::vector<std::vector<val_t> > next(1 << n, std::vector<val_t>(1 << n));
				
				for (int j = 0; j < 1 << n; j++) if (ok(j)) for (int k = 0; k < 1 << n; k++) if (ok(k)) for (int l = 0; l < 1 << n; l++) if (ok(l)) {
					if ((~j & ~k & l) & ((1 << n) - 1)) continue;
					if ((j & k & ~l) & ((1 << n) - 1)) continue;
					if (~j & ~(k >> 1) & (l >> 2) & 1) continue;
					if (j & (k >> 1) & ~(l >> 2) & 1) continue;
					next[k][l] += dp[j][k];
				}
				std::swap(dp, next);
				
				val_t cur_res = 0;
				for (auto &j : dp) cur_res += std::accumulate(j.begin(), j.end(), (val_t) 0);
				res.push_back(cur_res);
			}
			return res;
		}
		
		std::vector<std::vector<std::array<val_t, 2> > > dp(n, std::vector<std::array<val_t, 2> >(n, {2, 2}));
		std::vector<std::vector<std::array<val_t, 2> > > next_acc2(n, std::vector<std::array<val_t, 2> >(n));
		for (int i = 2; i < m_max; i++) {
			for (int l = 0; l < n; l++) {
				next_acc2[l].assign(n, {});
				for (int k = 0; k < n; k++) {
					{
						const auto &cur_val = dp[k][l][0];
						auto go_range = [&] (int next_j, int l_l, int l_r) {
							assert(0 <= l_l && l_l <= l_r && l_r <= n);
							next_acc2[l][l_l][next_j] += cur_val;
							if (l_l + 1 < n) next_acc2[l][l_l + 1][next_j] += cur_val;
							if (l_r < n) next_acc2[l][l_r][next_j] -= cur_val;
							if (l_r + 1 < n) next_acc2[l][l_r + 1][next_j] -= cur_val;
						};
						auto go = [&] (int next_j, int next_l) {
							assert(0 <= next_l && next_l < n);
							next_acc2[l][next_l][next_j] += cur_val;
							if (next_l + 2 < n) next_acc2[l][next_l + 2][next_j] -= cur_val;
						};
						auto go_range2 = [&] (int next_j, int l_l, int l_r) {
							assert(0 <= l_l && l_l <= l_r && l_r <= n + 1);
							assert(!((l_l ^ l_r) & 1));
							next_acc2[l][l_l][next_j] += cur_val;
							if (l_r < n) next_acc2[l][l_r][next_j] -= cur_val;
						};
						if ((k ^ l) & 1) {
							if (k + 3 <= l) go_range2(0, k, k + 4);
							else if (k + 1 == l) go_range(0, k, n);
							else if (k == l + 1) {
								if (k == n - 1) go_range(0, l, k + 1);
								else go(0, l), go_range(0, k + 1, n);
							} else {
								go_range2(0, l, k + 1);
								if (k == n - 1) go(0, k);
								else go_range(0, k + 1, n);
							}
						} else {
							if (k < l) go_range2(0, k, k + 4);
							else go_range2(0, l, k + 2);
						}
					}
					{
						const auto &cur_val = dp[k][l][1];
						auto go_range = [&] (int next_j, int l_l, int l_r) {
							assert(0 <= l_l && l_l <= l_r && l_r <= n);
							next_acc2[l][l_l][next_j] += cur_val;
							if (l_l + 1 < n) next_acc2[l][l_l + 1][next_j] += cur_val;
							if (l_r < n) next_acc2[l][l_r][next_j] -= cur_val;
							if (l_r + 1 < n) next_acc2[l][l_r + 1][next_j] -= cur_val;
						};
						auto go = [&] (int next_j, int next_l) {
							assert(0 <= next_l && next_l < n);
							next_acc2[l][next_l][next_j] += cur_val;
							if (next_l + 2 < n) next_acc2[l][next_l + 2][next_j] -= cur_val;
						};
						auto go_range2 = [&] (int next_j, int l_l, int l_r) {
							assert(0 <= l_l && l_l <= l_r && l_r <= n + 1);
							assert(!((l_l ^ l_r) & 1));
							next_acc2[l][l_l][next_j] += cur_val;
							if (l_r < n) next_acc2[l][l_r][next_j] -= cur_val;
						};
						if (k < l) {
							if (k == 0) go(0, 1), go(1, 0);
						} else if (k == l) {
							if (k == 0) go_range(0, 0, n), go_range(1, 0, n);
							else if (k == 1) go(1, k - 1), go_range(1, k + 1, n), go(0, k);
							else if (k == n - 1) go(1, k);
							else go_range(1, k + 1, n);
						} else {
							if ((k ^ l) & 1) {
								if (l < 2) {
									go(0, l);
									go_range2(1, k & 1, k + 2);
								} else go_range2(1, l + 1, k + 2);
							} else {
								if (l < 2) {
									go(0, l);
									go_range2(1, (k - 1) & 1, k + 1);
								} else go_range2(1, l + 1, k + 1);
								if (k == n - 1) go(1, k);
								else go_range(1, k + 1, n);
							}
						}
					}
				}
				for (int k = 0; k + 2 < n; k++) {
					next_acc2[l][k + 2][0] += next_acc2[l][k][0];
					next_acc2[l][k + 2][1] += next_acc2[l][k][1];
				}
			}
			std::swap(dp, next_acc2);
			
			val_t cur_res = 0;
			for (auto &i : dp) for (auto &j : i) cur_res += j[0] + j[1];
			res.push_back(cur_res);
		}
		
		return res;
	}
	template<typename val_t> std::vector<val_t> calc_all(int n) {
		std::vector<std::vector<int> > num;
		int x = 0;
		int y = 0;
		for (int i = 0; i < n; i++) {
			if (x >= (int) num.size()) num.resize(x + 1);
			num[x].push_back(i + 1);
			
			if (y) x++, y--;
			else y = x + y + 1, x = 0;
		}
		
		std::vector<val_t> res(n + 1);
		for (int i = 0; i < (int) num.size(); i++) {
			auto cur_res = calc_2d_all<val_t>(3 + i, 3 + num[i].size());
			for (int j = 0; j < (int) num[i].size(); j++) {
				res[num[i][j]] = cur_res[3 + j];
			}
		}
		return res;
	}
};
struct A202439 {
	template<typename val_t> val_t calc(int n) {
		return A202447().calc_2d_all<val_t>(2 + n, 2 + n).back();
	}
	
};
