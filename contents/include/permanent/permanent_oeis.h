#include "permanent.h"
#include "comb.h"

// A204233-

/*
	a(n) : permanent of matrix a[i][j] = n - abs(i - j) :
		[n    n-1  n-2  n-3  ...  1]
		[n-1  n    n-1  n-2  ...  2]
		[n-2  n-1  n    n-1  ...  3]
		[n-3  n-2  n-1  n    ...  4]
		[|    |    |    | \       |]
		[|    |    |    |    \    |]
		[|    |    |    |       \ |]
		[1    2    3    4    ...  n]
*/
struct A307783 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { return col_cnt * (n - i) + col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { return col_cnt * (n + i) - col_sum; }
		);
	}
};


/*
	a(n) : permanent of matrix a[i][j] = 1 + abs(i - j) :
		[1    2    3    4    ...  n  ]
		[2    1    2    3    ...  n-1]
		[3    2    1    2    ...  n-2]
		[4    3    2    1    ...  n-3]
		[|    |    |    | \       |  ]
		[|    |    |    |    \    |  ]
		[|    |    |    |       \ |  ]
		[n    n-1  n-2  n-3  ...  1  ]
*/
struct A204235 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (i + 1) - col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (1 - i) + col_sum; }
		);
	}
};
/*
	a(n) : permanent of matrix a[i][j] = abs(i - j)
*/
struct A085807 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * i  - col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * -i + col_sum; }
		);
	}
};



/*
	a(n) : permanent of matrix a[i][j] = max(2i - j, 2j - i) + 1 :
		[1    3    5    7    ...  2n-1]
		[3    2    4    6    ...  2n-2]
		[5    4    3    5    ...  2n-3]
		[7    6    5    4    ...  2n-4]
		[|    |    |    | \       |   ]
		[|    |    |    |    \    |   ]
		[|    |    |    |       \ |   ]
		[2n-1 2n-2 2n-3 2n-4  ... n   ]
*/
struct A204236 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (2 * i + 1) - col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (1 - i) + 2 * col_sum; }
		);
	}
};



/*
	a(n) : permanent of matrix a[i][j] = max(3i - j, 3j - i) + 2 :
		[2    5    8    11   ...  3n-1]
		[5    4    7    10   ...  3n-2]
		[8    7    6    9    ...  3n-3]
		[11   10   9    8    ...  3n-4]
		[|    |    |    | \       |   ]
		[|    |    |    |    \    |   ]
		[|    |    |    |       \ |   ]
		[3n-1 3n-2 3n-3 3n-4  ... 2n  ]
*/
struct A204239 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (3 * i + 2) - col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (2 - i) + 3 * col_sum; }
		);
	}
};


/*
	a(n) : permanent of matrix a[i][j] = max(3i - 2j, 3j - 2i) + 1 :
		[1    4    7    10   ...  3n-2]
		[4    2    5    8    ...  3n-4]
		[7    5    3    6    ...  3n-6]
		[10   8    6    4    ...  3n-8]
		[|    |    |    | \       |   ]
		[|    |    |    |    \    |   ]
		[|    |    |    |       \ |   ]
		[3n-2 3n-4 3n-6 3n-8  ... n   ]
*/
struct A204241 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (3 * i + 1) - 2 * col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (1 - 2 * i) + 3 * col_sum; }
		);
	}
};


/*
	a(n) : permanent of matrix a[i][j] = 2 + i + j :
		[2    3    4    5    ...  n+1]
		[3    4    5    6    ...  n+2]
		[4    5    6    7    ...  n+3]
		[5    6    7    8    ...  n+4]
		[|    |    |    | \       |  ]
		[|    |    |    |    \    |  ]
		[|    |    |    |       \ |  ]
		[n+1  n+2  n+3  n+4   ... 2n ]
*/
struct A204249 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {1, 2 + i}; },
			[] (int j) { return j; }
		);
	}
};


/*
	a(n) : permanent of matrix a[i][j] = i * j + 2i + 2j + 1 :
		[1    3    5    7    ...  2n-1]
		[3    6    9    12   ...  3n  ]
		[5    9    13   17   ...  4n+1]
		[7    12   17   22   ...  5n+2]
		[|    |    |    | \       |   ]
		[|    |    |    |    \    |   ]
		[|    |    |    |       \ |   ]
		[2n-1 3n   4n+1 5n+2  ... n^2+2n-2]
*/
struct A204251 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {(2 + i), 2*i + 1}; },
			[] (int j) { return j; }
		);
	}
};

/*
	a(n) : permanent of matrix a[i][j] = 1 + (i - j + n) % n
		[1    2    3    4    ...  n  ]
		[n    1    2    3    ...  n-1]
		[n-1  n    1    2    ...  n-2]
		[n-2  n-1  n    1    ...  n-3]
		[|    |    |    | \       |  ]
		[|    |    |    |    \    |  ]
		[|    |    |    |       \ |  ]
		[2    3    4    5    ...  1  ]
*/
struct A085719 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) {           return col_cnt * (n + 1 - i) + col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (1 - i) + col_sum; }
		);
	}
};

/*
	a(n) : (equivalent to) the permanent of matrix a[i][j] = (i - j + n) % n
		[0    1    2    3    ...  n-1]
		[n-1  0    1    2    ...  n-2]
		[n-2  n-1  0    1    ...  n-3]
		[n-3  n-2  n-1  0    ...  n-4]
		[|    |    |    | \       |  ]
		[|    |    |    |    \    |  ]
		[|    |    |    |       \ |  ]
		[1    2    3    4    ...  0  ]
*/
struct A086759 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) {           return col_cnt * (n - i) + col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * -i + col_sum; }
		);
	}
};



/*
	a(n) : permanent of matrix a[i][j] = 1 + (i + j) % 2 :
		[1    2    1    2   ...  ]
		[2    1    2    1   ...  ]
		[1    2    1    2   ...  ]
		[2    1    2    1   ...  ]
		[|    |    |    | \      ]
		[|    |    |    |    \   ]
		[|    |    |    |       \]
*/
struct A204252 {
	template<typename T> std::vector<T> calc_all(int n) {
		constexpr int MOD = 2;
		int n0 = (n + 1) / MOD;
		std::vector<T> dp(n0 + 1);
		dp[0] = 1;
		
		T cur_coef = 1;
		std::vector<T> res = { 1 };
		for (int i = 0; i < n; i++) {
			int mod = i % MOD;
			for (int j = n0; j >= 0; j--) {
				if (j < n0) dp[j + 1] += dp[j] * (1 + (0 + mod + 2) % MOD);
				dp[j] *= (1 + (1 + mod + 2) % MOD);
			}
			int cur_n = i + 1;
			auto cur_res = dp[(cur_n + 1) / MOD];
			cur_coef *= (i / MOD) + 1;
			res.push_back(cur_res * cur_coef);
		}
		return res;
	}
};


/*
	Mod-dependent matrices
*/

struct A179079 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (i + j + 1) % 3; });
	}
};
struct A204254 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return 1 + (i + j + 2) % 3; });
	}
};

struct A204256 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return 1 + (i + j + 2) % 4; });
	}
};
struct A204258 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return 1 + (i + 2 * j) % 3; });
	}
};
struct A204265 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (i + j + 2) % 3; });
	}
};
struct A204268 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (i + j) % 3; });
	}
};
/*
	a(n) : permanent of matrix a[i][j] = ((i + j) % 4) < 2 :
		[1 1 0 0 ..]
		[1 0 0 1 ..]
		[0 0 1 1 ..]
		[0 1 1 0 ..]
		[| | | | \ ]
*/
struct A204422 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 < 2; });
	}
};

struct A204424 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (2 * i + j) % 3; });
	}
};
struct A204426 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (2 * i + j + 1) % 3; });
	}
};
struct A204428 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (2 * i + j + 2) % 3; });
	}
};
struct A204430 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (2 * i + 2 * j + 1) % 3; });
	}
};
struct A204432 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (2 * i + 2 * j + 2) % 3; });
	}
};
struct A204434 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (2 * i + 2 * j) % 3; });
	}
};
struct A204436 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (i + j) % 3 != 1; });
	}
};
struct A204438 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (i + j) % 3 != 0; });
	}
};
struct A204440 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod3_dependent<T>(n, [] (int i, int j) { return (i + j) % 3 != 2; });
	}
};
struct A204442 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 != 3; });
	}
};
struct A204444 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 != 2; });
	}
};
struct A204446 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 != 1; });
	}
};
struct A204448 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 != 0; });
	}
};


/*
	a(n) : permanent of matrix a[i][j] = (i * n + j + 1) if i % 2 == 0, ((i + 1) * n - j) if i % 2 == 1
		[1    2    3    4    ...  n   ]
		[2n   2n-1 2n-2 2n-3 ...  n+1 ]
		[2n+1 2n+2 2n+3 2n+4 ...  3n  ]
		[4n   4n-1 4n-2 4n-3 ...  3n+1]
		[|    |    |    |   \     |     ]
		[|    |    |    |     \   |     ]
		[|    |    |    |       \ |     ]
*/
struct A322277 {
	template<typename T> T calc(int n) {
		return permanent_of_rowwise_linear_matrix<T>(n,
			[] (int n, int i) -> std::pair<int, int> {
				if (i & 1) return {-1, (i + 1) * n};
				else return {1, i * n + 1};
			},
			[] (int, int j) { return j; }
		);
	}
};
/* a(n) : permanent of matrix a[i][j] = (i + 1)^2 + (j + 1)^2 */
struct A278847 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {1, (i + 1) * (i + 1)}; },
			[] (int j) { return (j + 1) * (j + 1); }
		);
	}
};
/* a(n) : permanent of matrix a[i][j] = (i + 1)^3 + (j + 1)^3 */
struct A278925 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, s64> { return {1, (s64) (i + 1) * (i + 1) * (i + 1)}; },
			[] (int j) { return (s64) (j + 1) * (j + 1) * (j + 1); }
		);
	}
};
/* a(n) : permanent of matrix a[i][j] = (i + 1)^4 + (j + 1)^4 */
struct A278926 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, s64> { return {1, (s64) (i + 1) * (i + 1) * (i + 1) * (i + 1)}; },
			[] (int j) { return (s64) (j + 1) * (j + 1) * (j + 1) * (j + 1); }
		);
	}
};
/* a(n) : permanent of 2nx2n matrix a[i][j] = i - j */
struct A346934 {
	template<typename T> std::vector<T> calc_all(int n) {
		auto tmp = permanents_of_rowwise_linear_n_independent_matrix<T>(2 * n,
			[] (int i) -> std::pair<int, int> { return {-1, i}; },
			[] (int j) { return j; }
		);
		std::vector<T> res;
		for (int i = 0; i <= 2 * n; i += 2) res.push_back(tmp[i]);
		return res;
	}
};
/* a(n) : permanent of 2nx2n matrix a[i][j] = (i + j)^2 */
struct A278845 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_quadratic_n_independent_matrix<T>(n,
			[] (int i) { return std::make_tuple((i + 1) * (i + 1), 2 * (i + 1), 1); },
			[] (int j) { return j + 1; }
		);
	}
};
/* a(n) : permanent of 2nx2n matrix a[i][j] = (i - j)^2 */
struct A278857 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_quadratic_n_independent_matrix<T>(n,
			[] (int i) { return std::make_tuple((i + 1) * (i + 1), -2 * (i + 1), 1); },
			[] (int j) { return j + 1; }
		);
	}
};


/*
	a(n) : permanent of matrix a[i][j] = 2*i + j + 3
*/
struct A278927 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {1, 2 * i + 3}; },
			[] (int j) { return j; }
		);
	}
};

/*
	a(n) : permanent of matrix a[i][j] = (i / 2 + 1) * (j + 1)
*/
struct A203264 {
	template<typename T> std::vector<T> calc_all(int n) {
		std::vector<T> res = { 1 };
		T factorial = 1;
		T prod_half = 1;
		for (int i = 1; i <= n; i++) {
			factorial *= i;
			prod_half *= (i + 1) / 2;
			res.push_back(factorial * factorial * prod_half);
		}
		return res;
	}
};


// Diagonal/1 matrix

/*
	a(n) : permanent of matrix a[i][i] = (i + 1)^2, a[i][j] = 1 (i != j)
*/
struct A303000 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_diagonal_plus1_matrix<cpp_int>(n, [] (int i) { return (i + 1) * (i + 1); });
	}
};

/*
	a(n) : permanent of matrix a[i][i] = (i + 1) * (i + 2) / 2, a[i][j] = 1 (i != j)
*/
struct A303001 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_diagonal_plus1_matrix<cpp_int>(n, [] (int i) { return (i + 1) * (i + 2) / 2; });
	}
};


/*
	a(n) : permanent of matrix a[i][j] = (1 + j - i) if j >= i, (n + i - j) if j < i  :
		[1    2    3    4  ...  n  ]
		[n+1  1    2    3  ...  n-1]
		[n+2  n+1  1    2  ...  n-2]
		[n+3  n+2  n+1  1    ...n-3]
		[|    |    |    | \     |  ]
		[|    |    |    |   \   |  ]
		[|    |    |    |     \ |  ]
		[2n-1 2n-2 2n-3 2n-4    1  ]
*/
struct A322909 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) {           return col_cnt * (n + i) - col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { (void) n; return col_cnt * (1 - i) + col_sum; }
		);
	}
};
/*
	a(n) : permanent of (transpose of) the matrix a[i][j] = (2 * n - 1 - (j - i)) if j >= i, (n - (i - j)) if j < i  :
		[2n-1 2n-2 2n-3 2n-4 ... n  ]
		[n-1  2n-1 2n-2 2n-3 ... n+1]
		[n-2  n-1  2n-1 2n-2 ... n+2]
		[n-3  n-2  n-1  2n-1 ... n+3]
		[|    |    |    | \      |  ]
		[|    |    |    |   \    |  ]
		[|    |    |    |     \  |  ]
		[1    2    3    4       2n-1]
*/
struct A323255 {
	template<typename T> T calc(int n) {
		return permanent_of_half_linear_matrix<T>(n,
			[] (int n, int i, int col_cnt, int col_sum) { return col_cnt * (n - i) + col_sum; },
			[] (int n, int i, int col_cnt, int col_sum) { return col_cnt * (2 * n - 1 + i) - col_sum; }
		);
	}
};



// A204431 : wrong 

/* NOPE : A204261, A204266 */
/* FASTER : A204262, A204264 */

/* TODO(more) : A330287 A332566 A278858 A349108 A349107 */
