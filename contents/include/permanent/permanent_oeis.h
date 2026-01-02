#include "permanent.h"
#include "comb.h"
#include "prime.h"

// function returning one
#define CONSTANT(x) [] (int) { return x; }

// used when passing two lambda expression
#define PROD_WITH_CAPTURE(exp_of_i, exp_of_j, capture) [capture] (int i) { return exp_of_i; }, [capture] (int j) { return exp_of_j; }
#define PROD(exp_of_i, exp_of_j) PROD_WITH_CAPTURE(exp_of_i, exp_of_j, )


// A204233-

// A307783(n) : permanent of matrix a[i][j] = n - abs(i - j) :
struct A307783 {
	template<typename T> T calc(int n) {
		// one of the functions depends on n, so can't calc_all
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[n] (int i) { return n + i; }, CONSTANT(1), // n + i for upper
				[n] (int ij) { return n - ij; },
				[n] (int i) { return n - i; }, CONSTANT(1) // n - i for lower and diagonal
			},
			{
				CONSTANT(1), [] (int j) { return -j; }, // -j for upper
				[] (int ij) { return ij; },
				CONSTANT(1), [] (int j) { return j; } // j for lower and diagonal
			}
		).back();
	}
};

// A204235(n) : permanent of matrix a[i][j] = 1 + abs(i - j) :
struct A204235 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -i + 1; }, CONSTANT(1), // -i + 1 for upper
				CONSTANT(1),
				[] (int i) { return i + 1; }, CONSTANT(1) // i + 1 for lower
			},
			{
				CONSTANT(1), [] (int j) { return j; }, // j for upper
				CONSTANT(0),
				CONSTANT(1), [] (int j) { return -j; } // -j for lower
			}
		);
	}
};
// A085807(n) : permanent of matrix a[i][j] = abs(i - j)
struct A085807 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -i; }, CONSTANT(1), // -i for upper
				CONSTANT(0),
				[] (int i) { return i; }, CONSTANT(1) // i for lower
			},
			{
				CONSTANT(1), [] (int j) { return j; }, // j for upper
				CONSTANT(0),
				CONSTANT(1), [] (int j) { return -j; } // -j for lower
			}
		);
	}
};



// A204236(n) : permanent of matrix a[i][j] = max(2i - j, 2j - i) + 1 :
struct A204236 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -i + 1; }, CONSTANT(1), // -i + 1 for upper
				[] (int ij) { return 2*ij + 1; }, // 2i + 1 for lower and diagonal
				[] (int i) { return 2*i + 1; }, CONSTANT(1)
			},
			{
				CONSTANT(1), [] (int j) { return 2*j; }, // 2j for upper
				[] (int ij) { return -ij; }, // -j for lower and diagonal
				CONSTANT(1), [] (int j) { return -j; }
			}
		);
	}
};



// A204239(n) : permanent of matrix a[i][j] = max(3i - j, 3j - i) + 2 :
struct A204239 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -2*i + 2; }, CONSTANT(1), // -2i + 2 for upper
				[] (int ij) { return 3*ij + 2; }, // 3i + 2 for lower and diagonal
				[] (int i) { return 3*i + 2; }, CONSTANT(1)
			},
			{
				CONSTANT(1), [] (int j) { return 3*j; }, // 3j for upper
				[] (int ij) { return -2*ij; }, // -2j for lower and diagonal
				CONSTANT(1), [] (int j) { return -2*j; }
			}
		);
	}
};


// A204241(n) : permanent of matrix a[i][j] = max(3i - 2j, 3j - 2i) + 1 :
struct A204241 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -2*i + 1; }, CONSTANT(1), // -2i + 1 for upper
				[] (int ij) { return 3*ij + 1; }, // 3i + 1 for lower and diagonal
				[] (int i) { return 3*i + 1; }, CONSTANT(1)
			},
			{
				CONSTANT(1), [] (int j) { return 3*j; }, // 3j for upper
				[] (int ij) { return -2*ij; }, // -2j for lower and diagonal
				CONSTANT(1), [] (int j) { return -2*j; }
			}
		);
	}
};



// A085719(n) : permanent of matrix a[i][j] = 1 + (i + j) % n
//    equivalent to the case when a[i][j] = 1 + (i - j + n) % n (reversed columns)
struct A085719 {
	template<typename T> T calc(int n) {
		// one of the functions depends on n, so can't calc_all
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[n] (int i) { return 1 + i + n; }, CONSTANT(1), // 1 + i + n for upper
				[] (int ij) { return 1 + ij; },
				[] (int i) { return 1 + i; }, CONSTANT(1), // 1 + i for lower and diagonal
			},
			{
				CONSTANT(1), [] (int j) { return -j; }, // -j everywhere
				[] (int ij) { return -ij; },
				CONSTANT(1), [] (int j) { return -j; }
			}
		).back();
	}
};


// A086759(n) : (equivalent to) the permanent of matrix a[i][j] = (i - j + n) % n
struct A086759 {
	template<typename T> T calc(int n) {
		// one of the functions depends on n, so can't calc_all
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[n] (int i) { return i + n; }, CONSTANT(1), // i + n for upper
				[] (int ij) { return ij; },
				[] (int i) { return i; }, CONSTANT(1), // i for lower and diagonal
			},
			{
				CONSTANT(1), [] (int j) { return -j; }, // -j everywhere
				[] (int ij) { return -ij; },
				CONSTANT(1), [] (int j) { return -j; }
			}
		).back();
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
		// one of the functions depends on n, so can't calc_all
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return 1 - i; }, CONSTANT(1), // 1 - i for upper and diagonal
				[] (int ij) { return 1 - ij; },
				[n] (int i) { return n + i; }, CONSTANT(1), // n + i for lower
			},
			{
				CONSTANT(1), [] (int j) { return j; }, // j for upper and diagonal
				[] (int ij) { return ij; },
				CONSTANT(1), [] (int j) { return -j; } // -j for lower
			}
		).back();
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
		// one of the functions depends on n, so can't calc_all
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[n] (int i) { return 2 * n - 1 + i; }, CONSTANT(1), // 2n - 1 + i for upper and diagonal
				[n] (int ij) { return 2 * n - 1 + ij; },
				[n] (int i) { return n - i; }, CONSTANT(1), // n - i for lower
			},
			{
				CONSTANT(1), [] (int j) { return -j; }, // -j for upper and diagonal
				[] (int ij) { return -ij; },
				CONSTANT(1), [] (int j) { return j; } // j for lower
			}
		).back();
	}
};

// A204262(n): permanent of matrix a[i][j] = 1 + min(i, j)
struct A204262 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return i + 1; }, CONSTANT(1),
				[] (int ij) { return ij + 1; },
				CONSTANT(1), [] (int j) { return j + 1; }
			}
		);
	}
};

// A204264(n): permanent of matrix a[i][j] = 1 + max(i, j)
struct A204264 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_multilinear_per_triangular<T>(n, 
			{
				CONSTANT(1), [] (int j) { return j + 1; },
				[] (int ij) { return ij + 1; },
				[] (int i) { return i + 1; }, CONSTANT(1)
			}
		);
	}
};

// A204234(n): permanent of matrix a[i][j] = (i + j + 2) * (1 + min(i, j))
struct A204234 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return (T) (1 + i) * (1 + i); }, CONSTANT(1), // (1 + i)^2 for upper
				[] (int ij) { return (T) (1 + ij) * (1 + ij); },
				CONSTANT(1), [] (int j) { return (1 + j) * (1 + j); } // (1 + j)^2 for lower
			},
			{
				[] (int i) { return (T) (1 + i); }, // (1 + i)(1 + j) everywhere
				[] (int j) { return (T) (1 + j); },
				[] (int ij) { return (1 + ij) * (1 + ij); },
				[] (int i) { return (T) (1 + i); },
				[] (int j) { return (T) (1 + j); }
			}
		);
	}
};

// A278858(n): permanent of matrix a[i][j] = |(i+1)^2 - (j+1)^2|
struct A278858 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -(T) (1 + i) * (1 + i); }, CONSTANT(1), // - (1 + i)^2 for upper
				[] (int ij) { return -(T) (1 + ij) * (1 + ij); },
				[] (int i) { return (T) (1 + i) * (1 + i); }, CONSTANT(1), // (1 + i)^2 for lower
			},
			{
				CONSTANT(1), [] (int j) { return (T) (1 + j) * (1 + j); }, // (1 + j)^2 for upper
				[] (int ij) { return (T) (1 + ij) * (1 + ij); },
				CONSTANT(1), [] (int j) { return -(T) (1 + j) * (1 + j); }, // -(1 + j)^2 for upper
			}
		);
	}
};

// A347768(n): permanent of matrix a[i][j] = |i - j + 1|
struct A347768 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_sum_of_multilinear_per_triangular<T>(n, 
			{
				[] (int i) { return -i - 1; }, CONSTANT(1), // -i - 1 for upper
				[] (int ij) { return ij + 1; }, // i + 1 for diagonal and lower
				[] (int i) { return i + 1; }, CONSTANT(1),
			},
			{
				CONSTANT(1), [] (int j) { return j; }, // j for upper
				[] (int ij) { return -ij; }, // -j for diagonal and lower
				CONSTANT(1), [] (int j) { return -j; }, // -(1 + j)^2 for upper
			}
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
struct A204546 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 == 0 || (i + j) % 4 == 3; });
	}
};
struct A204548 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 >= 2; });
	}
};
struct A204550 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanent_of_mod4_dependent<T>(n, [] (int i, int j) { return (i + j) % 4 == 1 || (i + j) % 4 == 2; });
	}
};









// ----- k = 1 -----
// A203264(n) : permanent of matrix a[i][j] = (i / 2 + 1) * (j + 1)
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
// A330087(n) : permanent of matrix a[i][j] = (i+1) * prime(j)
struct A330087 {
	template<typename T> std::vector<T> calc_all(int n) {
		auto primes = get_first_n_primes<int>(n);
		std::vector<T> res = { 1 };
		T factorial = 1;
		T prod_half = 1;
		for (int i = 1; i <= n; i++) {
			factorial *= i;
			prod_half *= primes[i - 1];
			res.push_back(factorial * factorial * prod_half);
		}
		return res;
	}
};



// ----- k = 2 -----
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

// A204248(n) : permanent of matrix a[i][j] = 1 + i + j :
struct A204248 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {1, 1 + i}; },
			[] (int j) { return j; }
		);
	}
};

// A204249(n) : permanent of matrix a[i][j] = 2 + i + j :
struct A204249 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {1, 2 + i}; },
			[] (int j) { return j; }
		);
	}
};


// A204251(n) : permanent of matrix a[i][j] = i * j + 2i + 2j + 1 :
struct A204251 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {(2 + i), 2*i + 1}; },
			[] (int j) { return j; }
		);
	}
};
/* A278847(n) : permanent of matrix a[i][j] = (i + 1)^2 + (j + 1)^2 */
struct A278847 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, int> { return {1, (i + 1) * (i + 1)}; },
			[] (int j) { return (j + 1) * (j + 1); }
		);
	}
};
/* A278925(n) : permanent of matrix a[i][j] = (i + 1)^3 + (j + 1)^3 */
struct A278925 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, s64> { return {1, (s64) (i + 1) * (i + 1) * (i + 1)}; },
			[] (int j) { return (s64) (j + 1) * (j + 1) * (j + 1); }
		);
	}
};
/* A278926(n) : permanent of matrix a[i][j] = (i + 1)^4 + (j + 1)^4 */
struct A278926 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_linear_n_independent_matrix<T>(n,
			[] (int i) -> std::pair<int, s64> { return {1, (s64) (i + 1) * (i + 1) * (i + 1) * (i + 1)}; },
			[] (int j) { return (s64) (j + 1) * (j + 1) * (j + 1) * (j + 1); }
		);
	}
};

/*
	A322277(n) : permanent of matrix a[i][j] = (i * n + j + 1) if i % 2 == 0, ((i + 1) * n - j) if i % 2 == 1
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





// ----- k = 3 -----
/* A278845(n) : permanent of 2nx2n matrix a[i][j] = (i + j + 2)^2 */
struct A278845 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_quadratic_n_independent_matrix<T>(n,
			[] (int i) { return std::make_tuple((i + 1) * (i + 1), 2 * (i + 1), 1); },
			[] (int j) { return j + 1; }
		);
	}
};
/* A278857(n) : permanent of 2nx2n matrix a[i][j] = (i - j)^2 */
struct A278857 {
	template<typename T> std::vector<T> calc_all(int n) {
		return permanents_of_rowwise_quadratic_n_independent_matrix<T>(n,
			[] (int i) { return std::make_tuple((i + 1) * (i + 1), -2 * (i + 1), 1); },
			[] (int j) { return j + 1; }
		);
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


// A204431 : wrong 

/* NOPE : A204261, A204266 */
/* FASTER : A204262, A204264 */

/* TODO(more) : A330287 A332566 A278858 A349108 A349107 */
