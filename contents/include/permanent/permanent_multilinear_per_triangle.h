#include <vector>
#include <functional>
#include <utility>

/*
	Represents a multilinear matrix per triangular part.
	 - m[i][j] = upper_row_func(i) * upper_col_func(j) (j > i)
	 - m[i][i] = diagonal_func(i)
	 - m[i][j] = lower_row_func(i) * lower_col_func(j) (j < i)
*/
template<typename T> struct triangular_outprod {
	std::function<T (int)> upper_row_func;
	std::function<T (int)> upper_col_func;
	std::function<T (int)> diagonal_func;
	std::function<T (int)> lower_row_func;
	std::function<T (int)> lower_col_func;
};

/*
	Returns the permanents of all the principal submatrices of nxn matrix m that is multilinear per triangular part.
	See comments on triangular_outprod for how m should be given.
	
	return[i] (1 <= i <= n) is the permanent for the top-left ixi submatrix of m
*/
template<typename T> std::vector<T> permanents_of_multilinear_per_triangular(
	int n,
	triangular_outprod<T> m
) {
	/*
		dp[j]: number of ways to decide
		 - which of rows 0...i-1 and columns 0...i-1 will be used for upper element
			Let us call rows that we decided to use for upper element "upper rows".
			Similarly, we define "lower rows", "upper columns", "lower columns".
		 - for each "upper column", which "upper row" to pair with it (the row must be also in 0...i-1 because the element at the position must be in the upper triangle)
		 - for each "lower row", which "lower column" to pair with it (same)
		so that
		 - there are j more "upper rows" than "upper columns" (j can't be negative)
		 - there are j more "lower columns" than "lower rows" (this excess is always the same as the above excess)
		but each choice is weighed by the product of coefficients for row/column 0...i-1 (based on above decision)
	*/
	std::vector<T> dp = { 1 };
	std::vector<T> res = { 1 };
	
	for (int i = 0; i < n; i++) {
		std::vector<T> dp_next(i + 2);
		
		T diag = m.diagonal_func(i);
		T upper_row = m.upper_row_func(i);
		T upper_col = m.upper_col_func(i);
		T lower_row = m.lower_row_func(i);
		T lower_col = m.lower_col_func(i);
		for (int j = 0; j < (int) dp.size(); j++) if (dp[j] != 0) {
			dp_next[j] += dp[j] * diag; // use row i and column i for diagonal element
			dp_next[j] += dp[j] * upper_row * upper_col * j; // use row for upper element, col for upper
			dp_next[j + 1] += dp[j] * upper_row * lower_col; // use row for upper element, col for lower
			if (j) dp_next[j - 1] += dp[j] * lower_row * upper_col * j * j; // use (row, col) = (lower, upper)
			dp_next[j] += dp[j] * lower_row * lower_col * j; // use (row, col) = (lower, lower)
		}
		
		std::swap(dp, dp_next);
		
		res.push_back(dp[0]);
	}
	
	return res;
}


/*
	Returns the permanents of all the principal submatrices of nxn matrix m which can be represented as:
		m = m0 + m1
	where m0 and m1 are multilinear per triangular part.
	See comments on triangular_outprod for how m0 and m1 should be given.
	
	return[i] (1 <= i <= n) is the permanent for the top-left ixi submatrix of m
*/
template<typename T> std::vector<T> permanents_of_sum_of_multilinear_per_triangular(
	int n,
	triangular_outprod<T> m0,
	triangular_outprod<T> m1
) {
	/*
		dp[m0_upper_excess][m1_upper_excess][m0_lower_excess]:
			similar to permanents_of_multilinear_per_triangular version but
			 - m0/m1_upper_excess refers to row's excess over columns regarding those which use m0/m1's upper component
			 - m0_lower_excess refers to column's excess over rows regarding those which use m0's lower component
			 - "m1_lower_excess" is always (m0_upper_excess + m1_upper_excess - m0_lower_excess)
	*/
	std::vector<std::vector<std::vector<T> > > dp = { { { 1 } } };
	std::vector<T> res = { 1 };
	
	for (int i = 0; i < n; i++) {
		std::vector<std::vector<std::vector<T> > > dp_next(i + 2);
		for (int j = 0; j <= i + 1; j++) {
			// m0_upper_excess + m1_upper_excess <= i + 1
			// m0_lower_excess <= i + 1
			dp_next[j].resize(i + 2 - j, std::vector<T>(i + 2));
		}
		
		T m0_diag = m0.diagonal_func(i);
		T m0_upper_row = m0.upper_row_func(i);
		T m0_upper_col = m0.upper_col_func(i);
		T m0_lower_row = m0.lower_row_func(i);
		T m0_lower_col = m0.lower_col_func(i);
		T m1_diag = m1.diagonal_func(i);
		T m1_upper_row = m1.upper_row_func(i);
		T m1_upper_col = m1.upper_col_func(i);
		T m1_lower_row = m1.lower_row_func(i);
		T m1_lower_col = m1.lower_col_func(i);
		for (int m0_upper_excess = 0; m0_upper_excess < (int) dp.size(); m0_upper_excess++) {
			for (int m1_upper_excess = 0; m1_upper_excess < (int) dp[m0_upper_excess].size(); m1_upper_excess++) {
				for (int m0_lower_excess = 0; m0_lower_excess < (int) dp[m0_upper_excess][m1_upper_excess].size();
					m0_lower_excess++) {
					
					int m1_lower_excess = m0_upper_excess + m1_upper_excess - m0_lower_excess;
					
					// ignore if it is impossible to reach dp[0][0] in the following n - i updates
					if (std::max({m0_upper_excess + m1_upper_excess, m0_lower_excess + m1_lower_excess}) > n - i) break;
					
					auto &cur_val = dp[m0_upper_excess][m1_upper_excess][m0_lower_excess];
					if (cur_val == 0) continue;
					
					
					auto add = [&] (int m0_upper_excess_offset, int m1_upper_excess_offset, int m0_lower_excess_offset, T val) {
						if (m0_upper_excess + m0_upper_excess_offset < 0) return;
						if (m1_upper_excess + m1_upper_excess_offset < 0) return;
						if (m0_lower_excess + m0_lower_excess_offset < 0) return;
						dp_next[m0_upper_excess + m0_upper_excess_offset][m1_upper_excess + m1_upper_excess_offset]
							[m0_lower_excess + m0_lower_excess_offset] += val;
					};
					
					// use diagonal
					add(0, 0, 0, cur_val * (m0_diag + m1_diag));
					
					// (row, col) = (m0_upper, *)
					add(0, 0, 0, cur_val * m0_upper_row * m0_upper_col * m0_upper_excess);
					add(1, -1, 0, cur_val * m0_upper_row * m1_upper_col * m1_upper_excess);
					add(1, 0, 1, cur_val * m0_upper_row * m0_lower_col);
					add(1, 0, 0, cur_val * m0_upper_row * m1_lower_col);
					// (row, col) = (m1_upper, *)
					add(-1, 1, 0, cur_val * m1_upper_row * m0_upper_col * m0_upper_excess);
					add(0, 0, 0, cur_val * m1_upper_row * m1_upper_col * m1_upper_excess);
					add(0, 1, 1, cur_val * m1_upper_row * m0_lower_col);
					add(0, 1, 0, cur_val * m1_upper_row * m1_lower_col);
					
					// (row, col) = (m0_lower, *)
					add(-1, 0, -1, cur_val * m0_lower_row * m0_upper_col * m0_lower_excess * m0_upper_excess);
					add(0, -1, -1, cur_val * m0_lower_row * m1_upper_col * m0_lower_excess * m1_upper_excess);
					add(0, 0, 0, cur_val * m0_lower_row * m0_lower_col * m0_lower_excess);
					add(0, 0, -1, cur_val * m0_lower_row * m1_lower_col * m0_lower_excess);
					// (row, col) = (m1_lower, *)
					add(-1, 0, 0, cur_val * m1_lower_row * m0_upper_col * m1_lower_excess * m0_upper_excess);
					add(0, -1, 0, cur_val * m1_lower_row * m1_upper_col * m1_lower_excess * m1_upper_excess);
					add(0, 0, 1, cur_val * m1_lower_row * m0_lower_col * m1_lower_excess);
					add(0, 0, 0, cur_val * m1_lower_row * m1_lower_col * m1_lower_excess);
				}
			}
		}
		std::swap(dp, dp_next);
		res.push_back(dp[0][0][0]);
	}
	
	return res;
}

