#pragma once
#include <vector>

// using the Euler's pentagonal number theorem
template<typename T> std::vector<T> get_partition_sequence(int n) {
	std::vector<T> res(n + 1);
	res[0] = 1;
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j * (3 * j - 1) / 2 <= i; j++) {
			if (j & 1) {
				res[i] += res[i - j * (3 * j - 1) / 2];
				if (-j * (3 * -j - 1) / 2 <= i) res[i] += res[i - -j * (3 * -j - 1) / 2];
			} else {
				res[i] -= res[i - j * (3 * j - 1) / 2];
				if (-j * (3 * -j - 1) / 2 <= i) res[i] -= res[i - -j * (3 * -j - 1) / 2];
			}
		}
	}
	return res;
}
