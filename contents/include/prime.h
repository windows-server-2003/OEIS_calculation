#pragma once
#include <vector>

template<typename T> std::vector<int> get_is_prime(T max) {
	std::vector<bool> is_prime(max + 1, true);
	is_prime[0] = is_prime[1] = false;
	for (T i = 2; i <= max; i++) if (is_prime[i]) 
		for (T j = i + i; j <= max; j += i) is_prime[j] = false;
	return is_prime;
}
template<typename T> std::vector<bool> get_is_prime_until_nth_prime(T n) {
	T n_max6 = std::max(6, n);
	T max = n_max6 * (std::log(n_max6) + std::log(std::log(n_max6)));
	
	
	std::vector<bool> is_prime(max + 1, true);
	is_prime[0] = is_prime[1] = false;
	for (T i = 2; i <= max; i++) if (is_prime[i]) 
		for (T j = i + i; j <= max; j += i) is_prime[j] = false;
	return is_prime;
}
