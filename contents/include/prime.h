#pragma once
#include <vector>

template<typename T> std::vector<bool> get_is_prime(T max) {
	std::vector<bool> is_prime(max + 1, true);
	is_prime[0] = is_prime[1] = false;
	for (T i = 2; i <= max; i++) if (is_prime[i]) 
		for (T j = i + i; j <= max; j += i) is_prime[j] = false;
	return is_prime;
}
std::vector<bool> get_is_prime_until_nth_prime(size_t n) {
	size_t n_max6 = std::max<size_t>(6, n);
	size_t max = n_max6 * (std::log(n_max6) + std::log(std::log(n_max6)));
	
	return get_is_prime(max);
}

template<typename T> std::vector<T> get_first_n_primes(size_t n) {
	auto is_prime = get_is_prime_until_nth_prime(n);
	std::vector<T> res;
	for (size_t i = 0; i < is_prime.size(); i++) if (is_prime[i]) res.push_back(i);
	assert(res.size() >= n);
	res.resize(n);
	return res;
}
