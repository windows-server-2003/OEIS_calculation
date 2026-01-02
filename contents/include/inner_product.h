#include <vector>
#include <cassert>

// Inner product of (first `len` elements of, if specified) a and b, calculated in type T
template<typename T, typename A, typename B> T inner_product(const std::vector<A> &a, const std::vector<B> &b, size_t len = -1) {
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size());
	T res = 0;
	for (size_t i = 0; i < len; i++) res += (T) a[i] * b[i];
	return res;
}
// sum_i a[i]b[i]c[i], calculated in type T
template<typename T, typename A, typename B, typename C> T inner_product(
	const std::vector<A> &a, const std::vector<B> &b, const std::vector<C> &c, size_t len = -1) {
	
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size() && b.size() == c.size());
	T res = 0;
	for (size_t i = 0; i < len; i++) res += (T) a[i] * b[i] * c[i];
	return res;
}

// returns {dot(a, b), dot(b, b)} (reduced memory access compared to two calls on inner_product)
template<typename T, typename A, typename B> std::pair<T, T> inner_product_ab_bb(
	const std::vector<A> &a, const std::vector<B> &b, size_t len = -1) {
	
	if (len == (size_t) -1) len = a.size(); 
	assert(a.size() == b.size());
	T res_ab = 0;
	T res_bb = 0;
	for (size_t i = 0; i < len; i++) {
		res_ab += (T) a[i] * b[i];
		res_bb += (T) b[i] * b[i];
	}
	return {res_ab, res_bb};
}
