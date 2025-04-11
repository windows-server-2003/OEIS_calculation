#pragma once
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <iostream>
#include "types.h"

std::string bin(u64 x, int n) {
	std::string res;
	for (int i = 0; i < n; i++) res.push_back((x >> i & 1) + '0');
	return res;
}


template<typename A, typename B> std::ostream & operator << (std::ostream &stream, const std::pair<A, B> &a) {
	stream << "(" << a.first << "," << a.second << ")";
	return stream;
}
template<typename T> void write_(std::ostream &stream, T begin, T end) {
	for (auto i = begin; i != end; i++) stream << (i != begin ? " " : "") << *i;
}
template<typename T> std::ostream& operator << (std::ostream &stream, const std::vector<T> &a) {
	stream << "[";
	write_(stream, a.begin(), a.end());
	stream << "]";
	return stream;
}
template<typename T> std::ostream& operator << (std::ostream &stream, const std::set<T> &a) {
	stream << "{";
	write_(stream, a.begin(), a.end());
	stream << "}";
	return stream;
}
template<typename T> std::ostream& operator << (std::ostream &stream, const std::multiset<T> &a) {
	stream << "{";
	write_(stream, a.begin(), a.end());
	stream << "}";
	return stream;
}
