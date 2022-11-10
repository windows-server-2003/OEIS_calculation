#pragma once
#include <vector>
#include "types.h"

template<typename T, typename S, typename U> void type_stepped_calculate(int stage, const T &init, const S &transit_if_needed, const U &main_func) {
	enum Stage {
		U64,
		U128,
		U256,
		U512,
		U1024,
	};
	
	std::vector<u64> dp_u64;
	std::vector<u128> dp_u128;
	std::vector<u256> dp_u256;
	std::vector<u512> dp_u512;
	std::vector<u1024> dp_u1024;
	
	init(dp_u64);
	
	Stage type_stage = U64;
	for (int i = 0; i < stage; i++) {
		if (type_stage == U64 && transit_if_needed(dp_u64, dp_u128)) type_stage = U128;
		if (type_stage == U128 && transit_if_needed(dp_u128, dp_u256)) type_stage = U256;
		if (type_stage == U256 && transit_if_needed(dp_u256, dp_u512)) type_stage = U512;
		if (type_stage == U512 && transit_if_needed(dp_u512, dp_u1024)) type_stage = U1024;
		
		if (type_stage == U64)       main_func(i, dp_u64);
		else if (type_stage == U128) main_func(i, dp_u128);
		else if (type_stage == U256) main_func(i, dp_u256);
		else if (type_stage == U512) main_func(i, dp_u512);
		else                         main_func(i, dp_u1024);
	}
}

