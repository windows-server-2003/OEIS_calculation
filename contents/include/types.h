#pragma once
#include <cstdint>
#include <cinttypes>
#include <boost/multiprecision/cpp_int.hpp>

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using u128 = boost::multiprecision::uint128_t;
using u256 = boost::multiprecision::uint256_t;
using u512 = boost::multiprecision::uint512_t;
using u1024 = boost::multiprecision::uint1024_t;
using cu128 = boost::multiprecision::checked_uint128_t;
using cu256 = boost::multiprecision::checked_uint256_t;
using cu512 = boost::multiprecision::checked_uint512_t;
using cu1024 = boost::multiprecision::checked_uint1024_t;

using s8 = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;
using s128 = boost::multiprecision::int128_t;
using s256 = boost::multiprecision::int256_t;
using s512 = boost::multiprecision::int512_t;
using s1024 = boost::multiprecision::int1024_t;
using cs128 = boost::multiprecision::checked_int128_t;
using cs256 = boost::multiprecision::checked_int256_t;
using cs512 = boost::multiprecision::checked_int512_t;
using cs1024 = boost::multiprecision::checked_int1024_t;

using cpp_int = boost::multiprecision::cpp_int;

