#include "types.h"


/*
	* This is an internal function *
	Let f be a concave function.
	is_below_(s64 x, s64 y) must return whether (x, y) is below the graph of f
		guaranteed:
		 - x in [l, r]
	is_above_tangent(s64 x, s64 dx, s64 dy) must return whether (dx, dy) points more upward than the tangent vector of f at x
		guaranteed:
		 - x in [l, r]
		 - dx > 0
	
	This function returns sum_{x = l...r} (the largest integer y satisfying is_below_(x, y)).
*/
template<typename T, typename S> s128 curve_lattice_point_(s64 l, s64 r, const T &is_below_, const S &is_above_tangent_) {
	auto is_below = [&] (s64 x, s64 y) {
		if (x < l || x > r) return false;
		return is_below_(x, y);
	};
	auto is_above_tangent = [&] (s64 x, s64 dx, s64 dy) {
		if (x < l || x > r) return true;
		if (dx == 0) return false;
		return is_above_tangent_(x, dx, dy);
	};
	
	// find starting point at l
	s64 x = l;
	s64 y;
	{
		s64 low = -1000000000000000000;
		s64 high = 1000000000000000000;
		assert(is_below(l, low));
		assert(!is_below(l, high));
		while (high - low > 1) {
			s64 mid = (low + high) / 2;
			// std::cerr << "! " << l << " " << mid << " " << is_below(l, mid) << std::endl;
			if (is_below(l, mid)) low = mid;
			else high = mid;
		}
		y = low;
	}
	// Invariant I1: is_below(x, y) && !is_below(x, y + 1)
	
	
	s128 res = y; // lattice points on x = `x` is precounted here
	if (l == r) return res;
	
	// Invariant I2: every vector in stack keeps (x, y) inside the area, more downward one on the bottom of the stack
	std::stack<std::pair<s64, s64> > stack;
	// Invariant I3: stack.top() is the most upward vector (dx, dy) satisfying is_below(x + dx, y + dy)
	// tighten_stack ensures the above conditions, given another vector (dx2, dy2) which drives (x, y) outside the area.
	auto tighten_stack = [&] () {
		s64 dx2, dy2;
		std::tie(dx2, dy2) = stack.top();
		stack.pop();
		assert(stack.size());
		
		s64 dx1 = 0, dy1 = 0; // these are dummy values to suppress compiler warnings; consider them uninitialized
		// Invariant: (dx1, dy1) always keeps (x, y) inside the area and equals stack.top()
		// Invariant: (dx2, dy2) always keeps (x, y) outside the area
		
		// search for (dx1, dy1)
		while (!stack.empty()) {
			std::tie(dx1, dy1) = stack.top();
			if (is_below(x + dx1, y + dy1)) break;
			stack.pop();
			dx2 = dx1, dy2 = dy1;
		}
		assert(stack.size());
		// std::cerr << "  found (dx1, dy1)=("<< dx1 << "," << dy1 << ")" << std::endl;
		// std::cerr << "  found (dx2, dy2)=("<< dx2 << "," << dy2 << ")" << std::endl;
		
		while (1) {
			// std::cerr << "     (dx1, dy1)=("<< dx1 << "," << dy1 << ")" << std::endl;
			// std::cerr << "     (dx2, dy2)=("<< dx2 << "," << dy2 << ")" << std::endl;
			s64 dx12 = dx1 + dx2, dy12 = dy1 + dy2;
			if (is_below(x + dx12, y + dy12)) {
				stack.emplace(dx1 = dx12, dy1 = dy12);
			} else {
				// if (dx1, dy1) points more upward than tangent at (x + dx12, y + dy12),
				// addition of a (dx12, dy12) and any number of (dx1, dy1) will still drive (x, y) outside the area.
				// therefore, (x1, y1) is the most upward vector that keeps (x, y) inside the area.
				if (is_above_tangent(x + dx12, dx1, dy1)) break;
				dx2 = dx12, dy2 = dy12;
			}
		}
	};
	stack.push({0, -1});
	stack.push({1, 0});
	stack.push({0, 1});
	tighten_stack();
	
	while (1) {
		s64 dx1, dy1; std::tie(dx1, dy1) = stack.top();
		// Move (x, y) along (dx1, dy1) and add the number of lattice points below the straight line
		// We cannot miss lattice points here, because if there were such point P, vector (x, y) -> P must
		//  be less vertical than (dx1, dy1) contradicting I3
		while (is_below(x + dx1, y + dy1)) {
			// std::cerr << "Move: from (" << x << "," << y << ")  d=(" << dx1 << ", " << dy1 << ")" << std::endl;
			res += (s128) dx1 * (y - 1) + (((s128) (dx1 + 1) * (dy1 + 1)) >> 1);
			// std::cerr << "res += " << (s128) dx1 * (y - 1)  << " + " << (((s128) (dx1 + 1) * (dy1 + 1)) >> 1) << std::endl;
			x += dx1; y += dy1;
		}
		
		if (x == r) break;
		assert(x < r);

		tighten_stack();
	}
	
	return res;
}



/*
	For the following functions:
	
	Let f be a convex/concave function(depending on the function name)
	compare_f(s64 x, s64 y) must return -1, 0, 1 if y < f(x), y = f(x), y > f(x), respectively.
		guaranteed:
		 - x in [l, r]
	compare_tangent(s64 x, s64 dx, s64 dy) must return -1, 0, 1 if dy/dx < f'(x), dy/dx = f'(x), dy/dx > f'(x), respectively.
		guaranteed:
		 - x in [l, r]
		 - dx > 0
	
	The function returns sum_{x = l...r} {floor|ceil}(f(x))    (depending on the function name)
*/
template<typename T, typename S> s128 concave_floor_sum(s64 l, s64 r, const T &compare_f, const S &compare_tangent) {
	return curve_lattice_point_(l, r,
		[&] (s64 x, s64 y) { return compare_f(x, y) <= 0; },
		[&] (s64 x, s64 dx, s64 dy) { return compare_tangent(x, dx, dy) >= 0; }); 
}
template<typename T, typename S> s128 concave_ceil_sum(s64 l, s64 r, const T &compare_f, const S &compare_tangent) {
	return curve_lattice_point_(l, r,
		[&] (s64 x, s64 y) { return compare_f(x, y - 1) < 0; }, // y <= ceil(f(x))  <=> y-1 < f(x)
		[&] (s64 x, s64 dx, s64 dy) { return compare_tangent(x, dx, dy) >= 0; }); 
}

// if f is convex, floor sum of f is -(ceil sum of -f)
template<typename T, typename S> s128 convex_floor_sum(s64 l, s64 r, const T &compare_f, const S &compare_tangent) {
	return -curve_lattice_point_(l, r,
		[&] (s64 x, s64 y) { return compare_f(x, -(y - 1)) > 0; },
		[&] (s64 x, s64 dx, s64 dy) { return compare_tangent(x, dx, -dy) <= 0; }); 
}
// ceil sum of f is -(floor sum of -f)
template<typename T, typename S> s128 convex_ceil_sum(s64 l, s64 r, const T &compare_f, const S &compare_tangent) {
	return -curve_lattice_point_(l, r,
		[&] (s64 x, s64 y) { return compare_f(x, -y) >= 0; },
		[&] (s64 x, s64 dx, s64 dy) { return compare_tangent(x, dx, -dy) <= 0; }); 
}

template<typename T> int spaceship_operator(const T &x, const T &y) {
	return x < y ? -1 : x == y ? 0 : 1;
}
