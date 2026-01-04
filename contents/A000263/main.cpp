#include <cstdio>
#include <iostream>
#include <thread>
#include <stack>
#include <mutex>
#include <deque>
#include "utils.h"
#include "debug.h"
#include "types.h"
#include "curve_lattice_point.h"
#include "oeis.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}


/*
	Legacy code; not used anymore
	This has been generalized in curve_lattice_point.h and that version is currently used.
*/
struct A000263_old {
	/*
		Count \sum_{x = 1}^{n^2} ceil(2n sqrt(x))
	*/
	template<typename val_t> val_t solve_subroutine(int n) {
		// whether (x, y) is below y < ceil(2n sqrt(x)) + 1
		auto inside = [&] (s64 x, s64 y) {
			return (val_t) (y - 1) * (y - 1) < (val_t) 4 * n * n * x;
		};
		// whether (-dx1, -dy1) points less downward than tangent of the parabola at (x, y)
		auto cut = [&] (s64 x, s64 y, s64 dx1, s64 dy1) {
			(void) y;
			return (val_t) n * n * dx1 * dx1 >= (val_t) x * dy1 * dy1;
		};

		// Invariant I1: (x, y) is always inside the area, but (x, y + 1) is outside
		s64 x = (s64) n * n;
		s64 y = 2LL * n * n;

		// Invariant I2: more vertical vector on bottom
		// Invariant I3: stack.top() is the least vertical vector (dx, dy) that still keeps (x - dx, y, dy) inside the parabola.
		// In this case, the tangent of the parabola at x=`x` has direction (1, 1)
		auto stack = std::stack<std::pair<s64, s64> >({{0, 1}, {1, 1}});
		
		val_t res = 2LL * n * n; // lattic points on x = `x` is precounted here
		while (1) {
			s64 dx1, dy1; std::tie(dx1, dy1) = stack.top();
			stack.pop();
			// Move (x, y) along (dx1, dy1) and add the number of lattice points below the straight line
			// We cannot miss lattice points here, because if there were such point P, vector (x, y) -> P must
			//  be less vertical than (dx1, dy1) contradicting I3
			while (inside(x - dx1, y - dy1)) {
				x -= dx1; y -= dy1;
				res += (val_t) dx1 * (y - 1) + (((val_t) (dx1 + 1) * (dy1 + 1)) >> 1) - dy1;
			}

			// search for two consecutive vectors in stack, where one drives (x, y) outside the parabola and the other inside.
			s64 dx2 = dx1, dy2 = dy1;
			while (!stack.empty()) {
				std::tie(dx1, dy1) = stack.top();
				if (inside(x - dx1, y - dy1)) break;
				stack.pop();
				dx2 = dx1, dy2 = dy1;
			}
			if (stack.empty()) break;
			
			// Invariant I4:
			// (x - dx1, y - dy1): inside the parabola; equals stack.top()
			// (x - dx2, y - dy2): outside the parabola
			
			// binary search between (dx1, dy1) and (dx2, dy2)
			while (1) {
				s64 dx12 = dx1 + dx2, dy12 = dy1 + dy2;
				if (inside(x - dx12, y - dy12)) {
					stack.emplace(dx1 = dx12, dy1 = dy12);
				} else {
					// if (-dx1, -dy1) points less downward than tangent at (x - dx12, y - dy12),
					// addition of a (-dx12, -dy12) and any number of (-dx1, -dy1) will still drive (x, y) outside the area.
					// therefore, (x1, y1) is the least vertical vector that keeps (x, y) inside the area.
					if (cut(x - dx12, y - dy12, dx1, dy1)) break;
					dx2 = dx12, dy2 = dy12;
				}
			}
		}
		
		return res;
	}
	
	template<typename val_t> val_t solve(int n) {
		val_t res = (val_t) n * n * n * n + (val_t) ((u64) n * n) * ((u64) n * n + 1) / 2;
		res -= solve_subroutine<val_t>(n);
		val_t one = (u64) n * n / 4;
		assert(!((res - one) & 1));
		return (res - one) / 2;
	}
	template<typename val_t> val_t calc(int n) {
		return solve<s128>(n);
	}
};


struct A000263 {
	template<typename val_t> val_t calc(int n) {
		s64 l = 1;
		s64 r = n*n;
		/*
			Compare sqrt(x) + sqrt(y) against n without using floating-point.
			sqrt(x) + sqrt(y) < n
			<=> sqrt(y) < n - sqrt(x)
			<=> y < n^2 - 2n sqrt(x) + x  (n - sqrt(x) is nonnegative for x in range [l, r])
			<=> 2nsqrt(x) < n^2 + x - y
			<=> 4n^2 x < (n^2 + x - y)^2  &&  (n^2 + x - y) >= 0
		*/
		auto compare_f = [&] (s64 x, s64 y) {
			if (y <= 0) return -1;
			s64 t = n*n + x - y;
			if (t <= 0) return 1;
			return spaceship_operator((s128) 4 * n * n * x, (s128) t * t);
		};
		/*
			From above transformation, the curve is
				f(x) = n^2 - 2nsqrt(x) + x
			So
				f'(x) = 1 - n / sqrt(x)
			We want to compare dy/dx against this without using floating-point.
				dy/dx < f'(x)
				<=> n / sqrt(x) < 1 - dy/dx
				<=> n^2 / x < (1 - dy/dx)^2   &&  dy/dx < 1
				<=> n^2 dx^2 < x (dx - dy)^2   &&  dy/dx < 1
			
			* dx is guaranteed to be positive.
		*/
		auto compare_tangnent = [&] (s64 x, s64 dx, s64 dy) {
			if (dy >= dx) return 1;
			s128 dxn = (s128) dx * n;
			
			return spaceship_operator((s256) dxn * dxn, (s256) x * (dx - dy) * (dx - dy));
		};
		auto res = convex_floor_sum(l, r, compare_f, compare_tangnent);
		auto diagonal = (n * n / 4); // number of solutions on the x = y line
		return (res - diagonal) / 2;
	}
};

int main(int argc, char **argv) {
	std::cout << "Max n > ";
	int n_max = ri();
	std::cout << "# of threads > ";
	int n_threads = ri();
	
	OEISGenerator gen;
	int seq_num = 263;
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	opt.n_min = 1;
	opt.n_max = n_max;
	opt.thread_num = n_threads;
	
	gen.solve<A000263>(opt);
	
	return 0;
}
