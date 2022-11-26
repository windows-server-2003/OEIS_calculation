#include <cstdio>
#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <map>
#include "oeis.h"
#include "debug.h"
#include "types.h"

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}

struct A004204 {
	template<typename val_t> val_t calc_(int n) {
		std::vector<std::vector<val_t> > dp(1 << n, std::vector<val_t>(n));
		dp[0][0] = 1;
		for (int i = 0; i + 1 < 1 << n; i++) {
			int coef = __builtin_popcount(i);
			for (int j = 0; j < n; j++) if (!(i >> j & 1)) {
				int add = coef * j % n;
				for (int k = 0; k < n - add; k++) dp[i | 1 << j][k + add] += dp[i][k];
				for (int k = n - add; k < n; k++) dp[i | 1 << j][k + add - n] += dp[i][k];
			}
		}
		return dp.back()[0];
	}
	template<typename val_t> val_t calc(int n) {
		// val_t is not used internally
		return calc_<__uint128_t>(n);
	}
};

template<typename T, size_t B> struct HashSet {
	std::array<T, B> table;
	std::vector<size_t> table_changes;
	
	HashSet () { std::fill(table.begin(), table.end(), (T) -1); }
	u32 get_hash(T x) {
		x *= 0xbf58476d1ce4e5b9;
		x = x ^ (x >> 31);
		x = x ^ (x >> 13);
		return x & (B - 1);
	}
	bool set_1(T x) {
		auto hash = get_hash(x);
		while (table[hash] != (T) -1) {
			if (table[hash] == x) return true;
			hash = (hash + 1) & (B - 1);
		}
		table[hash] = x;
		table_changes.push_back(hash);
		return false;
	}
	void reset() {
		for (auto i : table_changes) table[i] = (T) -1;
		table_changes.clear();
	}
};

template<typename val_t> struct A004204_fast {
	int n;
	struct Task {
		u64 mask;
		u64 cnt[2];
	};
	std::deque<std::vector<Task> > que;
	std::vector<val_t> global_res;
	std::mutex lock;
	bool finish;
	
	static constexpr size_t PACKET_SIZE = 100000;
	
	template<typename T> T gcd(T a, T b) {
		while (a && b) {
			if (a > b) a %= b;
			else b %= a;
		}
		return a + b;
	}
	void process_thread_func() {
		while (1) {
			std::vector<Task> packet;
			{
				std::lock_guard<std::mutex> locking(lock);
				if (que.size()) packet = que.front(), que.pop_front();
			}
			if (!packet.size()) {
				if (finish) break;
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			
			std::cerr << "start packet " << packet.back().mask << std::endl;
			
			std::vector<val_t> res(n);
			for (auto task : packet) {
				u64 i = task.mask;
				auto &cnt = task.cnt;
				
				std::vector<val_t> dp(n);
				std::vector<val_t> next;
				dp[0] = 1;
				for (int j = 0; j < n; j++) {
					next.assign(n, 0);
					for (int k = 0; k < n; k++) if (i >> k & 1) {
						int add = k * j % n;
						for (int l = 0; l < n - add; l++) next[l + add] += dp[l];
						for (int l = n - add; l < n; l++) next[l + add - n] += dp[l];
					}
					std::swap(dp, next);
				}
				
				bool minus = (n - __builtin_popcountll(i)) & 1;
				for (int j = 0; j < n; j++) {
					val_t cur_res;
					if (n & 1) cur_res = dp[j] * (cnt[0] + cnt[1]);
					else {
						int index = j + n / 2;
						if (index >= n) index -= n;
						cur_res = dp[j] * cnt[0] + dp[index] * cnt[1];
					}
					if (minus) res[j] -= cur_res;
					else res[j] += cur_res;
				}
			}
			
			{
				std::lock_guard<std::mutex> locking(lock);
				for (int i = 0; i < n; i++) global_res[i] += res[i];
				std::cerr << "finish packet " << packet.back().mask << std::endl;
			}
		}
	}
	void main_thread_func() {
		std::vector<int> coprime;
		for (int i = 1; i < n; i++) if (gcd(n, i) == 1) coprime.push_back(i);
		
		constexpr int K = 8;
		constexpr int L = 5;
		std::vector<std::vector<std::vector<u64> > > mul_table(n, std::vector<std::vector<u64> >(L, std::vector<u64>(1 << K)));
		for (int i = 0; i < n; i++) for (int j = 0; j < L; j++) for (u64 k = 0; k < 1 << K; k++) {
			auto bit_mul_trans = [] (u64 x, int i, int n) {
				u64 res = 0;
				for (int j = 0; j < n; j++) if (x >> j & 1) res |= 1ULL << (i * j % n);
				return res;
			};
			mul_table[i][j][k] = bit_mul_trans(k << (j * K), i, n);
		}
		HashSet<u64, 1U << 18> used;
		
		// auto t0 = Timer::get();
		std::vector<Task> packet;
		packet.reserve(PACKET_SIZE);
		auto push = [&] () {
			while (1) {
				std::lock_guard<std::mutex> locking(lock);
				if (que.size() >= 50) std::this_thread::sleep_for(std::chrono::milliseconds(50));
				else {
					std::cerr << "push packet " << packet.back().mask << "   size:" << que.size() << std::endl;
					que.push_back(packet);
					break;
				}
			}
			packet.clear();
			// std::cerr << Timer::diff_ms(t0, Timer::get()) << std::endl;
		};
		for (u64 i = 0; i < 1ULL << n; i++) {
			if (i && !(i & 1)) continue;
			if ((i >> (n - 1) & 1) && i != ((1ULL << n) - 1)) continue;
			
			auto bit_rotate = [] (u64 x, int i, int n) { return x >> i | (x & ((1ULL << i) - 1)) << (n - i); };
			u64 cnt[2] = { 0 };
			bool is_min = true;
			for (auto j : coprime) {
				u64 cur = 0;
				for (int k = 0; k < L; k++) cur |= mul_table[j][k][i >> (k * K) & ((1 << K) - 1)];
				for (int k = 0; k < n; k++) {
					u64 t = bit_rotate(cur, k, n);
					if (t < i) {
						is_min = false;
						break;
					}
					if (!used.set_1(t)) cnt[k & 1]++;
				}
				if (!is_min) break;
			}
			used.reset();
			if (!is_min) continue;
			
			packet.push_back({i, {cnt[0], cnt[1]}});
			if (packet.size() >= PACKET_SIZE) push();
		}
		if (packet.size()) push();
		std::cerr << "===== push finish =====" << std::endl;
		finish = true;
	}

	static constexpr int THREAD_NUM = 10;
	val_t calc_(int n) {
		if (n == 1) {
			std::cout << std::vector<val_t>{1} << std::endl;
			return 1;
		}
		this->n = n;
		this->finish = false;
		this->global_res.assign(n, 0);
		std::thread main_thread(&A004204_fast::main_thread_func, this);
		std::vector<std::thread> process_threads;
		for (int i = 0; i < THREAD_NUM; i++) process_threads.push_back(std::thread(&A004204_fast::process_thread_func, this));
		for (auto &thread : process_threads) thread.join();
		main_thread.join();
		
		assert(que.empty());
		
		std::cout << global_res << std::endl;
		
		return global_res[0];
	}
	template<typename res_t> res_t calc(int n) { return calc_(n); }
};

int main(int argc, char **argv) {
	OEISGenerator gen;
	int seq_num = gen.get_sequence_number(argc, argv);
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	if (seq_num == -1) {
		std::cerr << "No valid sequence number found" << std::endl;
		return 1;
	}
	
	opt.n_min = 1;
	opt.n_max = ri();
	opt.thread_num = 1;
	
	switch (seq_num) {
		case   4204 : gen.solve<A004204_fast<u256> >(opt); break;
	}
	
	return 0;
}
