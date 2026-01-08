#include "utils.h"
#include "debug.h"
#include "types.h"
#include "oeis.h"

#include "prime.h"

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

int ri() {
	int n;
	scanf("%d", &n);
	return n;
}


/*
	Subexponential solution without any advanced assumption of the property of optimal solution.
	Used only when the attempt to prove "no two small primes" statement fails.
*/
struct A186736_slow {
	int n;
	
	std::vector<bool> is_prime;
	std::vector<int> primes;
	int sq; // number of small primes
	
	void init_primes(int max) {
		is_prime.assign(max + 1, true);
		primes.clear();
		sq = 0;
		
		is_prime[0] = is_prime[1] = false;
		for (int i = 2; i <= max; i++) if (is_prime[i]) {
			primes.push_back(i);
			for (int j = i + i; j <= max; j += i) is_prime[j] = false;
		}
		for (; sq + 1 < (int) primes.size() && primes[sq] * primes[sq + 1] <= n; sq++);
	}
	
	
	struct Info {
		u64 mask;
		int max_prod;
	};
	// prod_info[p]: list of (mask: usage of small primes, maximum single integer <= n consisting solely of p and small primes in `mask`)
	// p is the largest prime factor of .max_prod
	std::vector<std::vector<Info> > prod_info;
	
	void init_prod_info() {
		prod_info.assign(primes.size(), {});
		
		for (int i = 2; i <= n; i++) {
			u64 mask = 0;
			int first = -1;
			for (int j = primes.size(); j--; ) if (i % primes[j] == 0) {
				if (first == -1) first = j; // largest prime factor
				else assert(j < sq), mask |= 1ULL << j; // other prime factors
			}
			assert(first != -1);
			[&] () {
				for (auto &j : prod_info[first]) {
					if (j.mask == mask) {
						j.max_prod = i;
						return;
					}
				}
				prod_info[first].push_back({mask, i});
			}();
			// remove unnecessary one(if there is another combination of primes that is a subset of current one, with greater max_prod)
			for (auto &info : prod_info) {
				for (auto itr = info.begin(); itr != info.end(); ) {
					bool unnecessary = false;
					for (auto j : info) if (j.mask != itr->mask && (j.mask & itr->mask) == j.mask && j.max_prod >= itr->max_prod) {
						unnecessary = true;
						break;
					}
					if (unnecessary) itr = info.erase(itr);
					else itr++;
				}
			}
		}
		/*
		for (int i = 0; i < (int) prod_info.size(); i++) {
			std::cerr << primes[i] << " : " << std::endl;
			for (auto j : prod_info[i]) {
				std::cerr << "  " << j.mask << " : " << j.max_prod << std::endl;
			}
		}*/
	}
	
	u32 calc(int n) {
		/* val_t is not used internally */
		if (n == 1) return 1;
		
		this->n = n;
		
		init_primes(n);
		
		assert(n < (1 << 15)); // to avoid overflow
		assert(sq <= 63);
		init_prod_info();
		
		// std::cerr << "prod info inited" << std::endl;
		
		// dp[i]: maximum sum of integers with small factor usage i
		std::vector<u32> dp = { 0 };
		for (int i = primes.size(); i--; ) {
			// use one of integers in prod_info[i] (= integers whose largest prime factor is primes[i])
			int next_bits = 0;
			while (next_bits < (int) primes.size() && next_bits < i && primes[next_bits] * primes[i] <= n) next_bits++;
			size_t next_size = 1ULL << next_bits;
			
			if (i >= sq) { // dp table expanding
				assert(next_size >= dp.size());
				size_t dp_size_org = dp.size();
				dp.resize(next_size);
				for (size_t j = dp_size_org; j--; ) {
					for (int k = prod_info[i].size(); k--; ) if (!(j & prod_info[i][k].mask)) {
						auto &target = dp[j | prod_info[i][k].mask];
						target = std::max(target, dp[j] + prod_info[i][k].max_prod);
					}
				}
			} else { // dp table shrinking
				assert(next_size == dp.size() / 2);
				for (size_t j = next_size; j--; ) {
					if (!(j >> i & 1)) {
						for (int k = prod_info[i].size(); k--; ) if (!(j & prod_info[i][k].mask)) {
							auto &target = dp[j | prod_info[i][k].mask];
							target = std::max(target, dp[j] + prod_info[i][k].max_prod);
						}
					}
				}
				for (size_t j = next_size; j < dp.size(); j++) dp[j - next_size] = std::max(dp[j - next_size], dp[j]);
				dp.resize(next_size);
			}
			
			if (dp.size() >= (1 << 20)) {
				std::cerr << i << ": " << put_SI_prefix(dp.size() * 4) << "B " << prod_info[i].size() << std::endl;
			}
		}
		
		return dp.back() + 1; // 1 can be used
	}
};






struct A186736 {
	int n_min;
	std::vector<int> divisor; // divisor[n]: a prime factor of ns
	std::vector<int> primes;

	void init_primes(int max) {
		divisor.resize(max + 1, 0);
		for (int i = 2; i <= max; i++) if (divisor[i] == 0) {
			primes.push_back(i);
			for (int j = i; j <= max; j += i) divisor[j] = i;
		}
	}
	
	A186736 (int n_min, int n_max) : n_min(n_min) {
		init_primes(n_max);
	}
	
	/*
		Attempt to prove there is an optimal solution without any group with
			more than one small primes(primes <= sqrt(n)).
		Return:
		 - true if it succeeds in proving
		 - false if it fails. This does not imply the above statement is false.
	*/
	bool test_no_group_with_multiple_small_primes(int n) const {
		// m: number of primes <= ns
		// s: number of small primes
		int m = std::upper_bound(primes.begin(), primes.end(), n) - primes.begin();
		int s = std::upper_bound(primes.begin(), primes.end(), (int) std::sqrt(n)) - primes.begin();
		
		if (s >= 1) assert((s64) primes[s - 1] * primes[s - 1] <= n);
		if (s != m) assert((s64) primes[s] * primes[s] > n);
		
		/*
			Property of profit[i] :
				Let p = primes[i]. For *any* valid set S, one can create another valid set S' by removing from S the integer divisible by p and
					adding an integer either of form
					 - singleton p^l (l >= 1)
					 - pair p^l q (l >= 1), where q is a free large prime in S
					that has profit of profit[i].
				Moreover, one can do above thing simultaneously for multiple small primes: let p_1 = primes[i_1], ..., p_k = primes[i_k] distinct small primes.
				For *any* valid set S, one can create another valid set S' by removing from S all intergers divisible by at least one of p_i's and
					for each p_j adding an integer either of form
					 - singleton p_j^l (l >= 1)
					 - pair p_j^l q_j (l >= 1), where q_j is a free large prime in S, distinct across i
					that has profit of profit[i_j]
				
				This means, if an integer x in S_opt has small prime factors p_1, ..., p_k, f(x) >= sum profit[p_i] must hold because
					otherwise removing x and adding the above numbers would make S_opt better, contradicting optimality.
		*/
		std::vector<int> profit(s);
		for (int i = 0; i < s; i++) {
			int p = primes[i];
			
			int p_score = 0;
			for (s64 pp = p; pp <= n; pp *= p) {
				p_score = std::max<s64>(p_score, pp - p); // p^k
				
				// enumerate large prime q such that pp * q <= n, in descending order
				int t = std::upper_bound(primes.begin(), primes.end(), n / pp) - primes.begin();
				for (int i = t - 1; i >= s; i--) {
					int q = primes[i];
					assert((s64) q * q > n && (s64) pp * q <= n); // q is large && pp * q <= n
					
					int n_primes_can_take_q = std::upper_bound(primes.begin(), primes.end(), n / q) - primes.begin();
					assert(n_primes_can_take_q <= s);
					// Under S, each of q = primes[i] <, primes[i+1] < ... < primes[t - 1] is either free or grouped with (at least) one of n_primes_can_take_q small primes.
					// If n_primes_can_take_q <= s, after the group containing p is erased, there is at least one q' in primes[i...t-1] that can be paired with pp.
					if (n_primes_can_take_q <= t - i) { // 
						// we only know q' is in primes[i...t-1] but q' >= q holds so f(pp * q') >= f(pp * q) holds
						p_score = std::max<s64>(p_score, pp * q - p - q);
						break;
					}
				}
			}
			profit[i] = p_score;
			// std::cerr << p << " : " << profit[i] << std::endl;
		}
		
		/*
		for (int i = 1; i <= n; i++) {
			// consider using i
			s64 consumed_profit = 0;
			int n_small_primes_used = 0;
			for (int j = 0; j < s; j++) if (i % primes[j] == 0) {
				consumed_profit += profit[j];
				n_small_primes_used++;
			}
			s64 current_profit = i; // profit of i
			for (int j = 0; j < m; j++) if (i % primes[j] == 0) current_profit -= primes[j];
			
			if (n_small_primes_used >= 2 && consumed_profit < current_profit) {
				std::cerr << "possible: " << i << std::endl;
				return false;
			}
		}
		*/
		
		for (int i = 0; i < s; i++) for (int j = 0; j < i; j++) {
			// check if the optimal solution can contain a group including both primes[i] and primes[j]
			// profit of the group (which does not exceed n) must be larger than profit[i] & profit[j]
			if (profit[i] + profit[j] > n) continue; // so in this case, it cannot be optimal
			
			const int p = primes[i];
			const int q = primes[j];
			assert(q < p);
			assert(p * q <= n);
			// when using third small prime r other than p or q
			{
				int r_upperbound = n / q / p;
				for (int k = 0; k < s; k++) {
					int r = primes[k];
					if (r > r_upperbound) break;
					if (r == p || r == q) continue;
					// this condition is somewhat weaker because it bounds the profit of resulting integer(which is a multiple of pqr) by n,
					// but it turns out this condition is never satisfied anyway(at least for n <= 10^6)
					if (profit[i] + profit[j] + profit[k] <= n) {
						// std::cerr << "  (p, q, r) = (" << p << "," << q << "," << r << ") (third small prime) possible" << std::endl;
						return false;
					}
				}
			}
			for (s64 pp = p; pp <= n; pp *= p) {
				for (s64 qq = q; qq <= n / pp; qq *= q) { // range for which pp * qq <= n
					s64 ppqq = pp * qq;
					
					// when using only p and q
					if (ppqq - p - q > profit[i] + profit[j]) {
						// std::cerr << "  " << ppqq << " (ppqq) possible" << std::endl;
						return false;
					}
					// p, q, and a large prime r
					int r_upperbound = n / ppqq;
					for (int k = s; k < m; k++) {
						int r = primes[k];
						if (r > r_upperbound) break;
						int t = ppqq * r; // the value; its profit is t - p - q - r
						if (t - p - q - r > profit[i] + profit[j]) {
							// std::cerr << "  " << t << " (large) possible" << std::endl;
							return false;
						}
					}
				}
			}
			/*
				here, it is not possible that the optimal solution has a group with p & q because the code above confirmed
				1. the group cannot contain any small prime other than p and q
				2. the group cannot consist solely of (one or more) p and q
				3. the group cannot consist of (one or more) p and q and a large prime r
			*/
		}
		return true;
	}
	template<typename T> T solve_for_single_small_prime_per_group(int n) {
		int m = std::upper_bound(primes.begin(), primes.end(), n) - primes.begin();
		int s = std::upper_bound(primes.begin(), primes.end(), (int) std::sqrt(n)) - primes.begin();
		
		/*
			We can now reduce the problem into maximum weighted matching (in bipartite graph).
			That is
			 - Left vertices are the small primes
			 - Right vertices are the large primes
			 - For each small prime p, let base_profit(p) be its maximum profit when it is a singleton(power of p).
			   For each large prime q, we create edge between p and q of weight (maximum profit of integers of form p^l q) - base_profit(p)
			Note that choosing singleton of a small prime p corresponds to not using vertex p in the matching.
			
			However, that is too many edges (I think), so we prune edges by carefully examining the property of S_opt.
			Roughly, for a fixed small prime p, p wants to be paired with one of
			 - the largest large prime below n/pp where pp is of form p^l (1 <= l)
			Preference of multiple pp (for different p) may conflict, but one can say that p is paired in S_opt with
			 - the top-few largest prime below n/pp, depending on how much the preference is congested around those large primes.
		*/
		
		struct opt_pair_t {
			int small_id; // p = primes[small_id]
			int pp; // pp = p^l
			int preferred_large_id; // q = primes[preferred_large_id]
		};
		std::vector<opt_pair_t> opt_pairs;
		std::vector<int> base_profit(s); // profit when p is used as a singleton
		for (int i = 0; i < s; i++) {
			int p = primes[i];
			for (s64 pp = p; pp * pp <= n; pp *= p) {
				int t = std::upper_bound(primes.begin(), primes.end(), n/pp) - primes.begin();
				if (t <= s) break;
				// primes[t-1]: best(largest compatible) large prime with pp
				opt_pairs.push_back({i, (int) pp, t - 1});
			}
			s64 pp = p;
			while (pp * p <= n) pp *= p;
			base_profit[i] = pp - p;
		}
		
		std::sort(opt_pairs.begin(), opt_pairs.end(), [] (auto &i, auto &j) {
			return i.pp > j.pp;
		});
		// in descending order of pp, simulate the worst case: every opt_pair is used(which is impossible due to multiplicity of pp)
		// used_large_ids: set of (index of) large primes that may be used in S_opt
		std::set<int> used_large_ids;
		struct edge_t {
			int small_id;
			int large_id;
			int weight;
		};
		std::vector<edge_t> edges;
		for (auto [small_id, pp, preferred_large_id] : opt_pairs) {
			/*
			  if (in S_opt) primes[preferred_large_id] is used by another opt_pair with larger pp, search for smaller large prime.
			  if we find a free one, we can assert the following statement.
			    Let primes[preferred_large_id - t] be the largest prime below primes[preferred_large_id] that is not used in S_opt by pp's *larger* than current pp.
			    Then, it is optimal to pair pp with primes[preferred_large_id - t] if we are to use current pp (not another power of p or a singleton of p)
			  This is because otherwise
			   - pp is paired with yet smaller large prime primes[preferred_large_id - t - s]
			   - primes[preferred_large_id - t] must be used by a pp smaller than current pp
			   - but then swapping the pairing of the smaller pp and current pp will improve the profit
			  So below we simulate the worst case: assuming every pp is paired with a large prime
			*/
			for (int j = preferred_large_id; j >= s; j--) {
				int cur_profit = primes[j] * pp - primes[j] - primes[small_id];
				edges.push_back({small_id, j, cur_profit - base_profit[small_id]});
				if (!used_large_ids.count(j)) { // found a large prime not used by larger pp's even in the worst case
					used_large_ids.insert(j);
					break; // current pp must be paired with one of the large primes enumerated so far in this loop
				}
			}
		}
		
		using graph_t = boost::adjacency_list<
			boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property,
			boost::property<boost::edge_weight_t, int>
		>;
		// only a small number of large primes is used in the graph, so filter the right vertices
		std::vector<int> used_large_ids_vec(used_large_ids.begin(), used_large_ids.end());
		int t = used_large_ids_vec.size();
		// s small primes and t large primes
		graph_t graph(s + t);
		
		std::cerr << s << " : " << t << std::endl;
		
		// return the vertex id in `graph` of the vertex corresponding to primes[large_id]
		auto get_large_id_in_graph = [&] (int large_id) {
			int res = std::lower_bound(used_large_ids_vec.begin(), used_large_ids_vec.end(), large_id) - used_large_ids_vec.begin();
			assert(used_large_ids_vec[res] == large_id);
			return s + res; // vertices for large primes are after the s small prime vertices
		};
		
		auto weights = get(boost::edge_weight, graph);
		auto vindex = get(boost::vertex_index, graph);
		
		for (auto [small_id, large_id, weight] : edges) {
			// std::cerr << "edge small#" << small_id << ":" << primes[small_id] << "  large_id=" << large_id << ":" << primes[large_id] << " w=" << weight << std::endl;
			if (weight <= 0) continue;
			auto e = add_edge(small_id, get_large_id_in_graph(large_id), graph).first;
			weights[e] = weight;
		}
		
		std::vector<graph_t::vertex_descriptor> mate(s + t);
		boost::maximum_weighted_matching(graph, mate.data(), vindex, weights);
		
		s64 opt_profit = 0;
		// calculate weight sum of the maximum matching: enumerate over left vertices
		for (int i = 0; i < s; i++) if (mate[i] != graph_t::null_vertex()) { // vertex #i is matched
			auto [e, ok] = edge(i, mate[i], graph);
			assert(ok);
			opt_profit += weights[e];
		}
		// Enabling this will show information of the optimal solution
		if (0) {
			std::cerr << "  ---- Optimal Solution (small prime usage) ----" << std::endl;
			for (int i = 0; i < s; i++) {
				int p = primes[i];
				if (mate[i] == graph_t::null_vertex()) {
					int val = 1;
					int p_power = 0;
					while ((s64) val * p <= n) val *= p, p_power++;
					fprintf(stderr, "  %3d: %6d = %3d^%d (singleton)\n", p, val, p, p_power);
				} else {
					int q = primes[used_large_ids_vec[mate[i] - s]];
					int val = q;
					int p_power = 0;
					while ((s64) val * p <= n) val *= p, p_power++;
					fprintf(stderr, "  %3d: %6d = %3d^%d x %-5d (paired with large prime q=%d)\n", p, val, p, p_power, q, q);
				}
			}
			std::cerr << "  All large primes not present above is used as pure singleton" << std::endl;
		}
		// Add base_profit(p) for every p
		opt_profit += std::accumulate(base_profit.begin(), base_profit.end(), (s64) 0);
		// (sum of elements in S) = (sum of all primes) + (profit)
		for (int i = 0; i < m; i++) opt_profit += primes[i];
		
		return opt_profit + 1; // 1 is always in S_opt
	}
	
	template<typename T> T calc(int n) {
		// n is prime: no need for calculation. a(n) = a(n - 1) + n
		if (n != n_min && divisor[n] == n) return -1; // -1 means a(n) = a(n - 1) + n
		int n_small_primes = 0;
		for (int t = n; t != 1; ) {
			int p = divisor[t];
			while (t % p == 0) t /= p;
			if ((s64) p * p <= n) n_small_primes++;
		}
		
		if (test_no_group_with_multiple_small_primes(n)) {
			if (n != n_min && n_small_primes >= 2) return -2; // -2 means a(n) = a(n - 1)
			return solve_for_single_small_prime_per_group<T>(n);
		} else {
			return A186736_slow().calc(n);
		}
	}
	template<typename val_t> void finalize(int n_min, int n_max, std::vector<val_t> &result) {
		for (int n = n_min; n <= n_max; n++) {
			if (result[n] == -1) result[n] = result[n - 1] + n;
			else if (result[n] == -2) result[n] = result[n - 1];
		}
	}
};

int main(int argc, char **argv) {
	std::cerr << "Max n > ";
	int n_max = ri();
	std::cerr << "# of threads > ";
	int n_threads = ri();
	
	
	OEISGenerator gen;
	int seq_num = 186736;
	OEISGenerator::GenOption opt = gen.get_cmdline_options(argc, argv, seq_num);
	
	opt.n_min = 1;
	opt.n_max = n_max;
	
	opt.thread_num = n_threads;
	gen.solve<A186736>(opt);
	
	return 0;
}
