## Number of permutation of a specific spread
*Spread* of a permutation $p$ is defined by the remainder of $\displaystyle \sum_{i=0}^{n-1} ip_i$ divided by $n$.  
[A147679](https://oeis.org/A147679) gives the full table $T_{n, k}$ : number of permutation of $\[0...(n - 1)\]$ with spread $k$.  
Also, [A004204](https://oeis.org/A004204), [A004205](https://oeis.org/A004205), [A004206](https://oeis.org/A004206), [A004246](https://oeis.org/A004246) gives the number of permutation of $\[0...(n - 1)\]$ with spread $0, 1, 2, 3$, respectively.  

This document gives an $O(\frac{2^nn^2}{\phi(n)})$ algorithm to find the value of $T_{n, k}$ for all $k$, where $\phi(n)$ is the Euler's totient function.  
A part of the calculation of complexity $O(2^n)$ needs to be executed in a single thread, and the rest of the part can be parallelized.  

Implementation is available in [main.cpp](main.cpp), which is used to calculate up to $T(40, \*)$ spending 30 hours on an Intel Core i5-10600k machine.  
You can compile it with gcc with `-I../include -Wl,--stack=0x1000000`(for larger stack) and run ```[executable name] > stdout.txt``` with an integer from stdin which indicates the upper bound of $n$, then ```process.py``` to generate the b-files for the five sequences.  

## Algorithm
Inclusion-exclusion principle can be applied. Let $S$ be the set $\\{ 0, 1, 2, \dots, n - 1 \\}$. Then  
$\displaystyle T(n, k) = \sum_{S' \subset S} (-1)^{n - |S'|}|\\{a = (a_0, a_1, a_2, \dots, a_{n - 1}) \mid a_i \in S'(0 \le i \lt n), \sum_{i = 0}^{n - 1} ia_i = k \hspace{5pt}\mathrm{mod} \hspace{5pt} n\\}|$  
holds.  
Given a specific $S'$, the expression inside the sigma for all $k$ can be calculated using a simple $O(n^3)$ dynamic programming. This immediately gives an $O(2^NN^3)$ algorithm.  
The advantage of this over direct dynamic programming(saving which numbers are already used as the state) solutions without inclusion-exclusion is that its memory consumption is polynomial to $n$ and it can be easily and almost arbitrarily parallelized.  

## Symmetry 1 : adding a constant to every element in $S'$
Now, we use some symmetry to reduce the number of $O(N^3)$ calculation.  

Let $\displaystyle g(S', k) = |\\{a = (a_0, a_1, a_2, \dots, a_{n - 1}) \mid a_i \in S'(0 \le i \lt n), \sum_{i = 0}^{n - 1} ia_i = k \hspace{5pt}\mathrm{mod} \hspace{5pt} n\\}|$.  
Also, let $g(S', i) = \\{ (j + i) \hspace{5pt} \mathrm{mod} \hspace{5pt} n \mid j \in S' \\}$.  

We can see that $f(g(S', i), k)$ is equal to 

 - $f(S', k)$ if $n$ is odd or $i$ is even
 - $f(S', (k + \frac{n}{2}) \hspace{5pt}\mathrm{mod}\hspace{5pt} n)$ otherwise

Thus, it is sufficent to calculate $f(S', k)$ for all $k$ to determine the values of $f(T, k)$ for every $k$ and $T = g(S, 1), g(S, 2), g(S, 3), \dots, g(S, n - 1)$  
This reduces the complexity to $O(2^nn^2)$. (There are cases where there are duplicates in the $n$ sets, but it can be proven that the number of such cases is sufficiently small and does not affect the complexity) 

## Symmetry 2 : multiplying every element in $S'$ by a constant relatively prime to $n$

Let $h(S', i) = \\{ ji \hspace{5pt} \mathrm{mod} \hspace{5pt} n \mid j \in S' \\}$.  
If $\gcd(i, N) = 1$, $|h(S', i)| = |S'|$ and $f(h(S', i), k) = f(S', ki^{-1} \hspace{5pt} \mathrm{mod} \hspace{5pt} n)$ where $i^{-1}$ is the modulo inverse of $i$ which exists because $\gcd(i, n) = 1$.  
Combined with the previous improvement, this reduces the complexity to $O(\frac{2^nn^2}{\phi(n)})$.

