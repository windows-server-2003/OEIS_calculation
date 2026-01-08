# A186736(n) - Maximum sum of pairwise coprime integers no larger than n
We call primes smaller than or equal to $\sqrt{n}$ *small* and call primes exceeding $\sqrt{n}$ *large*.  
Let $S_{\mathrm{opt}}$ be the optimal solution (set of pairwise coprime integers).  
Intuition:  
 - There are much more large primes than small primes
 - Large primes cannot appear more than once in an integer no larger than $n$.
 - $\rightarrow$ Most large primes are "free": it is in $S_{\mathrm{opt}}$
 - Therefore, it is not worth grouping more than one small primes into a single integer in $S_{\mathrm{opt}}$  
   It is better to split the group and pair the small primes with the abundant large primes.

For specific $n$, we try to prove the following statement:  
$S_{\mathrm{opt}}$ does not contain an integer consisting of two or more different small primes.  

The proof details are in the code.  
The strategy I took in the code succeeds in proving the above statement for $n \le 10^6$ **except** for around 1600 values of $n$, all below 9000.  
It probably succeeds in proving for further $n$ with no exception anymore, though we have to actually run the program for those $n$'s to know whether it does.  

If the above statement is proved, the problem is reduced to a maximum weight matching in a bipartite graph (small primes on the left vertices, large primes on the right vertices).  
If it fails to prove the above statement, we use a less efficient sub-exponential algorithm which does not assume the above statements.  

It turns out (by comparison of the results of the two algorithms) that the statement is actually true for every $n \le 10^6$.

### Program
Compile with
```
g++ main.cpp -o main -O2 -Wall -Wextra -I../include -Wno-stringop-overflow
```

Run with `./main` and input
 - maximum $n$ for which the program computes terms
 - number of threads it uses

Running with `./main -f` will write the results to `b186736.txt` (Warning: this overwrites previous results).  

`b186736_slow.txt` is the results of the second sub-exponential algorithm up to $n = 19000$.  
