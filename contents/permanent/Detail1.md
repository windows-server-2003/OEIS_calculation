### All 1 except the diagonal
Let $m$ be an $n \times n$ matrix of form $m[i][i] = f(i)$, $m[i][j] = 1 \\, (i \neq j)$.  
By definition $\mathrm{per}(m) = \sum_{\sigma \in S(n)} m[i][\sigma(i)]$ where $S(n)$ is the set of permutation of $[0, 1, 2, \dots, n - 1]$.
Classifying $\sigma$ by the number of $i$ such that $\sigma(i) = i$, we get $\mathrm{perm}(m) = \sum_{i = 0}^n ([x^i]\prod_{j=0}^{n-1} (1+f(j)x))M_{n-i}$
where $M_i$ is the number of [derangements](https://en.wikipedia.org/wiki/Derangement) of $i$ items.
This can be calculated in $O(n^2)$. In addition, it is possible to calculate the permanents of all the principal submatrices of $m$ without changing the overall complexity.

