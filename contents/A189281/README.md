## A189281
This is the explanation of an $O(n^4)$ algorithm for calculating the $n$-th term of [A189281](https://oeis.org/A189281), which is the number of permutation $p$ of $[1, 2, 3, \dots, n]$ such that $p(i) \neq p(i + 2)$ holds for every $1 \le i \le n-2$.  

### Algorithm
We denote by $S(n)$ the set of all permutation of $[1, 2, 3, \dots, n]$.  
First, consider replacing both $i$ and $p_i$ by the rule $(1, 2, 3, \dots, n) \rightarrow (1, 1+\lfloor \frac{n+1}{2} \rfloor, 2, 2+\lfloor \frac{n+1}{2} \rfloor, 3, 3+\lfloor \frac{n+1}{2} \rfloor, \dots)$.  
For example, if $n = 8$, $[3, 5, 6, 2, 8, 7, 1, 4]$ becomes $[2, 7, 8, 1, 3, 5, 4, 6]$ by this conversion.  
This is a bijection from $S(n)$ to itself, and the original condition is equivalent to the following condition in the transformed permutation :

 - $p(i) + 1 \neq p(i + 1) (1 \le i \le n - 1)$ except that $p(i) + 1 = p(i+1)$ is allowed if $i = \lfloor \frac{n+1}{2} \rfloor$ or $p(i) = \lfloor \frac{n+1}{2} \rfloor$.

We will use a DP(dynamic programming) approach from here.  
In the figure below, which expresses the permutation $[4, 3, 1, 2]$, we call each of the red points a "gap" and number them $0, 1, 2, \dots$ from left to right.  

![](gap.png)

"Violating" gap is a gap between two numbers whose left number plus 1 equals the right number and the left number is not $\lfloor \frac{n+1}{2} \rfloor$ which is one of the two exceptions to the rule.  
Other gaps, including the gaps at the two ends of the sequence, are "nonviolating" gaps.  
In the example above, only gap $3$ is a violating gap.  

The main idea is composing $p$ by inserting $1, 2, 3, \dots, n$ to an empty sequence in order while maintaining the gap that eventually become the gap $\lfloor \frac{n+1}{2} \rfloor$ where $p(i) + 1 = p(i+1)$ is allowed.  

Let $\mathrm{dp}\_{i, j, k, l, m}$ the number of pairs of a permutation $p$ of $[1, 2, 3, \dots, i]$ and its gap $g$ such that :

 - $g$ is the gap $j$ of $p$
 - defining the gaps on the left of $g$ "left gaps" and the gaps on the right of $g$ "right gaps", 
 - there are $k$ violating left gaps
 - there are $l$ violating right gaps
 - the number $i$ is at the specified place depending on $m$ : 
    - if $m = 2$ : $i$ must be the left element of gap $j$
    - if $m = 1$ : $i$ must be on the left of the left element of gap $j$
    - if $m = 3$ : $i$ must be on the right of the left element of gap $j$

We can calculate $\mathrm{dp}\_{i+1, *, *, *, *}$ from $\mathrm{dp}\_{i, *, *, *, *}$ by considering inserting $i+1$ to an existing permutation of $[1, 2, 3, \dots, i]$.  
From $k, l$, we know how many violating and nonviolating, left and right gaps there are. From $m$, we (roughly) know where $i$ is.  
If we decide to insert $i + 1$ into a currently violating left gap, the number of violating left gap decreases by 1. (Note that a new violating gap cannot be created by the insertion in this case)  
If we decide to insert $i + 1$ into a currently nonviolating left gap, excluding the gap immediately before $i$(if it's a left gap), the number of violating left gap does not change.  
If we decide to insert $i + 1$ into a currently nonviolating left gap whose right element is $i$(if $i$ exists on the left of the left element of gap $j$), the number of violating left gap increases by 1. (However, when inserting $\lfloor \frac{n+1}{2} \rfloor + 1$, the number of violation does not change)  
Similar for inserting into a right gap.  
If we decide to insert $i + 1$ into gap $j$ when $m = 2$, the gap $j$ is split into one violating gap and one nonviolating gap. We try both of the gaps as the new $g$.  
If we decide to insert $i + 1$ into gap $j$ when $m \neq 1$, the gap $j$ is split into two nonviolating gaps. We try both of the gaps as the new $j$-th gap.  



Finally, the sum of $\mathrm{dp}\_{n, \lfloor \frac{n+1}{2} \rfloor, 0, 0, m}$ over all $m$ is the answer.  

Notes :
When transiting from $\mathrm{dp}\_{i, *, *, *, *}$ to $\mathrm{dp}\_{i+1, *, *, *, *}$, for any pair of $p$ and $g$, the previous $p$ and $g$ are uniquely determined, so there's no double-counting.  
Note that whether the gap $j$ is violating or not was not included in the information. This is because gap $i$ will eventually be the gap $\lfloor \frac{n+1}{2} \rfloor$ where violation is allowed, and when inserting into gap $j$, whether the two created gaps are violating or not does not depend on whether the original gap $j$ was violating.  

### Implementation
```main.cpp``` in this directory is a working C++ implementation of this algorithm using ```cpp_int``` in [boost](https://www.boost.org/) as the multiprecision integer type.  
After compiling with gcc, simply run it from command line. It will read $N$ from stdin and calculate $a(n)$ for $n = 0, 2, 3, \dots, N$.  

