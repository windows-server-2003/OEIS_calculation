### Permanent of matrix represented by sum of outer products
### Overview
Let $A$ be a $n \times n$ matrix of rank at most $k$.  
This document explains an $O(n^k)$ algorithm that calculates the permanents of all the principal submatrices of $A$ (including $A$ itself).  
Here, the time complexity is measured in the number of arithmetic operations.  

The method is due to **Barvinok, Alexander I. (1996). Two Algorithmic Results for the Traveling Salesman Problem. Mathematics of Operations Research, Vol. 21, No. 1, pp. 65–84.**.  

### Method($k=2$ case)
This section just explains the same thing as "General Case" in the case of $k=2$ purely for better understanding.  

$A$ can be represented as $A_{i, j} = a_ib_j + c_id_j$.  
The permanent is equal to the sum calculated as follows
 - choose permutation $\sigma$ of $[0, 1, \dots, n-1]$, and for each possible choice:
 - add $\prod_{i=0}^{n-1} (a_ib_{\sigma(i)} + c_id_{\sigma(i)})$

This can be rephrased as the sum calculated as follows
 - choose permutation $\sigma$ of $[0, 1, \dots, n-1]$, and for each possible choice:
 - for each $i$, choose which of $a_ib_{\sigma(i)}$ and $c_id_{\sigma(i)}$ contributes to the product, and for each possible choices:
 - add the product of the $n$ values chosen above

Let's flip the order of choice: first we choose the following 
 - for each row $i$ choose whether $a_i$ or $c_i$ contributes to the product
 - for each column $j$ choose whether $b_j$ or $d_j$ contributes to the product

then we choose the $\sigma$ that can realize the above choice. We see that
 - the number of $a_i$ choice and the number of $b_j$ choice must match(so does $c_i$ and $d_j$)
 - when matches, the rows with $a_i$ choice must map to the columns with $b_j$ by $\sigma$. Therefore, if the number of $a_i$ choice is $m$, there are $m! (n-m)!$ possible choice of $\sigma$.

Consequently, $\mathrm{per}(A)$ is the sum of the following over $m = 0, 1, \dots, n$:
$m!(n-m)! \times ([x^m] \prod_{i=0}^{n-1} (a_ix + c_i)) \times ([x^m] \prod_{i=0}^{n-1} (b_ix + d_i))$  
where $[x^m] P(x)$ is the $x^m$ coefficient of the polynomial $P(x)$.  

By iteratively calculating the partial products (e.g. $\prod_{i=0}^{n'-1} (a_ix + c_i) \hspace{7pt} (n' \le n)$), we can calculate the solutions for all $n' \le n$ in $O(n^2)$.

### Method(General case)
$\mathrm{per}(A)$  
$= \sum_{\sigma \in S_n} \prod_{i=0}^{n-1} \sum_{l = 0}^{k-1} a_{l,i} b_{l,\sigma(i)}$  
$= \sum_{\sigma \in S_n} \sum_{\lambda \in \\{0, \dots, k-1\\}^n} \prod_{i=0}^{n-1} a_{\lambda_i,i} b_{\lambda_i,\sigma(i)}$  
$= \sum_{\sigma \in S_n} \sum_{\lambda \in \\{0, \dots, k-1\\}^n} \prod_{i=0}^{n-1} a_{\lambda_i,i} b_{\lambda_{\sigma(i)},i}$  
$= \sum_{\lambda \in \\{0, \dots, k-1\\}^n} \sum_{\sigma \in S_n}  \prod_{i=0}^{n-1} a_{\lambda_i,i} b_{\lambda_{\sigma(i)},i}$  

Let $f(c_0, \dots, c_{k-1})$ be the set of $\lambda \in \\{0, \dots, k-1\\}^n$ such that $l$ appears exactly $c_l$ times in $\lambda$ (for every $l$). 
We can see that the above expression is equal to

$\sum_{c_0, \dots, c_{k-2}, c_{k-1}} c_0!c_1!\cdots c_{k-1}!\sum_{\lambda, \mu: f(c_0, \dots, c_{k-1})} \prod_{i=0}^{n-1} a_{\lambda_i,i} b_{\mu_i,i}$  
$=\sum_{c_0, \dots, c_{k-2}, c_{k-1}} c_0!c_1!\cdots c_{k-1}! \cdot \left(\sum_{\lambda: f(c_0, \dots, c_{k-1})} \prod_{i=0}^{n-1} a_{\lambda_i,i}\right) \cdot  \left(\sum_{\mu: f(c_0, \dots, c_{k-1})} \prod_{i=0}^{n-1} b_{\mu_i,i}\right)$  

because there are exactly $c_0!c_1!\cdots c_{k-1}! \hspace{10pt} \sigma$'s that map a specific $\lambda$ to a specific $\mu$ (if $\lambda$ and $\mu$ are both in $f(c_0, \dots, c_{k-1})$ ).  
Now, $\left(\sum_{\lambda: f(c_0, \dots, c_{k-1})} \prod_{i=0}^{n-1} a_{\lambda_i,i}\right)$ (or the equivalent $b$ version) can be calculated in $O(n^k)$ for every $c_0, \dots, c_{k-1}$ because it is
given by the following coefficient of polynomial and it can be calculated by multiplying polynomials in order $i = 0, 1, \dots$.  
$[x_1^{c_1} x_2^{c_2} \dots x_{k-1}^{c_{k-1}}] \prod_{i=0}^{n-1} (a_{0,i} + a_{1,i}x_1 + a_{2,i}x_2 + \dots + a_{k-1,i}x_{k-1})$  

In the process of calculating the above product of polynomials, we can get all the partial product (i.e. $\prod_{i=0}^{n'-1} (a_{0,i} + a_{1,i}x_1 + a_{2,i}x_2 + \dots + a_{k-1,i}x_{k-1})$ for all $n'$).
Therefore we can get the answers for all principal submatrices without changing the asymptotic time complexity.

## Reference
Alexander I. Barvinok, (1996) Two Algorithmic Results for the Traveling Salesman Problem. Mathematics of Operations Research 21(1):65-84. 

