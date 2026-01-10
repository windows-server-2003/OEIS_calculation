## Specific matrices whose permanent can be calculated efficiently
We will denote the order of the matrix by $N$, row index by $i$, and column index by $j$.  
"Complexity" means the complexity of the number of operations(addition, multiplication, etc.) on multiprecision integers and does not take in account the order of the growth of the numbers appearing in the calculation.  

### Low-rank matrix
 - Applicable condition : Matrix $A$ is of rank $k$ or less  
 - Complexity : $O(N^k)$ for calculating permanents of all the principal submatrices  
 - Method : See [this](Detail0.md)
 - List of sequences on OEIS to which this method applies :  

|                              OEIS                          |     $A[i][j]$       |   k   |  comment |
| ---------------------------------------------------------- | ------------------- | ----- | -------- |
| [A203264](https://oeis.org/A203264)                        | $\left(\left\lfloor \frac{i}{2} \right\rfloor + 1\right)(j+1)$ | 1 | |
| [A330087](https://oeis.org/A330087)                        | $(i+1) \times \mathrm{prime}(j)$ | 1 | |
| [A204248](https://oeis.org/A204248)                        |       $1+i+j$       |   2   | simple mathmatic formula already exists|
| [A204249](https://oeis.org/A204249)                        |       $2+i+j$       |   2   | |
| [A278927](https://oeis.org/A278927)                        |    $2i + j + 3$     |   2   | |
| [A346934](https://oeis.org/A346934)                        |      $i - j$        |   2   | |
| [A204251](https://oeis.org/A204251)                        |   $ij+2i+2j+1$      |   2   | |
| [A278847](https://oeis.org/A278847)                        | $(i+1)^2+(j+1)^2$   |   2   | |
| [A278925](https://oeis.org/A278925)                        | $(i+1)^3+(j+1)^3$   |   2   | |
| [A278926](https://oeis.org/A278926)                        | $(i+1)^4+(j+1)^4$   |   2   | |
| [A322277](https://oeis.org/A322277)                        |      Zigzag         |   2   | |
| [A278845](https://oeis.org/A278845)                        | $((i+1) + (j+1))^2$ |   3   | |
| [A278857](https://oeis.org/A278857)                        | $((i+1) - (j+1))^2$ |   3   | |


### "Low-rank per triangle" matrix
 - Applicable condition($k=2$) : Matrix $A$ formed by the lower triangular part of some rank-$2$ matrix $L$, upper triangular part of some rank-$2$ matrix $U$, and a diagonal $D$.  
 - Complexity : $O(N^4)$ for calculating permanents of all the principal submatrices  
 - Note: More generally, it can be solved in $O(N^{2k})$ if $L$ and $U$ are of rank-$k$, though it $k>2$ case does not seem to appear on OEIS and is not implemented either
 - List of sequences on OEIS to which this method applies :  
 
|                  OEIS               |              $A[i][j]$               |
| ----------------------------------- | ------------------------------------ |
| [A307783](https://oeis.org/A307783) |             $n-\|i-j\|$              |
| [A204235](https://oeis.org/A204235) |             $1+\|i-j\|$              |
| [A085807](https://oeis.org/A085807) |              $\|i-j\|$               |
| [A204236](https://oeis.org/A204236) |        $\max(2i-j, 2j-i)+1$          |
| [A204239](https://oeis.org/A204239) |        $\max(3i-j, 3j-i)+2$          |
| [A204241](https://oeis.org/A204241) |        $\max(3i-2j, 3j-2i)+1$        |
| [A085719](https://oeis.org/A085719) | $1 + ((i-j+n) \\,\mathrm{mod}\\, n)$ |
| [A086759](https://oeis.org/A086759) |   $(i-j+n) \\,\mathrm{mod}\\, n$   |
| [A322909](https://oeis.org/A322909) |                omitted               |
| [A323255](https://oeis.org/A323255) |                omitted               |
| [A204262](https://oeis.org/A204262) |          $1 + \min(i, j)$            |
| [A204264](https://oeis.org/A204264) |          $1 + \max(i, j)$            |
| [A204234](https://oeis.org/A204234) |    (i + j + 2)(1 + \min(i, j))       |
| [A347768](https://oeis.org/A347768) |            $|i - j + 1|$             |
| [A278858](https://oeis.org/A278858) |       $|(i+1)^2 - (j+1)^2|$          |


### All 1 except the diagonal
 - Complexity : $O(N^2)$ for calculating permanents of all the principal submatrices  
 - Applicable condition : Matrix $A$ of form: $A[i][i] = f(i)$, $a[i][j] = 1 (i \neq j)$
 - Method : See [this](Detail1.md)
 - Note: Special case of "low-rank per triangle" ($k=1$)
 - List of sequences on OEIS to which this method applies :  

|                  OEIS               |          $f(i)$        |
| ----------------------------------- | ---------------------- |
| [A303000](https://oeis.org/A278845) |        $(i+1)^2$       |
| [A303001](https://oeis.org/A278857) | $\frac{(i+1)(i+2)}{2}$ |


### Mod-dependent
 - Applicable condition : Matrix $A$ of form: $A[i][j] = f(i \\,\mathrm{mod}\\, k, j \\,\mathrm{mod}\\, k)$ for some $k$  
 - Complexity : $O(N^k)$ for calculating permanents of all the principal submatrices
 - Note: this is just a special case of low-rank matrices, but the implementation is easier and slightly faster
 - List of sequences on OEIS to which this method applies :  

|                 $k=3$               |                 $k=4$               |
| ----------------------------------- | ----------------------------------- |
| [A179079](https://oeis.org/A179079) | [A204256](https://oeis.org/A204256) |
| [A204254](https://oeis.org/A204254) | [A204422](https://oeis.org/A204422) |
| [A204258](https://oeis.org/A204258) | [A204442](https://oeis.org/A204442) |
| [A204265](https://oeis.org/A204265) | [A204444](https://oeis.org/A204444) |
| [A204268](https://oeis.org/A204268) | [A204446](https://oeis.org/A204446) |
| [A204424](https://oeis.org/A204424) | [A204448](https://oeis.org/A204448) |
| [A204426](https://oeis.org/A204426) |
| [A204428](https://oeis.org/A204428) |
| [A204430](https://oeis.org/A204430) |
| [A204432](https://oeis.org/A204432) |
| [A204434](https://oeis.org/A204434) |
| [A204436](https://oeis.org/A204436) |
| [A204438](https://oeis.org/A204438) |
| [A204440](https://oeis.org/A204440) |
