## Polynomial-time algorithm for the calculation of OEIS A137432
The sequence [A137432](http://oeis.org/A137432) represents the number of ways to place $n^2$ non-attacking kings on a $2n \times 2n$ cylindrical chessboard.  
Previously, it had the `hard` tag, indicating that the calculation of the sequence is hard, and only the first 31 terms were calculated.  
This document explains an algorithm that calculates the $n$-th term of the sequence in $O(n^4)$ multiprecision operations.  

## Implementation
An implementation of the algorithm in C++ is provided as `main.cpp` in the same directory. libboost is needed to compile it.

Performance(i5-10600k, Windows 10) :

|  n   |  Time   |
| ---- | ------- |
|  30  | 19.1 ms |
|  50  |  271 ms |
| 100  |  7.13 s |
| 200  |  159  s |

## Preliminaries
Let $(i, j) (0 \le i, j \lt n)$ denote the square at the row $i$ and column $j$.  
By "cylindrical chessboard", we consider the squares $(i, n - 1)$ and $(i, 0)$ are adjacent sharing an edge for every $0 \le i \lt n$. Consequently, two kings at $(1, n - 1)$ and $(2, 0)$ are considered mutually-attacking, and ones at $(0, 0)$ and $(n - 1, 0)$ are considered non-attacking, for example.  
Also, we denote by $a \hspace{5pt}\mathrm{mod}\hspace{5pt} b$ the remainder of $a$ divided by $b$.

## 1. Each king has only 4 possible locations
If we divide the $2n \times 2n$ board into $2 \times 2$ cells,  

 - Each cell can contain only one king
 - There are only $n^2$ cells

This leads to the fact that each cell must contain exactly one king. Let king $(i, j)$ denote the king at the cell containing squares $(2i + k, 2j + l) (0 \le k, l \le 1)$.
Also, let $s_{i, j}, t_{i, j}$ denote the row, column offset of king $(i, j)$, respectively. That is, king $(i, j)$ is at $(2i + s_{i, j}, 2j + t_{i, j})$.

## 2. Transform to a sequence counting problem
Considering the horizontal adjacency, we get $s_{i, j_1} = s_{i, j_2}$ for any $i, j_1, j_2$. (Remember the board is connected cyclically along the horizontal axis)  
We denote $b_i = s_{i, 0} = s_{i, 1} = \cdots = s_{i, n - 1}$. $b_i$ is either $0$ or $1$.  
Considering the vertical adjacency, we get $(0 \le ) \\, t_{0, j} \le t_{1, j} \le \cdots \le t_{n - 1, j} \\, (\le 1)$ for every $j$.  
We denote $c_j = \\#\\{ i \vert t_{i, j} = 0\\}$. $0 \le c_j \le n$ holds for every $j$, and $t_{i, j}$ equals $0$ if $i \lt c_j$ and $1$ if $i \ge c_j$.  

Considering the orthogonal adjacency, we get the following conditions on $b_i, c_j$ as the necessary and sufficcient condition.  

 - If $c_j < c_{(j + 1) \hspace{4pt}\mathrm{mod}\hspace{4pt} n}$, $b_i \le b_{i + 1}$ for every $c_j \le i \lt c_{(j + 1) \hspace{4pt}\mathrm{mod}\hspace{4pt} n} - 1$
 - If $c_j > c_{(j + 1) \hspace{4pt}\mathrm{mod}\hspace{4pt} n}$, $b_i \ge b_{i + 1}$ for every $c_{(j + 1) \hspace{4pt}\mathrm{mod}\hspace{4pt} n} \le i \lt c_j - 1$

Now, all we need is to count the number of $b$ and $c$ satisfying the above conditions.

## 3. Dynamic programming
The below figure represents an example of $b$ and $c$. The binary sequence at the top is $b$, and the location of vertices represents $c$ from the top to bottom.  
In this example, $b$ is $[0, 1, 1, 0, 0, 0, 0, 1, 0]$, and $c$ is $[3, 3, 8, 7, 5, 7, 8, 7, 5, 1]$.(The last vertex at the bottom is just to make the vertices cyclically-connected)  
The inequalities along the lines represent the condition on $b$ imposed by the two (cyclically) adjacent elements in $c$ at the end of the line.  

![](https://github.com/windows-server-2003/OEIS_calculation/blob/master/contents/A137432/0.png)

We can see that 

 - $(b_i, b_{i + 1})$ can be $(0, 1)$ iff every horizontal line crossing the $(i + 1)$-th vertical line from the right to left has one or more vertices on the vertical line.
 - $(b_i, b_{i + 1})$ can be $(1, 0)$ iff every horizontal line crossing the $(i + 1)$-th vertical line from the left to right has one or more vertices on the vertical line.
 - $(b_i, b_{i + 1})$ can always be $(0, 0)$ and $(1, 1)$

From this, we can do a DP(dynamic programming), scanning the vertical lines from the left to right and managing the number of lines crossing the current vertical line.
We have to hold the following values as the state of DP:

 - $i$, indicating we are currently at the vertical space between $i$-th and $(i + 1)$-th vertical lines(counted from left)
 - The number of line crossing the space per direction(the right to left or the left to right, which are always the same)
 - The number of vertices created(we have to create exactly $n - 1$ vertices until we reach the right end, excluding the start points at the top and the bottom)
 - The value of $b_i$.
 - Whether we have already passed the start point(3 in the example). (We can freely decide the start point, so we can at any point flip this value from 0 to 1)

A transition corresponds to advancing to the next right vertical space, passing a vertical line.
In a transition we do the following, all of which we can do by multiplying some binomial coefficient:

 - "Close" zero or more pairs of right-to-left line and left-to-right line  
     The left-to-right lines and the right-to-left lines are always sequenced alternately from the top to bottom, and we can decide which one is at the top from the state "Whether we have already passed the start point"
     We should create one or more vertices for this(=increase "The number of vertices created" in the dp key).
 - Add vertices to zero or more lines crossing the vertical line(we call this "knotting")
     We also allow transition to a different $b_i$ value, multiplying the number of ways to provide at least one vertices for every left-to-right/right-to-left(depending on the previous $b_i$) crossing lines.
 - "Open" zero or more pairs of left-to-right line and right-to-left line  
     Just like "closing", we can decide the number of places we can insert new pairs from "The number of line crossing" and "Whether we have already passed the start point" information.
     Again, we should create one or more vertices for this.

When we decide to pass the starting point, we have to try $2 \times 2$ times for which direciton(left or right) it goes first from the starting point / from which direction it returns to the starting point at the bottom, and manipulate "the number of line crossing per direction" by -1, 0, or 1.

This DP assumes that we eventually go either left or right after departing the start point, which misses one specific case: $c = [i, i, \dots, i]$ for some $0 \le i \le n$.
For the corresponding $b$, any binary sequence is allowed, so we have to add $(n + 1) \times 2^n$ to the answer.  
Precalculating the binomial coefficients, transitions from one DP state runs in $O(n^2)$ and the whole DP runs in $O(n^5)$. (ignoring the complexity of the multiprecision integer).  

If we do not allow adding redundant vertices(that is, fobid two vertices directly connected on the same vertical line or in other words forbid same values in a row in $c$) and
instead, after the DP finishes, multiply the number of ways to distribute the remaining vertices to the existing vertices(which is $n - 1$ minus the number of remaining vertices plus two because the starting points can also have redundant points) at the end of the DP,
transitions from one state runs in $O(n)$ and the total complexity will be $O(n^4)$. (Again, assuming $O(1)$ multiprecision operations)

Also, there is a DP key "the value of $b_i$", but it can be proved that regardless of this key the DP value is the same.
Therefore this key can be removed, resulting in 2x speedup. (`main.cpp` implements this)

With a little more observation, it is possible to get the answer for all $n' \le n$ with the same complexity, which is implemented in `main.cpp`.

