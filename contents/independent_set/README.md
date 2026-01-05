# Number of independent sets(in grid graphs)
There are many configurations on this problem.

Graph types:
 - square grid ($P_n \times P_m$)
 - square grid with a diagonal
 - square grid with a diagonal, alternating by rows
 - king graph(square grid with both diagonals)
 - etc...

Topology:  
 - Plain
 - Cylinder (a pair of edge sticked together)
 - Torus (two pairs of edge sticked together)

Counting:
 - Total number of independent sets
 - Number of maximum/maximal independent sets
 - Number of independent sets of specific size $k$


Below are configurations I worked on.

## Square grid
#### Plain
A006506(n): Number of independent sets in $n \times n$ square grid graph ($P_n \times P_n$).  
A089980(n, m): Number of independent sets in $n \times m$ square grid graph ($P_n \times P_m$).  

Status: Computation ongoing for $n \le 55, m \le 200$  

#### Cylinder
A212270(n): Number of independent sets in $n \times n$ cylindrical square grid graph $C_n \times P_n$.  
A286513(n, m): Number of independent sets in cylindrical square grid graph $C_n \times P_m$.  

#### Torus
A027683(n): Number of independent vertex sets of the n X n torus grid graph.  

Status: Primary computation completed for $n \le 35, m \le 100$.  
Implementing verification code(Part of results are unreliable due to RAM errors on one of the machies)  

**specific size**  
A232833(n, k): of size $k$ on $n \time n$ plain grid  
A201511(n): of size $n$ on $n \times n$ plain grid  
A201626(n): of size $n$ on $n \times n$ torus  
