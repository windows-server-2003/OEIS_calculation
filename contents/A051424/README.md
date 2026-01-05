## A051424(n) - Number of partition of n into pairwise coprime (positive) integers
Overview of algorithm is in the program comments.

The bottleneck is the memory usage.

Room for improvement:  
 - u16(or maybe u8) modint to reduce memory usage (currently u256)
 - Regarding the vectors in `dp`, only allocate memory for elements beyond first non-zero element  
   Experiment suggested this makes 4x memory saving. There may be zeros after first non-zero element but not that many.

