# SmithWaterman
Smith Waterman Algorithm for multiple sequence alignment
 This algorithm works only with nucleotide sequences. Only A,G,C,T characters are allowed.
 It directly returns on terminal the best alignments between the two inserted sequences.
 Optionally it can save on a text file all the other alignments between the two sequences. (see the options below)
 -----------------------
 Technical Info : 
  - This procedure implements a simple scoring schema, the user has the possibility to modify the values related to 
  		+ match score
 		+ mismatch score
 		+ gap penalty
  - Gap Score : the gap score (W()) is calulated using a simple formula, it is directly proportional to the number of gaps : gap_penalty * n_gaps 
  - This scoring schema is then used to build the scoring matrix, by using the flag -d the user can print it on the terminal.
  - Each cell of the scoring matrix is filled up with the max value between :
 		+ H[i-1,j-1]+S(a[i],b[j])
 		+ max {H[i,j-l]-W(l)}, l>=1
 		+ max {H[i-k,j]-W(k)}, k>=1
 		+ 0 					   
    Where W() is the gap score, and H is the scoring matrix that the procedure is filling up
  - Once the score matrix is created, the traceback procedure is used to find all the alignments
  - If this procedure is launched in standard mode (without flags), it returns the best alignment(s) with the score, the length,
     the number of gaps, the number of mismatches and the number of matches.
  - In the printed alignment, | stands for mismatch, * stands for match and ' ' stands for a gap.
  - Formatted results equal to  http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Smith-Waterman
  
============================================================================================================================,
