
#!/usr/bin/env python3

''' Walter Endrizzi

	Python implementation of Smith-Waterman algorithm for local alignment
	testing https://gtuckerkellogg.github.io/pairwise/demo/	
'''

import argparse    	
from argparse import RawTextHelpFormatter			#for creating a decent user-friendly menu
import os
import sys
from copy import copy, deepcopy 					#used to preserve the score_matrix during the traceback procedure
import re                      						#regex



def main():
	
	
	seq_r, seq_c, bool_override_tb, match, mismatch, gap, bool_save, bool_display_mat, min_length, min_score, bool_min_gap = parseCmdLine()
	print("")
	print(" Options ")
	print("------ ")
	if (bool_display_mat) : print("+ Score matrix display")
	if (bool_override_tb) : print("+ Allowed overlaps in backtrace procedure")
	if (bool_save) : print("+ Export results to result.txt")
	if (bool_min_gap) : print("+ Export only the alignments with the minimum number of gaps")
	if min_length > 1 : print("+ Export only the alignments with length >= "+str(min_length))
	if min_score > 0 : print("+ Export only the alignment with score >= "+str(min_score))
	print("------ ")
	print("")
	print(" Analyzed Sequences ")
	print("------ ")
	print(" seq_c: "+seq_c + " --length: "+ str(len(seq_c)))
	print(" seq_r: "+seq_r + " --length: "+ str(len(seq_r)))
	print("------ ")
	print("")
	print(" Scoring Scheme ")
	print("------ ")
	print(" Match Score    : " + str(match))
	print(" Mismatch Score :" + str(mismatch))
	print(" Gap Penalty    : " + str(gap))
	print("------ ")
	print("")
	print(" Generating the Scoring Matrix..")
	print("")

	score_matrix, max_score, max_pos  = generateScoreMatrix(seq_r,seq_c,match,mismatch,gap)

	if (bool_display_mat):
		print(" Score Matrix:")
		print("   ="+((len(seq_c)*4))*"=")
		#print('    '+'   '+'   '.join(str(seq_c)))

		print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in score_matrix]))
		print("   ="+((len(seq_c)*4))*"=")
		print("")
	
	
	''' traceback part: from the highest score 
	follow back the track that gives you the highest score
	until you reach a score of zero'''

	'''Here i iterate until there is not a maximum in the score_matrix
	at each iteration the traceback function is 
	run in a recursive way '''
	
	all_traces = []
	trace = []
	


	while max_pos:
		trace.append(max_pos)
		trace, max_pos = traceBack(score_matrix,max_pos, trace, gap, match, mismatch, seq_c, seq_r)
		if trace:
			all_traces.append([trace, max_score])
			trace = []
			max_pos, max_score = find_new_pos(all_traces,score_matrix,bool_override_tb)
		else:
			max_pos = None



	'''Now given all the traces in the list all_traces I can return the 
	real alignments with the characters and the relative scores:
	main_list = [["ACTAG","A-TAC",score],["ACGT","-ACG",score],...["","",score],..]
	
	'''
	main_list = getAlignments(all_traces,seq_r,seq_c)

	

	'''Filtering the results'''
	filtered_main_list = applyFilters(main_list, bool_min_gap, min_length, min_score, len(seq_r)) if bool_min_gap or min_length > 1 or min_score > 0 else main_list
	


	'''multiple best alignments ? '''
	i = 0
	for al in filtered_main_list:
		if filtered_main_list[0][2] == al[2]: 
			i = i +1 
		else: 
			break
	


	
	if i > 1 : 
		print(" Multiple alignments found with the same maximum score : \n")
		j = 0
		while j < i:
			print(" Alignment "+str(j+1))
			printInfo(filtered_main_list[j],filtered_main_list[j][2])
			j+=1
	else:
		print(" Best Alignment ")
		printInfo(filtered_main_list[0],filtered_main_list[0][2])

	'''Formatted print to result.txt '''
	if bool_save:
		export_res(filtered_main_list)
		print(" Exporting results to 'result.txt'..")
		print(" Use 'cat result.txt' to display all the alignments. ")
	else:
		print(" Results have not been exported.")
		if i > 1 : print(" Please use the flag '-f' to save the result and see all the other alignments.\n")
		

	print(" Done.")
	

def parseCmdLine():
	
	'''Parse the command line arguments.
    Create a help menu, take input from the command line, and validate the input
    '''


	parser = argparse.ArgumentParser(
			description=
			" ===================================== Smith Waterman Algorithm Implementation ==============================================\n"
			" This algorithm works only with nucleotide sequences. Only A,G,C,T characters are allowed.\n"
			" It directly returns on terminal the best alignments between the two inserted sequences.\n"
			" Optionally it can save on a text file all the other alignments between the two sequences. (see the options below)\n"
			" -----------------------\n"
			" Technical Info : \n"
			"  - This procedure implements a simple scoring schema, the user has the possibility to modify the values related to \n"
			"  		+ match score\n"
			" 		+ mismatch score\n"
			" 		+ gap penalty\n"
			"  - Gap Score : the gap score (W()) is calulated using a simple formula, it is directly proportional to the number of gaps : gap_penalty * n_gaps \n"
			"  - This scoring schema is then used to build the scoring matrix, by using the flag -d the user can print it on the terminal.\n"
			"  - Each cell of the scoring matrix is filled up with the max value between :\n"
			" 		+ H[i-1,j-1]+S(a[i],b[j])\n"
			" 		+ max {H[i,j-l]-W(l)}, l>=1\n "
			" 		+ max {H[i-k,j]-W(k)}, k>=1\n "
			" 		+ 0 					   \n "
			"    Where W() is the gap score, and H is the scoring matrix that the procedure is filling up\n"
			"  - Once the score matrix is created, the traceback procedure is used to find all the alignments\n"
			"  - If this procedure is launched in standard mode (without flags), it returns the best alignment(s) with the score, the length,\n"
			"     the number of gaps, the number of mismatches and the number of matches.\n"
			"  - In the printed alignment, | stands for mismatch, * stands for match and ' ' stands for a gap.\n"
			"  - Formatted results equal to  http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Smith-Waterman\n "
			"============================================================================================================================\n",
			formatter_class=RawTextHelpFormatter)


	parser.add_argument("seq_c1", help=" Input sequence 1 (on the columns)")
	parser.add_argument("seq_r2", help=" Input sequence 2 (on the rows)")
	parser.add_argument("-f", "--save-to-file", action="store_true", help= "Export all the alignments in a text file named 'result.txt'. Alignments are reported in an ordered way with respect to the score.\n"
																		  "Note that if result.txt already exists it will be overwritten. Use this flag if you want information about all the alignments (not only the best).\n"
																			  "   ")
	parser.add_argument("-d", "--display-scoring-matrix", action="store_true", help="Display the scoring matrix calculated by the algorithm using the scoring schema\n"
																			  "   ")
	parser.add_argument("-o", "--override-in-traceback", action="store_true", help="Allow overlaps in trace-back paths, by default they are not allowed.\n"
																			 	   "   ")
	parser.add_argument("-m", "--match" ,	 type = int, choices=[0,1,2,3,4,5,6,7,8,9], 		default = 3 ,help="Optional. set up the match score. (default = 3) \n"
																			  									  "   ")
	parser.add_argument("-s", "--mismatch" , type = int, choices=[-9,-8,-7,-6,-5,-4,-3,-2,-1,0], default =-1, help="Optional. set up the mismatch score. (default = -1) \n"
																			  										"   ")
	parser.add_argument("-g", "--gap" ,		 type = int, choices=[0,1,2,3,4,5,6,7,8,9], 		default = 1, help="Optional. set up the gap score. (default = 1) \n"
																			  									  "   ")
	parser.add_argument("-ml", "--minimum-length", type = int, default = 1, help="Optional. If -ml is specified, only the alignments with (length > ml) will be printed (or exported if -f is used)\n"
																			  "   ")
	parser.add_argument("-ms", "--minimum-score", type=int, default = 0, help="Optional. If -ms is specified, only the alignments with (score > ms) will be printed (and exported to the file result.txt if -f is used)\n"
																			  "   ")
	parser.add_argument("-mg", "--seq-with-min-number-of-gaps", action="store_true", help="Optional. If this flag is used, only the alignments with the minimum number of gaps are selected and printed (or exported if -f is used). \n"
																				" ")
	

	'''
	return all the possible
	local alignments, ordered by decreasing score, with length > 5, score > 4
	and gaps > 0. Gaps are included in the calculation of alignment
	length.
	'''
	min_length = 1 
	min_score = 1

	args = parser.parse_args()
	
	seq_c = args.seq_c1
	seq_r = args.seq_r2
	match = args.match
	mismatch = args.mismatch
	gap = args.gap
	bool_save = args.save_to_file
	bool_display_mat = args.display_scoring_matrix
	bool_override_tb = args.override_in_traceback
	min_length = args.minimum_length
	min_score = args.minimum_score
	bool_min_gap = args.seq_with_min_number_of_gaps

	#check string validity
	if(len(seq_c) <= 1 or len(seq_r) <= 1):
		print(' Sequences are too short')
		sys.exit()

	if any(c not in 'AGCT' for c in seq_c):  
		print(' Wrong character in sequence 1, see the help with -h.')
		sys.exit()
	

	if any(c not in 'AGCT' for c in seq_r):  
		print(' Wrong character in sequence 2, see the help with -h.')
		sys.exit()

	return seq_r, seq_c, bool_override_tb, match, mismatch, gap, bool_save, bool_display_mat, min_length, min_score, bool_min_gap
	

def generateScoreMatrix(seq_r,seq_c,match,mismatch,gap):
	
	cols, rows, score_matrix = initScoreMatrix(seq_r,seq_c)

	max_score = 0
	max_pos = 0
	for i in range(1,rows):
		for j in range(1,cols):
			
			score = getMaxScore(score_matrix,i,j,match,mismatch,gap,seq_r,seq_c)
			score_matrix[i][j]=score

			if (score>=max_score):
				max_score = score
				max_pos = (i,j)


	return score_matrix, max_score, max_pos


def initScoreMatrix(seq_r,seq_c):

	'''The scoring matrix has dimension [1+length_seqa][1+length_seqb]
		all the element in the first row and columns are 0'''


	cols = len(seq_c)+1
	rows = len(seq_r)+1


	score_matrix = [[0 for row in range(cols)] for col in range(rows)]


	return cols, rows , score_matrix


def getMaxScore(score_matrix, row, col,match,mismatch,gap,seq_r,seq_c):
	''' Function used to fill the score_matrix it checks which is the
	best score to insert in a certain cell of the matrix
	
		- the diagonal : H[i-1,j-1]+S(a[i],b[j])
		- right shift : max{H[i,j-l]-W(l)}, l>=1
		- vertical shift : max{H[i-k,j]-W(k)}, k>=1
		- 0
		
		W(k) in our case is linear with the number of gaps=> gap_penalty*n_gaps
		n_gaps simply is the distance (in number of cell orizontally or verically)
		from the current cell

		Similarity uses the value of the variables match and mismatch to score the current cell.
	'''

	supp_oriz = []
	supp_vert = []
	
	#here i iterate from m = 1 until m = col-1.
	for m in range(1,col):
		supp_oriz.append(score_matrix[row][m]-gap*(col-m))
	for n in range(1,row):
		supp_vert.append(score_matrix[n][col]-gap*(row-n))

	similarity =  match if seq_c[col-1] == seq_r[row-1] else mismatch

	
	diag = score_matrix[row-1][col-1]+similarity
	left = 0 if supp_oriz == [] else max(supp_oriz)
	up = 0 if supp_vert == [] else max(supp_vert)
	
	#simpler scoring schema, here I don't check if there is a cell in all the row or in all the column, i simply check the first after the current cell.
	'''
	similarity =  match if seq_c[col-1] == seq_r[row-1] else mismatch


	diag = score_matrix[row-1][col-1]+similarity
	left = score_matrix[row][col-1]-gap
	up = score_matrix[row-1][col]-gap
	'''
	return max(diag,left,up,0)


def traceBack(score_matrix,max_pos,trace,gap, match, mismatch, seq_c, seq_r):

	'''function that implements the traceback, it returns the 
	alignment given a starting position given by max_pos
	
	   WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2

    In trace we have the path of the current alignment
    at each recursive call we add one matrix position to trace,
    until we reach a score = 0 in the matrix.
    I append the first position of the trace before the call to traceback
    Here below I define i,j that take the position of the current cell from which i will find the next best position
    '''
	
	i,j=max_pos
		
	#here I have to add the gap penality to move laterally or vertically?????? if y, same result as expected
	similarity =  match if seq_c[j-1] == seq_r[i-1] else mismatch
	diag = score_matrix[i-1][j-1] + similarity 
	up   = score_matrix[i-1][j] -gap
	left = score_matrix[i][j-1] -gap
	
	
	mymax = max(diag,up,left)

	#computing the new position in the matrix that will be saved 
	if (mymax==diag):
		max_pos = (i-1,j-1)
	elif (mymax==left):
		max_pos = (i,j-1)
	else:
		max_pos = (i-1,j)

	i,j=max_pos
		
	#we discard a trace if it contains only one cell
	#If I mymax==0 zero, then I am at the beginning of the alignment, I have to stop the recursion.
	if score_matrix[i][j] == 0:
		if len(trace) <= 1:
			trace = None
		return trace, max_pos
	else:
		trace.append(max_pos)
		return traceBack(score_matrix,max_pos,trace, gap, match, mismatch, seq_c, seq_r)


def find_new_pos(all_traces, score_matrix, bool_override_tb):

	'''given the score_matrix find the new maximum 
	score omitting the ones that we have already found'''

	'''here I set to 0 all the values that are already present
	in an alignment, then I search for a maximum'''

	'''this function returns only the index of the new max
	found in the score_matrix
	'''

	score_matrix_sup = deepcopy(score_matrix)

	new_max = 0
	new_pos = None

	'''allow override in traceback, I remove from the scoring matrix 
	only the starting points of the paths already found.
	
	If not allowed I remove from the scoring matrix also all the traces
	already found.
	'''

	if bool_override_tb :
		for l in all_traces:
			tup = l[0][0]
			m,n= tup
			score_matrix_sup[m][n]=0
	else:
		for l in all_traces:
			for tup in l[0]:
				m,n = tup
				score_matrix_sup[m][n]=0

			


	for i in range(len(score_matrix_sup)):
		for j in range(len(score_matrix_sup[i])):
			if score_matrix_sup[i][j] >= new_max:
				new_max = score_matrix_sup[i][j]
				new_pos = (i,j)


	return new_pos, new_max


def getAlignments(all_traces,seq_r,seq_c):
	'''
	The main idea here is that, if both the indexes change, then there's
	a match or a mismatch.
	If only the index of the rows changes then, I have to insert a gap in the 
	sequence2 (index row)
	If only the index of the columns changes, then I have to insert a gap in 
	the sequence 1 (seq_c) (index col)

	'''
	main_list = []
	for l in all_traces:
		supp_list_1 = ""
		supp_list_2 = ""
		row_old = 1000
		col_old = 1000

		#print(l)
		for tup in (reversed(l[0])):

			row,col=tup
			#print(str(row)+"  old: "+str(row_old))
			#print(str(col)+ " old: "+str(col_old))
			if (row == row_old):
				supp_list_2= supp_list_2 + '-' 
			else:
				supp_list_2= supp_list_2 + seq_r[row-1]

			if (col == col_old):
				supp_list_1=  supp_list_1 + '-'
			else:
				supp_list_1=  supp_list_1 +seq_c[col-1] 

			#print(supp_list_2)
			#print(supp_list_1)
			row_old = row
			col_old = col
		
		main_list.append(([supp_list_1,supp_list_2,l[1]])) 
	
	return main_list


def applyFilters(main_list, bool_min_gap, min_length, min_score, max_gap):
	print(" Filtering the alignments..")
	print(" ")
	filtered_main_list = []
	mingap = max_gap

	if bool_min_gap: #if specified by the user, I save only the alignment which have the minimum gap. (it can be = 0)
		for alignment in main_list:
			mingap = alignment[0].count('-')+alignment[1].count('-') if (alignment[0].count('-')+alignment[1].count('-')) < mingap else mingap


	for alignment in main_list:
		if (alignment[2] < min_score):	#save only the alignment which has a score greater than '-ms' (specified by user)
			continue
		if (len(alignment[1]) < min_length):	#save only alignment which has a length grater than 'ml' (specified by user)
			continue
		if (bool_min_gap and alignment[0].count('-')+alignment[1].count('-') > mingap):
			continue
		
		filtered_main_list.append(alignment)


	if (filtered_main_list == []):
		print(" No alignments found, please change the sequences (or some parameters)")
		sys.exit()

	return filtered_main_list


def printInfo(alignment, score_):
	''' 
	print the number of matches, the number of mismatches, 
	the number of gaps, the length of the alignment
	'''
	sequence_a = alignment[0]
	sequence_b = alignment[1]

	align_length = len(sequence_a)
	n_match = 0
	n_mismatch = 0
	n_gaps = 0

	'''
	string which make easier the visualization
	'''
	visual_string = ""

	for i in range(len(sequence_a)):
		if (sequence_a[i]==sequence_b[i]):
			n_match+=1
			visual_string+="*"
		elif (sequence_a[i] == "-" or sequence_b[i] == "-"):
			n_gaps+=1
			visual_string+=" "
		else:
			n_mismatch+=1
			visual_string+="|"

	print("--------------------------")
	print(" "+sequence_a)
	print(" "+visual_string)
	print(" "+sequence_b)
	print("--------------------------")
	print(" alignment score      : " + str(score_))
	print(" alignment length     : " + str(align_length))
	print(" number of gaps       : " + str(n_gaps))
	print(" number of matches    : " + str(n_match))
	print(" number of mismatches : " + str(n_mismatch))
	print("--------------------------")
	print("   ")
	print("   ")


def export_res(main_list):


	original_stdout = sys.stdout
	sys.stdout = open('result.txt', 'w')
	i = 0
	for alignment in main_list:
		print(" Alignment number "+str(i+1))
		printInfo(alignment, alignment[2])
		i=i+1

	sys.stdout = original_stdout



if __name__== '__main__':
	sys.exit(main())
	sys.stdout.close()
