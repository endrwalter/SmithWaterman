#!/usr/bin/env python3

'''A Python implementation of the Smith-Waterman algorithm for local alignment
of nucleotide sequences.
'''


import argparse
import os
import sys


# simple scoring scheme 
match    =  2
mismatch = -1
gap      = -1
seq1     = None
seq2     = None



def main():
	try:
		parseCmdLine()
		print("seq1: "+seq1 + " --length: "+ str(len(seq1)))
		print("seq2: "+seq2 + " --length: "+ str(len(seq2)))
	except :
		print("error in parsing the command line..")


	score_matrix, max_score, max_pos  = generateScoreMatrix()

	print("Maximum score : "+str(max_score)+" in position "+str(max_pos))
	print("Score Matrix:")

	print('\n'.join([''.join(['{:4}'.format(item) for item in row]) 
      for row in score_matrix]))


	''' traceback part: from the highest score 
	follow back the track that gives you the highest score
	until you reach a score of zero'''

	'''Here i iterate until there is not a maximum in the score_matrix
	at each iteration the traceback function getBestAlignment is 
	run recursively '''
	
	all_alignments = []
	alignment = []

	while max_pos:
		alignment, max_pos = getBestAlignment(score_matrix,max_pos, alignment)
		if alignment:
			all_alignments.append(alignment)
			alignment = []
			max_pos = find_new_pos(score_matrix,all_alignments)
		else:
			max_pos=None


	print(all_alignments)

	#### I have to remove te alignments that contain only one element..


	'''Now given all the traces in the list all_alignments I can return the 
	real alignments with the characters'''

	''' 
	list of lists composed eah one by two aligned sequences

	main_list[["ACTAG","A-TAC"],["ACGT","-ACG"],...["",""],..]
	'''

	main_list=[] 

	for l in all_alignments:
		supp_list_a = ""
		supp_list_b = ""
		a_old = ''
		b_old = ''
		for tup in reversed(l):
			a,b=tup
			if (a == a_old):
				supp_list_a=supp_list_a+'-'
			else:
				supp_list_a= supp_list_a + seq1[a-1]

			if (b == b_old):
				supp_list_b= supp_list_b +'-'
			else:
				supp_list_b= supp_list_b + seq2[b-1]

			a_old = a
			b_old = b
		
		main_list.append(([supp_list_a,supp_list_b]))


	print(main_list)

	
	'''Formatted print for terminal '''





def parseCmdLine():

	'''Parse the command line arguments.
    Create a help menu, take input from the command line, and validate the
    input by ensuring it does not contain invalid characters (i.e. characters
    that aren't the bases A, C, G, or T).
    At the moment im using two constant sequences :
    '''

	global seq2
	global seq1
	seq1 = "ATATA"
	seq2 = "TATGT"


def generateScoreMatrix():

	'''Fill the scoring matrix iteratively using the schema (see getMaxScore):
		- the diagonal : H[i-1,j-1]+S(a[i],b[j])
		- right shift : max{H[i,j-l]-W(l)}, l>=1
		- vertical shift : max{H[i-k,j]-W(k)}, k>=1
		- 0
		This function uses the getMaxScore() function described below
	'''

	
	cols, rows, score_matrix = initScoreMatrix()


	max_score = 0
	for i in range(1,rows):
		for j in range(1,cols):
			
			score = getMaxScore(score_matrix,i,j)
			score_matrix[i][j]=score

			if (score>max_score):
				max_score = score
				max_pos = (i,j)


	return score_matrix, max_score, max_pos

def initScoreMatrix():
	'''function that initialize the score_matrix '''

	'''The scoring matrix has dimension [1+length_seqa][1+length_seqb]
		all the element in the first row and columns are 0'''


	rows = len(seq1)+1
	cols = len(seq2)+1


	score_matrix = [[0 for row in range(cols)] for col in range(rows)]


	return cols, rows , score_matrix

 


def getMaxScore(score_matrix, row, col):
	''' Function used to fill the score_matrix it checks which is the
	best score to insert in a certain cell of the matrix'''

	supp_oriz = []
	supp_vert = []

	for m in range(0,col):
		supp_oriz.append(score_matrix[row][m])
	for n in range(0,row):
		supp_vert.append(score_matrix[n][col])

	similarity =  match if seq1[row-1] == seq2[col-1] else mismatch

	
	diag = score_matrix[row-1][col-1]+similarity
	left = max(supp_oriz)+gap
	up = max(supp_vert)+gap

	#print("seq1 "+ str(seq1[row-1])+" seq2 "+str(seq2[col-1])+" "+str(max(diag,left,up,0)))

	return max(diag,left,up,0)


def getBestAlignment(score_matrix,max_pos,alignment):

	'''function that implements the traceback, it returns the 
	alignment given a starting position given by max_pos
	
	   WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2
    '''
	aligned_seq_a=""
	aligned_seq_b=""

	local_alignemnt=[]
	alignment.append(max_pos)

	i,j=max_pos
	if (i==j):
		aligned_seq_a=seq1[j-1]
		aligned_seq_b=seq2[i-1]
		
	#print(score_matrix[i][j])
	#print(alignment)
	diag = score_matrix[i-1][j-1]
	up = score_matrix[i-1][j]
	left = score_matrix[i][j-1]
	
	mymax = max(diag,up,left)
	
	if (mymax==diag):
		max_pos = (i-1,j-1)
		aligned_seq_a= seq1[j-1]+aligned_seq_a
		aligned_seq_b= seq2[i-1]+aligned_seq_b
	elif (mymax==left):
		aligned_seq_a= seq1[j-1]+aligned_seq_a
		aligned_seq_b= '-'+aligned_seq_b
		max_pos = (i,j-1)
		
	else:
		aligned_seq_a= '-'+aligned_seq_a
		aligned_seq_b= seq2[i-1]+aligned_seq_b
		max_pos = (i-1,j)
		
	
	if mymax==0:
		return alignment, max_pos
	else:
		return getBestAlignment(score_matrix,max_pos,alignment)

	
def find_new_pos(score_matrix,all_alignments):
	'''given the score_matrix find the new maximum 
	score omitting the ones that we have already found'''

	'''here I set to 0 all the values that are already present
	in an alignment, then I search for a maximum'''

	'''this function returns only the index of the new max
	found in the score_matrix
	'''

	new_max = 0
	new_pos = None
	for l in all_alignments:
		for tup in l:
			m,n = tup
			score_matrix[m][n]=0



	for i in range(len(score_matrix)):
		for j in range(len(score_matrix[i])):
			if score_matrix[i][j] > new_max:
				new_max = score_matrix[i][j]
				new_pos = (i,j)


	return new_pos




if __name__== '__main__':
	sys.exit(main())
	
