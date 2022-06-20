#This function uses a dynamic programming algorithm to calculate the optimal
#alignment score between two strings using gap, mismatch, and match scores
#Parameters: s1 and s2 are two strings (do not have to be the same length)
#            gap, mismatch, and match are all integers, gap and mismatch are
#            typically negative values, match is typically a positive value
#Returns: the optimal global alignment score, as we as a 2-D list of hints
#         to help determine the actual alignment
def NeedlemanWunsch(s1, s2, gap, mismatch, match):
    m = len(s1)
    n = len(s2)

    #Initialize 2 2-d lists of size (m+1)rows * (n+1)columns
    matrix = [[0 for i in range(n+1)] for j in range(m+1)]  #scores
    hints = [[0 for i in range(n+1)] for j in range(m+1)]   #backtracking hints
    
    #Below is the pseudocode for computing the alignment score table (matrix)
    #Once you get this working, add in code to keep track of the hints
    #You'll see in get_alignment, that it is expecting you to present hints as:
    # \\ (diagonal), | (from above), - (from left)
    #score1 is the diagonal, score2 is coming from the left, score3 if coming from above

    matrix[0][0] = 0
    for i in range(1, m+1):
        matrix[i][0] = matrix[i-1][0] + gap
        hints[i][0] = "|"

    for j in range(1, n+1):
        matrix[0][j] = matrix[0][j-1] + gap
        hints[0][j] = "-"

    for i in range(1, m+1):
        for j in range(1, n+1):
            if s1[i-1] == s2[j-1]:
                score1 = matrix[i-1][j-1] + match
            else:
                score1 = matrix[i-1][j-1] + mismatch

            score2 = matrix[i][j-1] + gap
            score3 = matrix[i-1][j] + gap
            matrix[i][j] = max(score1, score2, score3)
            if matrix[i][j] == score1:
                hints[i][j] = "\\"
            elif matrix[i][j] == score2:
                hints[i][j] = "-"
            elif matrix[i][j] == score3:
                hints[i][j] = "|"
            
    '''
    set matrix[0, 0] to 0
    for each i from 1 to m, inclusive
       matrix[i, 0] = matrix[i-1, 0] + gap

    for each j from 1 to n
       matrix[0, j] = matrix[0, j-1] + gap

    for each i from 1 to m, inclusive
        for each j from 1 to n, inclusive
            if s1[i-1] equals s2[j-1]
                score1 = matrix[i-1,j-1] + match
            else:
                score1 = matrix[i-1,j-1] + mismatch

            score2 = matrix[i][j-1] + gap
            score3 = matrix[i-1][j] + gap
            matrix[i,j] = max(score1, score2, score3)


    '''
    #for i in matrix:
        #print(i)
    return matrix[m][n], hints

#This function takes in 2 strings and a 2-D list of backtracking hints
#and returns the global alignment of the strings
#Parameters: A and B are two strings (do not have to be the same length)
#            hints is a 2-D lists containing \\, -, or | symbols
#Returns: two strings describing the global alignment between A and B

def get_alignment(A, B, hints):
    i = len(A)
    j = len(B)
    align1 = ''
    align2 = ''
    while i > 0 or j > 0:
        if hints[i][j] == '\\':
            align1 = A[i-1] + align1
            align2 = B[j-1] + align2
            i -= 1
            j -= 1
        elif hints[i][j] == '|':
            align1 = A[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = B[j-1] + align2
            j -= 1
            
    return align1, align2

            
def main():
    #feel free to hard code these while you're testing your code
    s1 = input("String 1: ")
    s2 = input("String 2: ")

    #you can get these are user inputs once the code is working
    gap_score = -2
    mismatch_score = -1
    match_score = 2
    
    score, hints = NeedlemanWunsch(s1, s2, gap_score, mismatch_score, match_score)
    print("Global alignment score is", score)
    
    align1, align2 = get_alignment(s1, s2, hints)
    print(align1)
    print(align2)

main()
