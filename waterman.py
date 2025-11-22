import numpy as np

def waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)

    #Create our empty matrix
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    max_score = 0
    max_position = (0,0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up_score = score_matrix[i-1][j] + gap
            
