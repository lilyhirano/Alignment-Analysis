import numpy as np
import pandas as pd

def read_fasta(filename):
    with open(filename, 'r') as f:
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"): #ignore first line with label
                pass
            else:
                sequence += line
    return sequence

class SmithWaterman:

    def __init__(self, seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=2):
        """
        Performs the Smith-Waterman local alignment algorithm.
        Returns the start and end indices of the alignment in seq1.
        """

        self.seq1 = seq1
        self.seq2 = seq2
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty


        self.rows = len(seq1) + 1
        self.cols = len(seq2) + 1

        self.max_score = 0
        self.max_pos = (0, 0)

        self.fill_matrix()

        self.score = self.find_score
        self.alignment = self.traceback()


    def initialize_matrix(self):
        # Initialize the scoring matrix with zeros
        self.matrix = np.zeros((self.rows, self.cols), dtype = int)
    
    def fill_matrix(self):

        self.initialize_matrix()

        # Fill the scoring matrix
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                score_diag = self.matrix[i - 1, j - 1] + (self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty)
                score_up = self.matrix[i - 1, j] - self.gap_penalty
                score_left = self.matrix[i, j - 1] - self.gap_penalty
            
                self.matrix[i, j] = max(0, score_diag, score_up, score_left)

    def find_score(self):
        max_score = 0
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                if self.matrix[i, j] > max_score:
                    max_score = self.matrix[i, j]
                    self.max_pos = (i, j)

        return max_score

    def traceback(self):
        i, j = self.max_pos
        aligned_seq1, aligned_seq2 = "", ""
    
        while self.matrix[i, j] != 0:
            current_score = self.matrix[i, j]
            score_diag = self.matrix[i - 1, j - 1]
            score_up = self.matrix[i - 1, j]
            score_left = self.matrix[i, j - 1]

            # Prioritize diagonal move if it's the source of the current score
            if current_score == score_diag + (self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_penalty):
                aligned_seq1 += self.seq1[i - 1]
                aligned_seq2 += self.seq2[j - 1]
                i -= 1
                j -= 1
            elif current_score == score_up - self.gap_penalty:
                aligned_seq1 += self.seq1[i - 1]
                aligned_seq2 += "-" # Gap in seq2
                i -= 1
            elif current_score == score_left - self.gap_penalty:
                aligned_seq1 += "-" # Gap in seq1
                aligned_seq2 += self.seq2[j - 1]
                j -= 1
            else:
                # Should not reach here if logic is correct
                break
        
        traceback_coords = []
        while i > 0:
            traceback_coords.append((i, j))
            aligned_seq1 = self.seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1

        while j > 0:
            traceback_coords.append((i, j))
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = self.seq2[j - 1] + aligned_seq2
            j -= 1

        self.traceback_coords = traceback_coords
        return aligned_seq1, aligned_seq2
    
    

example = SmithWaterman("AGGTA", "ACGT")
print(f"Score: {example.score}\n")
print(f"Alignment:\n{example.alignment[0]}\n{example.alignment[1]}")



            

