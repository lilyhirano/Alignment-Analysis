import numpy as np
import pandas as pd

class NeedlemanWunsch:
    def __init__(self, seq1, seq2, match = 1, mismatch = -1, gap = -1):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
    
        self.rows = len(seq1) + 1
        self.columns = len(seq2) + 1

        self.fill_matrix()
        self.score = self.find_score()

        self.alignment = self.traceback()

    def init_matrix(self):
        self.matrix = np.zeros((self.rows, self.columns), dtype = int)
        self.matrix[0, :] = np.arange(0, (self.columns) * self.gap, self.gap)
        self.matrix[:, 0] = np.arange(0, (self.rows) * self.gap, self.gap)

    def fill_matrix(self):

        self.init_matrix()

        for i in range(1, self.rows):
            for j in range(1, self.columns):
                if self.seq1[i - 1] == self.seq2[j - 1]:
                    match_score = self.match
                else:
                    match_score = self.mismatch
            
            up = self.matrix[i - 1, j] + self.gap
            left = self.matrix[i, j - 1] + self.gap
            diagonal = self.matrix[i - 1, j - 1] + match_score

            self.matrix[i, j] = max(diagonal, up, left) #save the best score

    def find_score(self):
        #score is the bottom right cell of the matrix
        return self.matrix[self.rows - 1 , self.columns - 1]
    
    def traceback(self):
        i = self.rows - 1
        j = self.columns - 1
        align1, align2 = "", ""

        while i > 0 and j > 0:
            current = self.matrix[i, j]
            x1 = self.seq1[i - 1]
            x2 = self.seq2[j - 1]


            if x1 == x2:
                match_score = self.match
            else:
                match_score = self.mismatch

            if current == self.matrix[i - 1, j - 1] + match_score:
                align1 = x1 + align1
                align2 = x2 + align2
                i -= 1
                j -= 1

            elif current == self.matrix[i - 1, j] + self.gap:
                align1 = x1 + align1
                align2 = "-" + align2
                i -= 1

            else:
                align1 = "-" + align1
                align2 = x2 + align2
                j -= 1

        #fill rest with gaps if one side already aligned
        while i > 0:
            align1 = self.seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1

        while j > 0:
            align1 = "-" + align1
            align2 = self.seq2[j - 1] + align2
            j -= 1

        return align1, align2
    
    def print_matrix(self):
        row_labels = ["-"] + list(self.seq1)
        column_labels = ["-"] + list(self.seq2)

        NMWdf = pd.DataFrame(self.matrix, index = row_labels, columns= column_labels)

        print(NMWdf)