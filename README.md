# Protein Alignment-Analysis

This project implements and compares four major methods used in protein sequence alignment: 

PAM250, BLOSUM62, Needleman–Wunsch, and Smith–Waterman

These dynammic programming algorithms and scoring systems have been used for the foundation of many modern bioinformatics tools like BLAST, Clustal Omega, and FASTA.

Overview/Implementations:

1. Needleman–Wunsch Algorithm (Global Alignment)
* Dynamic programming classical algorithm that computes an optimal end-to-end alignment of 2 FULL sequences
* Constant gap penalty
* Basis for many alignment scoring systems
* Aligns entire sequences end-to-end to maximize overall similarity

2. Smith–Waterman Algorithm (Local Alignment)
* Dynamic programming local algorithm that identifies the highest-scoring matching subsections between 2 sequences
* Returns the optimal local subsequence match
* Different traceback and scoring scheme than Needleman-Wunsch
* Ideal for detecting conserved regions within otherwise divergent sequences

3. BLOSUM62 Alignment (Global)
* Uses the Needleman–Wunsch algorithm withthe BLOSUM62 amino-acid substitution matrix as the scoring model
* Designed from conserved blocks in real proteins
* Best for aligning moderately diverged proteins (about 62% identity or lower)

4. PAM250 Alignment (Global)
* Uses the Needleman–Wunsch algorithm with the PAM250 evolutionary substitution matrix as the scoring model
* Based on evolutionary mutation probabilities
* Best suited for aligning more distantly related proteins ( about 250 PAM units)



| Method                         | Time Complexity | Space Complexity | Explanation                      |
| ------------------------------ | --------------- | ---------------- | -------------------------------- |
| **Needleman–Wunsch**           | O(nm)           | O(nm)            | n = len(seq1), m = len(seq2)     |
| **Smith–Waterman**             | O(nm)           | O(nm)            | Same DP structure but local      |
| **PAM250 alignment**           | O(nm)           | O(nm)            | NW + constant-time PAM lookup    |
| **Matrix lookup (PAM/BLOSUM)** | O(1)            | O(1)             | Dict/hashtable lookup            |
| **BLOSUM62 alignment**         | O(nm)           | O(nm)            | NW + constant-time BLOSUM lookup |


## Set Up

1. Clone or download the repository
2. Set up virtual environment
```bash
python3 -m venv venv
venv\Scripts\activate    # (Windows)
source venv/bin/activate   # (Mac / Linux)
```
3. Install Dependencies
    - This repository uses NumPy and Pandas
```bash
pip install -r requirements.txt
```
4. Run each cell in alignment-analysis.ipynb
