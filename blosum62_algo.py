from blosum62 import blosum62

GAP = -5    # gap penalty

def score(a, b):
    """Return BLOSUM62 score for characters a and b."""
    if (a, b) in blosum62:
        return blosum62[(a, b)]
    if (b, a) in blosum62:
        return blosum62[(b, a)]
    raise KeyError(f"No score for pair({a}, {b})")

def read_fasta(path):
    """
    Read a FASTA file containing one or more sequences.
    Returns a list of (header, sequence) tuples.
    """
    sequences = []
    header = None
    seq_lines = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # save previous sequence
                if header is not None:
                    sequences.append((header, "".join(seq_lines).upper()))
                    seq_lines = []
                header = line[1:]  # drop ">"
            else:
                seq_lines.append(line)

        # last sequence
        if header is not None:
            sequences.append((header, "".join(seq_lines).upper()))

    return sequences

def blosum_align(seq1, seq2):
    """Needleman-Wunsch global alignment using BLOSUM62."""
    n = len(seq1)
    m = len(seq2)

    dp = [[0] * (m+1) for _ in range(n+1)]
    ptr = [[None] * (m+1) for _ in range(n+1)]

    for i in range(1, n+1):
        dp[i][0] = i * GAP
        ptr[i][0] = 'U'
    for j in range(1, m+1):
        dp[0][j] = j * GAP
        ptr[0][j] = 'L'

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = dp[i-1][j-1] + score(seq1[i-1], seq2[j-1])
            delete = dp[i-1][j] + GAP
            insert = dp[i][j-1] + GAP

            dp[i][j] = max(match, delete, insert)

            if dp[i][j] == match:
                ptr[i][j] = 'D'
            elif dp[i][j] == delete:
                ptr[i][j] = 'U'
            else:
                ptr[i][j] = 'L'

    aligned1 = []
    aligned2 = []
    i, j = n, m

    while i > 0 or j > 0:
        move = ptr[i][j]
        if move == 'D':
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif move == 'U':
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        elif move == 'L':
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1

    return dp[n][m], "".join(reversed(aligned1)), "".join(reversed(aligned2))

def blosum_align_fasta(file1, file2):
    """
    Align the FIRST sequence in each FASTA file.
    Returns score, aligned_seq1, aligned_seq2.
    """
    seqs1 = read_fasta(file1)
    seqs2 = read_fasta(file2)

    if not seqs1 or not seqs2:
        raise ValueError("One of the FASTA files contains no sequences.")

    # Use first entry by default
    seq1 = seqs1[0][1]
    seq2 = seqs2[0][1]

    return blosum_align(seq1, seq2)

def main():
    files = [
        "data/HBB_COLLI.fasta",
        "data/HBB_HUMAN.fasta",
        "data/HBB1_MOUSE.fasta"
    ]

    # Load sequences
    seqs = {}
    for f in files:
        entries = read_fasta(f)
        if not entries:
            raise ValueError(f"{f} contains no sequences.")
        header, sequence = entries[0]
        seqs[f] = sequence

    # List of file names for pairwise alignment
    file_list = list(seqs.keys())

    print("Pairwise BLOSUM62 Global Alignments:\n")

    # Compare all three pairwise
    for i in range(len(file_list)):
        for j in range(i+1, len(file_list)):
            f1 = file_list[i]
            f2 = file_list[j]

            seq1 = seqs[f1]
            seq2 = seqs[f2]

            alignment_score, aln1, aln2 = blosum_align(seq1, seq2)

            print(f"=== {f1}  vs  {f2} ===")
            print(f"Alignment Score: {alignment_score}")
            print(aln1)
            print(aln2)
            print("\n")

if __name__ == "__main__":
    main()
