from blosum62 import blosum62

GAP = -5    # gap penalty, we can chnge this

def score(a, b):
    """Return BLOSUM62 score for characters a and b."""
    if (a, b) in blosum62:
        return blosum62[(a, b)]
    if (b, a) in blosum62:  # matrix is symmetric
        return blosum62[(b, a)]
    raise KeyError(f"No score for pair({a}, {b})")

def blosum_align(seq1, seq2):
    """Needleman-Wuncsch global alignment using BLOSUM62."""
    n = len(seq1)
    m = len(seq2)

    # DP matrices
    dp = [[0]*(m+1) for _ in range(n+1)]
    ptr = [[None]*(m+1) for _ in range(n+1)]    # traceback

    # initialize boundaries
    for i in range(1, n+1):
        dp[i][0] = i * GAP
        ptr[i][0] = 'U'
    for j in range(1, m+1):
        dp[0][j] = j * GAP
        ptr[0][j] = 'L'

    # fill DP
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

    # traceback
    aligned1 = []
    aligned2 = []
    i, j = n, m

    while i > 0 or j > 0:
        move = ptr[i][j]

        if move == 'D':
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i-=1
            j-=1
        elif move == 'U':
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i-=1
        elif move == 'L':
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j-=1
    
    return dp[n][m], ''.join(reversed(aligned1)), "".join(reversed(aligned2))