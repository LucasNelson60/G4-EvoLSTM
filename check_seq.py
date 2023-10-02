import numpy as np
from analyze_seq import analyze_seq

def check_seq(seq, GQs, scores, L):
    seq = (np.array(list(seq)) == 'G').astype(int)
    seqs = splitSeq(np.array(seq), L)
    Gregs = []
    x = np.zeros(len(seq))
    for i in range(len(seqs)):
        if sum(seqs[i][0]) >= 12:
            multi, Greg = analyze_seq(seqs[i][0], GQs, scores, seqs[i][1])
            for j in range(len(multi)):
                x[seqs[i][1]+j] = multi[j]
            if len(Greg) != 0:
                for element in Greg:
                    Gregs.append(element)
    return x, Gregs

def splitSeq(seq, L):
    # The sequence would be [0,...,0,seq,0,...,0,] for the move mean
    protected_seq = np.concatenate([np.zeros(L+1), seq, np.zeros(L+1)])
    A = np.convolve(protected_seq, np.ones(L+1), 'valid')/(L+1)
    seqs = []
    tag = 0
    for j in range(len(A)):
        # A[0] = 0
        # so the start index need to be subtracted by 1
        if A[j] > 0 and tag == 0:
            tag = 1
            startI = j-1
        elif A[j] == 0 and tag == 1:
            tag = 0
            endI = j-L
            if endI - startI + 1 >= 15 and sum(seq[startI:endI]) >= 12:
                seqs.append([seq[startI:endI], startI, endI-1])
                # store the start index and end index in the seqs
    return seqs