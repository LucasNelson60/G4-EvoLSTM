import numpy as np

def analyze_seq(seq, GQs, scores, startI):
    # to record the start indicies of all possible GQ forming squence
    start_ind = []
    # to record the end indicies of all possible GQ forming sequence
    end_ind = []
    # to record the estimated melting temperature
    score_GQs = []
    # to calculate the multiplicity of each nucleotide in GQ forming sequence
    # initialized with zero arrays with length of the input sequence
    multi = [0]*len(seq)

    # We sequentially go through the GQ matrix (GQs) by the length
    # and check if each G4 motif is present
    for i in range(len(GQs)):

        # if the size of the G4 matrix is less than the size of the sequence
        # we would use a sliding window to go through the full sequence
        if len(GQs[i][0]) < len(seq):
            longer = len(seq) - len(GQs[i][0])+1
            for j in range(longer):
                # to cut the subsequence with length GQs[i][0]
                # dot the matrix with same length and vector and compare
                # if the answer equals 12
                
                ind = np.dot(np.array(GQs[i]), seq[j:len(GQs[i][0])+j]) == 12
                # find where the ind equals True and store it in the index
                ind = np.where(ind)[0]
                if len(ind) != 0:
                    for each_ind in ind:
                        # for finding the score (estimated melting temperature)
                        score_GQs.append(scores[i][each_ind])

                        # store start index (start with 1)
                        start_ind.append(j+1)

                        # store end index (excluded in the seq)
                        end_ind.append(j+len(GQs[i][each_ind]))

                        # update multiplicity
                        for k in range(len(GQs[i][each_ind])):
                            if GQs[i][each_ind][k] == 1:
                                multi[k+j] += 1
        elif len(GQs[i][0]) == len(seq):
            # to cut the subsequence with length GQs[i][0]
            # dot the matrix with same length and vector and compare
            # if the answer equals 12
            ind = np.dot(np.array(GQs[i]), seq) == 12
            # find where the ind equals True and store it in the index
            ind = np.where(ind)[0]
            if len(ind) != 0:
                for each_ind in ind:
                    # for finding the score (estimated melting temperature)
                    score_GQs.append(scores[i][each_ind])

                    # store start index (start with 1)
                    start_ind.append(1)

                    # store end index (excluded in the seq)
                    end_ind.append(len(GQs[i][each_ind]))

                    # update multiplicity
                    for k in range(len(GQs[i][each_ind])):
                        if GQs[i][each_ind][k] == 1:
                            multi[k] += 1
        else:
            # sequence shorter than reference sequence length, break
            break
    
    # determine how many G4s can fold in tandem
    Gregs = []
    if len(start_ind) !=0:
        temp = sorted(zip(start_ind, end_ind, score_GQs))
        start_ind, end_ind, score_GQs= map(list, zip(*temp))

        # since python list starts at 0, the position need to -1 to get index
        sI = start_ind[0]-1
        eI = end_ind[0]
        for i in range(1, len(start_ind)):
            if start_ind[i] > eI:
                element = [multi[sI:eI], sI, eI]
                Gregs.append(element)
                sI = start_ind[i] - 1
                eI = end_ind[i]
            elif start_ind[i] <= eI and end_ind[i] > eI:
                eI = end_ind[i]
        element = [multi[sI:eI], sI, eI]
        Gregs.append(element)

    # Go through and find G4 containing regions, defined as overlapping G4s
    for i in range(len(Gregs)):
        index_start = np.array(start_ind) >= Gregs[i][1]
        index_end = np.array(end_ind) <= Gregs[i][2]
        score = []
        ind = []
        for j in range(len(index_start)):
            if index_start[j]==True and index_end[j]==True:
                score.append(score_GQs[j])
                ind.append(j)
        Gregs[i].append(score)
        
        e = np.array(end_ind)[ind]
        #print("e:", e)
        s = np.array(start_ind)[ind]
        temp = sorted(zip(e, s))
        e,s= map(list, zip(*temp))
        num_tan = 0
        while sum(ind) > 0:
            v = min(e)
            ind = s > v
            s = np.array(s)[ind]
            e = np.array(e)[ind]
            num_tan += 1
        if num_tan == 0:
            num_tan = 1
        Gregs[i].append(num_tan)
    
    # update the starting and ending index:
    for i in range(len(Gregs)):
        Gregs[i][1] = Gregs[i][1] + startI
        Gregs[i][2] = Gregs[i][2] + startI
    return multi, Gregs