import math

def get_GQs(max_loop, max_bulge, min_Tm):
    # GQ_dict store the sequence and corresponding melting temperature
    GQ_list = []
    Tm_list = []

    # to construct the first G-track with or without bulges
    # consider 3 cases
    # if T1 = 0, no bulges at all in first G tract
    # if T1 = 1:max_bulge+1, the bulge is between the first and second G in the first G tract
    # GNGG, GNNGG, GNNNGG, ..., GNN...NGG
    # if T1 = max_bulge+1:max_bulge*2+1, the bulge is between the second and thrid G in the first G tract
    # GGNG, GGNNG, GGNNNG, ..., GGNN...NG
    # It's the same for T2, T3 and T4
    
    for T1 in range(0, max_bulge*2+1):
        # to keep track the bulges and number of tracts in each loop and add them together
        # instead of calculating at the end of the loop could save a lot of time
        total_bulge_1 = 0
        num_tracts_1 = 0
        if T1 == 0:
            t1 = [1,1,1]
        elif 0 < T1 < max_bulge + 1:
            t1 = [1] + [13]*T1 + [1, 1]
            total_bulge_1 = T1
            num_tracts_1 = 1
        else:
            t1 = [1, 1] + [13]*(T1-max_bulge) + [1]
            total_bulge_1 = T1 - max_bulge
            num_tracts_1 = 1
        for L1 in range(1, max_loop + 1):
            # construct the first loop (length from 1 to max_loop_length)
            l1 = [0]*L1
            for T2 in range(0, max_bulge*2+1):
                total_bulge_2 = 0
                num_tracts_2 = 0
                if T2 == 0:
                    t2 = [1,1,1]
                elif 0< T2 < max_bulge + 1:
                    t2 = [1] + [13]*T2 + [1, 1]
                    total_bulge_2 = T2
                    num_tracts_2 = 1
                else:
                    t2 = [1, 1] + [13]*(T2-max_bulge) + [1]
                    total_bulge_2 = T2 - max_bulge
                    num_tracts_2 = 1
                for L2 in range(1, max_loop + 1):
                    l2 = [0]*L2
                    for T3 in range(0, max_bulge*2+1):
                        total_bulge_3 = 0
                        num_tracts_3 = 0
                        if T3 == 0:
                            t3 = [1, 1, 1]
                        elif 0 < T3 < max_bulge + 1:
                            t3 = [1] + [13]*T3 + [1, 1]
                            total_bulge_3 = T3
                            num_tracts_3 = 1
                        else:
                            t3 = [1, 1] + [13]*(T3-max_bulge) + [1]
                            total_bulge_3 = T3 - max_bulge
                            num_tracts_3 = 1
                        for L3 in range(1, max_loop + 1):
                            l3 = [0]*L3
                            for T4 in range(0, max_bulge*2+1):
                                total_bulge_4 = 0
                                num_tracts_4 = 0
                                if T4 == 0:
                                    t4 = [1, 1, 1]
                                elif 0 < T4 < max_bulge + 1:
                                    t4 = [1] + [13]*T4 + [1, 1]
                                    total_bulge_4 = T4
                                    num_tracts_4 = 1
                                else:
                                    t4 = [1, 1] + [13]*(T4-max_bulge) + [1]
                                    total_bulge_4 = T4 - max_bulge
                                    num_tracts_4 = 1
                                
                                total_bulge = total_bulge_1 + total_bulge_2 + total_bulge_3 + total_bulge_4
                                each_seq = t1 + l1 + t2 + l2 + t3 + l3 + t4

                                num_tracts = num_tracts_1 + num_tracts_2 + num_tracts_3 + num_tracts_4
                                # find the number of bulges
                                #total_bulge = each_seq.count(13.0)
                                #print(total_bulge)
                                # 13.0 in t1: return 1 if ture, 0 otherwise
                                # check how many G tracks in the sequence contain bulge
                                #num_tracts = (13.0 in t1) + (13.0 in t2) + (13.0 in t3) + (13.0 in t4)
                                L = math.log10(len(l1)*len(l2)*len(l3))
                                # The equation for calculating the melting temperature is:
                                # Tm_est = a - b * ln(l1*l2*l3) - d * Nb - f * (Lb - Nb)
                                Tm_est = 89.9 - 19.2*L - 20*num_tracts - 8.5*(total_bulge-num_tracts)
                                
                                # Check if the estimated melting temperature
                                # passes the threshold
                                if Tm_est >= min_Tm:
                                    GQ_list.append(each_seq)
                                    Tm_list.append(Tm_est)
    # to sort the GQ list and Tm list at the same time
    # sorting with respect to the length of each GQ length
    temp = sorted(zip(GQ_list, Tm_list), key=lambda x: len(x[0]))
    sorted_GQ, sorted_Tm = map(list, zip(*temp))

    # To make a list of lists to a list of matrix structure
    # The outside list has length max_len_sequence - min_len_sequence elements
    # inner two lists are used as a matrix structure
    # The middle list is a matrix that group all sequence with the same length together
    # The each inner list represents a possible GQ forming sequence
    # It helps with matrix multiplication in analyze seq function later
    new_GQ = []
    new_Tm = []
    curr_len = len(sorted_GQ[0])
    GQ_matrix = []
    Tm_matrix = []
    index = 0
    for element in sorted_GQ:
        if curr_len == len(element):
            GQ_matrix.append(element)
            Tm_matrix.append(sorted_Tm[index])
            index = index + 1
        else:
            new_GQ.append(GQ_matrix)
            new_Tm.append(Tm_matrix)
            GQ_matrix = [element]
            Tm_matrix = [sorted_Tm[index]]
            curr_len = len(element)
            index = index + 1
    new_GQ.append(GQ_matrix)
    new_Tm.append(Tm_matrix)
    return new_GQ, new_Tm