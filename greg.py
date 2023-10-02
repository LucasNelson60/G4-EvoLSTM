from get_GQs import get_GQs
from check_seq import check_seq
from statistics import median

def main(seq, max_loop, GQs, scores):

    G4CR_fields = ""
    _, GRegs = check_seq(seq, GQs, scores, max_loop)

    max_ntot = 0
    sum_ntots = 0
    max_ntand = 0
    sum_ntands = 0
    num_G4CRs = 0
    max_Tm_max = 0
    max_Tm_median = 0
    max_Tm_min = 0
    sum_G4CR_length = 0
    max_G4CR_length = 0
    avg_g_percentage = 0
    max_g_percentage = 0
    counter = 0
    
    for index, info in enumerate(GRegs):
        sequence = seq[info[1]:info[2]]
        each_multi = info[0]
        start_index = str(info[1]+1)
        end_index = str(info[2])
        length = str(info[2]-info[1])
        G_number = 0
        for nucleotide in sequence:
            if nucleotide == "G":
                G_number = G_number + 1
        G_content = str(round((G_number / int(length) * 100)))
        Ntot = str(int(sum(info[0])/12))
        Ntand = str(info[4])
        Tm_min = str(round(min(info[3]),1))
        Tm_median = str(round(median(info[3]), 1))
        Tm_max = str(round(max(info[3]), 1))

        if index != 0:
            this_string = "G4CR" + "," + "".join(sequence) + "," + ";".join(str(x) for x in each_multi) + "," + start_index + "," + end_index + "," + length + "," + G_content + "," + Ntot + "," + Ntand + "," + Tm_max + "," + Tm_median + "," + Tm_min + "\n"
        else:
            this_string = "G4CR" + "," + "".join(sequence) + "," + ";".join(str(x) for x in each_multi) + "," + start_index + "," + end_index + "," + length + "," + G_content + "," + Ntot + "," + Ntand + "," + Tm_max + "," + Tm_median + "," + Tm_min + "," + seq + "\n"

        G4CR_fields += this_string
        
        max_ntot = Ntot if int(Ntot) > int(max_ntot) else max_ntot
        max_ntand = Ntand if int(Ntand) > int(max_ntand) else max_ntand
        sum_ntots += int(Ntot)
        sum_ntands += int(Ntand)
        max_Tm_max = Tm_max if int(Tm_max) > int(max_Tm_max) else max_Tm_max
        max_Tm_median = Tm_median if int(Tm_median) > int(max_Tm_median) else max_Tm_median
        max_Tm_min = Tm_min if int(Tm_min) > int(max_Tm_min) else max_Tm_min
        num_G4CRs += 1
        sum_G4CR_length += int(length)
        max_G4CR_length = length if int(length) > int(max_G4CR_length) else max_G4CR_length
        avg_g_percentage += int(G_content)
        max_g_percentage = G_content if G_content > max_g_percentage else max_g_percentage
        counter += 1

    outputs = [G4CR_fields, num_G4CRs, max_ntot, sum_ntots, max_ntand, sum_ntands, max_Tm_max, max_Tm_median, max_Tm_min, sum_G4CR_length, max_G4CR_length, avg_g_percentage/counter, max_g_percentage]
    
    return (str(i) for i in outputs)
