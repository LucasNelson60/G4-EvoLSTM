# AUTHORS: Amos Zhang and Christopher Hennecker
from check_seq import check_seq
from statistics import median

def main(seq, seq_name, max_loop, GQs, scores):

    lines = []
    _, GRegs = check_seq(seq, GQs, scores, max_loop)

    for info in GRegs:
        
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
        
        line = [seq_name, "".join(sequence), ";".join(str(x) for x in each_multi), start_index, end_index, length, G_content, Ntot, Ntand, Tm_max, Tm_median, Tm_min]
        lines.append(line)
        
    return lines
