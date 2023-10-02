from tensorflow.keras.models import load_model
from tensorflow.keras.backend import clear_session
import numpy as np
from joblib import load
import os
from tqdm import tqdm
import gc

def decode_sequence(input_seq, encoder_model, decoder_model, encode_dimension, context_length, onehot_encoder, label_encoder):
    length = len(input_seq)
    nucleotide = label_encoder.inverse_transform(range(encode_dimension))

    initial_context = np.expand_dims(input_seq[0: context_length], axis=0)
    states_value = encoder_model.predict(initial_context)
    target_seq = np.zeros((1, 1, encode_dimension*2))
    target_seq[0][0]= np.hstack((input_seq[0], onehot_encoder.transform(np.ones(1).reshape(-1,1))[0]))
    decoded_seq = ''

    for i in tqdm(range(1, length)):
        if i%context_length == 0 :
            context = np.expand_dims(input_seq[i: i+context_length], axis=0)
            states_value = encoder_model.predict(context)
        output_tokens, h, c = decoder_model.predict([target_seq] + states_value)
        sampled_token_index = np.random.choice(encode_dimension, 1, p=output_tokens[0, -1, :])[0]
        sampled_nucleotide = nucleotide[sampled_token_index]

        if (sampled_nucleotide == '') or (not sampled_nucleotide.isprintable()) or (sampled_nucleotide.isspace()):
            decoded_seq += '0'
        else :
            decoded_seq += sampled_nucleotide

        # Update the target sequence (of length 1).
        target_seq = np.zeros((1, 1, encode_dimension*2))
        temp = np.zeros((encode_dimension))
        temp[sampled_token_index] = 1
        target_seq[0][0]= np.hstack((input_seq[i], temp))

        if i == length -1:
            output_tokens, h, c = decoder_model.predict(
                [target_seq] + states_value)
            sampled_token_index = np.random.choice(encode_dimension, 1, p=output_tokens[0, -1, :])[0]
            sampled_nucleotide = nucleotide[sampled_token_index]

            decoded_seq += sampled_nucleotide


        # Update states
        states_value = [h, c]
        
    del output_tokens, h, c, target_seq
    return decoded_seq

def simulate(gpu, anc):

    os.environ["CUDA_VISIBLE_DEVICES"]=gpu
    context_length = 15

    des = "TCGTACGATACATTTAAAGGAAGTGTAAATGAAGCTGGAATAACGTGTGAGCAAAACCATAATTTTGAAAATATTTTACCAGGAAAAGGAAGTGTGAAAAAAAGTCACCTTAGAGTTTACGGGGAGGAGTGATAGTTTGTAGACCGATAACAAGTGGTCGGGTTAAAAAAAGGGTCGTGAAAGGTTTCCGGCGAAACCAAAATTTTATTTTATTATTATTCGTTTTGGTTTATGGGGAGTGGAAGGTGGCACGGAGGGGAAGGGGTTATTTAGGTCACAGAACGAAATTTTGAACCACCAATCTTCTACTACCCAAAGATTCTGCACACCCGATTTCGAACAAAACGACAAATCCCAAAACAACCTTTAAAAAAAGCACACACATGAACAATTAATAAAGTGCAAAGGGTATTTGCCAAGAGGTATCCCACTACAAGTAACAGTCATCACTACCCAATTAAAAGTGGTGGCGAATACGCCAACTTATCGGCGGAGACTTGGTGAAAAAGGAGGTCATTAAGGAGAAGGAAGGCGGGGAGACGTCAGCCGGACTTCCTTAGTGTTTCTCCACCGACCTTTGAGCAAAATTCCTCGGCGGGGAGGAAGGGGCGACCTTTGGGGCGCGAGGCCTGCGGGGACGGGGGTGGGCTGGGGGCGGGAGCAACCGTAGGCCTGCGCCGCTAGAGCCGGCGGTCATCTCCCGTGTGAATGAAATGAAAGCGTTTGGGCTCGCGCCCACGTCGGGCCCCTCCCCCGCCCCTTTCTGCGAAACGTCGTTTTGGCTCGGATCCCTAACCAGCGAGGGGCGCAAACGCCGTTTCCGGACCTCCGACCTCATTAAACGTTAGGAATTTCGACTTAACACGGCACGCGGACTAAACCTTCAATGATGAGTAAATTGTTAACTTGCGACTCGACGTTTGAGTTGCCCATTATTGAATAGAGCTTGTGTACGTGCGATATGTGCGGGGGGGGGGGGGGCTTAAAAAAAGAAAAAACCCCCCCCCCCCCCCCCTTTTCAAATGAATTTTACGGAAACCCACTCCCTGGTCCCTACCCTTCTAACGGAAACAAAAACTACGGGACCATATTACGTTTTATCGTTTAGGGCTCCCTTATACGTAATATATAATTTATATCTAAGTAAGTCCCTCGCTTGTTTAGTACACAGCCCAACCCGTTGTACGACCCAGCTGCACATGTACTTTACACATATGTACATACGCCTAATATATATTACGTGAAAGTGATCATAAGTCTTTTTTTAAACTCAGTCACTTGATCCTTTAATTACGGACTTTACGCCGGTTAAAAATTAATCGAGTTCTGTGGGGGGGGGGGGTTTTTTCCGTGCCTTCATTATGAGGCGAGGAGAGGAAACTAGTCTTAGCTACGTAAAAAACACGTACTAGTGCTAAGGTTATTATTTCCCCTTTCTCCTGGGCCTTTCCTTAATTTGCAGGCCAAACAGGCCCCTCCTTTCTCAATTGCCAAAAGAAGTGGTCCCAGAGACGACTGAGGGGTCCGTGCCAGGCGTTCGAGCGGAGAGGGGGGAAGACCCCTCAGGCCAGGGCGCCAAGCGCGGGGGGGGGGGGGGGGTTATAAGAGGGCAGATCGCGGAAACTAAAGGGGGTTTGGGCCCTCGGGCTCTGACCACGTTTGGCCGCGGTGCCCCGCGTTTCTCCTAAACAGAGAAGACTTTGGACCGACGCCTTAACCCTTGAGGCACACTCTCCGCGCCCCCACCCCGCCCCCACGTACCGTCTCGGGGGCGTTGGCGCGAGGAGGGGATGAGGTCGAGCCCTTGTGCGTCTTTCTAGAGTCCCTACTTGTCGGCGGAGGGCGCGCCCCGGGGCACGTTCAGACGCGCCCGCCCCGCCCCCCCCCCCACAGACGCGAAACCGTCGTTTAACCCCCTAAGTAAGACCCACCTTCCACGGGTTAGCTCTATCGACATGCGTGTATTACGTATTGCGTAGTGTGAGGGGTTGTTTACGTCACCCACAAATAAGTATTGCGCGAGAGGTTCATATGCACCGTTACGCAACGACTCAATAAAATTATTAAGGTCCGTAGTAAAAGGAGGGATATGGAGATAGTAAATAGGGGTAGATGTGATTGTGGGGTGTGTGTCTCGCGCGCGGGTATTTATGGGGAGAAAGGAGGAGAGAAGCACCCTGAATAAAAGTTTTGCGACGGGAAAGGGGTCGGAATCCCTCCGCGGGCGGGCGGGCCCTGCACGCACCGCACCGGCACCCATGCGCCACATAAGCGTCACAACTCCCGTCGACAAGGCGGGCGCTACTAAATTCGTGCGTCCTGTTCATACGCCAAACAGTTTGTCGCGACGCGGTCTCCTCGTCGTCTCTTTCCCTCTCCCAAACTCGCCCTCGTTTTCTTTTACCATCCGCGCGCGTCAATTAAGTACGGCGCGAGAATGAGACAAATGTAGGCTCTCGACCTCACGGTCCGACGGCCCGACTCAGAGGAGGGGGCGGGGAGGGAGGGGTCTTCGCGGGGAGGGCCCAAGGGTTTCGTCTCCCGCCCCCTCTTTCTTTTTTCTAGGAGAGAGCGCTTAGGGGCGGGTGGCGGGAAATATTCCGCTCCCAGACGCGCCGGCTCCTGGGGGCTCGACGCGAAGAGCGCCGGCGGCGGCGGCGCGGGGCCGGCGGGGACCGAGGGGAGGACGGAGCTCTTCCCGTCCCGAAGAGTCTCCGAACCGCCCTTTTTCTTGCCTCCCTCCCTAGCGCGGCGCATATTTTCGGCCAAAAGACCCGAAATAGACTGAGCGACATCATTAAGGTCGCTCTCCGTCTCCCTCGCTCGCCCGCCGGGCGCTCCCACCTTCTCGGCTCGCTCGTCTCGGCGCGAGGCCCGCGGGCCCCTTCCCTCTAGGCCTCGCTTTCCCCCCGAAGCGGAGGCCGGGTCGGGCGGGCGGCTGGGGACGGTCGCCAGGCGTTGGGGGCGGCGTAGGCGCTTTGAAACGGGTAACGTCGCCCGCCCGTGAAACGTGACCTTGAATGTTGTGGGCTCGCTCCTGCGCTGAGAGGCCTGCGCCCCTCCGATAAGGCGGATAAACCCCTCTGAAGGGGCGGCGACGGGCCTGGGCGAGGAGACTTTCCGCGAGGAGCGGCGACAAGCCTGCGACCTAAAGGAAGCCTATCACCTTTTGGGCCATTCGTGGGCCTAGATAAACAGAAAATTAAAAAAAGAGATGGCGAAATTACGGCGCTACTCAGCTTACGGCTTAGCCCCACAGAAAAGAGGGTAAGGACGCGATAACTGTGAAAAGAGTCTCATCAACACCATCCAGCCCCACCCCACCCCCCGCTAGGTCTTGACCTAGCCCCATCTCACTGAACAGTTCTACCCTCCCCCCTTCCGTCTCCCTTTTGCCCCCACCAAAAACTTTGATCGGAAATCTCTAAAAACAAACGGATACAACTGCGACTGGGGCCGGCCGGCCTGTAAGGACGAAATAACGTAATTAACGAAAAACCCAACCCCCCCCCACCCGATACGAAGCGCCACCCGTCTTTCGGGAAACGTAGGACTCGAGGAACCTCATTCCTAATGTATAACGCACACACTCGCTCTCTCGAGGCGTCGGCGACTGAAAAGGGGCAGAGGCCCTCCCGTAAATTTAAAGCCGAATGGCGTAAAGACTGTCGGCCTCTGTCTGTGACGCCGCGCAGGGCGGGCGGCTAGGGGCGCCGCAAAGGTTGAGCGGGACGAGGAAAATTCTTCAACCGTAAACCGAAAAATTTTTTATTATTATTTTAATTTTTAGAACCAGAGGTCTCCACAATCCTGAACCACAACCCCTCCGCGTCCCTCCCCCTTCCCCCCCGTCCCTACACAGGCTAAGAGGACCTTAGCAACTGAACCCTTTTGGTCCCGCTTAGAGGCGTGGGTCGGAACCGAGGGGGCGGCGGCGGGCGCCCACAGGGGCGCGGGCTCTACGCCTCCCTGCCGCTCCTCGCCCCGAGGCCCGAAAAGG"

    anc = np.array(list(anc+'0'))
    des = np.array(list(des+'0'))
    
    label_encoder = load('/home/mcb/users/dlim63/research/alignment/label_encoder.joblib') 
    onehot_encoder = load('/home/mcb/users/dlim63/research/alignment/onehot_encoder.joblib') 
    
    integer_encoded_des = label_encoder.transform(des)
    integer_encoded_anc = label_encoder.transform(anc)
    integer_des = integer_encoded_des.reshape(len(integer_encoded_des), 1)
    encoded_des =onehot_encoder.transform(integer_des)
    integer_anc = integer_encoded_anc.reshape(len(integer_encoded_anc), 1)
    encoded_anc = onehot_encoder.transform(integer_anc)

    encode_dimension= len(encoded_des[0])
    encoder_model = load_model("/home/mcb/users/dlim63/research/models/E_insert3__HPGPNRMPC_hg38_10.h5")
    decoder_model = load_model("/home/mcb/users/dlim63/research/models/D_insert3__HPGPNRMPC_hg38_10.h5")
    decoded_seq = decode_sequence(encoded_anc, encoder_model, decoder_model, encode_dimension, context_length, onehot_encoder, label_encoder)

    del encoder_model
    del decoder_model
    gc.collect()

    os.environ.pop("CUDA_VISIBLE_DEVICES", None)
    clear_session()
    
    return decoded_seq
