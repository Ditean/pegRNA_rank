import torch
import pandas as pd
from collections import Counter
import numpy as np
import re
import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate pegRNA sequence.")
    parser.add_argument("--amplicon", type=str, help="Target amplicon sequence")
    parser.add_argument("--insert", type=str, help="Desired insertion sequence (optional)")
    return parser.parse_args()

args = parse_arguments()

amplicon = args.amplicon.upper()

if args.insert:
    insert = args.insert.upper()
else:
    default_insert1="atgatcctgacgacggagaccgccgtcgtcgacaagcc"
    default_insert2="gatgatcctgacgacggagaccgccgtcgtcgacaagccg"
    default_insert3="cggatgatcctgacgacggagaccgccgtcgtcgacaagccggc"
    default_insert4="ccggatgatcctgacgacggagaccgccgtcgtcgacaagccggcc"

def rc(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(s))


def pegRNA_generation(amplicon, insert):
    amplicon = amplicon.upper()
    insert = insert.upper()
    amplicon_plus = amplicon
    amplicon_minus = rc(amplicon_plus)

    ngg_plus = [m.start() for m in re.finditer('GG', amplicon_plus)]
    ngg_minus = [m.start() for m in re.finditer('GG', amplicon_minus)]
    ngg_plus = [x for x in ngg_plus if x > 20]
    ngg_minus = [x for x in ngg_minus if x > 20]

    guide_plus = []
    guide_minus = []
    for element in ngg_plus:
        guide = amplicon_plus[(element - 21):(element - 1)]
        guide_plus.append(guide)
    for element in ngg_minus:
        guide = amplicon_minus[(element - 21):(element - 1)]
        guide_minus.append(guide)

    pbs_length = np.arange(15) + 3
    RT_length = np.arange(30) + 5

    pegRNA = []
    scaffold = "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccGAGTCGGTGC"
    scaffold = scaffold.lower()

    for guide in guide_plus:
        for pbs_len in pbs_length:
            for rt_len in RT_length:
                guide_17 = rc(guide[:17])
                pbs = guide_17[:pbs_len]
                guide_right_3 = guide[17:20]
                rt_temp_pos = amplicon_plus.find(guide) + 20
                if rt_temp_pos + rt_len < (len(amplicon_plus)):
                    rt_temp = amplicon_plus[rt_temp_pos:(rt_temp_pos + rt_len)]
                    whole_temp = "".join([insert, guide_right_3, rt_temp])
                    rt_to_add = rc(whole_temp)
                    seq = ''.join([guide, scaffold, rt_to_add, pbs])
                    if seq[0] == "G":
                        seq = seq
                    else:
                        seq = "".join(["G", seq])
                    pegRNA.append(seq)

    for guide in guide_minus:
        for pbs_len in pbs_length:
            for rt_len in RT_length:
                guide_17 = rc(guide[:17])
                pbs = guide_17[:pbs_len]
                guide_right_3 = guide[17:20]
                rt_temp_pos = amplicon_minus.find(guide) + 20
                if rt_temp_pos + rt_len < (len(amplicon_minus)):
                    rt_temp = amplicon_minus[rt_temp_pos:(rt_temp_pos + rt_len)]
                    whole_temp = "".join([rc(insert), guide_right_3, rt_temp])
                    rt_to_add = rc(whole_temp)
                    seq = ''.join([guide, scaffold, rt_to_add, pbs])
                    if seq[0] == "G":
                        seq = seq
                    else:
                        seq = "".join(["G", seq])
                    pegRNA.append(seq)
    return pegRNA

character = ['A', 'T', 'C', 'G', 'N']
full_dict = []
for a in character:
    for b in character:
        for c in character:
            mer = ''.join([a, b, c])
            full_dict.append(mer)


def splittokmer(seq, sliding, full_dict):
    seq = seq.upper()
    y = np.zeros(len(full_dict))
    kmer = []
    for index in range(len(seq) - (sliding - 1)):
        a = seq[index:index + sliding]
        kmer.append(a)
    new_dict = Counter(kmer)
    for word in range(len(full_dict)):
        name = full_dict[word]
        y[word] = new_dict[name]
    return y


def truncate_seq(seq):
    if len(seq) > 198:
        seq = seq[:198]
    elif len(seq) < 198:
        insert = 'N' * (198 - len(seq))
        new = str(''.join([seq, insert]))
        seq = new
    else:
        seq = seq
    return seq

if args.insert:
    insert = args.insert.upper()
    pegRNA_list=pegRNA_generation(amplicon, insert)
else:
    group1=pegRNA_generation(amplicon, default_insert1)
    group2=pegRNA_generation(amplicon, default_insert2)
    group3=pegRNA_generation(amplicon, default_insert3)
    group4=pegRNA_generation(amplicon, default_insert4)
    pegRNA_list=group1 + group2 + group3 + group4


kmer_store = np.zeros((len(pegRNA_list), 125))

for index in range(len(pegRNA_list)):
    new_seq = truncate_seq(pegRNA_list[index])
    input = splittokmer(new_seq, 3, full_dict)
    kmer_store[index, :] = input

model = torch.load('mlp.pt')
model.eval()

output_filename = 'Result.csv'
counter = 1

while os.path.exists(output_filename):
    output_filename = f'Result_{counter}.csv'
    counter += 1

input_tensor = torch.Tensor(kmer_store)
output_tensor = model(input_tensor)
res = output_tensor.detach().numpy()
prob = res[:, 1]
sequence = pd.Series(pegRNA_list, name='pegRNA Sequence')
prob = pd.Series(prob, name='Score')
result_df = pd.concat([sequence, prob], axis=1)
result_df = result_df.sort_values(by=['Score'], ascending=False)
result_df.to_csv(output_filename, index=None)

print(f"Results saved to {output_filename}")