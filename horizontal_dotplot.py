#!/usr/bin/env python
import gzip
import regex as re


def reverse_complement(seq):
    """
    return reverse complement of sequence
    """
    # reverse sequence
    seq = seq.upper()[::-1]
    # complement
    seq = seq.replace("A", "t")
    seq = seq.replace("T", "a")
    seq = seq.replace("C", "g")
    seq = seq.replace("G", "c")
    # convert to upper case
    seq = seq.upper()
    return seq

def read_fasta_to_dict(fasta_file):
    """
    read fasta file and return dictionary with contig name as key and sequence as value
    """
    # first arch sequence is list of lines
    fasta_dict = {}
    #  fasta could be plain text or in gzip format
    if fasta_file.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(fasta_file, "r") as f:
        for line in f:
            # convert to string
            line = line.decode("utf-8")
            if line[0] == ">":
                contig_name = line.strip()[1:]
                fasta_dict[contig_name] = []
                continue
            else:
                fasta_dict[contig_name].append(line.strip())
    # concatenate list to str
    fasta_dict = {k: "".join(v) for k, v in fasta_dict.items()}
    return fasta_dict

s = read_fasta_to_dict("/mnt/raid/454_data/original_sequencing_data/0_submissions_to_databases/221003_Psativum_CEN6_assembly/CEN6_ver_220406.fasta.gz")
# to test - take one sequence from dict
s_str = list(s.values())[0]
contig_name = list(s.keys())[0]

N = len(s_str)
kmer_size = 21
window_size = 50000
kmer_start_indexes = range(int(window_size/2) - int(kmer_size/2), N - int(
    window_size/2) - int(kmer_size/2), 1)
window_start_indexes = range(0, N - window_size, 1)

frag = 0
for index, _ in enumerate(kmer_start_indexes[:1000000]):
    ksi = kmer_start_indexes[index]
    kmer = s_str[ksi:ksi + kmer_size]
    wsi = window_start_indexes[index]
    window = s_str[wsi:wsi + window_size]
    # m = re.finditer("(" + kmer + ")" + "{e<=3}", window)
    m = re.finditer(kmer, window)
    positions = [i.start() for i in m]
    # there is alway hit to itself
    if len(positions) > 1:
        #print(F'index: {index}; ksi: {ksi} ; {kmer} : {positions}')
        # convert indexes to position is s_str
        positions = [i + wsi for i in positions]
        #print(positions)
        for pos in positions:
            if pos != ksi:
                out = sorted([pos, pos + 1, ksi, ksi + 1])
                frag = +2
                print(F'0\t{contig_name}\t{out[0]}\t{frag}\t'
                      F'0\t{contig_name}\t{out[2]}\t{frag+1}')

    # same for reverse complement
    kmer = reverse_complement(kmer)
    m = re.finditer(kmer, window)
    positions = [i.start() for i in m]
    # there ishould not be hit to itself
    if len(positions) > 0:
        #print(F'index: {index}; ksi: {ksi} ; {kmer} : {positions}')
        # convert indexes to position is s_str
        positions = [i + wsi for i in positions]
        #print(positions)
        for pos in positions:
            if pos != ksi:
                out = sorted([pos, pos + 1, ksi, ksi + 1])
                frag = +2
                print(
                    F'1\t{contig_name}\t{out[0]}\t{frag}\t'
                    F'0\t{contig_name}\t{out[2]}\t{frag + 1}')