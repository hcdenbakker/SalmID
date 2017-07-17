#!/usr/bin/env python3

from sys import argv
import gzip
import io
import pickle
import os


def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'K', 'R': 'Y', 'W': 'W',
                            'S': 'S', 'Y': 'R', 'K': 'M', 'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V'}
    return "".join(complement[base] for base in reversed(sequence))


def createHashTable(string, kmer):
    kmer_table = {}
    sequence = string.strip('\n')
    for i in range(len(sequence)-kmer+1):
        new_mer =sequence[i:i+kmer]
        new_mer_rc = reverse_complement(new_mer)
        if new_mer in kmer_table:
            kmer_table[new_mer.upper()] += 1
        else:
            kmer_table[new_mer.upper()] = 1
        if new_mer_rc in kmer_table:
            kmer_table[new_mer_rc.upper()] += 1
        else:
            kmer_table[new_mer_rc.upper()] = 1
    return kmer_table



def target_read_kmerizer(file, k, kmerDict):
    # adaptation of function in stringMLST; if kmer of target sequence is found in middle of read,
    # read is either safed to file or used for further analysis
    # todo: make this paired-end?
    '''
    maybe use this here to make it faster
    with gzip.open('input.gz','r') as fin:        
            for line in fin:        
                print('got line', line) 

    initial implementation makes it slower!!!!
    check out io.BufferedReader (http://aripollak.com/pythongzipbenchmarks/)
    '''
    # readFile = open('matching_reads.fastq', 'w')
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = []
    for line in io.BufferedReader(gzip.open(file)):
        start = int((len(line) - k) // 2)
        if i % 4 == 2:
            s1 = line[start:k + start].decode()
            if s1 in kmerDict:
                n_reads += 1
                total_coverage += len(line)
                target_mers.append(set([k for k in createHashTable(str(line.decode()), k)]))
        i += 1
        if total_coverage >= 10000:
            break
    if len(target_mers) == 0:
        return 0
    else:
        return set.union(*target_mers)

def main():
    ex_dir = os.path.dirname(os.path.realpath(__file__))
    query_fastq_gz = argv[1]
    f = open(ex_dir + "/invA_mers_dict", "rb")
    sets_dict = pickle.load(f)
    f.close()
    allmers = sets_dict['allmers']
    uniqmers_I = sets_dict['uniqmers_I']
    uniqmers_IIa = sets_dict['uniqmers_IIa']
    uniqmers_IIb = sets_dict['uniqmers_IIb']
    uniqmers_IIIa = sets_dict['uniqmers_IIIa']
    uniqmers_IIIb = sets_dict['uniqmers_IIIb']
    uniqmers_IV = sets_dict['uniqmers_IV']
    uniqmers_VI = sets_dict['uniqmers_VI']
    uniqmers_VII = sets_dict['uniqmers_VII']
    uniqmers_VIII = sets_dict['uniqmers_VIII']
    uniqmers_bongori = sets_dict['uniqmers_bongori']
    target_mers = target_read_kmerizer(query_fastq_gz, 27, allmers)
    if target_mers == 0:
        print('No reads found matching invA, no Salmonella in sample?')
    else:
        p_bongori = (1 - len(uniqmers_bongori - target_mers) / len(uniqmers_bongori)) * 100
        p_I = (1 - len(uniqmers_I - target_mers) / len(uniqmers_I)) * 100
        p_IIa = (1 - len(uniqmers_IIa - target_mers) / len(uniqmers_IIa)) * 100
        p_IIb = (1 - len(uniqmers_IIb - target_mers) / len(uniqmers_IIb)) * 100
        p_IIIa = (1 - len(uniqmers_IIIa - target_mers) / len(uniqmers_IIIa)) * 100
        p_IIIb = (1 - len(uniqmers_IIIb - target_mers) / len(uniqmers_IIIb)) * 100
        p_VI = (1 - len(uniqmers_VI - target_mers) / len(uniqmers_VI)) * 100
        p_IV = (1 - len(uniqmers_IV - target_mers) / len(uniqmers_IV)) * 100
        p_VII = (1 - len(uniqmers_VII - target_mers) / len(uniqmers_VII)) * 100
        p_VIII = (1 - len(uniqmers_VIII - target_mers) / len(uniqmers_VIII)) * 100
        print('S. bongori: ' + str(round(p_bongori, 1)) + ', Subsp. I: ' + str(round(p_I, 1)) +
             ', Subsp. IIIb: ' +str(round(p_IIIb, 1)) +', subsp. IIa: '+ str(round(p_IIa, 1)) +
             ', subsp. IIb: '+ str(round(p_IIb, 1)) +', subsp. IIIa: ' + str(round(p_IIIa, 1)) +
             ', Subsp. IV: ' + str(round(p_IV, 1)) + ', Subsp. VI: ' + str(round(p_VI, 1)) +
              ', Subsp. VII: ' + str(round(p_VII, 1)) + ', Subsp. VIII (prov.): ' + str(round(p_VIII, 1)))

if __name__ == '__main__':
    main()
