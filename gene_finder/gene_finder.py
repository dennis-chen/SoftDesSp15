# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Dennis Chen

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
from load import load_seq

def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))

### YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    complement_dict = {'A':'T','T':'A','C':'G','G':'C'}
    return complement_dict[nucleotide]

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    complement_list = map(get_complement,dna)
    complement_list.reverse()
    return ''.join(complement_list)

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop_codons = ['TAG','TAA','TGA']
    no_stop_found = True
    index = 0
    while no_stop_found and index < len(dna):
        if dna[index:index+3] in stop_codons:
            no_stop_found = False
        else:
            index +=3
    return dna[:index]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    all_ORFS = []
    index = 0
    while index < len(dna):
        if dna[index:index+3] == 'ATG':
            ORF = rest_of_ORF(dna[index:]) 
            all_ORFS.append(ORF)
            index += len(ORF)
        else:
            index += 3
    return all_ORFS

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFS = []
    for i in range(3):
        all_ORFS += find_all_ORFs_oneframe(dna[i:])
    return all_ORFS

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_complement = get_reverse_complement(dna)
    return find_all_ORFs(dna) + find_all_ORFs(reverse_complement)  

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    return max(find_all_ORFs_both_strands(dna), key=len)

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_ORF_len = 0
    for i in range(num_trials):
        shuffled = shuffle_string(dna)
        ORF_len = len(longest_ORF(shuffled))
        if  ORF_len > longest_ORF_len:
            longest_ORF_len = ORF_len
    return longest_ORF_len

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    AA = ''
    for i in range(0,len(dna) - len(dna)/3,3):
        codon = dna[i:i+3]
        AA += aa_table[codon]
    return AA

def gene_finder(dna):
<<<<<<< HEAD
    """ Returns the amino acid sequences coded by all genes that have an ORF
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold =  longest_ORF_noncoding(dna,1500)
    ORFS = find_all_ORFs_both_strands(dna)
    long_ORFS = [i for i in ORFS if len(i) > threshold]
    AA = []
    for strand in coding_strands:
        AA.append(coding_strand_to_AA(strand))
    return AA

def get_threshold():
    """Returns a conservative threshold to use to get ORFS. Prints 789"""
    dna = load_seq('./data/X73525.fa')
    return longest_ORF_noncoding(dna,1500)

def run_gene_finder():
    """Loads gene and returns long_ORFS"""
    dna = load_seq('./data/X73525.fa')
    amino_acids = gene_finder(dna)
    return amino_acids
    
if __name__ == "__main__":
    #print get_threshold()
    print run_gene_finder()
    #import doctest
    #doctest.testmod()
