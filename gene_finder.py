# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: ISA BLANCETT

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    # A-->T, C-->G and vice versa
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'T':
        return 'A'


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
    rev_dna = "" # Initialize string
    for nucleotide in dna:
        rev_dna = rev_dna + get_complement(nucleotide) # Replace with complement one by one
    return rev_dna[::-1] # Return the reverse of complement string


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop_codons = ['TAG','TAA','TGA']

    for c in range(0,len(dna),3):
        if c+3 > len(dna): # Make sure enough codons left to make stop codon
            return dna # If not, return DNA string
        else:
            frame = dna[c:c+3] #Establish codon frame

        if frame not in stop_codons:
            return dna[:c] # if codon is stop codon, return everything before it

    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    all_ORFs = [] # Initialize ORF list
    c = 0 # Initialize DNA marker

    while c < len(dna):
        unread_dna = dna[c:] # What's left to read
        if unread_dna.startswith('ATG'):
            current_ORF = rest_of_ORF(unread_dna) # Get first ORF in sequence
            all_ORFs.append(current_ORF)
            c = c + len(current_ORF) # Move marker past ORF
        else:
            c = c + 3 # Move marker to next frame

    return all_ORFs



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []

    # Call find_all_ORFs_oneframe for all three frames
    for c in range(3):
        all_ORFs.extend(find_all_ORFs_oneframe(dna[c:]))

    # Return all frames together if ORF exists
    return [ORF for ORF in all_ORFs if ORF]

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    forward = find_all_ORFs(dna)
    backward = find_all_ORFs(get_reverse_complement(dna))
    return forward + backward # Return combined list of regular ORFs and reverse complement ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    if find_all_ORFs_both_strands(dna): # Make sure there are ORFs
        return max(find_all_ORFs_both_strands(dna), key=len)
    else:
        return " " # No sequence found


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    # Get list of ORFs, length of num_trials
    all_longest = [longest_ORF(shuffle_string(dna)) for i in range(num_trials)]

    if any(all_longest) != " ":
        return len(max(all_longest, key=len))  # If at least one ORF found, find the length of longest
    else:
        return 0 # No longest found



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
    # For every triple, look up codons in aa_table, join list into string
    codon_list = [aa_table[dna[c:c+3]] for c in range(0,len(dna),3) if c+3 <= len(dna)]
    return ''.join(codon_list)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500) # Set threshold
    all_ORFs = find_all_ORFs_both_strands(dna) # find ORFs of both strands
    # Return protein encoded ORFs longer than threshold
    return [coding_strand_to_AA(dna) for ORF in all_ORFs if len(all_ORFs)>threshold]


if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
