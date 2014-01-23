#! /usr/bin/python


"""
Demonstrate how to run 'reduce similarity' algorithm.

INPUT: a list of peptides of same length.
OUTPUT: a list of peptides with similar peptides removed.
"""

import seqsim
from seqsim.test_similarity import get_data_sample
from seqsim.generate_similarity_reduced import get_peptide_list_reduced

sim_cutoff = 0.80 # Two peptides are similar if at least 80% sequence identity.

peptide_list    = get_data_sample()
peptide_list_sr = get_peptide_list_reduced(peptide_list, sim_cutoff=sim_cutoff)

print 'sim_cutoff =', sim_cutoff
print 'INPUT =', peptide_list
print 'OUTPUT =', peptide_list_sr
