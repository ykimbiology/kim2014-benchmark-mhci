#! /usr/bin/python

"""
TODO
  Rename to util_similarity.py
"""

def get_distance_seq(pepa, pepb):
    """
    [0,1], where 1 indicates largest distnace and 0 means the two are teh same.
    """
    distance = 1.0
    if len(pepa) == len(pepb):
        xid = seq_identity(pepa,pepb)
        distance = 1.0 - xid  # To get distance.
    return distance
    
def seq_identity(pepa, pepb):
    """
    Two peptides will be considered similar if 80% sequence identity is observed.
    """
    assert len(pepa) == len(pepb)
    num_identity = sum([1 for (a,b) in zip(pepa,pepb) if a==b])
    fraction_identity = float(num_identity)/len(pepa)
    #print 'fraction_identity', fraction_identity, pepa, pepb
    return fraction_identity

def is_similar(pepa, pepb, sim_cutoff=0.80):
    """
    previously named is_homologous.
    """
    is_sim = False
    s = seq_identity(pepa,pepb)
    if s >= sim_cutoff: is_sim =  True
    #print 'is_similar', is_sim, s, pepa, pepb
    return is_sim

def count_neighbour(seq_list, sim_cutoff=0.80):
    """
    This version should be fast by a factor of 2. Checked this. See its test. use binary version of peptides? Only calcualte half of matrix?
    For a given peptide, how many similar peptides are there? Exclude self count.
    """
    d_count_neighbour = {}  # [sequence] = number of neighbours.
    seq_list_len = len(seq_list)
    for i in range(seq_list_len):
        pepi = seq_list[i]
        for j in range(i, seq_list_len):
            pepj = seq_list[j]
            count = 0
            if (is_similar(pepi, pepj, sim_cutoff=sim_cutoff) == True) and (i != j): count = 1
            d_count_neighbour[pepi] = d_count_neighbour.get(pepi, 0) + count
            d_count_neighbour[pepj] = d_count_neighbour.get(pepj, 0) + count
    return d_count_neighbour

def check_seq_list(seq_list, sim_cutoff):
    """
    Makes sure there are no similar peptides for each group: (a) binders, (b) nonbinders.
    """
    is_error=False
    seq_list_len = len(seq_list)
    for i in range(seq_list_len):
        pepi = seq_list[i]
        for j in range(i, seq_list_len):
            pepj = seq_list[j]
            if (i != j) and (is_similar(pepi, pepj, sim_cutoff=sim_cutoff)==True): is_error = True; break
    return is_error

def get_d_bindingdata(data_bd):
    """
    Returns a dictionary d[peptide] = a row of (species, mhc, length, cv, peptide, inequality, meas).
    """
    dic = {}
    for data in data_bd:
        (species, mhc, length, cv, peptide, inequality, meas) = data
        dic[peptide] = data
    return dic
