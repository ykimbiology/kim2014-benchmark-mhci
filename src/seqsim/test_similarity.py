#! /usr/bin/python

"""
To run tests, issue the following command:

    nosetests test_similarity.py

"""

def get_data_sample():
    """
    There are two clusters of peptides. Within each cluster, all peptide pairs have at least 80% seq identity.
    """
    plist = []
    #Cluster 1:
    plist.append('HEART')
    plist.append('AEART')
    plist.append('CEART')
    plist.append('DEART')

    #Cluster 1:
    plist.append('MILKY')
    plist.append('MELKY')
    plist.append('MALKY')
    return plist


def test_count_neighbour():
    from util_similarity import count_neighbour
    plist  = get_data_sample()
    dic = count_neighbour(plist)
    print dic
    assert dic['HEART'] == 3
    assert dic['MILKY'] == 2
    #assert False

def test_check_bindingdata():
    from util_similarity import check_seq_list
    sim_cutoff = 0.80

    plist = []
    plist.append('MILKY')
    plist.append('MILKY')
    #plist.append('AELKY')
    is_error_a = check_seq_list(plist, sim_cutoff=sim_cutoff)

    print 'is_error_a', is_error_a, plist
    plist = []
    plist.append('MILKY')
    #plist.append('MILKY')
    plist.append('AELKY')
    is_error_b = check_seq_list(plist, sim_cutoff=sim_cutoff)
    print 'is_error_b', is_error_b, plist

    assert is_error_a == True
    assert is_error_b == False
    #assert False

def test_check_similarity():
    from util_similarity import seq_identity
    pepa = 'AAAVVVLLLG' # 10-mer.
    pepb = 'AAAVVVLLPP'
    pepc = 'AAAVVVLHHH'

    s1 = seq_identity(pepa,pepa)  #
    s2 = seq_identity(pepa,pepb)
    s3 = seq_identity(pepa,pepc)
    print s1,s2,s3
    assert s1 == 1.0
    assert s2 >= 0.80
    assert s3 < 0.80


def test_is_homologous():
    from util_similarity import is_similar
    pepa = 'AAAVVVLLLG' # 10-mer.
    pepb = 'AAAVVVLLPP'
    pepc = 'AAAVVVLHHH'

    s1 = is_similar(pepa,pepa)  #
    s2 = is_similar(pepa,pepb)
    s3 = is_similar(pepa,pepc)
    print s1,s2,s3
    assert s1 == True
    assert s2 == True
    assert s3 == False
    #assert False
