import sys
from random import *
import math

'''The code was developed by Peng Wang 2010 for MHC-II binding data sets.
   The code will also be applied to MHC-I.'''

def check_similarity(pepa, pepb):
	''''''
	find = -1
	end = len(pepa) - 8
	for i in range(end):
        	temp_pep = pepa[i:i+9]
#		print "===", temp_pep, pepa, pepb
	        if pepb.find(temp_pep) > -1:
#       		    print "===", temp_pep, pepa, pepb
        	    find = 1
        	    return find
	idaa = 0
	for i in range(len(pepa) - 9):
		sub_pepa = pepa[i:]
		temp_count = 0
		if len(sub_pepa) <= len(pepb):
			end = len(sub_pepa)
		else:
			end = len(pepb)
		for j in range(end):
			if sub_pepa[j] == pepb[j]:
				temp_count += 1
		if temp_count > idaa:
			idaa = temp_count
	for i in range(len(pepb) - 9):
		sub_pepb = pepb[i:]
		temp_count = 0
		if len(sub_pepb) <= len(pepa):
			end = len(sub_pepb)
		else:
			end = len(pepa)
		for j in range(end):
			if sub_pepb[j] == pepa[j]:
				temp_count += 1
		if temp_count > idaa:
			idaa = temp_count
	if len(pepb) <= len(pepa):
		t_len = len(pepb)
	else:
		t_len = len(pepa)
#	print idaa, t_len
	perc = float(idaa)/float(t_len)
#	print perc
	if perc >= 0.80:
		find = 1
	return find

def reduce_seq(all_seq, seq_dict):
	''''''
	count = 0
	temp_seq = []
	ref_seq = []
	sorted_seq = []
	forward_seq = []
	for i in all_seq:
		temp_seq.append(i)
		ref_seq.append(i)
	sorted_list = sorted(seq_dict.iteritems(), key=lambda (k,v):(v,k))
	for (k, v) in sorted_list:
#		print "sorted ", k, v
		sorted_seq.append(k)
	forward_seq.append(sorted_seq[0])
	for i in sorted_seq:
		match = -1
		for j in forward_seq:
			dd = check_similarity(i, j)
			if dd == 1:
				match = 1
				break;
		if match == -1:
			forward_seq.append(i)
	for i in forward_seq:
		for j in forward_seq:
			dd = check_similarity(i, j)
			if dd == 1 and i != j:
				print "error: ", i, j
	temp_dict = {}
	for i in forward_seq:
		temp_dict[i] = 1
	final_seq = temp_dict.keys()
        return final_seq

def main():
#    a = "EKKYFAAEQFEPLAA"
#    b = "EKKYFAATIFEPLAA"
#    dd = check_similarity(a, b)
#    print dd, a, b
  for kk in range(10):
    cv_num = 5
    list_file = "class_II_list.txt"
#    list_file = "temp.txt"
    # do random split
    lines = open(list_file, "r").read().strip().split("\n")
    for line in lines:
        data_file = "class_II/" + line
        dataall = open(data_file, "r").read().strip().split("\n")
	shuffle(dataall)
        all_seq_binder = []
	temp_seq_binder = []
	data_dict_binder = {}
	seq_dict_binder = {}
        all_seq_nonbinder = []
	temp_seq_nonbinder = []
	data_dict_nonbinder = {}
	seq_dict_nonbinder = {}
	ref_seq = []
	data_dict = {}
        for data in dataall:
            elements = data.split("\t")
            allele = elements[1].strip()
	    if len(elements[4]) >= 15 and float(elements[6]) >= 1000:
	    	all_seq_nonbinder.append(elements[4].strip())
		data_dict_nonbinder[elements[4].strip()] = data
		data_dict[elements[4].strip()] = data
	    elif len(elements[4]) >= 15 and float(elements[6]) < 1000:
	    	all_seq_binder.append(elements[4].strip())
		data_dict_binder[elements[4].strip()] = data
		data_dict[elements[4].strip()] = data

	for i in all_seq_binder:
		match = 0
		for j in all_seq_binder:
			dd = check_similarity(i, j)
			if dd == 1 and i != j:
				match += 1
		seq_dict_binder[i] = match
	sorted_list_binder = sorted(seq_dict_binder.iteritems(), key=lambda (k,v):(v,k))
	for (k, v) in sorted_list_binder:
#		print "sorted ", k, v
		ref_seq.append(k)

	temp_seq_binder = reduce_seq(all_seq_binder, seq_dict_binder)
	for i in temp_seq_binder:
		try:
			ref_seq.remove(i)
		except ValueError:
			print "not exist ", allele, i

	for i in all_seq_nonbinder:
		match = 0
		for j in all_seq_nonbinder:
			dd = check_similarity(i, j)
			if dd == 1 and i != j:
				match += 1
		seq_dict_nonbinder[i] = match
	sorted_list_nonbinder = sorted(seq_dict_nonbinder.iteritems(), key=lambda (k,v):(v,k))
	for (k, v) in sorted_list_nonbinder:
#		print "sorted ", k, v
		ref_seq.append(k)

	temp_seq_nonbinder = reduce_seq(all_seq_nonbinder, seq_dict_nonbinder)
	for i in temp_seq_nonbinder:
		try:
			ref_seq.remove(i)
		except ValueError:
			print "not exist ", allele, i

	reduced_data_1000 = []
	for i in temp_seq_binder:
		reduced_data_1000.append(data_dict_binder[i])
	for i in temp_seq_nonbinder:
		reduced_data_1000.append(data_dict_nonbinder[i])

	for i in ref_seq:
		bstatus = 0
		nonbstatus = 0
		simb = []
		simnonb = []
		for dataj in reduced_data_1000:
			j = dataj.split("\t")[4].strip()
			if check_similarity(i, j) == 1 and i != j:
				ref_ic = float(data_dict[j].split("\t")[6])
				temp_ic = float(data_dict[i].split("\t")[6])
				if ref_ic >= 1000:
					nonbstatus = nonbstatus + 1
					tt = j + " " + str(ref_ic)
					simnonb.append(tt)
				elif ref_ic < 1000:
					bstatus = bstatus + 1
					tt = j + " " + str(ref_ic)
					simb.append(tt)
		if (bstatus == 0 and temp_ic <1000) or (nonbstatus == 0 and temp_ic >= 1000):
			reduced_data_1000.append(data_dict[i])
			print "error missed"

	for dataj in reduced_data_1000:
		j = dataj.split("\t")[4].strip()
		for datai in reduced_data_1000:
			i = datai.split("\t")[4].strip()
			if check_similarity(i, j) == 1 and i != j:
				ref_ic = float(data_dict[j].split("\t")[6])
				temp_ic = float(data_dict[i].split("\t")[6])
				if (ref_ic <1000 and temp_ic <1000) or (ref_ic >=1000 and temp_ic >= 1000):
					print "error", datai, dataj



	shuffle(reduced_data_1000)
	print allele, len(reduced_data_1000)
#        allele = elements[1]
#	temp_allele = allele.replace("/", "-")
#        for cv in range(cv_num):
#            train_set = []
#            test_set = []
#            for index in range(len(reduced_data_1000)):
#                if index % cv_num == cv:
#                    test_set.append(reduced_data_1000[index])
#                else:
#                    train_set.append(reduced_data_1000[index])
#            tt = len(train_set) + len(test_set)
#            print cv, tt
#            test_ofile = "class_II_similarity_reduced_maxsize_1000_5cv_sep/" + temp_allele + "_test_nooverlap_" + str(cv) + ".txt"
#            train_ofile = "class_II_similarity_reduced_maxsize_1000_5cv_sep/" + temp_allele + "_train_nooverlap_" + str(cv) + ".txt"
#            test_out = open(test_ofile, "w")
#            train_out = open(train_ofile, "w")
#            for item in test_set:
#                test_out.write(item)
#                test_out.write("\n")
#            for item in train_set:
#                train_out.write(item)
#                train_out.write("\n")
#            test_out.close()
#            train_out.close()

if __name__ == '__main__':
    main()
