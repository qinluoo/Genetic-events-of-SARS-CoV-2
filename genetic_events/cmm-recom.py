# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 16:02:43 2021

@author: qinluyao
"""


from sys import argv
##sequence file
fs = open(argv[1], "r")
##reference file reference.fasta
fr = open(argv[2], "r")
##co-mutations file co_mutations
fm = open(argv[3], "r") 
##group_mtations file group_mut
fx = open(argv[4], "r")
fw_recom = open(argv[5], "w")
fw_mut = open(argv[6], "w")
fw_group = open(argv[7], "w")
fw_g0 = open(argv[8], "w")


seqs = fs.readlines()
co_mutations = fm.readlines()
co_mutation = co_mutations[0].strip().split("\t")
co_mutation_mark = co_mutations[1].strip().split("\t")
group_mutations = fx.readlines()

recom_mut = []
recom_group = []
recom_name = []

s = 1
while s < len(seqs):
    iteration = 0
    mutations = []
    mutation_marks = []
    m = 0
    while m < len(co_mutation):
        if seqs[s][int(co_mutation[m].split("_")[0])-1].upper() == co_mutation[m].split("_")[1]:
            mutations.append(co_mutation[m])
            mutation_marks.append(co_mutation_mark[m])
        m = m + 1
    mutation_marks_unmodify = mutation_marks[:]
#check whether groups are mutually exclusive 
    while iteration == 0:
        t = 0
        while t < len(mutation_marks):
            count = 0
            length = len(mutation_marks[t].split("."))
            for mutation_mark in mutation_marks:
                if len(mutation_mark.split(".")) > length and mutation_marks[t] == ".".join(mutation_mark.split(".")[:length]):
                    count = 1
                    break
            if count == 1:
                break
            else:
                t = t + 1
                

#if not matually exclusive, change the group into a higher level                    
        if t < len(mutation_marks):
            mutation_marks_start = mutation_marks[:]
            mark = 0
            while mark < len(mutation_marks)-1:
                l1 = len(mutation_marks[mark].split("."))
                l2 = len(mutation_marks[mark+1].split("."))
                if l1 > l2 and mutation_marks[mark+1] == ".".join(mutation_marks[mark].split(".")[:l2]):
                    mutation_marks[mark+1] = mutation_marks[mark]
                elif l1 < l2 and mutation_marks[mark] == ".".join(mutation_marks[mark+1].split(".")[:l1]):
                    mutation_marks[mark] = mutation_marks[mark+1]
                mark = mark + 1
            if mutation_marks_start != mutation_marks:
                iteration = 0
            else:
                iteration = 1
        else:
            break
    if iteration == 1:
        for mutation_mark in mutation_marks:
            mark = 0
            while mark < len(mutation_marks):
                if mutation_mark in mutation_marks[mark] and mutation_mark != mutation_marks[mark]:
                    muation_mark = mutation_marks[mark]                        
                    break
                else:
                    mark = mark + 1

        
    mutations_original = mutations[:]
    mutation_marks_original = mutation_marks[:]
    
    if len(set(mutation_marks_original)) == 1:
        fw_group.write(seqs[s-1].strip() + "\t" + mutation_marks_original[0] + "\n")
        
    if len(set(mutation_marks_original)) == 0:
        fw_g0.write(seqs[s-1].strip() + "\t" + "G0" + "\n")
    
    if len(set(mutation_marks_original)) > 2:
        fw_mut.write(seqs[s-1].strip() + "\t" + max(mutation_marks_original, key = mutation_marks_original.count) + "\n")
        
    if len(set(mutation_marks_original)) == 2:
        if all(mutation_marks_unmodify.count(mutation_mark_original) > 1 for mutation_mark_original in list(set(mutation_marks_original))):
            if all(mutation_marks_original.count(mutation_mark_original) > 2 for mutation_mark_original in list(set(mutation_marks_original))):
                coti = 0
            else:
                fw_mut.write(seqs[s-1].strip() + "\t" + max(mutation_marks_original, key = mutation_marks_original.count) + "\n")
                coti = 2
            while coti == 0:
                nr_mutation_marks = list(set(mutation_marks))
                nr_mutation_marks.sort(key=mutation_marks.index)
                if len(set(mutation_marks)) == 2:
                    number = mutation_marks.index(nr_mutation_marks[1]) - mutation_marks.index(nr_mutation_marks[0])
                    if number > 1:
                        for group_mutation in group_mutations:
                            if group_mutation.strip().split("\t")[0] == nr_mutation_marks[0]:
                                number_all = (group_mutation.strip().split("\t")[1].split(",").index(mutations[mutation_marks.index(nr_mutation_marks[1])-1]) -
                                              group_mutation.strip().split("\t")[1].split(",").index(mutations[mutation_marks.index(nr_mutation_marks[0])]) + 1)
                                if number/number_all >= 0.9:
                                    del mutations[:mutation_marks.index(nr_mutation_marks[1])]
                                    del mutation_marks[:mutation_marks.index(nr_mutation_marks[1])]   
                                else:
                                    coti = 1
                                break
                    else:
                        coti = 1
                else:
                    break

            if coti == 0:
                for group_mutation in group_mutations:
                    if mutation_marks[0] == group_mutation.strip().split("\t")[0]:
                        number = len(mutation_marks)
                        number_all = group_mutation.strip().split("\t")[1].split(",").index(mutations[-1]) - group_mutation.strip().split("\t")[1].split(",").index(mutations[0]) + 1                    
                        if number > 1 and number/number_all >= 0.9:
                            fw_recom.write(seqs[s-1] + "\t".join(mutations_original) + "\n" + "\t".join(mutation_marks_original) + "\n")
                            break
                        else:
                            fw_mut.write(seqs[s-1].strip() + "\t" + max(mutation_marks_original, key = mutation_marks_original.count) + "\n")
                            break
            if coti == 1:
                fw_mut.write(seqs[s-1].strip() + "\t" + max(mutation_marks_original, key = mutation_marks_original.count) + "\n")
        else:
            nr_mutation_marks_original = list(set(mutation_marks_original))
            if mutation_marks_unmodify.count(nr_mutation_marks_original[0]) >= mutation_marks_unmodify.count(nr_mutation_marks_original[1]):
                fw_mut.write(seqs[s-1].strip() + "\t" + nr_mutation_marks_original[0] + "\n")
            else:
                fw_mut.write(seqs[s-1].strip() + "\t" + nr_mutation_marks_original[1] + "\n")
    s = s + 2

fw_recom.close()
fw_mut.close()
fw_group.close()
fw_g0.close()
               
            
                 
                                                                                                   
                                                                                                   
                     
                                                                                 
                                                                             
        
        
                        
                    
            
            
                