import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#This function help to design toehold complement
def displacement_miR(miRNA):
    length = len(miRNA)

    remain = 26 - length
    if remain < 0:
        return miRNA[:remain] , 0
    return miRNA, remain

#Modify secondary structure Toehold-Trigger prediction
def structure_define(miRNA,loop,linker):

    tail = displacement_miR(miRNA)[1]
    paired = 26 - tail
    Loop = len(loop)
    Final = len(linker) + 3 + 11

    if tail > 0:
      return  "(%d+.%d)%d.3(5.%d)5.%d" % (len(miRNA),tail,paired,Loop,Final)
    elif tail == 0:
      return "(%d+)%d.3(5.%d)5.%d" % (paired, paired, Loop, Final)


def rep_stop_codons(record, codon_stop_array):
    modify = False
    tempRecordSeq = list(record)
    for index in range(53, len(record), 3):
        codon = record[index:index + 3]
        if codon in codon_stop_array:
            tempRecordSeq[index:index + 3] = 'V', 'N', 'N'
            modify = True
    record = "".join(tempRecordSeq)
    return record, modify

#def aval_binding(structure):
