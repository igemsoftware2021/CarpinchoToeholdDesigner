import pandas as pd
import numpy as np
from nupack import *
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import RNA

#This function help to design toehold complement
def displacement_miR(miRNA,type):
    length = len(miRNA)
    if type == "c":
        remain = 26 - length
    elif type == "a2":
        remain = 30 - length
    elif type == "a1":
        remain = 30 - length

    if remain < 0:
        return miRNA[:remain] , 0
    return miRNA, remain

#Modify secondary structure Toehold-Trigger prediction
def structure_define(miRNA,loop,linker,reporter,type):

    tail = displacement_miR(miRNA, type)[1]* "N"
    lenLoop = len(loop)

    if type == "a2":
        paired = 30 - len(tail)
        Final = 6 + 9 + len(linker) + len(reporter)
        Toehold_structure = '.15(9.3(6.%d)6.3)9.6(3.4)3.%d' % (lenLoop, len (reporter)+ 5)
        mir_seq = displacement_miR(miRNA, type)[0]

        a_dom = ""
        b_dom = "WWW"
        c_dom = "WWWNNN"
        d_dom = "NNNNNNSSW"

        if len(tail) > 0:
            Trig_Toe_struct = "(%d+.%d)%d(3.%d)3.%d" % (len(miRNA), len(tail), paired, lenLoop, Final)


        elif len(tail) == 0:
            Trig_Toe_struct = "(%d+)%d.3(3.%d)3.%d" % (paired, paired, lenLoop, Final)

    if type == "c":

        paired = 26 - len(tail)
        Final = 3 + 11 + len(linker) + len(reporter)
        mir_seq = displacement_miR(miRNA, type)[0]
        Toehold_structure = '.15(11.3(5.%d)5.3)11.6(3.4)3.%d' % (lenLoop, len(reporter)+5)

        a_dom = "N3"
        b_dom = "WWWWW"
        c_dom = "WWWWW"
        d_dom = "N11"

        if len(tail) > 0:
            Trig_Toe_struct = "(%d+.%d)%d.3(5.%d)5.%d" % (len(mir_seq),len(tail),paired,lenLoop,Final)


        elif len(tail) == 0:
            Trig_Toe_struct = "(%d+)%d.3(5.%d)5.%d" % (paired, paired, lenLoop, Final)




    free_tail = Domain(tail, name="freetail")
    complement = Domain((mir_seq), name="complement")

    a = Domain(a_dom, name="a")
    b = Domain(b_dom, name="b")
    Loop = Domain(loop, name="Loop")
    c = Domain(c_dom, name="c")
    start = Domain("AUG", name = "start")
    d = Domain(d_dom, name = "d")
    Linker_dom = Domain(linker, name="Linker_dom")
    Reporter_dom = Domain(reporter, name="Reporter_dom")

    Domains = [free_tail, ~complement, a, b, Loop, c, start, d , Linker_dom, Reporter_dom]

    return Toehold_structure,Trig_Toe_struct, Domains, len(d)+len(linker)+len(reporter)



def rep_stop_codons(record, dist, codon_stop_array):
    modify = False
    tempRecordSeq = list(record)
    print(dist)
    n = len(record) - dist
    for index in range(n, len(record), 3):
        codon = record[index:index + 3]
        if codon in codon_stop_array:
            tempRecordSeq[index:index + 3] = 'V', 'N', 'N'
            modify = True
    record = "".join(tempRecordSeq)
    print(record)
    return record, modify


def save_file(CodonOpt, path_versions, seq, structure, prob, struct_theory,  format, model1):

    file_struct = (path_versions + "/%s_%s_struct.png") % (format, CodonOpt)
    file_prob = (path_versions + "/%s_%s_probab.png") % (format, CodonOpt)
    file_prob_100 = (path_versions + "/%s_%s_probab_theory.png") % (format, CodonOpt)
    name_fasta = (path_versions + "/%s_%s.fx") % (format,CodonOpt)



    output_file = ">%s\n%s\n%s" % (CodonOpt, seq, structure)
    with open(name_fasta, 'w') as f:
        f.write(''.join(output_file))



    if format == 'Toehold':
        plt.figure(figsize=(20, 20))
        cg = forgi.load_rna(name_fasta, allow_many=False)
        fvm.plot_rna(cg, text_kwargs={"fontweight": "black"}, lighten=0.7,
                 backbone_kwargs={"linewidth": 3})

        fig = plt.gcf()

        fig.savefig(file_struct, dpi=100)
        fig.clf()

    plt.imshow(prob.to_array())
    plt.xlabel('Base index')
    plt.ylabel('Base index')
    plt.title('Pair probabilities')
    plt.colorbar()
    plt.clim(0, 1)
    plt.savefig(file_prob, dpi = 100)
    plt.clf()

    theory_struct = des(structure=struct_theory, model=model1)
    prob_theory = pairs(strands=theory_struct, model = model1)
    plt.imshow(prob_theory.to_array())
    plt.xlabel('Base index')
    plt.ylabel('Base index')
    plt.title('Pair probabilities')
    plt.colorbar()
    plt.clim(0, 1)
    plt.savefig(file_prob_100, dpi=100)
    plt.clf()