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



#Modify secondary structure Toehold-Trigger prediction
def structure_define(miRNA,loop,linker,reporter,paired, unpaired):

    Final = 3 + 11 + len(linker)
    lenLoop = len(loop)
    Toehold_structure = '.3.%d(%d.3(5.%d)5.3)%d.6(3.4)3.5' % (unpaired, paired, lenLoop, paired)
    Trig_Toe_struct = "(%d+.3)%d.3(5.%d)5.%d" % (len(miRNA), paired+unpaired, lenLoop, Final)

    a_dom = "N3"
    b_dom = "W5"
    c_dom = "W5"
    d_dom = "N%d" % (paired)

    aaa = Domain("AAA", name = "aaaa")
    complement = Domain(miRNA, name = "complement")
    a = Domain(a_dom, name="a")
    b = Domain(b_dom, name="b")
    Loop = Domain(loop, name="Loop")
    c = Domain(c_dom, name="c")
    start = Domain("AUG", name = "start")
    d = Domain(d_dom, name = "d")
    Linker_dom = Domain(linker, name="Linker_dom")
    #Reporter_dom = Domain(reporter, name="Reporter_dom")

    Domains = [aaa, ~complement, a, b, Loop, c, start, d , Linker_dom]
    #Prevent1 = Pattern(["UAG", "UGA", "UAA"], scope = [d])
    Prevent2 = Pattern(['A4', 'C4', 'G4', 'U4', 'M6', 'K6', 'W6', 'S6', 'R6', 'Y6'])
    soft_cons = [Prevent2]

    return Toehold_structure,Trig_Toe_struct, Domains, len(d)+len(linker), soft_cons


def rep_stop_codons(record, dist, codon_stop_array):
    modify = False
    tempRecordSeq = list(record)

    n = len(record) - dist
    for index in range(n, len(record), 3):
        codon = record[index:index + 3]
        if codon in codon_stop_array:
            tempRecordSeq[index:index + 3] = 'V', 'N', 'N'
            modify = True
    record = "".join(tempRecordSeq)

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


 #   def score()