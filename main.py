

#in terminal
#python3 -m pip install -U nupack -f ~/Downloads/nupack-4.0.0.27/package
from Bio.Seq import Seq
from nupack import *
from test import *
from helper import *
import numpy as np
from random import randint
import pandas as pd
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import RNA
import os

path = "/home/samuel/Downloads/ResultsToehold"
os.makedirs(path, exist_ok=True)
'''''''''
def random(s):
    length = len(s)
    seq = list(s)

    nucl = "AUGC"
    lengthnucl = len(nucl)

    position_orgseq = randint(0, length - 1)
    position_nucl = randint(0, lengthnucl - 1)
    while seq[position_orgseq] == nucl[position_nucl]:
        position_orgseq = randint(0, length - 1)
        position_nucl = randint(0, lengthnucl - 1)

    seq[position_orgseq] = nucl[position_nucl]
    final = "".join(seq)

    return s, final
'''''
#Conditions
my_options = DesignOptions(f_stop=0.01,
                           seed=0,
                           wobble_mutations= Wobble)

my_soft_const = [Pattern(Prevent)]


model1 = Model(material = Material,
               celsius= T,
               sodium= Sodium,
               magnesium= Magnesium)

#Domains
mir_seq = displacement_miR(miRNA)[0]
tail = displacement_miR(miRNA)[1]*"N"

trigger = Domain(miRNA, name = "miRNA")
free_tail = Domain(tail, name = "freetail")
complement = Domain((mir_seq), name = "complement")

c = Domain("N3", name = "c")
d = Domain("W5", name = "d")
loop = Domain(Loop , name = "loop")
g = Domain("W5" , name = "g")
h = Domain("AUG" , name = "h")
i = Domain("N11", name= "i")
linker = Domain(Linker, name= "linker")
#reporter = Domain(Reporter, name = "reporter")

#Toehold Switch Complex - Define a estrutura do toehold

toe_strand = TargetStrand([free_tail,~complement,c,d,loop,g,h,i,linker], name='Toehold')

C1 = TargetComplex([toe_strand], Toehold_structure, name="C1")
toe_complex = complex_design(complexes = [C1],
                             model=model1,
                             soft_constraints=my_soft_const,
                             options=my_options)


toehold_complex = toe_complex.run(trials= Trials)




results_list = []

for trials in range(Trials):

    Trigger = TargetStrand([trigger], name="Trigger")
    toeholdcomplex = toehold_complex[trials].to_analysis(C1)

    codon = rep_stop_codons(str(toeholdcomplex),StopCodons)

    if codon[1]:
        toeholdcomplex = codon[0]
        #CodonOpt = codon[1]





    toehold_dom = Domain(str(toeholdcomplex), name = "toehold_dom")
    toehold_strand = TargetStrand([toehold_dom], name = "toehold_strand")

    #Toehold-Trigger duplex
    trig_complex = TargetComplex([Trigger], 'U21', name="Trig_complex")

    Trig_Toe_struct = structure_define(miRNA, Loop, Linker)

    C2 = TargetComplex([Trigger, toehold_strand], Trig_Toe_struct, name = "C2")

    t_toe = TargetTube(on_targets={C2: 1e-06},
                       off_targets=SetSpec(max_size=2),
                       name='t_toe')

    my_tube_design = tube_design(tubes=[t_toe],
                                 soft_constraints=my_soft_const,
                                 model=model1, options=my_options)

    tubetrials = 1

    if codon[1]:
        tubetrials = 5
    results= my_tube_design.run(trials=tubetrials)

    path_trial = (path + "/%s") % (trials)
    os.makedirs(path_trial, exist_ok=True)

    for tubetri in range(tubetrials):

        switchfinal = results[tubetri].to_analysis(toehold_strand)

        if rep_stop_codons(str(switchfinal), StopCodons)[1]:
            print("Tube Design Version Failed - Stop Codon detected")
            continue

        binding = results[tubetri].to_analysis(C2)
        triggerfinal = results[tubetri].to_analysis(Trigger)
        my_mfe_trigswitch = mfe(strands=binding, model=model1)
        my_mfe_switch = mfe(strands = switchfinal, model = model1)
        my_mfe_trigger = mfe(strands=triggerfinal, model=model1)
        my_mfe_rbslinker = mfe(strands=str(switchfinal)[29:], model=model1)
        diff_mfe = my_mfe_trigswitch[0].energy - (my_mfe_trigger[0].energy + my_mfe_switch[0].energy)
        print(my_mfe_switch[0].structure)

        if not "."*15 in str(my_mfe_switch[0].structure)[:15]:
            continue
        print(my_mfe_trigswitch[0].structure)

        print((my_mfe_trigswitch[0].structure).pairlist())

        CodonOpt = "%d.%d" % (trials,tubetri)
        print(my_mfe_switch[0].structure)
        print(switchfinal)
        print(CodonOpt)


        path_versions = (path_trial + "/%s") % (CodonOpt)

        os.makedirs(path_versions, exist_ok=True)

        name_file = (path_versions + "/Toehold_%s.png") % (CodonOpt)
        name_fasta = (path_versions + "/%s.fx") % (CodonOpt)

        output_file = ">%s\n%s\n%s" % (CodonOpt, str(switchfinal), str(my_mfe_switch[0].structure))


        with open(name_fasta, 'w') as f:
            f.write(''.join(output_file))



        cg = forgi.load_rna(name_fasta, allow_many=False)
        fvm.plot_rna(cg, text_kwargs={"fontweight": "black"}, lighten=0.7,
                     backbone_kwargs={"linewidth": 3})

        fig = plt.gcf()
        fig.savefig(name_file, dpi=100)
        fig.clf()
        #Append results
        appends = [CodonOpt, miRNA,str(Seq(mir_seq).reverse_complement()), str(switchfinal),
               str(my_mfe_switch[0].structure), str(my_mfe_trigger[0].structure),
               str(my_mfe_trigswitch[0].structure),str(my_mfe_trigger[0].energy),
               str(my_mfe_switch[0].energy), str(my_mfe_trigswitch[0].energy),
               str(my_mfe_rbslinker[0].energy), str(diff_mfe)#, CodonOpt
               ]
        results_list.append(appends)

results_tab = pd.DataFrame(results_list, columns=["Trial", "miRNA.sequence",
                                                  "miRNA.Toehold.Complement", "Full.Toehold.seq",
                                                  "Toehold.structure", "Trigger.structure",
                                                  "Full.TriggerToehold.structure",
                                                  "Trigger.energy", "Toehold.energy",
                                                  "BindingToeholdTrigger.energy", "DeltaG_RBSLinker",
                                                  "Diff_energy"]).set_index("Trial")
print(results_tab)

results_tab.to_csv(r'test_toehold1.csv',index = True, header = True)



