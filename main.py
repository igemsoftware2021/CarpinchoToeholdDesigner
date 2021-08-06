

#in terminal
#python3 -m pip install -U nupack -f ~/Downloads/nupack-4.0.0.27/package
from Bio.Seq import Seq
from nupack import *
from Input import *
from helper import *
import numpy as np
from random import randint
import pandas as pd
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import RNA
import os

path = (PATH + "/ResultsToehold")
os.makedirs(path, exist_ok=True)


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

mir_seq = displacement_miR(miRNA, structure_type)[0]
trigger = Domain(mir_seq, name= "trigger")
#Toehold Switch Complex - Define a estrutura do toehold
Structure_def = structure_define(miRNA, Loop, Linker, Reporter, structure_type)

toe_strand = TargetStrand(Structure_def[2], name='Toehold')

C1 = TargetComplex([toe_strand], Structure_def[0], name="C1")
toe_complex = complex_design(complexes = [C1],
                             model=model1,
                             soft_constraints=my_soft_const,
                             options=my_options)


toehold_complex = toe_complex.run(trials= Trials)




results_list = []

for trials in range(Trials):

    Trigger = TargetStrand([trigger], name="Trigger")
    toeholdcomplex = toehold_complex[trials].to_analysis(C1)

    codon = rep_stop_codons(str(toeholdcomplex), Structure_def[3],StopCodons)

    if codon[1]:
        toeholdcomplex = codon[0]
        #CodonOpt = codon[1]





    toehold_dom = Domain(str(toeholdcomplex), name = "toehold_dom")
    toehold_strand = TargetStrand([toehold_dom], name = "toehold_strand")

    #Toehold-Trigger duplex
    trig_complex = TargetComplex([Trigger], 'U21', name="Trig_complex")



    C2 = TargetComplex([Trigger, toehold_strand], Structure_def[1], name = "C2")

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

        if rep_stop_codons(str(switchfinal), Structure_def[3], StopCodons)[1]:
            print("Tube Design Version Failed - Stop Codon detected")
            continue

        binding = results[tubetri].to_analysis(C2)
        triggerfinal = results[tubetri].to_analysis(Trigger)
        my_mfe_trigswitch = mfe(strands=binding, model=model1)
        my_mfe_switch = mfe(strands = switchfinal, model = model1)
        my_mfe_trigger = mfe(strands=triggerfinal, model=model1)
        my_mfe_rbslinker = mfe(strands=str(switchfinal)[29:], model=model1)
        #my_complex_ensemble_defect = defect(strands=[switchfinal],
        #                                   structure= Structure(Structure_def[0]), model = model1)
        diff_mfe = my_mfe_trigswitch[0].energy - (my_mfe_trigger[0].energy + my_mfe_switch[0].energy)


        if not "."*15 in str(my_mfe_switch[0].structure)[:15]:
            continue

        CodonOpt = "%d.%d" % (trials,tubetri)

        path_versions = (path_trial + "/%s") % (CodonOpt)

        os.makedirs(path_versions, exist_ok=True)

        probability_matrix_toehold = pairs(strands = str(switchfinal),model = model1)
        probability_matrix_binding= pairs(strands=str(binding), model=model1)

        save_file(CodonOpt, path_versions, str(binding), str(my_mfe_trigswitch[0].structure),
                                probability_matrix_binding, Structure_def[1], "Binding_Toehold_Trigger", model1)

        save_file(CodonOpt, path_versions, str(switchfinal), str(my_mfe_switch[0].structure),
                                probability_matrix_toehold, Structure_def[0], "Toehold", model1)
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



