
########################### GETTING START########################################
######## Download the NUPACK package using the command:
######## python3 -m pip install -U nupack -f ~/Downloads/nupack-4.0.0.27/package
######## for more information go to: https://docs.nupack.org/start/
#################################################################################



########################## Import Libraries ##################################
from Bio.Seq import Seq                    #v1.79
from nupack import *                       #v4.0.027
import numpy as np                         #v1.21.1
from random import randint                 #-
import pandas as pd                        #v1.3.1
import matplotlib.pyplot as plt            #v3.4.2
import forgi.visual.mplotlib as fvm        #v2.0.2
import forgi                               #v2.0.2
import RNA                                 #v0.7.5
import os                                  #-

#Import secondary commands
from Input import *
from helper import *

#Create Results folder
path = (PATH + "/ResultsToehold")
os.makedirs(path, exist_ok=True)


#Design Conditions
#Experiment conditions from  input file

#Design option specify any non-default job option from previous Wolfe17 paper
# (see more here: https://docs.nupack.org/design/#job-options)
#The Wobble option allows you to set whether the design will follow WATSON-CRICK base pair rules

my_options = DesignOptions(f_stop=0.01,
                           seed=0,
                           wobble_mutations= Wobble)

#Setting soft constrainst is an alternative to hard constraint
# which prohibits certain pattern violations
# (such as equal sequences), in which case the constraints are "softer".
# For more information: https://docs.nupack.org/design/#specify-soft-constraints

my_soft_const = [Pattern(Prevent)]

model1 = Model(material = Material,
               celsius= T,
               sodium= Sodium,
               magnesium= Magnesium)


# Create Toehold Switch Complex
#First define the domains that are part of the toehold structure
# and then generate the re-designed toehold determined by the sequence
# and structure in TargetComplex.

Structure_def = structure_define(miRNA, Loop, Linker, Reporter, paired, unpaired)

toe_strand = TargetStrand(Structure_def[2], name='Toehold')

C1 = TargetComplex([toe_strand], Structure_def[0], name="C1")

toe_complex = complex_design(complexes = [C1],
                             model=model1,
                             soft_constraints=my_soft_const,
                             options=my_options)

toehold_complex = toe_complex.run(trials= Trials)


#Stop codons identification and replacing in tube design trials
results_list = []

for trials in range(Trials):


    trigger = Domain(miRNA, name="trigger")
    Trigger = TargetStrand([trigger], name="Trigger")
    toeholdcomplex = toehold_complex[trials].to_analysis(C1)
    codon = rep_stop_codons(str(toeholdcomplex), Structure_def[3],StopCodons)

    if codon[1]:               #replace stop codons by "VNN"
        toeholdcomplex = codon[0]

    #Create TargetStrand to run test tube design
    # and optimizate the stop codon replaced

    toehold_dom = Domain(str(toeholdcomplex), name = "toehold_dom")
    toehold_strand = TargetStrand([toehold_dom], name = "toehold_strand")

    # Trials for test tube design
    C2 = TargetComplex([Trigger, toehold_strand], Structure_def[1], name = "C2")

    t_toe = TargetTube(on_targets={C2: 1e-06},
                       off_targets=SetSpec(max_size=2),
                       name='t_toe')

    my_tube_design = tube_design(tubes=[t_toe],
                                 soft_constraints=my_soft_const,
                                 model=model1, options=my_options)

    #if have replaced stop codon, run more 10 times
    tubetrials = 1

    if codon[1]:
        tubetrials = Trials_stop

    results= my_tube_design.run(trials=tubetrials)

    #Define file name and directory
    path_trial = (path + "/%s") % (trials)
    os.makedirs(path_trial, exist_ok=True)

    for tubetri in range(tubetrials):

        switchfinal = results[tubetri].to_analysis(toehold_strand)

        #Re-evaluates the presence of stop codons
        if rep_stop_codons(str(switchfinal), Structure_def[3], StopCodons)[1]:
            print("Tube Design Version Failed - Stop Codon detected")
            continue

        binding = results[tubetri].to_analysis(C2)
        triggerfinal = results[tubetri].to_analysis(Trigger)

        #Get MFE values for binding trigger-toehold, toehold, trigger and RBS-linker
        my_mfe_trigswitch = mfe(strands=binding, model=model1)
        my_mfe_switch = mfe(strands = switchfinal, model = model1)
        my_mfe_trigger = mfe(strands=triggerfinal, model=model1)
        my_mfe_rbslinker = mfe(strands=str(switchfinal)[paired+unpaired+8:], model=model1)

        #my_complex_ensemble_defect = defect(strands=[switchfinal], structure= Structure_def[0], model = model1)

        #Difference between Binding MFE and the sum of unbound trigger and toehold energies
        diff_mfe = my_mfe_trigswitch[0].energy - (my_mfe_trigger[0].energy + my_mfe_switch[0].energy)


        if not "."*(unpaired+3) in str(my_mfe_switch[0].structure)[:unpaired+3]:
            continue

        CodonOpt = "%d.%d" % (trials,tubetri)

        sampled_structures = sample(strands=[str(switchfinal)], num_sample=1, model=model1)


        ensemble_defect_binding = defect(strands=[miRNA,str(switchfinal)[:51]], structure='(21+.3)21.3(5.11)5.3', model=model1)
        ensemble_defect_complex = defect(strands= [str(switchfinal)], structure=str(my_mfe_switch[0].structure), model=model1)
        ensemble_defect_toeholddom = defect(strands=[str(switchfinal)[3:13]], structure = '.10', model=model1)
        ensemble_defect_trigger = defect(strands= [str(miRNA)], structure=".21", model=model1)
        print("complex")
        print(ensemble_defect_complex)
        print("binding")
        print(ensemble_defect_binding)
        print("single strand complex")
        print(ensemble_defect_toeholddom)
        print("trigger")
        print(ensemble_defect_trigger)

        score = 5*ensemble_defect_trigger + 4*ensemble_defect_toeholddom + 3*ensemble_defect_complex + 1*ensemble_defect_binding
        path_versions = (path_trial + "/%s") % (CodonOpt)

        os.makedirs(path_versions, exist_ok=True)

        probability_matrix_toehold = pairs(strands = str(switchfinal),model = model1)
        probability_matrix_binding= pairs(strands=str(binding), model=model1)

        save_file(CodonOpt, path_versions, str(binding), str(my_mfe_trigswitch[0].structure),
                                probability_matrix_binding, Structure_def[1], "Binding_Toehold_Trigger", model1)

        save_file(CodonOpt, path_versions, str(switchfinal), str(my_mfe_switch[0].structure),
                                probability_matrix_toehold, Structure_def[0], "Toehold", model1)
        #Append results
        appends = [CodonOpt, miRNA,str(Seq(miRNA).reverse_complement()), str(switchfinal),
               str(my_mfe_switch[0].structure), str(my_mfe_trigger[0].structure),
               str(my_mfe_trigswitch[0].structure),str(my_mfe_trigger[0].energy),
               str(my_mfe_switch[0].energy), str(my_mfe_trigswitch[0].energy),
               str(my_mfe_rbslinker[0].energy), str(diff_mfe), ensemble_defect_complex,
                   ensemble_defect_binding, score #, CodonOpt
               ]
        results_list.append(appends)

results_tab = pd.DataFrame(results_list, columns=["Trial", "miRNA.sequence",
                                                  "miRNA.Toehold.Complement", "Full.Toehold.seq",
                                                  "Toehold.structure", "Trigger.structure",
                                                  "Full.TriggerToehold.structure",
                                                  "Trigger.energy", "Toehold.energy",
                                                  "BindingToeholdTrigger.energy", "DeltaG_RBSLinker",
                                                  "Diff_energy", "ensemble_defect_complex", "ensemble_defect_binding", "Score"]).set_index("Trial")
print(results_tab)

results_tab.to_csv(r'test_toehold1.csv',index = True, header = True)



