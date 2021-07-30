

#in terminal
#python3 -m pip install -U nupack -f ~/Downloads/nupack-VERSION/package
from Bio.Seq import Seq
from nupack import *
from test import *
from helper import *
import numpy as np
from random import randint
import pandas as pd

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
test = Domain("N21", name = "test")
c = Domain("N3", name = "c")
d = Domain("W5", name = "d")
loop = Domain(Loop , name = "loop")
g = Domain("W5" , name = "g")
h = Domain("AUG" , name = "h")
i = Domain("N11", name= "i")
linker = Domain(Linker, name= "linker")


#Toehold Switch Complex - Define a estrutura do toehold

toe_strand = TargetStrand([free_tail,~complement,c,d,loop,g,h,i,linker], name='Toehold')

C1 = TargetComplex([toe_strand], Toehold_structure, name="C1")
toe_complex = complex_design(complexes = [C1],
                             model=model1,
                             soft_constraints=my_soft_const,
                             options=my_options)

toehold_complex = toe_complex.run(trials= Trials)
'''''''''
    # Test for stop codons
stop_codons = ['UGA', 'UAA', 'UAG']
test_region = i + linker

    # Split the test region into codons to check one at a time
    test_region_codons = [test_region[i:i + 3] for i in range(0, len(test_region), 3)]

    # Loop through the list and check if there are any stop codons in this reading frame
    for codon in test_region_codons:
        # assert not codon in stop_codons, 'The generated switch contains a stop codon'
        if codon in stop_codons:
            print("Stop")
            break
    print(s)
    s = mut_seq
'''

print(complement.reverse_complement())
for trials in range(Trials):

    Trigger = TargetStrand([trigger], name="Trigger")
    toeholdcomplex = toehold_complex[trials].to_analysis(C1)
    print(toeholdcomplex)
    toehold_dom = Domain(str(toeholdcomplex),
                     name = "toehold_dom")
    toehold_strand = TargetStrand([toehold_dom], name = "toehold_strand")

    #Toehold-Trigger duplex
    trig_complex = TargetComplex([Trigger], 'U21', name="Trig_complex")

    Trig_Toe_struct = structure_define(miRNA, Loop, Linker)
    print(Trig_Toe_struct)
    C2 = TargetComplex([Trigger, toehold_strand],Trig_Toe_struct, name = "C2")
    print(C2)

    t_toe = TargetTube(on_targets={C2: 1e-06},
                       off_targets=SetSpec(max_size=2),
                       name='t_toe')

    my_tube_design = tube_design(tubes=[t_toe],
                                 soft_constraints=my_soft_const,
                                 model=model1, options=my_options)

    results= my_tube_design.run(trials=1)


    binding = results[0].to_analysis(C2)
    fullswitch_seq = results[0].to_analysis(toehold_strand)


    #MFE
    my_mfe_trigswitch = mfe(strands=binding, model=model1)
    my_mfe_switch = mfe(strands=toeholdcomplex, model=model1)
    my_mfe_trigger = mfe(strands=miRNA, model=model1)
    my_mfe_rbslinker = mfe(strands=str(toeholdcomplex)[29:],
                           model=model1)
    #Append results
    results_tab["miRNA.sequence"].append(miRNA)
    results_tab["miRNA.Toehold.Complement"].append(str(Seq(mir_seq).reverse_complement()))
    results_tab["Full.switch.seq"].append(str(fullswitch_seq))
    results_tab["Switch.structure"].append(str(my_mfe_switch[0].structure))
    results_tab["Trigger.structure"].append(str(my_mfe_trigger[0].structure))
    results_tab["Full.TriggerSwitch.structure"].append(str(my_mfe_trigswitch[0].structure))
    #result_tab["Trigger.mfe"]


    #print('Partition function =', my_pfunc)
    print(my_mfe_trigswitch[0].energy)
    print(my_mfe_rbslinker[0].energy)
    print(my_mfe_switch[0].energy)
    print(my_mfe_trigger[0].energy)
    print(my_mfe_trigswitch[0].energy - (my_mfe_trigger[0].energy + my_mfe_switch[0].energy))
    print(my_mfe_trigswitch[0].structure)
    print(str(Seq(miRNA).reverse_complement())[:-5])

print(results_tab)
