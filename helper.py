import pandas as pd
import numpy as np


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


# creating the Results File
results_tab = {"miRNA.sequence":[],
               "miRNA.Toehold.Complement": [],
               "Full.switch.seq": [],
               "Switch.structure": [],
               "Trigger.structure": [],
               "Full.TriggerSwitch.structure": []
               }