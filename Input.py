

structure_type = "c"     # a_series | b_series

paired = 11
unpaired = 10
#structure
#Toehold_structure = 'U15 D11(U3 D5 U11 U3) U6 D3 U4 U5'
#Trigger_Toehold_structure = "(21+.5)21.3(5.11)5.35"

#Pattern Prevention
Prevent = ['A4', 'C4', 'G4', 'U4', 'M6', 'K6', 'W6', 'S6', 'R6', 'Y6']

StopCodons = ["UAG", "UGA", "UAA"]
#Loop sequence
Loop = "UAAUAAGGAGG" #ref https://microbialcellfactories.biomedcentral.com/track/pdf/10.1186/s12934-020-01404-2.pdf

#Trials
Trials = 10

#Conditions/Model
Material = "rna95"
T = 37
Sodium = 1.0
Magnesium = 0.0
Wobble = True #True/False - do not follow WATSON-CRICK base pair rules

#miR
miRNA = "UAAAUAUCAGCUGGUAAUUCU"

#Linker and Reporter
Linker = "AACCUGGCGGCAGCGCAAAAG"
Reporter = "AUGAACAUCAAAAAGUUUGCA" #SacB

#Output Path
PATH = "/home/samuel/Downloads"
