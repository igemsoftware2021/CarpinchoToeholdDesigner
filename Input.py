
# Input information for the Toehold design


paired = 11                #Number of paired nucleotides (see readme for more information)
unpaired = 10              #Number of unpaired nucleotides



#Pattern Prevention
Prevent = ['A4', 'C4', 'G4', 'U4', 'M6', 'K6', 'W6', 'S6', 'R6', 'Y6']   #Avoiding Sequence patterns
StopCodons = ["UAG", "UGA", "UAA"] #Avoiding Stop Codons


#Spacers and RBS
Loop = "UAAUAAGGAGG" #ref RBS for bacillus subtilis https://microbialcellfactories.biomedcentral.com/track/pdf/10.1186/s12934-020-01404-2.pdf

#Trials
Trials = 5                #Trials to generate Toehold
Trials_stop = 2            #Trials to generate Toehold with tube design after stopcodon replacement

#Conditions/Model
Material = "rna95"
T = 37
Sodium = 1.0
Magnesium = 0.0
Wobble = True #True- do not follow WATSON-CRICK base pair rules

#miR to detect
miRNA = "UAAAUAUCAGCUGGUAAUUCU"

#Linker and Reporter
Linker = "AACCUGGCGGCAGCGCAAAAG"
Reporter = "AUGAACAUCAAAAAGUUUGCA" #SacB

#Output Path
PATH = "/home/samuel/Downloads"
