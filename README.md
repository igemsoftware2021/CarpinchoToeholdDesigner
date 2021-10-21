 <h1 align =  "center">iGEM 2021 project: Carpincho Toehold Designer</h1>

---
<h4 align="center"> 
	🚧  To be improved... (See limitations) 🚧
</h4>
This repository contains the code for Toehold Design based on the 
NUPACK package[1]. The program was built as part of the BioPank 
project for the iGEM 2021 Contest. The iGEM UNILA_LatAm 
implemented this code to design adapted Toeholds for 
_Bacillus subtilis_ in the detection of miRNAs acquired from 
the extracellular environment.

The project was uploaded on GitHub because of iGEM’s rule



## Toehold Template

Carpincho Tool segue estritamente o modelo Type II  descrito por 
Green et al.(2020) [2] e baseia seu funcinamento em NUPACK Design and Analysis algorithms.


![toehold_git.png](https://github.com/chagas98/CarpinchoToeholdDesigner/blob/master/toehold_git.png)







## Requirements

---

- NUPACK 4.0.0.27 (http://nupack.org/downloads)
- Biopython v1.79
- NumPy v1.21.1
- Pandas v1.3.1
- Matplotlib v3.4.2
- Forgi v2.0.2
- RNA v0.7.5

## Getting Started

---

- In the root directory:

mk dir build

cd build

cmake ../

install the libraries versions (pip3 install (library name)... == version )


 - Dowload Carpincho Toehold Designer (Git link)
and go to the root directory and type

python3 Carpincho.py

## Usage

---

You can fill the inputs for your Toehold Design preferences in th Input.py file. The variables are:

|    Variable    |    Description    |
| -------------- | ----------------- |
| Unpaired       |  Seed sequence for trigger binding that is not connected to hairpin|
| Paired         |  Lower stem nucleotides
| Prevent        | Sequence patterns to  prevent |
| StopCodons     | Stop codons to Prevent between RBS and Linker |
| Loop           | RBS sequence (with or without spacers) |
| Trials         | Number of iterations for Toehold generation |
| Trials_stop    | Number of iterations if stop codons are identified |
| Material       | Temperature-dependent RNA and DNA free energy parameter sets ([NUPACK User Guide](https://docs.nupack.org/model/)|
| T              | Temperature specified (default: celsius = 37) ([NUPACK User Guide](https://docs.nupack.org/model/)|
| Sodium         | Salt Conditions default = 1M ([NUPACK User Guide](https://docs.nupack.org/model/)|
| Magnesium      | Salt Conditions default = 1M ([NUPACK User Guide](https://docs.nupack.org/model/)|
| Wobble         | True- do not follow WATSON-CRICK base pair rules (G-U, for instance) |
| miRNA (trigger)| Sequence Trigger for miRNA target |
| Linker         | Linker sequence |
| Reporter       | Reporter (firstest 5 codons) |
| PATH           | Directory output path |


## RESULTS

---
After installation and running, Carpincho Tool saves each version in a folder that contains subversions if there are stop-codons (Trials_stop). 
These sub-versions are stored in *version.sub-version* folders.

###Output Files
 - **Binding_Toehold_Trigger.fx**: Trigger+Toehold sequence and Structure notation;
 - **Toehold_version.fx**: Toehold sequence and Structure notation;
 - **Toehold_version_struct.png**: Toehold conformation from ViennaRNA package[3];
 - **Binding_Toehold_Trigger_probab.png**: Pair binding probabilities
 - **Binding_Toehold_Trigger_probab_Theory.png**: Pair bindgin probabilities for toehold+trigger binding structure in theory (standard)
 - **Toehold_version_probab.png**: Pair probabilities for Toehold Structure
 - **Toehold_version_probab_Theory.png**: Pair probabilities for Toehold Structure in theory
 - **.csv file**: All thermodynamic parameters, ensemble defects, sequences and structures information

### Results files
ResultsToehold_Bs_RBS.zip shows an example of output files. This Toehold were designed for B. subtilis to detect miR-283


### Limitations
- Limited by the standardized structure at unpaired =11 and paired = 10
- Slow for large libraries
- The code does not contain the best toeholds score



[1] J. N. Zadeh, C. D. Steenberg, J. S. Bois, B. R. Wolfe, M. B. Pierce, A. R. Khan, R. M. Dirks, N. A. Pierce. NUPACK: analysis and design of nucleic acid systems. J Comput Chem, 32:170–173, 2011. 

[2] Green, A.A. et al. (2020). Complex cellular logic computation using ribocomputing devices. Nature, 548, p. 117-121. 
