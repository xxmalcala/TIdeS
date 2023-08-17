# TIdeS

**T**ranscript **Ide**ntification and **S**election (TIdeS) is a method to identify putative open reading frames (pORFs) from a given transcriptome and is able to aid in the bulk decontamination of sequences from "messy" transcriptomic data. Overall, TIdeS couples sequence composition with ML approaches to discern pORFs in the correct reading frame and to identify target sequences from contaminated datasets (e.g. kleptoplastic organisms). 

## Dependencies
+ [Python 3.7+](https://www.python.org/downloads/)
  - [BioPython](https://biopython.org/wiki/Download)
  - [Pandas](https://pandas.pydata.org/)
  - [scikit-learn](https://scikit-learn.org/stable/)
  - [Optuna](https://optuna.org/#installation)
+ [DIAMOND](https://github.com/bbuchfink/diamond)
+ [CD-HIT](https://github.com/weizhongli/cdhit)
+ [Barrnap](https://github.com/tseemann/barrnap)
+ Patience!

## Quick Start

### Installation

Python's pip can be used to install the necessary python version and related packages.

```
pip install biopython pandas scikit-learn optuna
```

Followed by downloading the precompiled executables for the remaining dependencies.

Alternatively, you can do this through conda (note this will be updated):
```
# Create a new environment for TIdeS
conda create -n tides-ml
conda activate tides-ml

# Install the necessary packages
conda install -c bioconda -c conda-forge diamond cd-hit barrnap
conda install biopython pandas scikit-learn optuna
```

Clone the repository.
```
git clone https://github.com/xxmalcala/TIdeS.git
```

### Prepare the reference protein database
Create a reference protein database for TIdeS (note you can use your own if you choose!).
This will generate a database from six diverse eukaryotes, representing a broad yet compact database for subsequent ORF-calling.

Note that this database will be prepared from whichever directory you call upon this script. This may be changed in subsequent updates.

```
TIdes/util/prep_tides_db.sh
```

### ORF-Calling and Assessment
**Inputs**
- FASTA formatted transcriptome assembly
- Taxon name (e.g., Homo sapiens, Op_me_Hsap)
- Protein database (can be prepared by "prep_tides_db.sh" in the **util** folder)

```
python3 tides.py --fin <transcriptome-assembly> --taxon <taxon-name> --db <protein-database>
```
##### Required Arguments
```
-f, --fin           Input file in FASTA format
-n, --taxon         Taxon-name
-d, --db            Protien database (FASTA or DIAMOND format)
```
##### Optional Arguments
```
-t, --threads       Number of CPU threads (default = 1)
-m, --model         Previously trained TIdeS model (".pkl" file)
-k, --kmer          kmer size for generating sequence features (default = 3)
-ov, --overlap      Permit overlapping kmers (see --kmer)
--step              Step-size for overlapping kmers (default is kmer-length/2)
-q, --quiet         No console output
-gz, --gzip         Tar and gzip TIdeS output
```
#### Example

```cmd
python3 tides.py -f example/TBD -n TBD -d TBD
```

### Decontamination
**Note that these options are currently being streamlined/changed...**
**Inputs**
- FASTA formatted transcriptome assembly
- Taxon name (e.g., Myrionecta rubrum, Sr_ci_Mrub)
- Table of sequence names annotated as 'target' or 'non-target'
```
python3 tides.py --fin <transcriptome-assembly> --taxon <taxon-name> --db <protein-database>
```
Kmer-size and overlap can have dramatic impact on the inference of target/non-target sequences. For contaminated taxa with highly similar composition (e.g. Monocystis agilis [parasite] and Helobdella robusta [host]), recommended to include ```-k 4 -ov``` in the command.

##### Required Arguments
```
-f, --fin           Input file in FASTA format
-n, --taxon         Taxon-name or PhyloToL taxon-code
-c, --contam        Table of sequences annotated as 'target' or 'non-target'
```
#### Example

```cmd
python3 tides.py -f example/TBD -n TBD -d TBD -c
```

### Planned Updates - 08-2023
- [ ] Conda and/or PyPi packaging
- [ ] Prepare basic examples (pORF-calls, contamination, EDS for contamination)
