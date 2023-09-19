# TIdeS

**T**ranscript **Ide**ntification and **S**election (TIdeS) is a method to identify putative open reading frames (pORFs) from a given transcriptome and is able to aid in the bulk decontamination of sequences from "messy" transcriptomic data. 

Overall, TIdeS couples sequence composition with ML approaches to discern pORFs in the correct reading frame with substantial improvement over other popular tools, while providing support for additional non-standard genetic codes. Additionally, TIdeS can be used to classify ORFs into several user-defined categories from highly contaminated datasets (e.g., parasite + host, kleptoplasts, big "dirty" protists).

## Dependencies
+ [Python 3.7+](https://www.python.org/downloads/)
  - [BioPython](https://biopython.org/wiki/Download)
  - [Pandas](https://pandas.pydata.org/)
  - [scikit-learn](https://scikit-learn.org/stable/)
  - [Optuna](https://optuna.org/#installation)
  - [seaborn](https://seaborn.pydata.org/installing.html)
+ [DIAMOND](https://github.com/bbuchfink/diamond)
+ [CD-HIT](https://github.com/weizhongli/cdhit)
+ [Barrnap](https://github.com/tseemann/barrnap)

## Installation
Note that TIdeS is only supported on UNIX systems (linux and MacOS).

Python's pip can be used to install the necessary python version and related packages.

```
pip install biopython pandas scikit-learn optuna seaborn
```

Followed by downloading the precompiled executables for the remaining dependencies.

Alternatively, you can do this through conda (note this will be updated):
```
# Create a new environment for TIdeS
conda create -n tides-ml
conda activate tides-ml

# Install the necessary packages
conda install -c bioconda -c conda-forge diamond cd-hit barrnap
conda install biopython pandas scikit-learn optuna seaborn
```

Clone the repository.
```
git clone https://github.com/xxmalcala/TIdeS.git
```

## Running TIdeS
The general syntax to run TIdeS is:

```
python3 tides.py --fin <transcriptome-assembly> --taxon <taxon-name> --db <protein-database>
```

Several example command lines and uses for TIdeS (i.e., ORF-calling and ORF classifying) are included in the examples folder. To run the examples, you need to be within the examples folder (e.g., `./orf_call_and_decontam.sh`)

### List of all options

|    Command                |  Comment  |
|---------------------------|-----------|
| `-h`, `--help`  | Print the help message |
| `-f`, `--fin <STRING>`  | FASTA formatted file. |
| `-n`,  `--taxon <STRING>`  | Name for your taxon, project, outputs. |
| `-t`, `--threads <INTEGER>`  | Number of available threads to use. Default value is `4`. |
| '-d`, `--db <STRING>`  | Path to FASTA or DIAMOND formatted proteome database. |
|`-p`, `--partials`  | Include partial ORFs for ORF calling. |
|`-id`, `--id <INTEGER>`  | Minimum % identity to remove redundant transcripts. Default value is `97`. |
|`-l`, `--min-orf <INTEGER>`  | Minimum transcript length (bp) for ORF calling. Default value is `300`. |
|`-ml`, `--max-orf <INTEGER>`  | Maximum transcript length (bp) for ORF calling. Default value is `10000`. |
|`-e`, `--evalue <REAL>`  | Maximum e-value to infer reference ORFs. Default value is `1e-30`. |
|`-gc`, `--gencode <STRING/INTEGER>`  | Genetic code to use to for ORF calling and translation. Default is `1`. |
|`-s`, `--strand <STRING>`  | Strands to call ORFs (both/minus/plus). Default value is `both`. |
|`-c`, `--contam <STRING>`  | Path to annotated sequence table. If unset, TIdeS will assume a prior model is provided as well. |
| `m`, `--model <STRING>`  | Path to a prior TIdeS run's model. These are the ".pkl" files. |
|`-k`, `--kmer <INTEGER>`   kmer length to use. Default value is `3`. |
|`-ov`, `--overlap`  | Permit overlapping kmers. |
|`--step <INTEGER>`  | Step-size for overlapping kmers. Default value is `kmer-length/2`. | 
|`--clean`  | Remove intermediate filter-step files. |
|`-gz`, `--gzip`  | Compress TIdeS outputs when finished. | 

### Prepare the reference protein database
Create a reference protein database for TIdeS (note you can use your own if you choose!).
This will generate a database from six diverse eukaryotes, representing a broad yet compact database for subsequent ORF-calling.

Note that this database will be prepared from whichever directory you call upon this script.

```
./TIdes/util/prep_tides_db.sh
```

## ORF-Calling and Assessment
**Inputs**
- FASTA formatted transcriptome assembly
- Taxon name (e.g., Homo sapiens, Op_me_Hsap)
- Protein database (can be prepared by "prep_tides_db.sh" in the **util** folder)

```
python3 tides.py -f examples/Durinskia_baltica.Example.Contam.fasta --taxon Durinskia_baltica -d tides_aa_db.dmnd
```

## Decontamination
**Inputs**
- FASTA formatted transcriptome assembly
- Taxon/project name (e.g., Durisnkia baltica, Dinotoms)
- Table of annotated sequence names (see examples folder)

```
python3 tides.py -f <transcriptome-assembly> -n <taxon-name> -c <annotated-seqs-table>
```
### Table of annotated sequences
The `<annotated-seqs-table>` should include sequence names and their category separated by tabs. Note that these sequences _should_ be present within the input FASTA file as well. Please aim to include at least 25 sequences for each category, although more (up to ~200) is great!

```
seq1  human
seq2  lunch
seq3  lunch
seq4  human
seq5  lunch
...
```

## Additional uses/approaches
More on how to run TIdeS and its uses can be found in the ```examples``` folder, including:
+ ORF-calling and sequence classifying with a previously trained model
+ Preparing a simple proteome database and ORF-calling
+ Naive approaches to inferring contamination
+ Example FASTA and `<annotated-seqs-table>` files
