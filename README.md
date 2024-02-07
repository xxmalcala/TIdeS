# TIdeS

**T**ranscript **Ide**ntification and **S**election (TIdeS) is a machine learning approach to discern pORFs in the correct reading frame with substantial improvement over other popular tools, while providing support for additional non-standard genetic codes. Additionally, TIdeS can be used to classify ORFs into several user-defined categories from highly contaminated datasets (e.g., parasite + host, kleptoplasts, big "dirty" protists) or broadly into "eukaryotic" _versus_ "non-eukaryotic" using the metagenomic classifier Kraken2.

# Installation
Note that TIdeS is only supported on UNIX systems (linux and MacOS).

## Install with [mamba](https://mamba.readthedocs.io/en/latest/index.html) (recommended)
```
mamba create -n tides-ml
mamba activate tides-ml
mamba install -c bioconda tides-ml
```

## Install with pip
```
pip install tides-ml
```

Afterwards, ensure that the following dependencies are installed and in your path:

**Dependencies**
+ [DIAMOND 2.0.13+](https://github.com/bbuchfink/diamond)
+ [CD-HIT 4.8.1+](https://github.com/weizhongli/cdhit)
+ [Barrnap 0.9](https://github.com/tseemann/barrnap)
+ [Kraken2 2.1.0+](https://github.com/DerrickWood/kraken2).


# ORF Prediction
### Prepare a reference protein database

Feel free to use your own if you choose!
Alternatively, we provide a bash script to create a database from six diverse eukaryotes, representing a broad yet compact database.

```
./TIdeS/util/prep_tides_db.sh
```
 
### ORF Prediction Inputs
- FASTA formatted transcriptome assembly
- Taxon name (e.g., Homo sapiens, Op_me_Hsap)
- Protein database (can be prepared by "prep_tides_db.sh" in the **util** folder)
- **Note**: examplar commands can be found in the [orf_call_and_decontam.sh](https://github.com/xxmalcala/TIdeS/blob/main/examples/orf_call_and_decontam.sh) script found in the examples folder.

```
tides -i <transcriptome-assembly> -o <taxon-name> -d <protein-database>
```

**TIdeS** does support several alternative genetic codes (i.e., reassigned stop-to-sense codons)
For example, using 'ciliate' genetic code (translation table 6; [NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi):

```
tides -i <transcriptome-assembly> -o <taxon-name> -d <protein-database> -g 6
```

# ORF Classification and Decontamination
**Inputs**
- FASTA formatted transcriptome assembly
- Taxon/project name (e.g., Durisnkia baltica, Dinotoms)
- Table of annotated sequence names (see examples folder) OR path to a formatted [Kraken2 database](https://benlangmead.github.io/aws-indexes/k2)
- Python scripts for generating composition plots (`orf_composition.py`) and selection of sequences based on composition metrics (`seqs_by_composition.py`) can be found in the **util** folder
- **Note**: examplar commands can be found in the `orf_call_and_decontam.sh` script found in the examples folder.

Using user-defined table of annotated sequences:
```
tides -i <Predicted-ORFs> -o <taxon-name> -c <annotated-seqs-table>
```

Using Kraken2 to identify non-eukaryotic sequences:
```
tides -i <Predicted-ORFs> -o <taxon-name> -c -k <kraken2-database>
```

Please note that there must be at least 25 annotated sequences for **each** class (this includes automatic classification with Kraken2).

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

# Deploy a previously trained TIdeS model
**Inputs**
- FASTA formatted transcriptome assembly
- Taxon/project name (e.g., Durisnkia baltica, Dinotoms)
- Path to a trained TIdeS model ('.pkl' file)
- **Note**: examplar commands for using a [prior ORF-prediction model](https://github.com/xxmalcala/TIdeS/blob/main/examples/prev_model_orf_call.sh) and [prior ORF classification model](https://github.com/xxmalcala/TIdeS/blob/main/examples/prev_model_decontam.sh) are provided
```
tides -i <transcriptome-assembly> -o <taxon-name> -m <TIdeS.pkl-file>
```

# List of all options

|    Command                |  Comment  |
|---------------------------|-----------|
| `-h`, `--help`  | Print the help message |
| `-i`, `--fin <STRING>`  | FASTA formatted file. |
| `-o`, `--taxon <STRING>`  | Name for your taxon, project, outputs. |
| `-t`, `--threads <INTEGER>`  | Number of available threads to use. Default value is `4`. |
| `-d`, `--db <STRING>`  | Path to FASTA or DIAMOND formatted proteome database. |
| `-p`, `--partials`  | Include partial ORFs for ORF calling. |
| `-id`, `--id <INTEGER>`  | Minimum % identity to remove redundant transcripts. Default value is `97`. |
| `-l`, `--min-orf <INTEGER>`  | Minimum transcript length (bp) for ORF calling. Default value is `300`. |
| `-ml`, `--max-orf <INTEGER>`  | Maximum transcript length (bp) for ORF calling. Default value is `10000`. |
| `-e`, `--evalue <REAL>`  | Maximum e-value to infer reference ORFs. Default value is `1e-30`. |
| `--memory <INTEGER>`  | memory limit (MB) for CD-HIT. Default value is `2000`, unlimited is `0`. |
| `-g`, `--gencode <STRING/INTEGER>`  | Genetic code to use to for ORF calling and translation. Default is `1`. |
| `-s`, `--strand <STRING>`  | Strands to call ORFs (both/minus/plus). Default value is `both`. |
| `-c`, `--contam <STRING>`  | Path to annotated sequence table. If unset, TIdeS will assume a prior model is provided as well. |
| `-k`, `--kraken <STRING>`  | kraken2 database to identify and filter non-eukaryotic sequences. |
| `--no-filter` | Skip the rRNA and transcript clustering steps. |
| `m`, `--model <STRING>`  | Path to a prior TIdeS run's model. These are the ".pkl" files. |
| `--kmer <INTEGER>`  | kmer length to use. Default value is `3`. |
| `--overlap`  | Permit overlapping kmers. |
| `--step <INTEGER>`  | Step-size for overlapping kmers. Default value is `kmer-length/2`. | 
| `--clean`  | Remove intermediate filter-step files. |
| `-gz`, `--gzip`  | Compress TIdeS outputs when finished. | 

## Additional uses/approaches
More on how to run TIdeS and its uses can be found in the ```examples``` folder, including:
+ ORF Prediction
+ ORF Classification
+ ORF prediction and classification with a previously trained model
+ Preparing a simple proteome database and ORF-calling
+ Naive approaches to inferring contamination
+ Example FASTA and `<annotated-seqs-table>` files
