# TIdeS

**T**ranscript **Ide**ntification and **S**election

## Motivation

Many projects simply use the "longest" complete ORF (using tools such as TransDecoder) prior to phylogenomic studies. This has many potential implications on the results, especially when identifying lineage-specific gene families (LSGFs) from transcriptomic data, given the partial nature of most RNA-seq approaches.

**TIdeS** uses a Random Forest Classifiers to select the most "likely" open reading frame(s) from a transcript. After aligning transcripts against a protein database, TIdeS accounts for the composition of alignable ORFs to infer **complete** putative ORFs (pORFs) from the transcriptome.

In practice, this may also extend towards extraction of transcripts from a targeted taxon in a contaminated transcriptome... Development for this is TBD.

**TIdeS** will (hopefully) help mitigate these issues, focusing on extracting the transcripts and ORFs most like a small robustly curated samples.

## Dependencies
+ [Python 3.6+](https://www.python.org/downloads/)
  - [BioPython](https://biopython.org/wiki/Download)
  - [Pandas](https://pandas.pydata.org/)
  - [scikit-learn](https://scikit-learn.org/stable/)
+ [DIAMOND](https://github.com/bbuchfink/diamond)
+ [CD-HIT](https://github.com/weizhongli/cdhit)
+ [Barrnap](https://github.com/tseemann/barrnap)
+ Patience!

## Quick Start

### ORF-Calling and Assessment

**Inputs**
- FASTA formatted transcriptome assembly
- Taxon name (e.g., Homo sapiens, Op_me_Hsap)
- Protein database (can be prepared by "prep_tides_db.sh" in the **util** folder)

```
python3 tides.py --fin <transcriptome-assembly> --taxon <taxon-name> --db <protein-database>
```

#### Arguments

##### Required

```
-f, --fin           Input file in FASTA format
-n, --taxon         Taxon-name or PhyloToL taxon-code
-d, --db            Protien database (FASTA or DIAMOND format)
```

##### Optional
-p, --threads       Number of CPU threads (default = 1)
-m, --model         Previously trained TIdeS model (".pkl" file)
-k, --kmer          kmer size for generating sequence features (default = 3)
-q, --quiet         No console output
-ov, --overlap      Permit overlapping kmers (see --kmer)
-gz, --gzip         Tar and gzip TIdeS output
```

##### Optional


To see all options:
    
    python3 tides.py --help

### Planned Updates - 05-2023
- [ ] Compare XGBoost to Sci-kit Learn
- [ ] Incorporate HyperOpt
- [ ] Conda and PyPi packaging
- [ ] Prepare basic examples (orientation/contamination)

