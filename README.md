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

You can run **TIdeS** on a _de novo_ assembled transcriptome:
  
    python3 tides.py --fin myTranscriptome.fasta --taxon myTaxon --db proteinDB

Additionally, several genetic codes/translation tables are supported ([NCBI translation tables: 1,5,6,10,12,26,29,30](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).

    python3 tides.py --fin myTranscriptome.fasta --taxon myTaxon --db proteinDB --genetic-code 6

To see all options:
    
    python3 tides.py --help

### Planned Updates - 03-2023
- [x] Auto-optimize Random Forests (GridSearchCV)
- [x] Determine most useful set of composition criteria on broad phylogenomic scale (euks for now)
- [X] Automate TIdeS pipeline
- [X] tarball packaging
- [X] Rethinking folder(s) structure and outputs...
- [X] Support translation tables by NAME too, not just number (e.g. "Universal", "Ciliate", etc)
- [ ] Conda and PyPi packaging
- [ ] Prepare basic examples (orientation/contamination)
