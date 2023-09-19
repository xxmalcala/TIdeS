# Note this script can only be run from the examples folder given its reliance on relative paths 
tar -zxvf Durinskia_baltica.Example.Contam.fasta.tar.gz

# Make a small but appropriate proteome database
../util/prep_tides_db.sh

# Make the initial ORF calls
python3 ../tides.py -d tides_aa_db.dmnd -t 1 -f Durinskia_baltica.Example.Contam.fasta -n Durinskia_baltica_ORF_Call -gz

# Generate the ORF-Composition plots -- quick way to assess obvious contamination
python3 ../util/orf_composition.py Durinskia_baltica_ORF_Call_TIdeS/Durinskia_baltica_ORF_Call.TIdeS.fasta

# Grab a subsample from each of the apparent clusters -- look at the ".png" in the ORF-composition folder
# Left-most cluster
python3 ../util/seqs_by_composition.py -f Durinskia_baltica_ORF_Call_TIdeS/Durinskia_baltica_ORF_Call.TIdeS.fasta \
-t ORF_Composition/Durinskia_baltica_ORF_Call.TIdeS.CompSummary.tsv \
-o Cluster1 --min_gc12 40 --max_gc12 55 --min_gc3 40 --max_gc3 55

# Right-most cluster
python3 ../util/seqs_by_composition.py -f Durinskia_baltica_ORF_Call_TIdeS/Durinskia_baltica_ORF_Call.TIdeS.fasta \
-t ORF_Composition/Durinskia_baltica_ORF_Call.TIdeS.CompSummary.tsv \
-o Cluster2 --min_gc12 40 --max_gc12 60 --min_gc3 70 --max_gc3 90

# Merge the cluster files together
cat Durinskia_baltica_ORF_Call_TIdeS/*.Cluster*txt > Durinskia_baltica.Cluster_Labels.txt

# Train TIdeS to group the sequences into separate clusters
python3 ../tides.py -t 1 -f Durinskia_baltica_ORF_Call_TIdeS/Durinskia_baltica_ORF_Call.TIdeS.fasta \
-n Durinskia_baltica_Decontam -c Durinskia_baltica.Cluster_Labels.txt -gz

