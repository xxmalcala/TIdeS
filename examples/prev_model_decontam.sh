# Decontaminate a file from a prior TIdeS run!

# Need the model ".pkl" file from a prior run and a FASTA file of existing ORFs
tar -zxvf Durinskia_baltica.ORFs.fasta.tar.gz
tar -zxvf Durinskia_baltica_Decontam.TIdeS.pkl.tar.gz

# Classify ORFs from a previously trained model
python3 ../tides.py -t 1 -c -m Durinskia_baltica_Decontam.TIdeS.pkl \
-f Durinskia_baltica.ORFs.fasta \
-n Durinskia_baltica_Classify_ORFs

