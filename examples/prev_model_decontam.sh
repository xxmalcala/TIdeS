# Decontaminate a file from a prior TIdeS run!

# Need the model ".pkl" file from a prior run and a FASTA file of existing ORFs
tar -zxvf Durinskia_baltica.ORFs.fasta.tar.gz
tar -zxvf Durinskia_baltica_Decontam.TIdeS.pkl.tar.gz

# Classify ORFs from a previously trained model
tides -t 1 -c -m Durinskia_baltica_Decontam.TIdeS.pkl \
-i Durinskia_baltica.ORFs.fasta \
-o Durinskia_baltica_Classify_ORFs

