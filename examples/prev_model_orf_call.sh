# Decontaminate a file from a prior TIdeS run!

# Need the model ".pkl" file from a prior run and a FASTA file of existing ORFs
tar -zxvf Durinskia_baltica.Example.Contam.fasta.tar.gz
tar -zxvf Durinskia_baltica_ORF_Call.TIdeS.pkl.tar.gz

# Classify ORFs from a previously trained model
tides -t 1 -m Durinskia_baltica_ORF_Call.TIdeS.pkl \
-i Durinskia_baltica.Example.Contam.fasta \
-o Durinskia_baltica_Call_ORFs

