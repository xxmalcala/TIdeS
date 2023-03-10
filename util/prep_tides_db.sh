#!/bin/bash

mkdir 'TIdeS_Prot_DB'

curl ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/955/GCF_000150955.2_ASM15095v2/GCF_000150955.2_ASM15095v2_protein.faa.gz --output TIdeS_Prot_DB/phaeodactylum.aa.fas.gz

curl ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz  --output TIdeS_Prot_DB/arabidopsis.aa.fas.gz

curl ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/985/GCF_000004985.1_V1.0/GCF_000004985.1_V1.0_protein.faa.gz --output TIdeS_Prot_DB/naegleria.aa.fas.gz

curl ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/203/815/GCF_000203815.1_DFas_2.0/GCF_000203815.1_DFas_2.0_protein.faa.gz --output TIdeS_Prot_DB/cavenderia.aa.fas.gz

curl ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/271/745/GCF_000271745.1_FO_FOSC_3_a_V1/GCF_000271745.1_FO_FOSC_3_a_V1_protein.faa.gz --output TIdeS_Prot_DB/fusarium.aa.fas.gz

curl ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/113/875/GCF_022113875.1_Hydra_105_v3/GCF_022113875.1_Hydra_105_v3_protein.faa.gz --output TIdeS_Prot_DB/hydra.aa.fas.gz

cat TIdeS_Prot_DB/*fas.gz > TIdeS_Prot_DB/tides.aa_db.fas.gz

rm TIdeS_Prot_DB/*aa.fas.gz

gunzip TIdeS_Prot_DB/tides.aa_db.fas.gz

cd-hit -c 0.8 -i TIdeS_Prot_DB/tides.aa_db.fas -o TIdeS_Prot_DB/tides_db.fas

rm TIdeS_Prot_DB/*aa_db.fas

diamond makedb --in TIdeS_Prot_DB/tides_db.fas -d tides_aa_db.dmnd

rm -rf TIdeS_Prot_DB
