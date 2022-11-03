def classify_contam(args):
    from bin import (codon_counts as ccnt,
                    parse_sister as ps,
                    rfc_classify as rfcc)
    import shutil

    # Parse tree-walking summary and mark known contam and "clean" for training.
    print('\nGathering marked contaminant and "clean" sequences for training TIdeS.')
    mlen_fas, ref_cntm_fas, taxon_dir = ps.bin_seqs(
                                args.fin,
                                args.taxon,
                                args.sister_table,
                                args.min_len)

    print('Generating codon counts for training data.')
    train_orf_tsv = ccnt.codon_counts_fasta(
                            ref_cntm_fas,
                            taxon_dir,
                            args.taxon,
                            True, False)

    query_orf_tsv = ccnt.codon_counts_fasta(
                            mlen_fas,
                            taxon_dir,
                            args.taxon,
                            False)

    # Use Random Forest Classifier to predict contaminant and "clean" genes.
    rfc_fas = rfcc.classify_seqs(
                train_orf_tsv,
                query_orf_tsv,
                mlen_fas,
                taxon_dir,
                args.taxon,
                1,
                False)

    shutil.copy2(rfc_fas, taxon_dir)
    shutil.copy2(rfc_fas.replace("fas","Contam.fas"), taxon_dir)
    return taxon_dir
