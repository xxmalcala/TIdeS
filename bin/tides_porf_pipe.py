def classify_pORFs(args):
    from bin import (filter_trans as ft,
            dmnd_search as dms,
            random_orientation as rando,
            codon_counts as ccnt,
            call_orfs as co,
            rfc_classify as rfcc,
            translate_porfs as tp)

    co.eval_gcode_ttable(f'{args.genetic_code}')

    print(f'\nIntial Filtering of FASTA-FILE: {args.fin.split("/")[-1]}')
    filt_fas, taxon_dir = ft.filter_txpts(
                            args.fin,
                            args.taxon,
                            args.min_len,
                            args.threads)

    print(f'Generating set of reference ORFs.')
    ref_orfs = dms.extract_orf_hits(
                    filt_fas,
                    args.taxon,
                    args.db,
                    taxon_dir,
                    args.threads,
                    args.evalue)

    print(f'Generating mixed orientation data for: {ref_orfs.split("/")[-1]}')
    ref_rand_orfs = rando.gen_rand_orientation(ref_orfs)

    print('Calling "complete" ORFs.')
    query_orf_fas = co.call_orfs(
                        filt_fas,
                        taxon_dir,
                        args.genetic_code,
                        True,
                        False,
                        args.min_len)

    print('Generating codon counts for training data.')
    train_orf_tsv = ccnt.codon_counts_fasta(
                            ref_rand_orfs,
                            taxon_dir,
                            args.taxon,
                            True)

    query_orf_tsv = ccnt.codon_counts_fasta(
                            query_orf_fas,
                            taxon_dir,
                            args.taxon,
                            False)

    # reorganize inputs for consistency...
    rfc_fas = rfcc.classify_seqs(
                train_orf_tsv,
                query_orf_tsv,
                query_orf_fas,
                taxon_dir,
                args.taxon,
                args.reps)


    rfc_fin_fas, rfc_fin_aa = tp.prep_final_pORFs(
                                rfc_fas,
                                args.genetic_code,
                                taxon_dir)

    return taxon_dir
