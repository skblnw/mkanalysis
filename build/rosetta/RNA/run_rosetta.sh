stepwise.default.linuxgccrelease -in:file:fasta step1_motifs.fasta -s combine_rebuild_4142.pdb -out:file:silent 4142_rebuild.out -motif_mode -score:weights stepwise/rna/rna_res_level_energy7beta.wts -cycles 20 -terminal_res 1 78