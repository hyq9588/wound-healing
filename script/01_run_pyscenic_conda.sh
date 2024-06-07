## inputs
f_loom_grn=output/00-2.mc_mat_for_step1.loom

## outputs
grn_output=output/01-step1_adj.tsv
ctx_output=output/01-step2_reg.tsv

## reference
f_tfs=cisTarget_db/hsa_hgnc_tfs.motifs-v10.txt
f_motif_path=cisTarget_db/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
f_db_500bp=cisTarget_db/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
f_db_10kb=cisTarget_db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


#### 1. build GRN
## ~60 min
time arboreto_with_multiprocessing.py \
$f_loom_grn \
$f_tfs \
--method grnboost2 \
--output $grn_output \
--num_workers 10 \
--seed 777


#### 2. cisTarget
## ~33 min
time pyscenic ctx \
$grn_output \
$f_db_500bp $f_db_10kb \
--annotations_fname $f_motif_path \
--expression_mtx_fname $f_loom_grn \
--output $ctx_output \
--num_workers 10
