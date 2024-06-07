## inputs
f_loom_grn=00-2.mc_mat_for_step1.loom

## outputs
grn_output=01-step1_adj.tsv
ctx_output=01-step2_reg.tsv

## 创建文件夹
data_path=`pwd`/../output

## reference
db_path=`pwd`/../cisTarget_db
f_tfs=hsa_hgnc_tfs.motifs-v10.txt
f_motif_path=motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
f_db_500bp=mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
f_db_10kb=mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


#### 1. build GRN
## ~17.5 mins
time singularity run --bind $data_path:/data,$db_path:/cisdb \
../aertslab-pyscenic-0.9.18.sif \
	pyscenic grn \
		--seed 777 \
		--num_workers 10 \
		--method grnboost2 \
        -o /data/$grn_output \
        /data/$f_loom_grn \
        /cisdb/$f_tfs

#### 2. cisTarget
## ~30 mins 内存占用小
time singularity run --bind $data_path:/data,$db_path:/cisdb \
../aertslab-pyscenic-0.12.1.sif \
	pyscenic ctx \
		/data/$grn_output \
		/cisdb/$f_db_500bp /cisdb/$f_db_10kb \
		--annotations_fname /cisdb/$f_motif_path \
		--expression_mtx_fname /data/$f_loom_grn \
		--output /data/$ctx_output \
		--num_workers 10

