
REC_FEATS=/proj/berzelius-2023-267/users/x_patbr/results/molecular_glues/VHL_KRAS/afm_scoring/VHL_KRAS/features.pkl
PEP_SEQS=../../data/VHL_KRAS/metrics.csv
PARAMDIR=/proj/berzelius-2023-267/users/x_patbr/software/AFProfile/data/
MODEL_PRESET='multimer'
NUM_RECYCLES=20 #Number of recycles
CYCLIC=1
OUTDIR=/proj/berzelius-2023-267/users/x_patbr/results/molecular_glues/VHL_KRAS/afm_scoring/VHL_KRAS/

#Run
python3 ./run_AFM.py --receptor_features $REC_FEATS \
--data_dir $PARAMDIR --model_preset $MODEL_PRESET \
--peptide_sequences $PEP_SEQS \
--num_recycles $NUM_RECYCLES \
--cyclic_offset $CYCLIC \
--outdir $OUTDIR
