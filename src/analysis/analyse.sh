VHL_KRAS_MET_DIR=../../data/VHL_KRAS/metrics/
VHL_BRD4_MET_DIR=../../data/VHL_BRD4/metrics/
VHL_KRAS_EXP=../../data/VHL_KRAS/degradation_results.csv
OUTDIR=../../data/plots/

python3 ./calcs.py --vhl_kras_metric_dir $VHL_KRAS_MET_DIR \
--vhl_brd4_metric_dir $VHL_BRD4_MET_DIR \
--vhl_kras_experimental_results $VHL_KRAS_EXP \
--outdir $OUTDIR
