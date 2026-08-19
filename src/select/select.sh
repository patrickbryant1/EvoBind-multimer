#Select from EvoBind
VHL_KRAS_MET_DIR=../../data/VHL_KRAS/metrics/
VHL_BRD4_MET_DIR=../../data/VHL_BRD4/metrics/
OUTDIR=../../data/

python3 select_evo.py --vhl_kras_metric_dir $VHL_KRAS_MET_DIR \
--vhl_brd4_metric_dir $VHL_BRD4_MET_DIR \
--outdir $OUTDIR
