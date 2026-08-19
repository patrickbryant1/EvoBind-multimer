#!/bin/bash -x

#############INPUT PARAMETERS#############
BASE=. #Change this depending on your local path
DATADIR=$BASE/data # Input files path
###IDS###
ID1=VHL   # E3 ligase protein
ID2=KRAS  # protein target to degrade
###Fasta sequences###
FASTA1=$DATADIR/$ID1'_'$ID2/$ID1'.fasta'
FASTA2=$DATADIR/$ID1'_'$ID2/$ID2'.fasta'
###Receptor interface residues - to be provided with --target_*_if_residues=$IFRES* \ if using
IFRES1='' # 0 index of first residue, comma separated list of residues
IFRES2='' # 0 index of first residue, comma separated list of residues
#Get length of the peptide molecular glue to be designed
PEPTIDELENGTH=8
NITER=1000 # edit

#########Step1: Create MSAs #########
MSA1=$DATADIR/$ID1'_'$ID2/$ID1'.a3m'
MSA2=$DATADIR/$ID1'_'$ID2/$ID2'.a3m'
#Write MSAs for the two protein sequences
HHBLITSDB=$BASE/data/uniclust30_2018_08/uniclust30_2018_08
if test -f $MSA1; then
	echo $MSA1 exists    
else
	$BASE/hh-suite/build/bin/hhblits -i $ID1 -d $HHBLITSDB -E 0.001 -all -n 2 -oa3m $MSA1
fi
if test -f $MSA2; then
	echo $MSA2 exists
else
	$BASE/hh-suite/build/bin/hhblits -i $ID2 -d $HHBLITSDB -E 0.001 -all -n 2 -oa3m $MSA2
fi
###Pair and block MSAs for the two input proteins###
if test -f $DATADIR/$ID1'_'$ID2/$ID1'_'$ID2'_blocked.a3m'; then
        echo $DATADIR/$ID1'_'$ID2/$ID1'_'$ID2'_blocked.a3m' exists
else
        echo Making MSAs for $ID1 and $ID2
        #Pair MSAs
        python3 $BASE/scr/pair_msas.py --a3m1 $MSA1 --a3m2 $MSA2 --max_gap_fraction 0.9 --outname $DATADIR/$ID1'_'$ID2/$ID1'_'$ID2'_paired.a3m'
        #Block diagonalize MSAs
        python3 $BASE/scr/block_msas.py --a3m1 $MSA1 --a3m2 $MSA2 --max_gap_fraction 0.9 --outname $DATADIR/$ID1'_'$ID2/$ID1'_'$ID2'_blocked.a3m'
fi
PAIRED_MSA=$(ls $DATADIR/$ID1'_'$ID2/$ID1'_'$ID2'_paired.a3m')
BLOCKED_MSA=$(ls $DATADIR/$ID1'_'$ID2/$ID1'_'$ID2'_blocked.a3m')
MSAS=${PAIRED_MSA},${BLOCKED_MSA} # Comma separated list of msa paths

#########Step3: Design molecular glue#########
##### AF2 CONFIGURATION ####
OUTDIR=$DATADIR/$ID1'_'$ID2/designs/$PEPTIDELENGTH
PARAM=$BASE'/src/AF2/'
PRESET='full_dbs' #Choose preset model configuration - no ensembling (full_dbs) and (reduced_dbs) or 8 model ensemblings (casp14).
MAX_RECYCLES=8 #max_recycles (default=3)
MODEL_NAME='model_1' #model_1_ptm

mkdir -p $OUTDIR

#Run the molecular glue design
python3 $BASE/src/mc_design.py \
        --target_1_fasta=$FASTA1 \
        --target_2_fasta=$FASTA2 \
        --target_1_if_residues=$IFRES1 \
        --target_2_if_residues=$IFRES2 \
        --msas=$MSAS \
        --peptide_length=$PEPTIDELENGTH \
        --model_names=$MODEL_NAME \
        --output_dir=$OUTDIR \
        --data_dir=$PARAM \
        --max_recycles=$MAX_RECYCLES \
        --num_iterations=$NITER \
        --predict_only=False \
        --cyclic_offset=1 # remove this if not using cyclic offset
