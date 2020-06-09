#!/bin/bash
anno=$1
prefix=$2
export PYTHONPATH=/zfssz4/B2C_RD_P2/PMO/fangzhonghai/software/pjg/:$PYTHONPATH
/zfssz4/B2C_RD_P2/PMO/fangzhonghai/software/anaconda3/envs/python37/bin/python autopvs1_in_bgianno.py -i $anno -o $prefix -p 1
