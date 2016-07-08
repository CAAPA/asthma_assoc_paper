#!/bin/bash

#Set parameters
if [ "$#" -eq  "0" ]
then
    echo "Usage: ${0##*/} <study_name> <cohort_name> <out_file_name> <err_p_file_name>"
    echo "Eg:"
    echo "bash run_qq.sh BRIDGE jhu_abr"
    exit
fi

study_name=$1
cohort_name=$2
in_file_name=/gpfs/barnes_share/caapa_metal/data/input/${cohort_name}.txt
out_file_name=../data/output/qq_plots/${cohort_name}.png
err_p_file_name=../data/output/err_p_vals/${cohort_name}.txt

#Create output directories in case they do not exist
mkdir ./data/
mkdir ./data/output/
mkdir ./data/output/qq_plots
mkdir ./data/output/err_p_vals

#Load R
module load R/3.2.5

#Run QQ plot
cat draw_qq.R | R --vanilla --args $study_name $in_file_name $out_file_name $err_p_file_name
