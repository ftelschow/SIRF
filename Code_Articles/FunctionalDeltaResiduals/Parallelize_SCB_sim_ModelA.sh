#!/bin/bash
for parallel_count in {1..2}
do
for model in "ModelA" "ModelB" "ModelC"
do
for transform in "skewness" "skewness (normality)" "kurtosis" "kurtosis (normality)"
do
	nohup Rscript --vanilla Sim_Script.R $model 500 $transform $parallel_count "2022_01_11" "~/test/" &
done
done
done
