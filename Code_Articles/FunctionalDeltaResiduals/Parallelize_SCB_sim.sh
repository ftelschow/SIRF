#!/bin/bash
for parallel_count in {1..20}
do
for model in "ModelA" "ModelB" "ModelC"
do
for transform in "skewness" "skewnessN" "kurtosis" "kurtosisN"
do
	nohup Rscript --vanilla Sim_Script.R $model 500 $transform $parallel_count "2022_01_11" "~/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/" &
done
done
done
