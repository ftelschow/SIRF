#!/bin/bash
for parallel_count in {21..40}
do
for model in "ModelA" "ModelB"
do
for transform in "cohensd" "skewness" "skewnessN"
do
for obs in 0.05 0.1
do
	nohup Rscript --vanilla Sim_Script_obs.R $model 250 $parallel_count $x $obs "2022_01_20" "~/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"  $transform &
done
done
done
done
