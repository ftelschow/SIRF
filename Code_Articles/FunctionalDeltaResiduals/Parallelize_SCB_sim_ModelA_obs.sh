#!/bin/bash
for parallel_count in {1..20}
do
for x in 50 100 175
do
for obs in 0.05 0.1
do
for transform in "cohensd" "skewness" "skewness (normality)"
do
	nohup Rscript --vanilla Sim_Script_obs.R "ModelA" 500 $parallel_count $x $obs "2022_01_20" "~/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/" $transform &
do
done
done
done
