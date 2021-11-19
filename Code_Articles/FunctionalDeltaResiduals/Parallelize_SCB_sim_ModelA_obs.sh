#!/bin/bash
for parallel_count in {1..20}
do
for x in 50 100 175
do
for obs in 0.05 0.1
do
	nohup Rscript --vanilla Sim_Script2.R "ModelA" 500 $parallel_count $x $obs "2021_11_19" "~/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/" &
done
done
done
