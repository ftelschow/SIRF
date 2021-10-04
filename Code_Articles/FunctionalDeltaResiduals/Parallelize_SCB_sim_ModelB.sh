#!/bin/bash
for parallel_count in {1..20}
do
  nohup Rscript --vanilla Sim_Script.R "ModelB" 500 $parallel_count "2021_05_20" "~/Rpackages/SCBfun/DeltaResiduals/" &
done
