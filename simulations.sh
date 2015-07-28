#!/bin/bash 

echo "SIMULATOR V0.1"

# create a directory to store our simulations
cd ~
mkdir ~/simulation_out
echo "simulation_out created"

DEFAULT_CORES=2
DEFAULT_TIMES=1
SIMULATION_TIMES="${1:-$DEFAULT_TIMES}"
# SIMULATION_CORES="${2:-$DEFAULT_CORES}"
SIMULATION_CORES="$DEFAULT_CORES"
SIMULATIONS_NAME="${@:2}"
echo "CORES $SIMULATION_CORES"
echo "TIMES $SIMULATION_TIMES"


# run the simulations
# cd CC3D_*
# echo "running python dispatch_simulations"
# #this file dispatches simulations to the core
# python dispatch_simulations.py $SIMULATIONS_NAME
# echo "completed dispatch_simulations"
# cd ~

echo "running python dispatch_simulations"
cd ~/summer15/metastasis
python dispatch_simulations.py $SIMULATION_TIMES $SIMULATION_CORES $SIMULATIONS_NAME > ~/simulation_out/sim_log.txt
echo "completed dispatch_simulations"
cd ~

# move all simulations to the metastasis folder
mv ~/simulation_out/ ~/summer15/metastasis/
echo 'all simulations were moved'

cd ~/summer15/metastasis
echo "now in metastasis directory"

SIMULATIONS_OUT=$(ls ./simulation_out)

for var in $SIMULATIONS_OUT; do
	./pipeline init ./simulation_out/$var 
done
echo 'initializations done'

python dispatch_analysis.py $SIMULATIONS_OUT > ./simulation_out/analysis_log.txt
echo 'dispatching completed'

echo "SIMULATIONS AND ANALYSIS COMPLETED $(date)" > x.txt

# git add x.txt
# git commit -m 'auto command #39'
# git push origin master

# git add simulation_out/*/pipe_out_*/*
# git commit -m 'results of simulation $(date) #39'
# git push origin master