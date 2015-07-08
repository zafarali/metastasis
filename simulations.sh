#!/bin/bash 

echo "SIMULATOR V0.1"

# create a directory to store our simulations
cd ~
mkdir ~/simulation_out
echo "simulation_out created"

SIMULATION_TIMES="$1"
SIMULATIONS_NAME="${@:2}"

# run the simulations
# cd CC3D_*
# echo "running python dispatch_simulations"
# #this file dispatches simulations to the core
# python dispatch_simulations.py $SIMULATIONS_NAME
# echo "completed dispatch_simulations"
# cd ~

echo "running python dispatch_simulations"
cd ~/summer15/metastasis
python dispatch_simulations.py $SIMULATION_TIMES $SIMULATIONS_NAME
echo "completed dispatch_simulations"
cd ~

# move all simulations to the metastasis folder
mv ~/simulation_out/ ~/summer15/metastasis/
echo 'all simulations were moved'

cd ~/summer15/metastasis
echo "now in metastasis directory"

SIMULATIONS_OUT=$(ls ./simulation_out)

for var in $SIMULATIONS_OUT; do
	python pipeline.py init ./simulation_out/$var
done
echo 'initializations done'

python dispatch_analysis.py $SIMULATIONS_OUT
echo 'dispatching now.'
	