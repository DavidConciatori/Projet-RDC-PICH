#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --account=def-conciato
#SBATCH --mem=4G
module load  intel/2018.3 openmpi/3.1.2 openfoam/7


blockMesh
scalarTransportFoam