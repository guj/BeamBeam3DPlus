#!/bin/bash
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -q debug
#SBATCH -C knl
#SBATCH -J bp
#SBATCH -e %x%j.err 
#SBATCH -o %x%j.out


module swap craype-haswell craype-mic-knl
module swap PrgEnv-intel PrgEnv-gnu
module load cmake/3.14.4
module load emacs 

export CRAYPE_LINK_TYPE=dynamic

export LD_LIBRARY_PATH=/global/cscratch1/sd/junmin/software/adios2/master/install/lib64/:${LD_LIBRARY_PATH}
#lfs setstripe -c 256  -S 16M beam3d.bp


#load python3
module load python
module load ffmpeg

export PYTHONPATH=/global/cscratch1/sd/junmin/software/adios2/master/install/lib/python3.7/site-packages/


srun -n 32 -N 1 ../../b2/bin/adios2.xmain 
srun -n 32 -N 1 ../../b2/bin/adios2.mod_tunefoot 

# in file plot.in, it specified figs/ as the output directory. Remove old contents before using it again.
rm -rf figs/
mkdir figs

srun -n 1 -N 1 python3 ../../test/adios2/tunePlot.py -c plot.in 





