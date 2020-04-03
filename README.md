# BeamBeam3DPlus
Extending BeamBeam3D code with I/O options

PREREQUEST:  module load cmake

==== For ADIOS2 ===
PREREQUEST:  install adios2:
      git clone https://github.com/ornladios/ADIOS2
      
      
=> How to build:

export ADIOS2_DIR=<YOUR/ADIOS2/INSTALL/>

mkdir build-adios2

cd build-adios2

cmake ../source/adios2

make 

=> To Test:   

cd ../test/adios2 
 
It contains a sample input files:  beam1.in beam2.in + adios2-config.xml

- Simulation:

      mpirun -n 4 ../../build-adios2/bin/adios2.xmain 

- Analysis:

      ../../build-adios2/bin/adios2.mod_tunefoot
  
  
===== FOR SENSEI ===

PREREQUEST: install sensei: 

      git clone https://gitlab.kitware.com/sensei/sensei
  
=> How to build:

export SENSEI_DIR=/YOUR/SENSEI/INSTALL/lib/cmake/

(NOTE: may need to export ADIOS2_DIR=<YOUR/ADIOS2/INSTALL> if you used ADIOS2 in SENSEI)

(NOTE: same goes for HDF5)

mkdir build-sensei

cd build-sensei

cmake ../source/sensei

make

=> To test

cd ../test/sensei

- A simple fortran test: (needs sensei.xml to be present) 

   ./build-sensei/bin/sensei.simpleTest 

- Simulation: (needs beam1.in beam2.in)

      mpirun -n 4 ../../build-sensei/bin/sensei.xmain
  
- Analysis:

      ../../build-sensei/bin/sensei.tune_foot



==== About sensei.xml ====

      - SUPOORTED transport type: hdf5 or adios2

      - samples: see sensei-h5.xml or sensei-a2.xml

      - users can customize the filename and reader options (bunch, particles range, etc)
