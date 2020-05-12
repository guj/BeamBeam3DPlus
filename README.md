# BeamBeam3DPlus
Extending BeamBeam3D code with I/O options

PREREQUEST:  module load cmake

================== For using with ADIOS2 =============

PREREQUEST:  install adios2:

      git clone https://github.com/ornladios/ADIOS2
      
=> How to build with adios2:

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
  
without argument, if looks for file tunefoot.in. which looks like this:

	more tunefoot.in

	0                       #bunch <0 or 1>

	0 1000                  #particleStart & count

	0 4096                  #turnStart & count

	1 0 1 0 1 0             #x #px #y #py #z #pz


=> To plot with tunePlot.py 

PREREQUEST:

	- Have ffmeg in your file path (to save to movie) 

	- Point pythonpath to ADIOS python library, e.g.      

      export PYTHONPATH=/Users/junmin/software/adios2/master/install/lib/python3.7/site-packages/

      - Command to use the python script to visualize:

      python3 tunePlot.py -i tunefoot.bp

      "Note: tunefoot.bp is the filename hardcoded in the analysis code to store every tunefoot created in the file"

      
Workflow can be 

      Simulation --- SST/BPFile --- Analysis ---- SST/BPFile --- Python 

where SST or BPFile mode are specified in the adios2_config.xml 
* Simulation -- Analysis uses "SimulationOutput" I/O 
* Analysis -- Python uses "TuneFoot" I/O



            
            
  

:pensive:  ================== FOR using with SENSEI ===============. :pensive:

PREREQUEST: install sensei: 

      git clone https://gitlab.kitware.com/sensei/sensei
  
=> How to build:

      export SENSEI_DIR=/YOUR/SENSEI/INSTALL/lib/cmake/

      (NOTE: need to export ADIOS2_DIR=<YOUR/ADIOS2/INSTALL> if you used ADIOS2 in SENSEI)

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
      
      mpirun -n 4 ../../build-sensei/bin/sensei.xmain



  
    
