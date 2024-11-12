***HOW TO COMPILE AND RUN***
This project is separated into four modes: mode1, mode2, mode3, and mode2_skewed, each with its own compilation and execution commands. 
Please compile each mode first before running any commands. Example commands are provided for each mode at the end.

g++ -std=c++11 -O3 -o mode1 mode1.cpp
./mode1 <cutoff_value>
./mode1 17500

g++ -std=c++11 -pthread -O3 -o mode2 mode2.cpp 
./mode2 <cutoff_value> <thread_#>
./mode2 17500 8

mpic++ -std=c++11 -pthread -O3 -o mode3 mode3.cpp
mpirun -np <process_#> ./mode3 <cutoff_value> <thread_#>
mpirun -np 2 ./mode3 17500 8

g++ -std=c++11 -pthread -O3 -o mode2_skewed mode2_skewed.cpp
./mode2_skewed <cutoff_value> <thread_#>
./mode2_skewed 26000 8



***DEPENDENCIES***
To compile and run this project, you need the following dependencies installed:

1. GNU Compiler Collection (GCC) - Specifically `g++` for C++ compilation. 
   - Install via: `sudo apt-get install g++` (on Debian/Ubuntu) or `brew install gcc` (on macOS with Homebrew).

2. MPI Library - Required for multi-process parallelism in `mode3`.
   - Install via: `sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev` (on Debian/Ubuntu) or `brew install open-mpi` (on macOS with Homebrew).



***HOW TO REPRODUCE ALL TEST CASES***
All experimental data are provided as screenshots in the folder named 'experiment_data'. 
You can follow the commands shown in each screenshot to reproduce the corresponding results.



***HOW TO REPRODUCE ALL FIGURES IN THE REPORT***
All scripts used to generate the figures in the report are saved in the 'plotting scripts' folder. 
These scripts are written in Python. To generate a figure, use the command 'python3 plotting_scripits/<filename>', replacing 
`<filename>` with the name of the specific script.
