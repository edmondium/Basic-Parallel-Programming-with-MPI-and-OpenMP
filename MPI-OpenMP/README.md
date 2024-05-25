# Basic Parallel Programming with MPI and OpenMP (2 days)

In this course we will be using the [Snellius supercomputer](https://servicedesk.surf.nl/wiki/display/WIKI/Snellius) hosted at SURF. You will need to connect to it through the ssh protocol (natively installed on Linux and Mac).

For Windows users, we recommend to connect via ssh using [MobaXterm](https://mobaxterm.mobatek.net/) or, on Windows 10, using the [native bash environment](https://msdn.microsoft.com/en-us/commandline/wsl/install_guide)

## Before you connect to Cartesius

Please check your email inbox for a message with your access credentials and follow the instructions to connect to the user portal, so that you can accept the usage agreement and change your password.

## Connect to Snellius

To connect to Snellius, please open a terminal and use the following command:

    ssh scurXXX@snellius.surf.nl

In case you cannot connect directly because of IP whitelisting, please use the alternative connection via our doornode. The same user name and password as for Snellius applies, and you need to write it twice after selecting the connection to Snellius in the drop-down menu. The command is as follows:

    ssh scurXXX@doornode.surfsara.nl

For both commands above "scurXXX" is the login name that you received via email. Type in the password you received together with your login name and press Enter (note that the cursor may not move while introducing the password, but the typing will succeed).

This access to Snellius is made available for duration of the course and at most one week after its finalization.

## OpenMP hands-on exercises
    
In order to execute the OpenMP tests, please load the following GCCcore module.

    module load 2023
    module load foss/2023a
    
You may compile OpenMP codes like this:

    gcc -fopenmp openmp_code.c          ### C compilation
    gfortran -fopenmp openmp_code.f     ### Fortran compilation

Any OpenMP codes that are low in computational demands may be executed interactively in the terminal where the interactive session is opened. However, if you wish to check performance data or make sure that the codes are executed without disturbance from other users, you MUST use the template script provided.

After compilation, OpenMP codes are eventually executed as follows:

    ./a.out
    
If you wish to indicate the number of threads for the execution using the environment variable OMP_NUM_THREADS, you may run the tests including this variable before the executable.

    OMP_NUM_THREADS=4 ./a.out


## MPI hands-on exercises

For the use of MPI, please purge the default environment and load the Intel compiler module as follows:

    module load 2023
    module load foss/2023a
    
You may compile MPI codes like this:

    mpicc mpi_code.c        ### C compilation
    mpif90 mpi_code.f       ### Fortran compilation
    
All MPI codes must be executed using the wrapper script provided. The correct way to execute MPI codes is as follows (e.g. for the use of 4 processes using an executable with name "a.out").

    srun -n 4 ./a.out


## Using the wrapper script

Please use the scripts mpi_script_template.sh or openmp_script_template in order to run your MPI or OpenMP tasks respectively. You should modify this template conveniently and save it in a file in your home directory of Snellius (e.g. script.sh)

In order to run your code, you should submit your script in the terminal like this:

    sbatch script.sh
    
Use this command to check that your job is running:

    squeue –u $(whoami)
   
After executing squeue you will see the JobID of the script ("NNNNNNN"), and if you do:

    ls slurm*

you may also see a new file called "slurm-NNNNNNN.out". This file contains the output of your code, so that you can check if the execution was successful.
