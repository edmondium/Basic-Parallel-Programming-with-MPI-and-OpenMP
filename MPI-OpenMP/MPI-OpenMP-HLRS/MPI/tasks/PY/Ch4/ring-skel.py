#!/usr/bin/env python3

################################################################# 
#                                                               # 
#  This file has been written as a sample solution to an        # 
#  exercise in a course given at the High Performance           # 
#  Computing Centre Stuttgart (HLRS).                           # 
#  The examples are based on the examples in the MPI course of  # 
#  the Edinburgh Parallel Computing Centre (EPCC).              # 
#  It is made freely available with the understanding that      # 
#  every copy of this file must include this header and that    # 
#  HLRS and EPCC take no responsibility for the use of the      # 
#  enclosed teaching material.                                  # 
#                                                               # 
#  Authors: Joel Malard, Alan Simpson,            (EPCC)        # 
#           Rolf Rabenseifner, Traugott Streicher,              #
#           Tobias Haas (HLRS)                                  # 
#                                                               # 
#  Contact: rabenseifner@hlrs.de                                # 
#                                                               # 
#  Purpose: A program to try MPI_Issend and MPI_Recv.           # 
#                                                               # 
#  Contents: Python code, buffer send version (comm.Send)       # 
#                                                               # 
#################################################################

from mpi4py import MPI
import numpy as np

rcv_buf = np.empty((),dtype=np.intc)
status = MPI.Status()

comm_world = MPI.COMM_WORLD
my_rank = comm_world.Get_rank()
size = comm_world.Get_size()
right = (my_rank+1)      % size
left  = (my_rank-1+size) % size

sum = 0
snd_buf = np.array(my_rank, dtype=np.intc)

for i in range(size):
   comm_world.Send((snd_buf, 1, MPI.INT), dest=right, tag=17)
   # WRONG program, because if MPI_Send is implemented with a
   # synchronous communication protocol then this program will deadlock!
   # Task of this exercise: substitute MPI_Send through a nonblocking call
   #                        as descibed in the slides
   comm_world.Recv((rcv_buf, 1, MPI.INT), source=left,  tag=17, status=status)
   
   np.copyto(snd_buf, rcv_buf) # We make a copy here. What happens if we write snd_buf = rcv_buf instead?
   sum += rcv_buf

print(f"PE{my_rank}:\tSum = {sum}")
