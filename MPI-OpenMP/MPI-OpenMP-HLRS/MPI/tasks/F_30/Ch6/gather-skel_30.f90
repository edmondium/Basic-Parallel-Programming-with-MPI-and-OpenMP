PROGRAM first_example

!==============================================================!
!                                                              !
! This file has been written as a sample solution to an        !
! exercise in a course given at the High Performance           !
! Computing Centre Stuttgart (HLRS).                           !
! The examples are based on the examples in the MPI course of  !
! the Edinburgh Parallel Computing Centre (EPCC).              !
! It is made freely available with the understanding that      !
! every copy of this file must include this header and that    !
! HLRS and EPCC take no responsibility for the use of the      !
! enclosed teaching material.                                  !
!                                                              !
! Authors: Joel Malard, Alan Simpson,            (EPCC)        !
!          Rolf Rabenseifner, Traugott Streicher (HLRS)        !
!                                                              !
! Contact: rabenseifner@hlrs.de                                !
!                                                              !
! Purpose: Gathering data from all processes                   !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE

  INTEGER :: n               ! application-related data
  DOUBLE PRECISION :: result ! application-related data
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: result_array
  INTEGER :: my_rank, num_procs, rank  ! MPI-related data

  CALL MPI_Init()

  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  CALL MPI_Comm_size(MPI_COMM_WORLD, num_procs)

  !  doing some application work in each process, e.g.:
  result = 100.0 + 1.0 * my_rank
  WRITE(*,'(A,I3,A,I3,A,I2,A,I5,A,F9.2)') &
   &        'I am process ', my_rank, ' out of ', num_procs, &
   &        ' handling the ', my_rank, 'th part of n=', n, ' elements, result=', result

  IF (my_rank == 0) THEN  
    ALLOCATE(result_array(0:num_procs-1))
  ENDIF

  ! ----------- the following gathering of the results should --------------
  ! -----------   be substituted by one call to MPI_Gather    --------------
  IF (my_rank /= 0) THEN  ! in all processes, except "root" process 0
    !  sending some results from all processes (except 0) to process 0:
    CALL MPI_Send(result, 1, MPI_DOUBLE_PRECISION, 0, 99, MPI_COMM_WORLD)
  ELSE  ! only in "root" process 0
    result_array(0) = result  ! process 0's own result
    !  receiving all these messages
    DO rank=1, num_procs-1  ! result of processes 1, 2, ...
      CALL MPI_Recv(result_array(rank), 1, MPI_DOUBLE_PRECISION, rank, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
    END DO
  ENDIF
  ! ------------------ end of the gathering algorithm ----------------------

  IF (my_rank == 0) THEN  
    DO rank=0, num_procs-1
      WRITE(*,'(A,I3,A,F9.2)') &
       &      'I''m proc 0: result of process ', rank, ' is ', result_array(rank) 
    END DO
  ENDIF

  CALL MPI_Finalize()

END PROGRAM
