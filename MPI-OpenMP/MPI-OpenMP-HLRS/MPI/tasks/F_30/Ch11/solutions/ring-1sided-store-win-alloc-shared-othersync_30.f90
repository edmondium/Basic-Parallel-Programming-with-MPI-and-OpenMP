PROGRAM ring

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
! Purpose: A program to try out one-sided communication        !
!          with window=rcv_buf and MPI_PUT to put              !
!          local snd_buf value into remote window (rcv_buf).   !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
  IMPLICIT NONE

  INTEGER :: my_rank_world, size_world
  INTEGER :: my_rank_sm,    size_sm
  TYPE(MPI_Comm) :: comm_sm
  INTEGER :: right, left
  INTEGER :: i, sum
  INTEGER :: snd_buf  ! no longer ASYNCHRONOUS, because no MPI_Put(snd_buf, ...)
  INTEGER, POINTER, ASYNCHRONOUS :: rcv_buf(:)  ! "(:)" because it is an array
  TYPE(C_PTR) :: ptr_rcv_buf
  TYPE(MPI_Win) :: win 
  INTEGER :: disp_unit
  INTEGER(KIND=MPI_ADDRESS_KIND) :: integer_size, lb
  INTEGER(KIND=MPI_ADDRESS_KIND) :: rcv_buf_size, target_disp

  TYPE(MPI_Status)  :: status
  TYPE(MPI_Request) :: request
  INTEGER :: snd_dummy, rcv_dummy

  CALL MPI_Init()
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank_world)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size_world)

  CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_sm)
  CALL MPI_Comm_rank(comm_sm, my_rank_sm) 
  CALL MPI_Comm_size(comm_sm, size_sm) 
  IF (my_rank_sm == 0) THEN
    IF (size_sm == size_world) THEN
      write (*,*) 'MPI_COMM_WORLD consists of only one shared memory region'
    ELSE
      write (*,*) 'MPI_COMM_WORLD is split into 2 or more shared memory islands'
    END IF
  END IF

  right = mod(my_rank_sm+1,         size_sm)
  left  = mod(my_rank_sm-1+size_sm, size_sm)

! ALLOCATE THE WINDOW.
  CALL MPI_Type_get_extent(MPI_INTEGER, lb, integer_size)
  rcv_buf_size = 1 * integer_size
  disp_unit = integer_size
  CALL MPI_Win_allocate_shared(rcv_buf_size, disp_unit, MPI_INFO_NULL, comm_sm, ptr_rcv_buf, win)
  CALL C_F_POINTER(ptr_rcv_buf, rcv_buf, (/1/)) ! if rcv_buf is an array
  rcv_buf(0:) => rcv_buf ! change lower bound to 0 (instead of default 1) ! if rcv_buf is an array
! CALL C_F_POINTER(ptr_rcv_buf, rcv_buf) ! if rcv_buf is not an array

  CALL MPI_Win_lock_all(MPI_MODE_NOCHECK, win)

  sum = 0
  snd_buf = my_rank_sm

  DO i = 1, size_sm

!    ... The compiler may move the read access to rcv_buf
!        in the previous loop iteration after the following 
!        1-sided MPI calls, because the compiler has no chance
!        to see, that rcv_buf will be modified by the following 
!        1-sided MPI calls.  Therefore a dummy routine must be 
!        called with rcv_buf as argument:
 
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
 
!    ... Now, the compiler expects that rcv_buf was modified,
!        because the compiler cannot see that MPI_F_SYNC_REG
!        did nothing. Therefore the compiler cannot move any
!        access to rcv_buf across this dummy call. 
 
     CALL MPI_Win_sync(win)
 
!      ... tag=17: posting to left that rcv_buf is exposed to left, i.e.,
!                  the left process is now allowed to store data into the local rcv_buf
     CALL MPI_Irecv(rcv_dummy, 0, MPI_INTEGER, right, 17, comm_sm, request)
     CALL MPI_Send (snd_dummy, 0, MPI_INTEGER, left,  17, comm_sm)
     CALL MPI_Wait(request, status)

     CALL MPI_Win_sync(win)

 
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
!    ... that the "rcv_buf(0+(right-my_rank_sm)) = snd_buf" statement
!        cannot be moved across then "CALL MPI_Win_fence" above  

!    target_disp = 0
!    CALL MPI_Put(snd_buf, 1, MPI_INTEGER, right, target_disp, 1, MPI_INTEGER, win)
     rcv_buf(0+(right-my_rank_sm)) = snd_buf
 
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
!    ... that the "rcv_buf(0+(right-my_rank_sm)) = snd_buf" statement
!        cannot be moved across then "CALL MPI_Win_fence" below

!    ... The following pair of syncs is needed that the write-read-rule is fulfilled.
     CALL MPI_Win_sync(win)
 
!    ... The following communication synchronizes the processors in the way
!        that the origin processor has finished the store
!        before the target processor starts to load the data.
!    ... tag=18: posting to right that rcv_buf was stored from left
     CALL MPI_Irecv(rcv_dummy, 0, MPI_INTEGER, left,  18, comm_sm, request)
     CALL MPI_Send (snd_dummy, 0, MPI_INTEGER, right, 18, comm_sm)
     CALL MPI_Wait(request, status)

     CALL MPI_Win_sync(win)
 
!    ... The compiler has no chance to see, that rcv_buf was
!        modified. Therefore a dummy routine must be called
!        with rcv_buf as argument:
 
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
 
!    ... Now, the compiler expects that rcv_buf was modified,
!        because the compiler cannot see that MPI_F_SYNC_REG
!        did nothing. Therefore the compiler will use the new
!        value on the memory, instead of some old value in a
!        register.

!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
!    ... no longer needed, because the access to snd_buf
!        no longer a nonblocking MPI call.
!        Now, it is a directly executed expression.   

     snd_buf = rcv_buf(0)
     sum = sum + rcv_buf(0)

  END DO

  WRITE(*,*) 'World:',  my_rank_world,' of ',size_world, &
   &         'comm_sm:',my_rank_sm,   ' of ',size_sm,   '; Sum =', sum

  CALL MPI_Win_unlock_all(win)

  CALL MPI_Win_free(win)

  CALL MPI_Finalize()

END PROGRAM
