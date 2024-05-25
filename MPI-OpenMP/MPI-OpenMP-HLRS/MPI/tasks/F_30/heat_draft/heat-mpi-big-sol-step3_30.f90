program heat
  use mpi_f08
  implicit none
! naming: i = 1st array coordinate = 1st process coordinate = y
!         k = 2nd array coordinate = 2nd process coordinate = x
!                                                    Caution: y,x and not x,y
  integer i, k, it

  logical, parameter :: prt=.true.,  prt_halo=.true.
! logical, parameter :: prt=.false., prt_halo=.false.
  integer, parameter :: imax=15, kmax=12
  integer, parameter :: istart=0, kstart=0, b1=1
  integer, parameter :: itmax=20000
  double precision, parameter :: eps=1.d-08
 
  double precision, allocatable :: phi(:,:), phin(:,:)
  double precision dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax

! Preparation for parallelization with domain decomposition:

!   integer is=istart, ie=imax, ks=kstart, ke=kmax ! now defined and calculated below

! The algortithm below is already formulated with lower and upper 
! index values (is:ie), (ks:ke) instead of (istart:imax), (kstart:kmax).
! With the MPI domain decomposition, (is:ie) and (ks:ke) will define the index rang 
! of the subdomain in a given MPI process.

! additional variables for parallelization
  integer numprocs, my_rank, right, left, upper, lower, ierror
  TYPE(MPI_Comm) :: comm  ! Cartesian communicator
  integer dims(1:2), coords(1:2) ! helper variables only for MPI_CART_CREATE and ..._COORDS
  logical period(1:2)            !                  only for MPI_CART_CREATE and ..._COORDS
  integer idim, kdim, icoord, kcoord
  integer isize, iinner1, in1,  ksize, kinner1, kn1 ! only for calculation of iinner, ...
  integer iinner, iouter, is, ie, kinner, kouter, ks, ke 
  integer ic, i_first, i_last,  kc, k_first, k_last ! only for printing
  integer, parameter :: stride=10 ! for reduced calculation of the abort criterion
  TYPE(MPI_Request) :: req(4)
  TYPE(MPI_Status) :: statuses(4)
  TYPE(MPI_Datatype) :: horizontal_border, vertical_border ! datatype
  double precision dphimaxpartial ! for local/global calculation of the abort criterion
  double precision start_time, end_time, comm_time, criterion_time

!  naming: originally: i=istart..imax,      now: i=is..ie,  icoord=0..idim-1
!          originally: k=kstart..kmax,      now: k=ks..ke,  kcoord=0..kdim-1
!                      with istart=kstart=0      s=start index,  e=end index

  call MPI_Init()
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs)

! the 2-dimensional domain decomposition:

  dims = (/0,0/)
  period = (/.false.,.false./)
  call MPI_Dims_create(numprocs, 2, dims)
  idim = dims(1)
  kdim = dims(2)
  call MPI_Cart_create(MPI_COMM_WORLD,2,dims,period,.true.,comm)
  call MPI_Comm_rank(comm, my_rank)
  call MPI_Cart_coords(comm, my_rank, 2, coords) 
  icoord = coords(1) 
  kcoord = coords(2) 
  call MPI_Cart_shift(comm, 0, 1, left, right)
   ! the ranks (left,right) represent the coordinates ((icoord-1,kcoord), (icoord+1,kcoord)) 
  call MPI_Cart_shift(comm, 1, 1, lower,upper)
   ! the ranks (lower,upper) represent the coordinates ((icoord,kcoord-1),(icoord,kcoord+1)) 

! whole y indecees   |------------- isize = 1+imax-istart ------------|
!    start/end index  ^-istart                                  imax-^
! 
! 1. interval        |--- iouter1---|   
!                    |--|--------|--|
!                     b1  iinner1 b1 
!    start/end index  ^-is      ie-^ 
! 2. interval                 |--|--------|--|   
! 3. interval                          |--|--------|--| 
! 4. interval                                   |--|-------|--| 
!                                                   iinner0 = inner1 - 1
! 5. interval = idim's interval                         |--|-------|--|
!
! In each iteration on each interval, the inner area is computed
! by using the values of the last iteration in the whole outer area. 
!
! icoord = number of the interval - 1
! 
! To fit exactly into isize, we use in1 intervals of with iinner1
! and (idim-in1) intervals of with (iinner1 - 1) 
!
! isize <= 2*b1 + idim * iinner1  
!
! ==>   iinner1 >= (isize - 2*b1) / idim
! ==>   smallest valid integer value:


! computing is:ie, ks:ke  ---  details see heat-mpi1-big.f 
  isize = 1 + imax - istart
  iinner1 = (isize - 2*b1 - 1) / idim + 1 
  in1 = isize - 2*b1 - idim * (iinner1 - 1) 
  if (icoord .lt. in1) then
    iinner = iinner1
    is = istart + icoord * iinner 
  else
    iinner = iinner1 - 1
    is = istart + in1 * iinner1 + (icoord-in1) * iinner 
  endif 
  iouter = iinner + 2*b1 
  ie = is + iouter - 1
! same for x coordinate: 
  ksize = 1 + kmax - kstart
  kinner1 = (ksize - 2*b1 - 1) / kdim + 1 
  kn1 = ksize - 2*b1 - kdim * (kinner1 - 1) 
  if (kcoord .lt. kn1) then
    kinner = kinner1
    ks = kstart + kcoord * kinner 
  else
    kinner = kinner1 - 1
    ks = kstart + kn1 * kinner1 + (kcoord-kn1) * kinner 
  endif 
  kouter = kinner + 2*b1 
  ke = ks + kouter - 1

  if (my_rank.eq.0) then
    write (*,'(/,2(1x,a,i3),2(1x,a,i3,1x,a,i3),2(1x,a,i3))') &
     &           'isize=',isize,'idim=',idim, &
     &                         '|| in1*iinner1=',in1,'*',iinner1,  '+ (idim-in1)*(iinner1-1)=',idim-in1,'*',iinner1-1, &
     &                         '+ 2*b1=2*',b1,'|| sum = ',in1*iinner1 + (idim-in1)*(iinner1-1)+2*b1
    write (*,'(/,2(1x,a,i3),2(1x,a,i3,1x,a,i3),2(1x,a,i3))') &
     &           'ksize=',ksize,'kdim=',kdim, &
     &                         '|| kn1*kinner1=',kn1,'*',kinner1,  '+ (kdim-kn1)*(kinner1-1)=',kdim-kn1,'*',kinner1-1, &
     &                         '+ 2*b1=2*',b1,'|| sum = ',kn1*kinner1 + (kdim-kn1)*(kinner1-1)+2*b1
  endif 
 
! It must be guaranteed that the array 'phin' is not an empty array,
! because otherwise, the halo communication would not work: 
  if(((isize - 2*b1) .lt. idim) .or. &
 &     ((ksize - 2*b1) .lt. kdim)) then
    if(my_rank.eq.0) then
      print *, 'phin is in some processes an empty array because ', &
 &             'isize-2*b1=',isize-2*b1,'< idim=',idim,' or ', &
 &             'ksize-2*b1=',ksize-2*b1,'< kdim=',kdim
    end if
    goto 500
  end if

  allocate(phi(is:ie,ks:ke))
  allocate(phin(is+b1:ie-b1,ks+b1:ke-b1))

! create a MPI Vector
  call MPI_Type_vector(kinner,b1,iouter, MPI_DOUBLE_PRECISION,vertical_border)
  call MPI_Type_commit(vertical_border)

  call MPI_Type_vector(b1,iinner,iouter, MPI_DOUBLE_PRECISION,horizontal_border)
  call MPI_Type_commit(horizontal_border)
!         
! naming: i = 1st array coordinate = 1st process coordinate = y
!         k = 2nd array coordinate = 2nd process coordinate = x
!                                                    Caution: y,x and not x,y
  dx=1.d0/(kmax-kstart)
  dy=1.d0/(imax-istart)
  dx2=dx*dx
  dy2=dy*dy
  dx2i=1.d0/dx2
  dy2i=1.d0/dy2
  dt=min(dx2,dy2)/4.d0
! start values 0.d0 
  do k=ks, min(ke,kmax-1)
    do i=max(1,is), min(ie,imax-1) ! do not overwrite the boundary condition
      phi(i,k)=0.d0
    enddo
  enddo
! start values 1.d0
 if (ke.eq.kmax) then 
  do i=is,ie
    phi(i,kmax)=1.d0
  enddo
 endif 
! start values dx
 if (is.eq.istart) then 
!       phi(0,0)=0.d0
!       do k=1,kmax-1 
!         phi(0,k)=phi(0,k-1)+dx 
!         phi(imax,k)=phi(imax,k-1)+dx 
!       enddo 
!       ... substitute algorithmus by a code, 
!           that can be computed locally: 
  do k=ks,min(ke,kmax-1)
    phi(istart,k) = (k-kstart)*dx 
  enddo
 endif 
 if (ie.eq.imax) then 
  do k=ks,min(ke,kmax-1)
    phi(imax,k) = (k-kstart)*dx 
  enddo
 endif 
 
! print details
 if (my_rank.eq.0) then
       write (*,'(/,a)') 'Heat Conduction 2d'
       write (*,'(/,4(a,1pg12.4))') 'dx =',dx,', dy =',dy,', dt =',dt,', eps =',eps
 endif 
  start_time = MPI_Wtime() 
  comm_time = 0
  criterion_time = 0 
 
! iteration
  do it=1,itmax
      dphimax=0.
      do k=ks+b1,ke-b1
       do i=is+b1,ie-b1
        dphi=(phi(i+1,k)+phi(i-1,k)-2.*phi(i,k))*dy2i &
 &             +(phi(i,k+1)+phi(i,k-1)-2.*phi(i,k))*dx2i
        dphi=dphi*dt
        dphimax=max(dphimax,dphi)
        phin(i,k)=phi(i,k)+dphi
       enddo
      enddo
! save values
      do k=ks+b1,ke-b1
       do i=is+b1,ie-b1
        phi(i,k)=phin(i,k)
       enddo
      enddo
 
! for optimization: allreduce only each stride's loop:
    criterion_time = criterion_time - MPI_Wtime() 
      if (mod(it,stride) .eq. 0) then 
        if (numprocs .gt. 1) then 
          dphimaxpartial = dphimax
          call MPI_Allreduce(dphimaxpartial, dphimax, 1,MPI_DOUBLE_PRECISION,MPI_MAX,comm)
        endif 
        if(dphimax.lt.eps) then
          criterion_time = criterion_time + MPI_Wtime() 
          goto 10 ! Finish the timestep loop "do it=…"
        endif 
      endif 
    criterion_time = criterion_time + MPI_Wtime() 
 
    comm_time = comm_time - MPI_Wtime() 
! send and receive to/from upper/lower
    if (kdim.gt.1) then ! otherwise in all processes both, lower and upper are MPI_PROC_NULL
      call MPI_Irecv(phi(is+b1,ks),   1,horizontal_border, lower,1,comm, req(1))
      call MPI_Irecv(phi(is+b1,ke),   1,horizontal_border, upper,2,comm, req(2))
      call MPI_Isend(phi(is+b1,ke-b1),1,horizontal_border, upper,1,comm, req(3))
      call MPI_Isend(phi(is+b1,ks+b1),1,horizontal_border, lower,2,comm, req(4))
      call MPI_Waitall(4, req, statuses)
    endif

! send and receive to/from left/right
    if (idim.gt.1) then ! otherwise in all processes both, left and right are MPI_PROC_NULL
      call MPI_Irecv(phi(is,ks+b1),   1,vertical_border, left, 3,comm, req(1))
      call MPI_Irecv(phi(ie,ks+b1),   1,vertical_border, right,4,comm, req(2))
      call MPI_Isend(phi(ie-b1,ks+b1),1,vertical_border, right,3,comm, req(3))
      call MPI_Isend(phi(is+b1,ks+b1),1,vertical_border, left, 4,comm, req(4))
      call MPI_Waitall(4, req, statuses)
    endif
    comm_time = comm_time + MPI_Wtime() 
 
  enddo
10      continue
  end_time = MPI_Wtime() 

  if (prt) then
    do ic = 0, idim-1 
      do kc = 0, kdim-1
        if ((ic.eq.icoord) .and. (kc.eq.kcoord)) then 
          i_first = is
          i_last  = ie
          k_first = ks
          k_last  = ke
          if (.not. prt_halo) then
            if (ic.gt.0) i_first = is + b1     ! do not print halo at the beginning
            if (ic.lt.idim-1) i_last = ie - b1 ! do not print halo at the end
            if (kc.gt.0) k_first = ks + b1     ! do not print halo at the beginning
            if (kc.lt.kdim-1) k_last = ke - b1 ! do not print halo at the end
          endif
          if (kc.eq.0) then
            write(*,'(/a,i3,a)') 'printing the ',ic,'th horizontal block'
            write(*,'(''                i='',1000(1x,i4))') &
                                             (i,i=i_first,i_last) 
          else
            if (prt_halo) write(*,*) ! additional empty line between the processes
          endif
          do k=k_first,k_last
            write(*,'(''ic='',i2,'' kc='',i2,'' k='',i3,'':'',1000(1x,f4.2))') &
                              ic,         kc,        k,       (phi(i,k),i=i_first,i_last) 
          enddo
        endif 
        call MPI_Barrier(comm) ! ito separate the printing of each block by different processes.
                                      ! Caution: This works in most cases, but does not guarantee that the
                                      !          output lines on the common stdout are in the expected sequence.
      enddo 
    enddo 
  endif 

  if (my_rank.eq.0) then
   write (*,'(/a/a/)') &
 &'!numprocs=idim    iter-  wall clock time communication part abort criterion', &
 &'!          x kdim ations    [seconds]    method [seconds]   meth. stride [seconds]'
   write(*, &
 &   '(''!''i6,'' ='',i3,'' x'',i3,1x,i6,2x,g12.4,4x,i2,2x,g12.4,2x,i2,2x,i6,2x,g12.4)') &
 &    numprocs,    idim,    kdim,   it, end_time - start_time, 1, comm_time, 1, stride, criterion_time 
  endif 
!

500 continue

  call MPI_Finalize()

end program heat
