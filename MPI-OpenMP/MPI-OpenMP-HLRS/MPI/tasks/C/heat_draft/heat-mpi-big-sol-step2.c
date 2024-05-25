#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))
#define idx(i,k)  (((i)-is)  *kouter + (k)-ks)
#define idxn(i,k) (((i)-is-1)*kinner + (k)-ks-1)
int main(int argc, char *argv[])
{
// naming: i = 1st array coordinate = 1st process coordinate = y
//         k = 2nd array coordinate = 2nd process coordinate = x
//                                                    Caution: y,x and not x,y
  int i, k, it;
  int const prt=1, prt_halo=1; // 1 or 0 for printing / not printing the result / plus halos
  int const imax=15, kmax=12, istart=0, kstart=0, b1=1, itmax=20000;
  double const eps=1.e-08;
  double *phi, *phin, dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax;

// Preparation for parallelization with domain decomposition:

//   int is=istart, ie=imax, ks=kstart, ke=kmax; // now defined and calculated below

// The algortithm below is already formulated with lower and upper 
// index values (is:ie), (ks:ke) instead of (istart:imax), (kstart:kmax).
// With the MPI domain decomposition, (is:ie) and (ks:ke) will define the index rang 
// of the subdomain in a given MPI process.

// additional variables for parallelization
  int numprocs, my_rank, right, left, upper, lower, ierror;
  MPI_Comm comm;  // Cartesian communicator
  int dims[2], coords[2]; // helper variables only for MPI_CART_CREATE and ..._COORDS
  int period[2];          //                  only for MPI_CART_CREATE and ..._COORDS
  int idim, kdim, icoord, kcoord;
  int isize, iinner1, in1,  ksize, kinner1, kn1; // only for calculation of iinner, ...
  int iinner, iouter, is, ie, kinner, kouter, ks, ke;
  int ic, i_first, i_last,  kc, k_first, k_last; // only for printing
  int const stride=10; // for reduced calculation of the abort criterion
  MPI_Request req[4];
  MPI_Status statuses[4];
  MPI_Datatype horizontal_border, vertical_border; // datatype
  double dphimaxpartial; // for local/global calculation of the abort criterion
  double start_time, end_time, comm_time, criterion_time;

//  naming: originally: i=istart..imax,      now: i=is..ie,  icoord=0..idim-1
//          originally: k=kstart..kmax,      now: k=ks..ke,  kcoord=0..kdim-1
//                      with istart=kstart=0      s=start index,  e=end index

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

// the 2-dimensional domain decomposition:

  dims[0] = 0; dims[1] = 0;
  period[0] = 0; period[1] = 0;
  MPI_Dims_create(numprocs, 2, dims);
  idim = dims[0];
  kdim = dims[1];
  MPI_Cart_create(MPI_COMM_WORLD,2,dims,period,1,&comm);
  MPI_Comm_rank(comm, &my_rank);
  MPI_Cart_coords(comm, my_rank, 2, coords);
  icoord = coords[0];
  kcoord = coords[1];
  MPI_Cart_shift(comm, 0, 1, &left, &right);
;   // the ranks (left,right) represent the coordinates ((icoord-1,kcoord), (icoord+1,kcoord)) 
  MPI_Cart_shift(comm, 1, 1, &lower,&upper);
;   // the ranks (lower,upper) represent the coordinates ((icoord,kcoord-1),(icoord,kcoord+1)) 

// whole y indecees   |------------- isize = 1+imax-istart ------------|
//    start/end index  ^-istart                                  imax-^
// 
// 1. interval        |--- iouter1---|   
//                    |--|--------|--|
//                     b1  iinner1 b1 
//    start/end index  ^-is      ie-^ 
// 2. interval                 |--|--------|--|   
// 3. interval                          |--|--------|--| 
// 4. interval                                   |--|-------|--| 
//                                                   iinner0 = inner1 - 1
// 5. interval = idim's interval                         |--|-------|--|
;
// In each iteration on each interval, the inner area is computed
// by using the values of the last iteration in the whole outer area. 
;
// icoord = number of the interval - 1
// 
// To fit exactly into isize, we use in1 intervals of with iinner1
// and (idim-in1) intervals of with (iinner1 - 1) 
;
// isize <= 2*b1 + idim * iinner1  
;
// ==>   iinner1 >= (isize - 2*b1) / idim
// ==>   smallest valid integer value:


// computing is:ie, ks:ke  ---  details see heat-mpi1-big.f 
  isize = 1 + imax - istart;
  iinner1 = (isize - 2*b1 - 1) / idim + 1;
  in1 = isize - 2*b1 - idim * (iinner1 - 1);
  if (icoord < in1) {
    iinner = iinner1;
    is = istart + icoord * iinner;
  } else {
    iinner = iinner1 - 1;
    is = istart + in1 * iinner1 + (icoord-in1) * iinner;
  }
  iouter = iinner + 2*b1;
  ie = is + iouter - 1;
// same for x coordinate: 
  ksize = 1 + kmax - kstart;
  kinner1 = (ksize - 2*b1 - 1) / kdim + 1;
  kn1 = ksize - 2*b1 - kdim * (kinner1 - 1);
  if (kcoord < kn1) {
    kinner = kinner1;
    ks = kstart + kcoord * kinner;
  } else {
    kinner = kinner1 - 1;
    ks = kstart + kn1 * kinner1 + (kcoord-kn1) * kinner;
  }
  kouter = kinner + 2*b1;
  ke = ks + kouter - 1;

  if (my_rank == 0) {
    printf("\n isize=%3d idim=%3d || in1*iinner1=%3d *%3d + (idim-in1)*(iinner1-1)=%3d *%3d + 2*b1=2*%1d || sum = %3d\n",
               isize,    idim,       in1,iinner1,            idim-in1,  iinner1-1,              b1,
                                                                             in1*iinner1 + (idim-in1)*(iinner1-1)+2*b1);
    printf("\n ksize=%3d kdim=%3d || kn1*kinner1=%3d *%3d + (kdim-kn1)*(kinner1-1)=%3d *%3d + 2*b1=2*%1d || sum = %3d\n",
               ksize,    kdim,       kn1,kinner1,            kdim-kn1,  kinner1-1,              b1,
                                                                             kn1*kinner1 + (kdim-kn1)*(kinner1-1)+2*b1);
  }



// It must be guaranteed that the array 'phin' is not an empty array,
// because otherwise, the halo communication would not work: 
  if(((isize - 2*b1) < idim) ||
     ((ksize - 2*b1) < kdim)) {
    if(my_rank == 0) {
        printf("phin is in some processes an empty array because isize-2*b1=%3d < idim=%3d or ksize-2*b1=%3d < kdim=%3d\n",
                                                                 isize-2*b1,      idim,       ksize-2*b1,      kdim);
    }
    MPI_Finalize();
    exit(0);
  }

  phi  = (double *) malloc(iouter*kouter*sizeof(double)); // index range: (is:ie,ks:ke)
  phin = (double *) malloc(iinner*kinner*sizeof(double)); // index range: (is+b1:ie-b1,ks+b1:ke-b1)

// create a MPI Vector
  MPI_Type_vector(b1,kinner,kouter, MPI_DOUBLE, &vertical_border);
  MPI_Type_commit(&vertical_border);

  MPI_Type_vector(iinner,b1,kouter, MPI_DOUBLE, &horizontal_border);
  MPI_Type_commit(&horizontal_border);
//         
// naming: i = 1st array coordinate = 1st process coordinate = y
//         k = 2nd array coordinate = 2nd process coordinate = x
//                                                    Caution: y,x and not x,y
  dx=1.e0/(kmax-kstart);
  dy=1.e0/(imax-istart);
  dx2=dx*dx;
  dy2=dy*dy;
  dx2i=1.e0/dx2;
  dy2i=1.e0/dy2;
  dt=min(dx2,dy2)/4.e0;
// start values 0.d0 
  for(i=max(1,is); i<=min(ie,imax-1); i++) { // do not overwrite the boundary condition
    for(k=ks; k<=min(ke,kmax-1); k++) {
      phi[idx(i,k)]=0.e0;
    } /*for*/
  } /*for*/
// start values 1.d0
 if (ke == kmax) { 
  for(i=is; i<=ie; i++) {
    phi[idx(i,kmax)]=1.e0;
  } /*for*/
 }
// start values dx
 if (is == istart) { 
//       phi[idx(0,0)]=0.d0
//       for(k=1; k<=kmax-1 ; k++) { 
//         phi[idx(0,k)]=phi[idx(0,k-1)]+dx 
//         phi[idx(imax,k)]=phi[idx(imax,k-1)]+dx 
//       } /*for*/ 
//       ... substitute algorithmus by a code, 
//           that can be computed locally: 
  for(k=ks; k<=min(ke,kmax-1); k++) {
    phi[idx(istart,k)] = (k-kstart)*dx;
  } /*for*/
 }
 if (ie == imax) { 
  for(k=ks; k<=min(ke,kmax-1); k++) {
    phi[idx(imax,k)] = (k-kstart)*dx;
  } /*for*/
 }

// print details
 if (my_rank == 0) {
    printf("\nHeat Conduction 2d\n");
    printf("\ndx =%12.4e, dy =%12.4e, dt=%12.4e, eps=%12.4e\n", dx, dy, dt, eps); 
 }
  start_time = MPI_Wtime();
  comm_time = 0;
  criterion_time = 0;

// iteration
  for(it=1; it<=itmax; it++) {
      dphimax=0.;
      for(i=is+b1; i<=ie-b1; i++) {
       for(k=ks+b1; k<=ke-b1; k++) {
        dphi=(phi[idx(i+1,k)]+phi[idx(i-1,k)]-2.*phi[idx(i,k)])*dy2i
               +(phi[idx(i,k+1)]+phi[idx(i,k-1)]-2.*phi[idx(i,k)])*dx2i;
        dphi=dphi*dt;
        dphimax=max(dphimax,dphi);
        phin[idxn(i,k)]=phi[idx(i,k)]+dphi;
       } /*for*/
      } /*for*/
// save values
      for(i=is+b1; i<=ie-b1; i++) {
       for(k=ks+b1; k<=ke-b1; k++) {
        phi[idx(i,k)]=phin[idxn(i,k)];
       } /*for*/
      } /*for*/

// for optimization: allreduce only each stride's loop:
    criterion_time = criterion_time - MPI_Wtime();
      if ((it % stride) == 0) { 
        if (numprocs > 1) { 
          dphimaxpartial = dphimax;
          MPI_Allreduce(&dphimaxpartial, &dphimax, 1, MPI_DOUBLE, MPI_MAX, comm);
        }
        if(dphimax < eps) {
          criterion_time = criterion_time + MPI_Wtime();
          goto endOfLoop; // Finish the timestep loop "do it=…"
        }
      }
    criterion_time = criterion_time + MPI_Wtime();

    comm_time = comm_time - MPI_Wtime();
// send and receive to/from upper/lower
    if (kdim > 1) { // otherwise in all processes both, lower and upper are MPI_PROC_NULL
//    MPI_...
//    MPI_...
//    MPI_...
//    MPI_...
//    MPI_...
    }

// send and receive to/from left/right
    if (idim > 1) { // otherwise in all processes both, left and right are MPI_PROC_NULL
//    MPI_...
//    MPI_...
//    MPI_...
//    MPI_...
//    MPI_...
    }
    comm_time = comm_time + MPI_Wtime();

  } /*for*/
endOfLoop:
  end_time = MPI_Wtime();

  if (prt) {
    for(ic = 0; ic <= idim-1; ic++) {
      for(kc = 0; kc <= kdim-1; kc++) {
        if ((ic == icoord) && (kc == kcoord)) { 
          i_first = is ;
          i_last  = ie ;
          k_first = ks ;
          k_last  = ke ;
          if (! prt_halo) {
            if (ic > 0) i_first = is + b1;     // do not print halo at the beginning
            if (ic < idim-1) i_last = ie - b1; // do not print halo at the end
            if (kc > 0) k_first = ks + b1;     // do not print halo at the beginning
            if (kc < kdim-1) k_last = ke - b1; // do not print halo at the end
          }
          if (kc == 0) {
              printf("\nprinting the %3dth horizontal block\n", ic);
              printf("               i=");
              for (i=i_first; i<=i_last; i++) printf(" %4d",i); printf("\n");
          } else {
            if (prt_halo)  printf("\n"); // additional empty line between the processes
          }
          for(k=k_first; k<=k_last; k++) {
              printf("ic=%2d kc=%2d k=%3d", ic, kc, k);
              for (i=i_first; i<=i_last; i++) printf(" %4.2f", phi[idx(i,k)]); printf("\n");
          } /*for*/
        }
        MPI_Barrier(comm); // ito separate the printing of each block by different processes.
                           // Caution: This works in most cases, but does not guarantee that the
                           //          output lines on the common stdout are in the expected sequence.
      } /*for*/
    } /*for*/
  } /*if prt*/

  if (my_rank == 0) {
   printf("\n!numprocs=idim    iter-   wall clock time  communication part  abort criterion\n");
   printf(  "!          x kdim ations     [seconds]     method [seconds]    meth. stride [seconds]\n");
   printf(  "!%6d =%3d x%3d %6d %12.4g      %2d %12.4g    %2d %6d %12.4g\n",
          numprocs,idim,kdim,it, end_time - start_time,
                                            1, comm_time, 1, stride, criterion_time);
  }





  MPI_Finalize();

}
