
#include <iostream>
#include <array>
#include <span>
#include <print>
#include <format>
#include <chrono>
#include <cmath>

//******************************************************
//? Constantes
constexpr int nx = 21;
constexpr int ny = 7;

constexpr int maxrow = 27*ny;
constexpr int nelemn = 2*(nx-1)*(ny-1);
constexpr int mx = 2*nx-1;
constexpr int my = 2*ny-1;
constexpr int np = mx*my;
constexpr int maxeqn = 2*mx*my+nx*ny;
constexpr int nnodes = 6;
constexpr int nquad = 3;


//******************************************************
void timestamp() {
  auto now = std::chrono::system_clock::now();

  std::println("Current time: {}", std::format("{:%Y-%m-%d %H:%M:%S}", now));
}

//******************************************************
void setgrd(
    std::array<std::array<int, 2>, np> &indx,
    std::array<int, np> &insc,
    std::array<int, nelemn> &isotri,
    int &iwrite, bool &_long, int &neqn,
    std::array<std::array<int, nnodes>, nelemn> &node
);

//******************************************************
void setban(
  std::array<std::array<int, 2>, np>& indx,
  std::array<int, np>& insc,
  int& nband,
  int& nlband,
  std::array<std::array<int, nnodes>, nelemn>& node,
  int& nrow
);

//******************************************************
void setlin(
  std::array<int, my>& iline,
  std::array<std::array<int, 2>, np>& indx,
  int& iwrite,
  bool& _long,
  int& nodex0,
  double& xlngth
);

//******************************************************
void setxy(
  int iwrite,
  bool _long,
  std::array<double, np> &xc,
  double &xlngth,
  std::array<double, np> &yc,
  double &ylngth
);

//******************************************************
void setqud(
  std::array<double, nelemn> &area,
  std::array<std::array<int, nnodes>, nelemn> &node,
  std::array<double, np> &xc,
  std::array<std::array<double, nquad>, nelemn> &xm,
  std::array<double, np> &yc,
  std::array<std::array<double, nquad>, nelemn> &ym
);

//******************************************************
void setbas(
  std::array<std::array<int, nnodes>, nelemn> &node,
  std::array<double, np> &xc,
  std::array<double, np> &yc,
  std::array<std::array<std::array<std::array<double, 3>, nnodes>, nquad>, nelemn> &phi,
  std::array<std::array<std::array<double, nnodes>, nquad>, nelemn> &psi,
  std::array<std::array<double, nquad>, nelemn> &xm,
  std::array<std::array<double, nquad>, nelemn> &ym
);

//******************************************************
double bsp(
  double xq,
  double yq,
  int it,
  int iq,
  int id,
  const std::array<std::array<int, nnodes>, nelemn> &node,
  const std::array<double, np> &xc,
  const std::array<double, np> &yc
);

//******************************************************
void qbf(
  double x,
  double y,
  int it,
  int in,
  double &bb,
  double &bx,
  double &by,
  const std::array<std::array<int, nnodes>, nelemn> &node,
  const std::array<double, np> &xc,
  const std::array<double, np> &yc
);

//******************************************************
void nstoke(
  std::array<std::array<double, maxeqn>, maxrow> &a,
  std::array<double, nelemn> &area,
  std::array<double, maxeqn> &f,
  std::array<double, maxeqn> &g,
  std::array<std::array<int, 2>, np> &indx,
  std::array<int, np> &insc,
  std::array<int, maxeqn> &ipivot,
  int iwrite,
  int &maxnew,
  int neqn, int nlband,
  std::array<std::array<int, nnodes>, nelemn> &node,
  int &nrow,
  int &numnew,
  double &para,
  std::array<std::array<std::array<std::array<double, 3>, nnodes>, nquad>, nelemn> &phi,
  std::array<std::array<std::array<double, nnodes>, nquad>, nelemn> &psi,
  double &reynld,
  double &tolnew,
  std::array<double, np>& yc
);

//******************************************************
void linsys(
  std::array<std::array<double, maxeqn>, maxrow>& a,
  std::array<double, nelemn>& area,
  std::array<double, maxeqn>& f,
  std::array<double, maxeqn>& g,
  std::array<std::array<int, 2>, np>& indx,
  std::array<int, np>& insc,
  std::array<int, maxeqn>& ipivot,
  int neqn,
  int nlband,
  std::array<std::array<int, nnodes>, nelemn>& node,
  int nrow,
  double para1,
  double para2,
  std::array<
    std::array<
      std::array<std::array<double, 3>, nnodes>,
      nquad
    >, nelemn
  >& phi,
  std::array<
    std::array<
      std::array<double, nnodes>, nquad
    >, nelemn
  >& psi,
  double reynld,
  std::array<double, np>& yc
);

//******************************************************
void uval(
  std::array<double, maxeqn> &g,
  std::array<std::array<int, 2>, np> &indx,
  int iquad,
  int it,
  int neqn,
  std::array<std::array<int, nnodes>, nelemn> &node,
  double para,
  std::array<
  std::array<
  std::array<std::array<double, 3>, nnodes>,
  nquad
  >, nelemn
  > &phi,
  std::array<double, 2> &un,    // ← tamaño fijo 2
  std::array<double, 2> &unx,   // ← tamaño fijo 2
  std::array<double, 2> &uny,   // ← tamaño fijo 2
  std::array<double, np> &yc
);
//******************************************************
double ubdry (double y, double para);

int main(void)
{
  std::array<std::array<double, maxeqn>, maxrow> a;
  double a2;
  double abound;
  double anew;
  double aold;
  std::array<double, nelemn> area;
  std::array<double, my> dcda;
  std::array<double, maxeqn> f;
  std::string fileg;
  std::string fileu;
  std::string filex;
  std::array<double, maxeqn> g;
  std::array<std::array<double, my>, my> gr;
  int i;
  std::array<int, my> iline;
  std::array<std::array<int, 2>, np> indx;
  std::array<int, np> insc;
  int iounit;
  std::array<int, maxeqn> ipivot;
  std::array<int, nelemn> isotri;
  int iter;
  int ivunit;
  int iwrite;
  int ixunit;
  int j;
  bool _long;
  int maxnew;
  int maxsec;
  int nband;
  int neqn;
  int nlband;
  std::array<std::array<int, nnodes>, nelemn> node;
  int nodex0;
  int npara;
  int nrow;
  int numnew;
  int numsec;
  double para;
  std::array<
    std::array<
      std::array<std::array<double, 3>, nnodes>,
      nquad
    >, nelemn
  > phi;
  std::array<
    std::array<
      std::array<double, nnodes>, nquad
    >, nelemn
  > psi;
  std::array<double, my> r;
  std::array<double, maxeqn> res;
  double reynld;
  double rjpnew;
  double rjpold;
  double temp;
  double test;
  double tolnew;
  double tolsec;
  std::array<double, my> ui;
  std::array<double, my> unew;
  std::array<double, np> xc;
  double xlngth;
  std::array<std::array<double, nquad>, nelemn> xm;
  std::array<double, np> yc;
  double ylngth;
  std::array<std::array<double, nquad>, nelemn> ym;

  timestamp();
  std::println(" ");
  std::println("channel():");
  std::println("  C++ version");
  std::println("  Channel flow control problem");
  std::println(" ");
  std::println("  Flow control problem:");
  std::println("    Inflow controlled by one parameter.");
  std::println("    Velocities measured along vertical line.");
  std::println("    Try to match specified velocity profile.");
//!
//!  Set input data
//!
  fileg = "display.txt";
  fileu = "uv.txt";
  filex = "xy.txt";
  iounit = 2;
  ivunit = 4;
  iwrite = 10;
  ixunit = 3;
  maxnew = 10;
  maxsec = 8;
  npara = 1;
  numnew = 0;
  numsec = 0;
  reynld = 1.0;
  rjpnew = 0.0;
  tolnew = 1.0E-04;
  tolsec = 1.0E-06;
  xlngth = 10.0;
  ylngth = 3.0;

  std::println("");
  std::println("NX = {}", nx);
  std::println("NY = {}", ny);
  std::println("Number of elements = {}", nelemn);
  std::println("Reynolds number = {}", reynld);
  std::println("Secant tolerance = {}", tolsec);
  std::println("Newton tolerance = {}", tolnew);
  std::println("");
//!
//!  SETGRD constructs grid, numbers unknowns, calculates areas,
//!  and points for midpoint quadrature rule.
//!
  setgrd(indx, insc, isotri, iwrite, _long, neqn, node);
//!
//!  Compute the bandwidth
//!
  setban(indx, insc, nband, nlband, node, nrow);
//!
//!  Record variable numbers along profile sampling line.
//!
  setlin(iline, indx, iwrite, _long, nodex0, xlngth);
//!
//!  Set the coordinates of grid points.
//!
  setxy(iwrite, _long, xc, xlngth, yc, ylngth);
//!
//!  Set quadrature points
//!
  setqud(area, node, xc, xm, yc, ym);
//!
//!  Evaluate basis functions at quadrature points
//!
  setbas(node, xc, yc, phi, psi, xm, ym);
//!
//!  NSTOKE now solves the Navier Stokes problem for an inflow
//!  parameter of 1.0.
//!
  para = 1.0;
  std::println(" ");
  std::println("Solve Navier Stokes problem with parameter = {}", para);
  std::println("for profile at x = {}", xc[nodex0-1]);
  for (int i = 0; i < neqn; ++i) {
    g[i] = 1.0;
  }

  nstoke(
    a, area, f, g, indx, insc, ipivot, iwrite,
    maxnew, neqn, nlband, node, nrow, numnew,
    para, phi, psi,
    reynld, tolnew, yc
  );
//!
//!  RESID computes the residual at the given solution
//!
//  if ( 1 <= iwrite ) then
//    call resid (area,f,indx,insc,iwrite,nelemn,neqn, &
//      nnodes,node,np,nquad,para,phi,psi,res,reynld,yc)
//  end if
//!
//!  GETG computes the internal velocity profile at X = XC(NODEX0), which will
//!  be used to measure the goodness-of-fit of the later solutions.
//!
//  call getg ( f, iline, my, neqn, ui )
//
//  if ( 1 <= iwrite ) then
//    std::println(" ")
//    std::println("U profile:")
//    std::println(" ")
//    write (*,'(5g14.6)') ui(1:my)
//  end if
//!
//!  GRAM generates the Gram matrix GR and the vector
//!  R = line integral of ui*phi
//!
//  call gram (gr,iline,indx,iwrite,my,nelemn,nnodes,node, &
//    nodex0,np,para,r,ui,xc,yc)
//!
//!  GDUMP dumps information for graphics display by DISPLAY.
//!
//  if ( .false. ) then
//    std::println("Writing graphics data to file '//fileg
//    call delete(fileg)
//    open (unit = iounit,file=fileg,form='formatted',status='new', &
//      err = 50)
//    rjpnew = 0.0
//    call gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, &
//      nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
//  end if
//!
//!  Write the XY data to a file.
//!
//  if ( .false. ) then
//    call delete(filex)
//    open(unit = ixunit,file=filex,form='formatted',status='new')
//    call xy_plot3d (ixunit,long,np,nx,ny,xc,yc)
//    close(unit = ixunit)
//  else
//    call delete ( filex )
//    open ( unit = ixunit, file = filex, form = 'formatted', &
//      status = 'new')
//    call xy_table ( ixunit, np, xc, yc )
//    close ( unit = ixunit )
//  end if
//!
//!  Write the velocity data to a file.
//!
//  if ( .false. ) then
//    call delete(fileu)
//    open(unit = ivunit,file=fileu,form='formatted',status='new')
//    call uv_plot3d (f,indx,insc,ivunit,long,mx,my, &
//      nelemn,neqn,nnodes,node,np,para,a,reynld,yc)
//    close(unit = ivunit)
//  else
//    call delete(fileu)
//    open(unit = ivunit,file=fileu,form='formatted',status='new')
//    call uv_table ( f, indx, ivunit, neqn, np, para, yc )
//    close(unit = ivunit)
//  end if
//!
//!  Destroy information about true solution
//!
//  f(1:neqn) = 0.0
//  g(1:neqn) = 0.0
//!
//!  Secant iteration loop
//!
//  aold = 0.0
//  rjpold = 0.0
//  anew = 0.1
//
//  do iter = 1, maxsec
//
//    numsec = numsec+1
//    std::println(" ")
//    std::println("Secant iteration ',iter
//!
//!  Solve for unew at new value of parameter anew
//!
//    std::println(" ")
//    std::println("Solving Navier Stokes problem for parameter = ',anew
//!
//!  Use solution F at previous value of parameter for starting point.
//!
//    call dcopy ( neqn, f, 1, g, 1 )
//    para = anew
//
//    call nstoke (a,area,f,g,indx,insc,ipivot,iwrite, &
//      maxnew,maxrow,nelemn,neqn,nlband,nnodes,node, &
//      np,nquad,nrow,numnew,para,phi,psi,reynld,tolnew,yc)
//!
//!  Get velocity profile
//!
//    call getg ( f, iline, my, neqn, unew )
//
//    if ( 1 <= iwrite ) then
//      std::println(" ")
//      std::println("Velocity profile:")
//      std::println(" ")
//      write (*,'(5g14.6)') unew(1:my)
//    end if
//!
//!  Solve linear system for du/da
//!
//    para = anew
//    abound = 1.0
//    call linsys (a,area,g,f,indx,insc,ipivot, &
//      maxrow,nelemn,neqn,nlband,nnodes,node, &
//      np,nquad,nrow,para,abound,phi,psi,reynld,yc)
//!
//!  Output in DCDA
//!
//    call getg ( g, iline, my, neqn, dcda )
//
//    if ( 2 <= iwrite ) then
//      std::println(" ")
//      std::println("Sensitivities:")
//      std::println(" ")
//      write (*,'(5g14.6)') dcda(1:my)
//    end if
//!
//!  Evaluate J prime at current value of parameter where J is
//!  functional to be minimized.
//!
//!  JPRIME = 2.0 * DCDA(I) * (GR(I,J)*UNEW(J)-R(I))
//!
//    rjpnew = 0.0
//    do i = 1, my
//      temp = -r(i)
//      do j = 1, my
//        temp = temp + gr(i,j) * unew(j)
//      end do
//      rjpnew = rjpnew + 2.0 * dcda(i) * temp
//    end do
//
//    std::println(" ")
//    std::println("Parameter  = ',anew,' J prime = ',rjpnew
//!
//!  Dump information for graphics
//!
//    if ( .false. ) then
//      para = anew
//      call gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, &
//        nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
//    end if
//!
//!  Update the estimate of the parameter using the secant step
//!
//    if (iter == 1) then
//      a2 = 0.5
//    else
//      a2 = aold-rjpold*(anew-aold)/(rjpnew-rjpold)
//    end if
//
//    aold = anew
//    anew = a2
//    rjpold = rjpnew
//    test = abs(anew-aold)/abs(anew)
//
//    std::println("New value of parameter = ',anew
//    std::println("Convergence test = ',test
//
//    if (abs(anew-aold) < abs(anew)*tolsec) then
//      std::println("Secant iteration converged.")
//      go to 40
//    end if
//
//  end do
//
//  std::println("Secant iteration failed to converge.")
//   40 continue
//
//  write (*,*)'Number of secant steps = ', numsec
//  write (*,*)'Number of Newton steps = ', numnew
//!
//!  Close graphics file
//!
//  if ( .false. ) then
//    close ( unit = iounit )
//  end if
//!
//!  Terminate.
//!
//  std::println( " ")
//  std::println( "CHANNEL:")
//  std::println( "  Normal end of execution.")
//
//  std::println( " ")
//  call timestamp ( )
//
//  stop
//!
//!  Error opening graphics file
//!
//   50 continue
//  std::println("CHANNEL could not open the graphics file!")
}

//******************************************************
double bsp(
  double xq,
  double yq,
  int it,
  int iq,
  int id,
  const std::array<std::array<int, nnodes>, nelemn> &node,
  const std::array<double, np> &xc,
  const std::array<double, np> &yc
) {
//? bsp() evaluates the linear basis functions associated with pressure.
  int iq1 = iq;
  int iq2 = (iq%3);
  int iq3 = ((iq+1)%3);
  int i1 = node[it][iq1]-1;
  int i2 = node[it][iq2]-1;
  int i3 = node[it][iq3]-1;

  double d = (xc[i2] - xc[i1])*(yc[i3] - yc[i1]) -
    (xc[i3] - xc[i1])*(yc[i2] - yc[i1]);

  if (id == 1) {
    return 1.0+(
      (yc[i2]-yc[i3])*(xq-xc[i1])+(xc[i3]-xc[i2])*(yq-yc[i1])
    )/d;
  } else if (id == 2) {
    return (yc[i2]-yc[i3])/d;
  } else if (id == 3) {
    return (xc[i3]-xc[i2])/d;
  } else {
    std::println("BSP - fatal error!");
    std::println("unknown value of id = {}", id);
    std::exit(1);
  }
}

//subroutine daxpy ( n, da, dx, incx, dy, incy )
//
//******************************************************
//!
//!! daxpy() computes constant times a vector plus a vector.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    Jack Dongarra
//!
//!  Parameters:
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) dx(*),dy(*),da
//  integer i,incx,incy,ix,iy,m,n
//!
//  if ( n <= 0)return
//  if (da  ==  0.0 ) return
//  if ( incx == 1.and.incy == 1)go to 20
//!
//!        code for unequal increments or equal increments
//!          not equal to 1
//!
//  ix = 1
//  iy = 1
//  if ( incx < 0)ix = (-n+1)*incx + 1
//  if ( incy < 0)iy = (-n+1)*incy + 1
//  do i = 1,n
//    dy(iy) = dy(iy) + da*dx(ix)
//    ix = ix + incx
//    iy = iy + incy
//  end do
//  return
//!
//!        code for both increments equal to 1
//!
//!
//!        clean-up loop
//!
//   20 m = mod(n,4)
//  if (  m  ==  0 ) go to 40
//  do 30 i = 1,m
//    dy(i) = dy(i) + da*dx(i)
//   30 continue
//  if (  n  <  4 ) return
//   40 continue
//
//  do i = m+1, n, 4
//    dy(i) = dy(i) + da*dx(i)
//    dy(i + 1) = dy(i + 1) + da*dx(i + 1)
//    dy(i + 2) = dy(i + 2) + da*dx(i + 2)
//    dy(i + 3) = dy(i + 3) + da*dx(i + 3)
//  end do
//
//  return
//end
//subroutine dcopy ( n, dx, incx, dy, incy )
//
//******************************************************
//!
//!! dcopy() copies a vector, x, to a vector, y.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    Jack Dongarra
//!
//!  Parameters:
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) dx(*),dy(*)
//  integer i,incx,incy,ix,iy,m,n
//
//  if ( n <= 0)return
//  if ( incx == 1.and.incy == 1)go to 20
//!
//!        code for unequal increments or equal increments
//!          not equal to 1
//!
//  ix = 1
//  iy = 1
//  if ( incx < 0)ix = (-n+1)*incx + 1
//  if ( incy < 0)iy = (-n+1)*incy + 1
//  do i = 1,n
//    dy(iy) = dx(ix)
//    ix = ix + incx
//    iy = iy + incy
//  end do
//  return
//!
//!        code for both increments equal to 1
//!
//!
//!        clean-up loop
//!
//   20 m = mod(n,7)
//  if (  m  ==  0 ) go to 40
//
//  dy(1:m) = dx(1:m)
//
//  if (  n  <  7 ) return
//   40 continue
//  do i = m+1, n ,7
//    dy(i) = dx(i)
//    dy(i + 1) = dx(i + 1)
//    dy(i + 2) = dx(i + 2)
//    dy(i + 3) = dx(i + 3)
//    dy(i + 4) = dx(i + 4)
//    dy(i + 5) = dx(i + 5)
//    dy(i + 6) = dx(i + 6)
//  end do
//
//  return
//end
//function ddot ( n, dx, incx, dy, incy )
//
//******************************************************
//!
//!! ddot() forms the dot product of two vectors.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    Jack Dongarra
//!
//!  Parameters:
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) ddot
//  real ( kind = rk8 ) dx(*),dy(*),dtemp
//  integer i,incx,incy,ix,iy,m,n
//
//  ddot = 0.0
//  dtemp = 0.0
//  if ( n <= 0)return
//  if ( incx == 1.and.incy == 1)go to 20
//!
//!        code for unequal increments or equal increments
//!          not equal to 1
//!
//  ix = 1
//  iy = 1
//  if ( incx < 0)ix = (-n+1)*incx + 1
//  if ( incy < 0)iy = (-n+1)*incy + 1
//  do i = 1,n
//    dtemp = dtemp + dx(ix)*dy(iy)
//    ix = ix + incx
//    iy = iy + incy
//  end do
//
//  ddot = dtemp
//  return
//!
//!        code for both increments equal to 1
//!
//!
//!        clean-up loop
//!
//   20 m = mod(n,5)
//  if (  m  ==  0 ) go to 40
//  do i = 1,m
//    dtemp = dtemp + dx(i)*dy(i)
//  end do
//
//  if (  n  <  5 ) go to 60
//   40 continue
//  do i = m+1, n, 5
//    dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
//        dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
//  end do
//   60 ddot = dtemp
//  return
//end
//subroutine delete ( filnam )
//
//******************************************************
//!
//!! delete() deletes a file.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  character ( len = * ) filnam
//  integer ios
//
//  open (unit = 99,file=filnam,status='old', iostat = ios )
//
//  if ( ios /= 0 ) then
//    return
//  end if
//
//  close (unit = 99,status='delete', iostat = ios)
//
//  return
//end
//subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )
//
//******************************************************
//!
//!! dgbfa() factors a band matrix by elimination.
//!
//!  Discussion:
//!
//!     dgbfa is usually called by dgbco, but it can be called
//!     directly with a saving in time if  rcond  is not needed.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    Cleve Moler
//!
//!  Parameters:
//!
//!     on entry
//!
//!        abd     real ( kind = rk8 )(lda, n)
//!                contains the matrix in band storage.  the columns
//!                of the matrix are stored in the columns of  abd  and
//!                the diagonals of the matrix are stored in rows
//!                ml+1 through 2*ml+mu+1 of  abd .
//!                see the comments below for details.
//!
//!        lda     integer
//!                the leading dimension of the array  abd .
//!                2*ml + mu + 1 <= LDA.
//!
//!        n       integer
//!                the order of the original matrix.
//!
//!        ml      integer
//!                number of diagonals below the main diagonal.
//!                0 <= ml < n .
//!
//!        mu      integer
//!                number of diagonals above the main diagonal.
//!                0 <= mu < n .
//!                more efficient if  ml <= mu .
//!     on return
//!
//!        abd     an upper triangular matrix in band storage and
//!                the multipliers which were used to obtain it.
//!                the factorization can be written  a = l*u  where
//!                l  is a product of permutation and unit lower
//!                triangular matrices and  u  is upper triangular.
//!
//!        ipvt    integer(n)
//!                an integer vector of pivot indices.
//!
//!        info    integer
//!                = 0  normal value.
//!                = k  if  u(k,k) == 0.0 .  this is not an error
//!                     condition for this subroutine, but it does
//!                     indicate that dgbsl will divide by zero if
//!                     called.  use  rcond  in dgbco for a reliable
//!                     indication of singularity.
//!
//!     band storage
//!
//!           if  a  is a band matrix, the following program segment
//!           will set up the input.
//!
//!                   ml = (band width below the diagonal)
//!                   mu = (band width above the diagonal)
//!                   m = ml + mu + 1
//!                   do j = 1, n
//!                      i1 = max ( 1, j-mu )
//!                      i2 = min ( n, j+ml )
//!                      do i = i1, i2
//!                         k = i - j + m
//!                         abd(k,j) = a(i,j)
//!                      end do
//!                   end do
//!
//!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
//!           in addition, the first  ml  rows in  abd  are used for
//!           elements generated during the triangularization.
//!           the total number of rows needed in  abd  is  2*ml+mu+1 .
//!           the  ml+mu by ml+mu  upper left triangle and the
//!           ml by ml  lower right triangle are not referenced.
//!
//!     linpack. this version dated 08/14/78 .
//!     cleve moler, university of new mexico, argonne national lab.
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer lda
//  integer n
//
//  real ( kind = rk8 ) abd(lda,n)
//  integer i
//  integer i0
//  integer info
//  integer ipvt(n)
//  integer idamax
//  integer j
//  integer j0
//  integer j1
//  integer ju
//  integer jz
//  integer k
//  integer l
//  integer lm
//  integer m
//  integer ml
//  integer mm
//  integer mu
//  real ( kind = rk8 ) t
//
//  m = ml + mu + 1
//  info = 0
//!
//!  Zero initial fill-in columns.
//!
//  j0 = mu + 2
//  j1 = min ( n, m ) - 1
//
//  do jz = j0, j1
//     i0 = m + 1 - jz
//     do i = i0, ml
//        abd(i,jz) = 0.0
//     end do
//  end do
//
//  jz = j1
//  ju = 0
//!
//!  Gaussian elimination with partial pivoting.
//!
//  do k = 1, n-1
//!
//!  Zero next fill-in column.
//!
//     jz = jz + 1
//     if ( jz <= n ) then
//       abd(1:ml,jz) = 0.0
//     end if
//!
//!  Find L = pivot index.
//!
//     lm = min ( ml, n-k )
//     l = idamax ( lm+1, abd(m,k), 1 ) + m - 1
//     ipvt(k) = l + k - m
//!
//!  Zero pivot implies this column already triangularized.
//!
//     if ( abd(l,k) == 0.0 ) then
//
//       info = k
//!
//!  Interchange if necessary.
//!
//     else
//
//        if ( l /= m ) then
//           t = abd(l,k)
//           abd(l,k) = abd(m,k)
//           abd(m,k) = t
//        end if
//!
//!  Compute multipliers.
//!
//        t = -1.0 / abd(m,k)
//        call dscal ( lm, t, abd(m+1,k), 1 )
//!
//!  Row elimination with column indexing.
//!
//        ju = min ( max ( ju, mu+ipvt(k) ), n )
//        mm = m
//
//        do j = k+1, ju
//           l = l - 1
//           mm = mm - 1
//           t = abd(l,j)
//           if ( l /= mm ) then
//              abd(l,j) = abd(mm,j)
//              abd(mm,j) = t
//           end if
//           call daxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
//        end do
//
//     end if
//
//  end do
//
//  ipvt(n) = n
//
//  if ( abd(m,n) == 0.0 ) then
//    info = n
//  end if
//
//  return
//end
//subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )
//
//******************************************************
//!
//!! dgbsl() solves a banded system factored by DGBFA.
//!
//!  Discussion:
//!
//!    SGBSL can solve either a * x = b  or  trans(a) * x = b.
//!
//!  Parameters:
//!
//!     on entry
//!
//!        abd     real ( kind = rk8 )(lda, n)
//!                the output from dgbco or dgbfa.
//!
//!        lda     integer
//!                the leading dimension of the array  abd .
//!
//!        n       integer
//!                the order of the original matrix.
//!
//!        ml      integer
//!                number of diagonals below the main diagonal.
//!
//!        mu      integer
//!                number of diagonals above the main diagonal.
//!
//!        ipvt    integer(n)
//!                the pivot vector from dgbco or dgbfa.
//!
//!        b       real ( kind = rk8 )(n)
//!                the right hand side vector.
//!
//!        job     integer
//!                = 0         to solve  a*x = b ,
//!                = nonzero   to solve  trans(a)*x = b , where
//!                            trans(a)  is the transpose.
//!
//!     on return
//!
//!        b       the solution vector  x .
//!
//!     error condition
//!
//!        a division by zero will occur if the input factor contains a
//!        zero on the diagonal.  technically this indicates singularity
//!        but it is often caused by improper arguments or improper
//!        setting of lda .  it will not occur if the subroutines are
//!        called correctly and if dgbco has set 0.0 < RCOND
//!        or dgbfa has set info == 0 .
//!
//!     to compute  inverse(a) * c  where  c  is a matrix
//!     with  p  columns
//!           call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
//!           if (rcond is too small) go to ...
//!           do j = 1, p
//!              call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
//!           end do
//!
//!     linpack. this version dated 08/14/78 .
//!     cleve moler, university of new mexico, argonne national lab.
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer lda
//  integer n
//
//  real ( kind = rk8 ) abd(lda,n)
//  real ( kind = rk8 ) b(n)
//  integer ipvt(n)
//  integer job
//  integer k
//  integer l
//  integer la
//  integer lb
//  integer lm
//  integer m
//  integer ml
//  integer mu
//  real ( kind = rk8 ) ddot
//  real ( kind = rk8 ) t
//
//  m = mu + ml + 1
//!
//!  JOB = 0, Solve  a * x = b.
//!
//!  First solve l*y = b.
//!
//  if ( job == 0 ) then
//
//     if ( 0 < ml ) then
//
//        do k = 1, n-1
//           lm = min ( ml, n-k )
//           l = ipvt(k)
//           t = b(l)
//           if ( l /= k ) then
//              b(l) = b(k)
//              b(k) = t
//           end if
//           call daxpy ( lm, t, abd(m+1,k), 1, b(k+1), 1 )
//        end do
//
//     end if
//!
//!  Now solve u*x = y.
//!
//     do k = n, 1, -1
//        b(k) = b(k) / abd(m,k)
//        lm = min ( k, m ) - 1
//        la = m - lm
//        lb = k - lm
//        t = -b(k)
//        call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
//     end do
//!
//!  JOB nonzero, solve  trans(a) * x = b.
//!
//!  First solve  trans(u)*y = b.
//!
//  else
//
//     do k = 1, n
//        lm = min ( k, m ) - 1
//        la = m - lm
//        lb = k - lm
//        t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
//        b(k) = ( b(k) - t ) / abd(m,k)
//     end do
//!
//!  Now solve trans(l)*x = y
//!
//     if ( 0 < ml ) then
//
//        do k = n-1, 1, -1
//           lm = min ( ml, n-k )
//           b(k) = b(k) + ddot ( lm, abd(m+1,k), 1, b(k+1), 1 )
//           l = ipvt(k)
//           if ( l /= k ) then
//              t = b(l)
//              b(l) = b(k)
//              b(k) = t
//           end if
//        end do
//
//      end if
//
//  end if
//
//  return
//end
//subroutine dscal ( n, da, dx, incx )
//
//******************************************************
//!
//!! dscal() scales a vector by a constant.
//!
//!     uses unrolled loops for increment equal to one.
//!     jack dongarra, linpack, 3/11/78.
//!     modified 3/93 to return if incx  <=  0.
//!     modified 12/3/93, array(1) declarations changed to array(*)
//!
//!  Parameters:
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) da,dx(*)
//  integer i,incx,m,n,nincx
//!
//  if (  n <= 0 .or. incx <= 0 )return
//  if ( incx == 1)go to 20
//!
//!        code for increment not equal to 1
//!
//  nincx = n*incx
//  do i = 1,nincx,incx
//    dx(i) = da*dx(i)
//  end do
//  return
//!
//!        code for increment equal to 1
//!
//!
//!        clean-up loop
//!
//   20 continue
//
//  m = mod(n,5)
//  if (  m  ==  0 ) go to 40
//  dx(1:m) = da*dx(1:m)
//  if (  n  <  5 ) return
//
//40 continue
//
//  do i = m+1,n,5
//    dx(i) = da*dx(i)
//    dx(i + 1) = da*dx(i + 1)
//    dx(i + 2) = da*dx(i + 2)
//    dx(i + 3) = da*dx(i + 3)
//    dx(i + 4) = da*dx(i + 4)
//  end do
//
//  return
//end
//subroutine gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, &
//  nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
//
//******************************************************
//!
//!! gdump() writes information to a file.
//!
//!  Discussion:
//!
//!    The information can be used to create
//!    graphics images.  In order to keep things simple, exactly one
//!    value, real or integer, is written per record.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//!    Input, integer NPARA, the number of parameters.  Fixed at 1
//!    for now.
//!
//!    Input, real ( kind = rk8 ) PARA(MAXPAR), the parameters.
//!
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer nelemn
//  integer neqn
//  integer nnodes
//  integer np
//
//  real ( kind = rk8 ) f(neqn)
//  real ( kind = rk8 ) fval
//  integer i
//  integer indx(np,2)
//  integer insc(np)
//  integer iounit
//  integer, save :: iset = 0
//  integer isotri(nelemn)
//  integer j
//  logical long
//  integer node(nelemn,nnodes)
//  integer npara
//  integer nx
//  integer ny
//  real ( kind = rk8 ) para
//  real ( kind = rk8 ) reynld
//  real ( kind = rk8 ) rjpnew
//  real ( kind = rk8 ) ubdry
//  real ( kind = rk8 ) xc(np)
//  real ( kind = rk8 ) yc(np)
//
//  iset = iset+1
//
//  write(iounit,*)long
//  write (iounit,*) nelemn
//  write (iounit,*) np
//  write (iounit,*) npara
//  write (iounit,*) nx
//  write (iounit,*) ny
//!
//!  Pressures
//!
//  do i = 1, np
//    j = insc(i)
//    if (j <= 0) then
//      fval = 0.0
//    else
//      fval = f(j)
//    end if
//    write (iounit,*) fval
//  end do
//!
//!  Horizontal velocities, U
//!
//  do i = 1, np
//    j = indx(i,1)
//    if (j == 0) then
//      fval = 0.0
//    else if (j < 0) then
//      fval = ubdry(yc(i),para)
//    else
//      fval = f(j)
//    end if
//    write (iounit,*) fval
//  end do
//!
//!  Vertical velocities, V
//!
//  do i = 1, np
//    j = indx(i,2)
//    if (j <= 0) then
//      fval = 0.0
//    else
//      fval = f(j)
//    end if
//    write (iounit,*) fval
//  end do
//
//  do i = 1, np
//    write (iounit,*) indx(i,1)
//    write (iounit,*) indx(i,2)
//  end do
//
//  do i = 1, np
//    write (iounit,*) insc(i)
//  end do
//
//  do i = 1, nelemn
//    write (iounit,*) isotri(i)
//  end do
//
//  do i = 1, nelemn
//    do j = 1, 6
//      write (iounit,*) node(i,j)
//    end do
//  end do
//
//  write (iounit,*) para
//  write (iounit,*) reynld
//  write (iounit,*) rjpnew
//
//  do i = 1, np
//    write (iounit,*) xc(i)
//  end do
//
//  do i = 1, np
//    write (iounit,*) yc(i)
//  end do
//
//  std::println("GDUMP wrote data set ',iset,' to file.")
//
//  return
//end
//subroutine getg ( f, iline, my, neqn, u )
//
//******************************************************
//!
//!! getg() outputs field values along the profile line X = XZERO.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer neqn
//  integer my
//
//  real ( kind = rk8 ) f(neqn)
//  integer iline(my)
//  integer j
//  integer k
//  real ( kind = rk8 ) u(my)
//
//  do j = 1, my
//    k = iline(j)
//    if ( 0 < k ) then
//      u(j) = f(k)
//    else
//      u(j) = 0.0
//    end if
//  end do
//
//  return
//end
//subroutine gram ( gr, iline, indx, iwrite, my, nelemn, nnodes, node, &
//  nodex0, np, para, r, ui, xc, yc )
//
//******************************************************
//!
//!! gram() computes the Gram matrix, GR(I,J) = INTEGRAL PHI(I)*PHI(J).
//!
//!  and the vector R(I) = INTEGRAL UI*PHI(I).
//!
//!  The integrals are computed along the line where the profile is
//!  specified.  The three point Gauss quadrature rule is used for the
//!  line integral.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer my
//  integer nelemn
//  integer nnodes
//  integer np
//
//  real ( kind = rk8 ) ar
//  real ( kind = rk8 ) bb
//  real ( kind = rk8 ) bbb
//  real ( kind = rk8 ) bbx
//  real ( kind = rk8 ) bby
//  real ( kind = rk8 ) bma2
//  real ( kind = rk8 ) bx
//  real ( kind = rk8 ) by
//  real ( kind = rk8 ) gr(my,my)
//  integer i
//  integer igetl
//  integer ii
//  integer iline(my)
//  integer indx(np,2)
//  integer ip
//  integer ipp
//  integer iq
//  integer iqq
//  integer iquad
//  integer it
//  integer iun
//  integer iwrite
//  integer j
//  integer jj
//  integer k
//  integer kk
//  integer node(nelemn,nnodes)
//  integer nodex0
//  real ( kind = rk8 ) para
//  real ( kind = rk8 ) r(my)
//  real ( kind = rk8 ) ubc
//  real ( kind = rk8 ) ubdry
//  real ( kind = rk8 ) ui(my)
//  real ( kind = rk8 ) uiqdpt
//  real ( kind = rk8 ) wt(3)
//  real ( kind = rk8 ) x
//  real ( kind = rk8 ) xc(np)
//  real ( kind = rk8 ) xzero
//  real ( kind = rk8 ) y
//  real ( kind = rk8 ) yq(3)
//  real ( kind = rk8 ) yc(np)
//!
//!  Input values for 3 point Gauss quadrature
//!
//  wt(1) = 5.0 / 9.0
//  wt(2) = 8.0 / 9.0
//  wt(3) = wt(1)
//  yq(1) = -0.7745966692
//  yq(2) = 0.0
//  yq(3) = -yq(1)
//!
//!  zero arrays
//!
//  r(1:my) = 0.0
//  gr(1:my,1:my) = 0.0
//!
//!  Compute line integral by looping over intervals along line
//!  using three point Gauss quadrature
//!
//  xzero = xc(nodex0)
//  do 70 it = 1, nelemn
//!
//!  Check to see if we are in a triangle with a side along line
//!  x = xzero.  If not, skip out
//!
//    k = node(it,1)
//    kk = node(it,2)
//
//    if ( 1.0D-04 < abs(xc(k)-xzero) ) then
//      cycle
//    end if
//
//    if ( 1.0D-04 < abs(xc(kk)-xzero) ) then
//      cycle
//    end if
//
//    do 60 iquad = 1, 3
//      bma2 = (yc(kk)-yc(k))/2.0
//      ar = bma2*wt(iquad)
//      x = xzero
//      y = yc(k)+bma2*(yq(iquad)+1.0 )
//!
//!  Compute u internal at quadrature points
//!
//      uiqdpt = 0
//      do 30 iq = 1, nnodes
//        if ( 4 < iq ) go to 30
//        if (iq == 3) go to 30
//        call qbf (x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
//        ip = node(it,iq)
//        iun = indx(ip,1)
//        if ( 0 < iun ) then
//          ii = igetl(iun,iline,my)
//          uiqdpt = uiqdpt+bb*ui(ii)
//        else if (iun < 0) then
//          ubc = ubdry(yc(ip),para)
//          uiqdpt = uiqdpt+bb*ubc
//        end if
//   30     continue
//!
//!  Only loop over nodes lying on line x = xzero
//!
//      do 50 iq = 1, nnodes
//        if ( iq == 1.or.iq == 2.or.iq == 4 ) then
//          ip = node(it,iq)
//          call qbf (x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
//          i = indx(ip,1)
//          if (i <= 0) go to 50
//          ii = igetl(i,iline,my)
//          r(ii) = r(ii)+bb*uiqdpt*ar
//
//          do iqq = 1, nnodes
//            if ( iqq == 1.or.iqq == 2.or.iqq == 4 ) then
//              ipp = node(it,iqq)
//              call qbf (x,y,it,iqq,bbb,bbx,bby,nelemn,nnodes, &
//                node,np,xc,yc)
//              j = indx(ipp,1)
//              if (j /= 0) then
//                jj = igetl(j,iline,my)
//                gr(ii,jj) = gr(ii,jj)+bb*bbb*ar
//              end if
//            end if
//          end do
//
//        end if
//   50     continue
//   60   continue
//   70 continue
//
//  if ( 2 <= iwrite ) then
//    std::println(" ")
//    std::println("Gram matrix:")
//    std::println(" ")
//    do i = 1,my
//      do j = 1,my
//        write(*,*)i,j,gr(i,j)
//      end do
//    end do
//    std::println(" ")
//    std::println("R vector:")
//    std::println(" ")
//    do i = 1,my
//      write(*,*)i,r(i)
//    end do
//  end if
//
//  return
//end

//******************************************************
int idamax(int n, std::span<const double> dx, int incx) {
//? idamax() finds the index of element having max. absolute value.

  // Si el tamaño es menor que 1 o el incremento no es positivo,
  // retorna -1 (equivalente a 0 en Fortran para indicar "no encontrado")
  if (n < 1 || incx <= 0) return -1;

  int idamax_val = 0; // Guardará el índice basado en 0
  if (n == 1) return idamax_val;

  if (incx == 1) {
    // Código para incremento igual a 1
    double dmax = std::abs(dx[0]);

    for (int i = 1; i < n; ++i) {
      if (dmax < std::abs(dx[i])) {
        idamax_val = i;
        dmax = std::abs(dx[i]);
      }
    }
  } else {
    // Código para incremento diferente de 1
    int ix = 0; // Índice basado en 0
    double dmax = std::abs(dx[ix]);
    ix += incx;

    for (int i = 1; i < n; ++i) {
      if (std::abs(dx[ix]) > dmax) {
        idamax_val = i;
        dmax = std::abs(dx[ix]);
      }
      ix += incx;
    }
  }

  return idamax_val;
}

//function igetl ( k, iline, my )
//
//******************************************************
//!
//!! igetl() gets the local unknown number along the profile line.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer my
//
//  integer igetl
//  integer iline(my)
//  integer j
//  integer k
//
//  do j = 1, my
//    if (iline(j) == k) then
//      igetl = j
//      return
//    end if
//  end do
//
//  write ( *, * ) ' ")
//  std::println("IGETL - fatal error!")
//  std::println("  Unable to get local unknown number for ")
//  std::println("  Global variable number ',k
//  igetl = 0
//  stop
//end

//******************************************************
void linsys(
  std::array<std::array<double, maxeqn>, maxrow>& a,
  std::array<double, nelemn>& area,
  std::array<double, maxeqn>& f,
  std::array<double, maxeqn>& g,
  std::array<std::array<int, 2>, np>& indx,
  std::array<int, np>& insc,
  std::array<int, maxeqn>& ipivot,
  int neqn,
  int nlband,
  std::array<std::array<int, nnodes>, nelemn>& node,
  int nrow,
  double para1,
  double para2,
  std::array<
    std::array<
      std::array<std::array<double, 3>, nnodes>,
      nquad
    >, nelemn
  >& phi,
  std::array<
    std::array<
      std::array<double, nnodes>, nquad
    >, nelemn
  >& psi,
  double reynld,
  std::array<double, np>& yc
) {
  //? linsys() sets up and solves the linear system.
  //?    The G array contains the previous solution.
  //?
  //?    The F array contains the right hand side initially and then the
  //?    current solution.

  int i;
  int info;
  int ip;
  double bb;
  double bx;
  double by;
  double bbl;
  int ihor;
  int iver;
  int iprs;
  int ipp;
  double bbb;
  double bbx;
  double bby;
  double bbbl;
  int ju;
  int jv;
  int jp;
  int iuse;
  std::array<double, 2> un;
  std::array<double, 2> unx;
  std::array<double, 2> uny;
  double ar;
  double uu;
  int ioff = nlband + nlband;
  double visc = 1.0 / reynld;
  for (int j = 0; j < neqn; j++) {
    f[j] = 0.0;
  }
  for (int k = 0; k < nrow; k++) {
    for (int j = 0; j < neqn; j++) {
      a[k][j] = 0.0;
    }
  }
  //!  For each element,
  for (int it = 0; it < nelemn; ++it) {

    ar = area[it] / 3.0;
    //!  and for each quadrature point in the element,
    for (int iquad = 0; iquad < nquad; ++iquad) {
      //!  Evaluate velocities at quadrature point
      uval(g, indx, iquad, it, neqn, node, para1, phi, un, unx, uny, yc);
      //!  For each basis function,
      for (int iq = 0; iq < nnodes; ++iq) {
        ip = node[it][iq]-1;
        bb = phi[it][iquad][iq][0];
        bx = phi[it][iquad][iq][1];
        by = phi[it][iquad][iq][2];
        bbl = psi[it][iquad][iq];
        ihor = indx[ip-1][0]-1;
        iver = indx[ip-1][1]-1;
        iprs = insc[ip-1]-1;

        if ( 0 < ihor+1 ) {
          f[ihor] += ar*bb*(un[0]*unx[0]+un[1]*uny[0]);
        }

        if ( 0 < iver+1 ) {
          f[iver] += ar*bb*(un[0]*unx[1]+un[1]*uny[1]);
        }
        //!  For another basis function,
        //do iqq = 1, nnodes
        for (int iqq = 0; iqq < nnodes; ++iqq) {
          ipp = node[it][iqq]-1;
          bbb = phi[it][iquad][iqq][0];
          bbx = phi[it][iquad][iqq][1];
          bby = phi[it][iquad][iqq][2];
          bbbl = psi[it][iquad][iqq];
          ju = indx[ipp-1][0]-1;
          jv = indx[ipp-1][1]-1;
          jp = insc[ipp-1]-1;
          //!  Horizontal velocity variable
          if ( 0 < ju+1 ) {
            if ( 0 < ihor+1 ) {
              iuse = ihor-ju+ioff;
              a[iuse][ju] += ar*(
                visc*(by*bby+bx*bbx)
                + bb*(bbb*unx[0]+bbx*un[0]+bby*un[1])
              );
            }

            if ( 0 < iver+1 ) {
              iuse = iver-ju+ioff;
              a[iuse][ju] += ar*bb*bbb*unx[1];
            }

            if ( 0 < iprs+1 ) {
              iuse = iprs-ju+ioff;
              a[iuse][ju] += ar*bbx*bbl;
            }

          } else if ( ju+1 < 0 ) {
            uu = ubdry(yc[ipp],para2);
            if ( 0 < ihor+1 ) {
              f[ihor] -= ar*uu*(
                visc*(by*bby+bx*bbx) +
                bb*(bbb*unx[0]+bbx*un[0]+bby*un[1])
              );
            }

            if ( 0 < iver+1 ) {
              f[iver] -= ar*uu*bb*bbb*unx[1];
            }

            if ( 0 < iprs+1 ) {
              f[iprs] -= ar*uu*bbx*bbl;
            }
          }
          //!  Vertical velocity variable
          if ( 0 < jv+1 ){

            if ( 0 < ihor+1 ) {
              iuse = ihor-jv+ioff;
              a[iuse][jv] += ar*bb*bbb*uny[0];
            }

            if ( 0 < iver+1 ) {
              iuse = iver-jv+ioff;
              a[iuse][jv] += ar*(
                visc*(by*bby+bx*bbx) +
                bb*(bbb*uny[1]+bby*un[1]+bbx*un[0])
              );
            }

            if ( 0 < iprs+1 ) {
              iuse = iprs-jv+ioff;
              a[iuse][jv] += ar*bby*bbl;
            }
          }
          //!  Pressure variable
          if ( 0 < jp+1 ) {
            if ( 0 < ihor+1 ) {
              iuse = ihor-jp+ioff;
              a[iuse][jp] -= ar*bx*bbbl;
            }
            if ( 0 < iver+1 ) {
              iuse = iver-jp+ioff;
              a[iuse][jp] -= ar*by*bbbl;
            }
          }
        }
      }
    }
  }
  //!  To avoid singularity of the pressure system, the last pressure
  //!  is simply assigned a value of 0.
  f[neqn-1] = 0.0;
  //do j = neqn-nlband, neqn-1
  for (int j = neqn-nlband; j < neqn; ++j) {
    i = neqn-j+ioff;
    a[i][j] = 0.0;
  }
  a[ioff][neqn-1] = 1.0;
  //!  Factor the matrix
  //call dgbfa ( a, maxrow, neqn, nlband, nlband, ipivot, info )

  if ( info != 0 ) {
    std::println(" ");
    std::println("LINSYS - fatal error!");
    std::println("DGBFA returns INFO = {}", info);
    std::exit(1);
  }
  //!  Solve the linear system
  int job = 0;
//  call dgbsl ( a, maxrow, neqn, nlband, nlband, ipivot, f, job )
}

//subroutine nstoke (
//  a,area,f,g,indx,insc,ipivot,iwrite,maxnew,maxrow, &
//  nelemn,neqn,nlband,nnodes,node,np,nquad,nrow,&
//  numnew,para,phi,psi,reynld,tolnew,yc)
//
//******************************************************
void nstoke(
  std::array<std::array<double, maxeqn>, maxrow> &a,
  std::array<double, nelemn> &area,
  std::array<double, maxeqn> &f,
  std::array<double, maxeqn> &g,
  std::array<std::array<int, 2>, np> &indx,
  std::array<int, np> &insc,
  std::array<int, maxeqn> &ipivot,
  int iwrite,
  int &maxnew,
  int neqn, int nlband,
  std::array<std::array<int, nnodes>, nelemn> &node,
  int &nrow,
  int &numnew,
  double &para,
  std::array<std::array<std::array<std::array<double, 3>, nnodes>, nquad>, nelemn> &phi,
  std::array<std::array<std::array<double, nnodes>, nquad>, nelemn> &psi,
  double &reynld,
  double &tolnew,
  std::array<double, np>& yc
) {
  //!! nstoke() solves the Navier Stokes equation using Taylor-Hood elements.
  //!
  //!  The G array contains the previous iterate.
  //!
  //!  The F array contains the right hand side initially and then the current iterate.
  for (int iter = 1; iter <= maxnew; ++iter) {
    numnew ++;

    linsys(
      a, area, f, g, indx, insc, ipivot,
      neqn, nlband, node,
      nrow, para, para, phi, psi, reynld, yc
    );
    //! Check for convergence
    for (int i = 0; i < neqn; ++i) {
      g[i] = g[i] - f[i];
    }
    double diff = std::abs( g[idamax(neqn,g,1)] );

    if ( 1 <= iwrite ) {
      std::println("NSTOKE iteration {} Mnorm = {}", iter, diff);
    }

    for (int i = 0; i < neqn; ++i) {
      g[i] = f[i];
    }

    if ( diff <= tolnew ) {
      std::println("Navier Stokes iteration converged in {} iterations.", iter);
      return;
    }

  }

  std::println("Navier Stokes solution did not converge!");
}

//subroutine pval (g,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)
//
//******************************************************
//!
//!! pval() computes a table of pressures.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer nelemn
//  integer neqn
//  integer nnodes
//  integer np
//
//  real ( kind = rk8 ) g(neqn)
//  integer i
//  integer insc(np)
//  integer ip
//  integer iq
//  integer it
//  integer ivar
//  integer j
//  logical long
//  integer mx
//  integer my
//  integer node(nelemn,nnodes)
//  real ( kind = rk8 ) press(mx,my)
//
//  press(1:mx,1:my) = 0.0
//!
//!  Read the pressures where they are computed.
//!  These are "(odd, odd)" points.
//!
//  do it = 1, nelemn
//    do iq = 1, 3
//
//      ip = node(it,iq)
//      ivar = insc(ip)
//
//      if ( long ) then
//        i = ((ip-1)/my)+1
//        j = mod(ip-1,my)+1
//      else
//        i = mod(ip-1,mx)+1
//        j = ((ip-1)/mx)+1
//      end if
//
//      if ( 0 < ivar ) then
//        press(i,j) = g(ivar)
//      else
//        press(i,j) = 0.0
//      end if
//
//    end do
//  end do
//!
//!  Interpolate the pressures at points (even, odd) and (odd, even).
//!
//  do i = 2,mx-1,2
//    do j = 1,my,2
//      press(i,j) = 0.5*(press(i-1,j)+press(i+1,j))
//    end do
//  end do
//
//  do j = 2,my-1,2
//    do i = 1,mx,2
//      press(i,j) = 0.5*(press(i,j-1)+press(i,j+1))
//    end do
//  end do
//!
//!  Interpolate the pressures at points (even,even).
//!
//  do j = 2,my-1,2
//    do i = 2,mx-1,2
//      press(i,j) = 0.5*(press(i-1,j-1)+press(i+1,j+1))
//    end do
//  end do
//
//  return
//end

//******************************************************
void qbf(
  double x,
  double y,
  int it,
  int in,
  double &bb,
  double &bx,
  double &by,
  const std::array<std::array<int, nnodes>, nelemn> &node,
  const std::array<double, np> &xc,
  const std::array<double, np> &yc
) {
  //? qbf() evaluates a quadratic basis function in a triangle.
  int in1;
  int in2;
  int in3;
  int i1;
  int i2;
  int i3;
  double d;
  double t;
  double c;
  double s;
  int inn;
  int j1;
  int j2;
  int j3;

  if (in < 3) {
    in1 = in;
    in2 = (in%3);
    in3 = (in+1)%3;
    i1 = node[it][in1];
    i2 = node[it][in2];
    i3 = node[it][in3];
    d = (xc[i2] - xc[i1])*(yc[i3] - yc[i1]) -
      (xc[i3] - xc[i1])*(yc[i2] - yc[i1]);
    t = 1.0+(
      (yc[i2]-yc[i3])*(x-xc[i1])+(xc[i3]-xc[i2])*(y-yc[i1])
    )/d;
    bb = t*(2.0*t-1.0);
    bx = (yc[i2]-yc[i3])*(4.0*t-1.0)/d;
    by = (xc[i3]-xc[i2])*(4.0*t-1.0)/d;
  } else {
    inn = in-3;
    in1 = inn;
    in2 = inn % 3;
    in3 = (inn+1) % 3;
    i1 = node[it][in1];
    i2 = node[it][in2];
    i3 = node[it][in3];
    j1 = i2;
    j2 = i3;
    j3 = i1;
    d = (xc[i2]-xc[i1])*(yc[i3]-yc[i1]) -
      (xc[i3]-xc[i1])*(yc[i2]-yc[i1]);
    c = (xc[j2]-xc[j1])*(yc[j3]-yc[j1]) -
      (xc[j3]-xc[j1])*(yc[j2]-yc[j1]);
    t = 1.0+(
      (yc[i2] - yc[i3])*(x - xc[i1]) +
      (xc[i3] - xc[i2])*(y - yc[i1])
    )/d;
    s = 1.0+(
      (yc[j2] - yc[j3])*(x - xc[j1]) +
      (xc[j3] - xc[j2])*(y - yc[j1])
    )/c;
    bb = 4.0*s*t;
    bx = 4.0*(t*(yc[j2] - yc[j3])/c +
      s*(yc[i2] - yc[i3])/d);
    by = 4.0*(t*(xc[j3] - xc[j2])/c +
      s*(xc[i3] - xc[i2])/d);
  }
}

//subroutine resid (area,g,indx,insc,iwrite,nelemn,neqn,nnodes, &
//  node,np,nquad,para,phi,psi,res,reynld,yc)
//
//******************************************************
//!
//!! resid() computes the residual.
//!
//!  Discussion:
//!
//!    The G array contains the current iterate.
//!
//!    The RES array will contain the value of the residual.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer nelemn
//  integer neqn
//  integer nnodes
//  integer np
//  integer nquad
//
//  real ( kind = rk8 ) aijpu
//  real ( kind = rk8 ) aijpv
//  real ( kind = rk8 ) aijup
//  real ( kind = rk8 ) aijuu
//  real ( kind = rk8 ) aijuv
//  real ( kind = rk8 ) aijvp
//  real ( kind = rk8 ) aijvu
//  real ( kind = rk8 ) aijvv
//  real ( kind = rk8 ) ar
//  real ( kind = rk8 ) area(nelemn)
//  real ( kind = rk8 ) bb
//  real ( kind = rk8 ) bbb
//  real ( kind = rk8 ) bbbl
//  real ( kind = rk8 ) bbl
//  real ( kind = rk8 ) bbx
//  real ( kind = rk8 ) bby
//  real ( kind = rk8 ) bx
//  real ( kind = rk8 ) by
//  real ( kind = rk8 ) g(neqn)
//  integer i
//  integer ibad
//  integer ihor
//  integer imax
//  integer indx(np,2)
//  integer insc(np)
//  integer ip
//  integer ipp
//  integer iprs
//  integer iq
//  integer iqq
//  integer iquad
//  integer it
//  integer iver
//  integer iwrite
//  integer j
//  integer jp
//  integer ju
//  integer jv
//  integer node(nelemn,nnodes)
//  real ( kind = rk8 ) para
//  real ( kind = rk8 ) phi(nelemn,nquad,nnodes,3)
//  real ( kind = rk8 ) psi(nelemn,nquad,nnodes)
//  real ( kind = rk8 ) res(neqn)
//  real ( kind = rk8 ) reynld
//  real ( kind = rk8 ) rmax
//  real ( kind = rk8 ) test
//  real ( kind = rk8 ) ubdry
//  real ( kind = rk8 ) un(2)
//  real ( kind = rk8 ) unx(2)
//  real ( kind = rk8 ) uny(2)
//  real ( kind = rk8 ) uu
//  real ( kind = rk8 ) visc
//  real ( kind = rk8 ) yc(np)
//
//  visc = 1.0 / reynld
//
//  res(1:neqn) = 0.0
//!
//!  For each element,
//!
//  do 90 it = 1, nelemn
//
//    ar = area(it) / 3.0
//!
//!  and for each quadrature point in the element,
//!
//    do 80 iquad = 1, nquad
//!
//!  Evaluate velocities at quadrature point
//!
//      call uval (g,indx,iquad,it,nelemn,neqn,nnodes,node, &
//        np,nquad,para,phi,un,unx,uny,yc)
//!
//!  For each basis function,
//!
//      do 70 iq = 1, nnodes
//        ip = node(it,iq)
//        bb = phi(it,iquad,iq,1)
//        bx = phi(it,iquad,iq,2)
//        by = phi(it,iquad,iq,3)
//        bbl = psi(it,iquad,iq)
//        iprs = insc(ip)
//        ihor = indx(ip,1)
//        iver = indx(ip,2)
//
//        if ( 0 < ihor ) then
//          res(ihor) = res(ihor)+(un(1)*unx(1)+un(2)*uny(1))*bb*ar
//        end if
//
//        if ( 0 < iver ) then
//          res(iver) = res(iver)+(un(1)*unx(2)+un(2)*uny(2))*bb*ar
//        end if
//!
//!  For another basis function,
//!
//        do iqq = 1, nnodes
//
//          ipp = node(it,iqq)
//          bbb = phi(it,iquad,iqq,1)
//          bbx = phi(it,iquad,iqq,2)
//          bby = phi(it,iquad,iqq,3)
//          bbbl = psi(it,iquad,iqq)
//          ju = indx(ipp,1)
//          jv = indx(ipp,2)
//          jp = insc(ipp)
//
//          if ( 0 < ju ) then
//            if ( 0 < ihor ) then
//              aijuu = visc*(by*bby+bx*bbx) &
//                + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2))
//              res(ihor) = res(ihor)+aijuu*ar*g(ju)
//            end if
//            if ( 0 < iver ) then
//              aijvu = bb*bbb*unx(2)
//              res(iver) = res(iver)+aijvu*ar*g(ju)
//            end if
//            if ( 0 < iprs ) then
//              aijpu = bbx*bbl
//              res(iprs) = res(iprs)+aijpu*ar*g(ju)
//            end if
//          else if ( ju < 0 ) then
//            uu = ubdry(yc(ipp),para)
//            if ( 0 < ihor ) then
//              aijuu = visc*(by*bby+bx*bbx) &
//                + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2))
//              res(ihor) = res(ihor)+ar*aijuu*uu
//            end if
//            if ( 0 < iver ) then
//              aijvu = bb*bbb*unx(2)
//              res(iver) = res(iver)+ar*aijvu*uu
//            end if
//            if ( 0 < iprs ) then
//              aijpu = bbx*bbl
//              res(iprs) = res(iprs)+ar*aijpu*uu
//            end if
//          end if
//
//          if ( 0 < jv ) then
//            if ( 0 < ihor ) then
//              aijuv = bb*bbb*uny(1)
//              res(ihor) = res(ihor)+aijuv*ar*g(jv)
//            end if
//            if ( 0 < iver ) then
//              aijvv = visc*(by*bby+bx*bbx) &
//                +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1))
//              res(iver) = res(iver)+aijvv*ar*g(jv)
//            end if
//            if ( 0  < iprs ) then
//              aijpv = bby*bbl
//              res(iprs) = res(iprs)+aijpv*ar*g(jv)
//            end if
//          end if
//
//          if ( 0 < jp ) then
//            if ( 0 < ihor ) then
//              aijup = -bx*bbbl
//              res(ihor) = res(ihor)+aijup*ar*g(jp)
//            end if
//            if ( 0 < iver ) then
//              aijvp = -by*bbbl
//              res(iver) = res(iver)+aijvp*ar*g(jp)
//            end if
//          end if
//
//        end do
//
//   70     continue
//   80   continue
//   90 continue
//
//  res(neqn) = g(neqn)
//
//  rmax = 0.0
//  imax = 0
//  ibad = 0
//
//  do i = 1,neqn
//
//    test = abs(res(i))
//
//    if ( rmax < test ) then
//      rmax = test
//      imax = i
//    end if
//
//    if ( 1.0D-03 < test ) then
//      ibad = ibad+1
//    end if
//
//  end do
//
//  if ( 1 <= iwrite ) then
//    std::println(" ")
//    std::println("RESIDUAL INFORMATION:")
//    std::println(" ")
//    std::println("Worst residual is number ',IMAX
//    std::println("of magnitude ',RMAX
//    std::println(" ")
//    std::println("Number of "bad" residuals is ',IBAD,' out of ',NEQN
//    std::println(" ")
//  end if
//
//  if ( 2 <= iwrite ) then
//    std::println("Raw residuals:")
//    std::println(" ")
//    i = 0
//
//    do j = 1,np
//
//      if ( 0 < indx(j,1) ) then
//        i = i+1
//        if ( abs(res(i)) <= 1.0D-03 ) then
//          write(*,'(1x,a1,2i5,g14.6)')'U',i,j,res(i)
//        else
//          write(*,'(a1,a1,2i5,g14.6)')'*','U',i,j,res(i)
//        end if
//      end if
//
//      if ( 0 < indx(j,2) ) then
//        i = i+1
//        if ( abs(res(i)) <= 1.0D-03 ) then
//          write(*,'(1x,a1,2i5,g14.6)')'V',i,j,res(i)
//        else
//          write(*,'(a1,a1,2i5,g14.6)')'*','V',i,j,res(i)
//        end if
//      end if
//
//      if ( 0 < insc(j) ) then
//        i = i+1
//        if ( abs(res(i)) <= 1.0D-03 ) then
//          write(*,'(1x,a1,2i5,g14.6)')'P',i,j,res(i)
//        else
//          write(*,'(a1,a1,2i5,g14.6)')'*','P',i,j,res(i)
//        end if
//      end if
//
//    end do
//
//  end if
//
//  return
//end

//******************************************************
void setban(
  std::array<std::array<int, 2>, np> &indx,
  std::array<int, np> &insc,
  int &nband,
  int &nlband,
  std::array<std::array<int, nnodes>, nelemn> &node,
  int &nrow
){
//? setban() computes the half band width.
  nlband = 0;
  int i, j;

  //do it = 1, nelemn
  for (int it = 1; it <= nelemn; ++it) {
    //do iq = 1, nnodes
    for (int iq = 1; iq <= nnodes; ++iq) {
      int ip = node[it-1][iq-1];
      for (int iuk = 1; iuk <= 3; ++iuk) {
        if (iuk == 3) {
          i = insc[ip-1];
        } else {
          i = indx[ip-1][iuk-1];
        }

        if ( 0 < i ){
          //do iqq = 1, nnodes
          for(int iqq = 1; iqq <= nnodes; ++iqq){
            int ipp = node[it-1][iqq-1];
            //do iukk = 1, 3
            for(int iukk = 1; iukk <= 3; ++iukk){
              if (iukk == 3) {
                j = insc[ipp-1];
              } else {
                j = indx[ipp-1][iukk-1];
              }
              nlband = std::max(nlband,j-i);
            }
          }
        }
      }
    }
  }

  nband = nlband + nlband + 1;
  nrow = nlband + nlband + nlband + 1;

  std::println("Lower bandwidth = {}", nlband);
  std::println("Total bandwidth = {}", nband);
  std::println("NROW  = {}", nrow);
  if (maxrow < nrow) {
    std::println("SETBAN - NROW is too large!");
    std::println("The maximum allowed is {}", maxrow);
  }
}

//******************************************************
void setbas(
  std::array<std::array<int, nnodes>, nelemn> &node,
  std::array<double, np> &xc,
  std::array<double, np> &yc,
  std::array<std::array<std::array<std::array<double, 3>, nnodes>, nquad>, nelemn> &phi,
  std::array<std::array<std::array<double, nnodes>, nquad>, nelemn> &psi,
  std::array<std::array<double, nquad>, nelemn> &xm,
  std::array<std::array<double, nquad>, nelemn> &ym
) {
//? setbas() computes the basis functions at each integration point.
//  integer nelemn
//  integer nnodes
//  integer np
//  integer nquad
//
  double bb;
  double bx;
  double by;

  double x, y;
  //do it = 1,nelemn
  for (int it = 0; it < nelemn; ++it) {
    //do j = 1,nquad
    for (int j = 0; j < nquad; ++j) {
      x = xm[it][j];
      y = ym[it][j];
      //do iq = 1,6
      for (int iq = 0; iq < 6; ++iq) {
        psi[it][j][iq] = bsp(x,y,it,iq,1,node,xc,yc);
        qbf(x,y,it,iq,bb,bx,by,node,xc,yc);
        phi[it][j][iq][0] = bb;
        phi[it][j][iq][1] = bx;
        phi[it][j][iq][2] = by;
      }
    }
  }
}

//******************************************************
void setgrd(
//******************************************************
    std::array<std::array<int, 2>, np> &indx,
    std::array<int, np> &insc,
    std::array<int, nelemn> &isotri,
    int &iwrite, bool &_long, int &neqn,
    std::array<std::array<int, nnodes>, nelemn> &node
){
//? setgrd() sets up the grid for the problem..

//!  Determine whether region is long or skinny.
//!  This will determine how we number the nodes and elements.

  if (ny < nx) {
    _long = true;
    std::println("Using vertical ordering.");
  } else {
    _long = false;
    std::println("Using horizontal ordering.");
  }
//!
//!  Set parameters for Taylor Hood element
//!
  std::println(" ");
  std::println("SETGRD: Taylor Hood element");
//!
//!  Construct grid coordinates, elements, and ordering of unknowns
//!
  neqn = 0;
  int ielemn = 0;

  for (int ip = 1; ip <= np; ++ip) {
    int ic, jc;

    if (_long) {
      ic = ((ip - 1) / my) + 1;
      jc = ((ip - 1) % my) + 1;
    } else {
      ic = ((ip - 1) % mx) + 1;
      jc = ((ip - 1) / mx) + 1;
    }

    int icnt = ic % 2;   // igual que mod(ic,2)
    int jcnt = jc % 2;   // igual que mod(jc,2)
//!
//!  If both the row count and the column count are odd,
//!  and we're not in the last row or top column,
//!  then we can define two new triangular elements based at the node.
//!
//!  For horizontal ordering,
//!  given the following arrangement of nodes, for instance:
//!
//!    21 22 23 24 25
//!    16 17 18 19 20
//!    11 12 13 14 15
//!    06 07 08 09 10
//!    01 02 03 04 05
//!
//!  when we arrive at node 13, we will define
//!
//!  element 7: (13, 23, 25, 18, 24, 19)
//!  element 8: (13, 25, 15, 19, 20, 14)
//!
//!
//!  For vertical ordering,
//!  given the following arrangement of nodes, for instance:
//!
//!    05 10 15 20 25
//!    04 09 14 19 24
//!    03 08 13 18 23
//!    02 07 12 17 22
//!    01 06 11 16 21
//!
//!  when we arrive at node 13, we will define
//!
//!  element 7: (13, 25, 23, 19, 24, 18)
//!  element 8: (13, 15, 25, 14, 20, 19)
//!
    if (icnt == 1 && jcnt == 1 && ic != mx && jc != my) {

      if (_long) {
        int ip1 = ip + my;
        int ip2 = ip + my + my;   // o 2 * my

        // Primer elemento
        ++ielemn;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2 + 2;
        node[ielemn - 1][2] = ip2;
        node[ielemn - 1][3] = ip1 + 1;
        node[ielemn - 1][4] = ip2 + 1;
        node[ielemn - 1][5] = ip1;
        isotri[ielemn - 1] = 0;

        // Segundo elemento
        ++ielemn;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip + 2;
        node[ielemn - 1][2] = ip2 + 2;
        node[ielemn - 1][3] = ip + 1;
        node[ielemn - 1][4] = ip1 + 2;
        node[ielemn - 1][5] = ip1 + 1;
        isotri[ielemn - 1] = 0;
      } else {   // .NOT. long
        int ip1 = ip + mx;
        int ip2 = ip + mx + mx;   // o 2 * mx

        // Primer elemento
        ++ielemn;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2;
        node[ielemn - 1][2] = ip2 + 2;
        node[ielemn - 1][3] = ip1;
        node[ielemn - 1][4] = ip2 + 1;
        node[ielemn - 1][5] = ip1 + 1;
        isotri[ielemn - 1] = 0;

        // Segundo elemento
        ++ielemn;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2 + 2;
        node[ielemn - 1][2] = ip + 2;
        node[ielemn - 1][3] = ip1 + 1;
        node[ielemn - 1][4] = ip1 + 2;
        node[ielemn - 1][5] = ip + 1;
        isotri[ielemn - 1] = 0;
      }
    }

//!
//!  Consider whether velocity unknowns should be associated with this node.
//!
//!  If we are in column 1, horizontal velocities are specified, and
//!  vertical velocities are zero.
//!
    if (ic == 1 && jc > 1 && jc < my) {
      indx[ip - 1][0] = -1;
      indx[ip - 1][1] = 0;
    }
//!
//!  If we are in column MX, horizontal velocities are unknown, and
//!  vertical velocities are zero.
//!
    else if (ic == mx && jc > 1 && jc < my) {
      ++neqn;
      indx[ip - 1][0] = neqn;
      indx[ip - 1][1] = 0;
    }
//!
//!  Otherwise, if we are in row 1 or row MY, both horizontal and
//!  vertical velocities are zero.
//!
    else if (jc == 1 || jc == my) {
      indx[ip - 1][0] = 0;
      indx[ip - 1][1] = 0;
    }
//!
//!  Otherwise, we are at an interior node
//!
    else {
      neqn += 2;
      indx[ip - 1][0] = neqn - 1;
      indx[ip - 1][1] = neqn;
    }
//!
//!  Consider whether a pressure unknown should be associated with this node.
//!  The answer is yes if both nodes are odd.
//!
    if (jcnt == 1 && icnt == 1) {
      ++neqn;
      insc[ip - 1] = neqn;
    } else {
      insc[ip - 1] = 0;
    }
  }
//!
//!  If debugging is requested, print out data.
//!
  if ( 2 <= iwrite ) {
    std::println(" ");
    std::println("    I      INDX 1 & 2, INSC");
    std::println(" ");
    //do i = 1,np
    for (int i = 1; i <= np; ++i) {
      //write (*,'(2xi6,2x,i6,2x,i6,2x,i6)') i,indx(i,1:2),insc(i)
      std::println(
        "{:6} {:6} {:6} {:6}",
        i, indx[i-1][0], indx[i-1][1], insc[i-1]
      );
    }
    std::println(" ");
    std::println("    IT    NODE(IT,1:6)");
    std::println(" ");
    for (int it = 1; it <= nelemn; ++it) {
      std::println(
        "{:6} {:6} {:6} {:6} {:6} {:6} {:6}",
        it, node[it-1][0], node[it-1][1], node[it-1][2], node[it-1][3], node[it-1][4], node[it-1][5]
      );
    }
  }

  std::println("Number of unknowns = {}", neqn);
  if ( maxeqn < neqn ) {
    std::println("SETGRD - Too many unknowns!");
    std::println("The maximum allowed is MAXEQN = {}", maxeqn);
    std::println("This problem requires NEQN = {}", neqn);
    std::exit(1);
  }
}

//subroutine setlin(iline,indx,iwrite,long,mx,my,nodex0,np, &
//  nx,ny,xlngth)
//******************************************************
void setlin(
  std::array<int, my> &iline,
  std::array<std::array<int, 2>, np> &indx,
  int &iwrite,
  bool &_long,
  int &nodex0,
  double &xlngth
){
//? setlin() gets the unknown indices along the profile line.
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer my
//  integer np
//
//  integer i
//  integer iline(my)
//  integer indx(np,2)
//  integer ip
//  integer itemp
//  integer iwrite
//  logical long
//  integer mx
//  integer nodex0
//  integer nx
//  integer ny
//  real ( kind = rk8 ) xlngth
//!
//!  Determine the number of a node on the profile line
//!
  int ip;
  int itemp = std::lround((2.0 * (nx - 1) * 9.0) / xlngth);
  if ( _long ) {
    nodex0 = itemp*(2*ny-1)+1;
  } else {
    nodex0 = itemp + 1;
  }

  //do i = 1, my
  for (int i = 1; i <= my; ++i) {
    if (_long) {
      ip = nodex0 + (i-1);
    } else {
      ip = nodex0 + mx*(i-1);
    }
    iline[i-1] = indx[ip-1][0];
  }

  if ( 1 <= iwrite ) {
    std::println(" ");
    std::println("SETLIN: unknown numbers along line:");
    std::println(" ");
    for (int i = 0; i < my; ++i) {
      std::print("{:>5}", iline[i]);
      if ((i + 1) % 15 == 0 && (i + 1) < my) {
          std::print("\n ");
      }
    }
    std::println(" ");
  }
}

//******************************************************
//subroutine setqud ( area, nelemn, nnodes, node, np, nquad, xc, xm, yc, ym )
void setqud(
  std::array<double, nelemn> &area,
  std::array<std::array<int, nnodes>, nelemn> &node,
  std::array<double, np> &xc,
  std::array<std::array<double, nquad>, nelemn> &xm,
  std::array<double, np> &yc,
  std::array<std::array<double, nquad>, nelemn> &ym
) {
//? setqud() sets midpoint quadrature rule information.
//  integer nelemn
//  integer nnodes
//  integer np
//  integer nquad
//
//  real ( kind = rk8 ) area(nelemn)
//  integer ip1
//  integer ip2
//  integer ip3
//  integer it
//  integer node(nelemn,nnodes)
//  real ( kind = rk8 ) x1
//  real ( kind = rk8 ) x2
//  real ( kind = rk8 ) x3
//  real ( kind = rk8 ) xc(np)
//  real ( kind = rk8 ) xm(nelemn,nquad)
//  real ( kind = rk8 ) y1
//  real ( kind = rk8 ) y2
//  real ( kind = rk8 ) y3
//  real ( kind = rk8 ) yc(np)
//  real ( kind = rk8 ) ym(nelemn,nquad)
//
  //do it = 1, nelemn
  int ip1, ip2, ip3;
  double x1, x2, x3, y1, y2, y3;
  for (int it = 0; it < nelemn; ++it) {
    ip1 = node[it][0];
    ip2 = node[it][2];
    ip3 = node[it][3];
    x1 = xc[ip1-1];
    x2 = xc[ip2-1];
    x3 = xc[ip3-1];
    y1 = yc[ip1-1];
    y2 = yc[ip2-1];
    y3 = yc[ip3-1];
    xm[it][0] = 0.5*(x1+x2);
    xm[it][1] = 0.5*(x2+x3);
    xm[it][2] = 0.5*(x3+x1);
    ym[it][0] = 0.5*(y1+y2);
    ym[it][1] = 0.5*(y2+y3);
    ym[it][2] = 0.5*(y3+y1);
    area[it] = 0.5*std::abs(
      (y1+y2)*(x2-x1)+(y2+y3)*(x3-x2)
      +(y3+y1)*(x1-x3)
    );
  }
}

//******************************************************
void setxy(
    int iwrite,
    bool _long,
    std::array<double, np> &xc,
    double &xlngth,
    std::array<double, np> &yc,
    double &ylngth
){
  //? setxy() sets the X, Y coordinates of grid points.
  int ic, jc;
  //!
  //!  Construct grid coordinates
  //!
  //do ip = 1,np
  for (int ip = 1; ip <= np; ++ip) {
    if (_long) {
      ic = ((ip-1)/my)+1;
      jc = ((ip-1)%my)+1;
    } else {
      ic = ((ip-1)%mx)+1;
      jc = ((ip-1)/mx)+1;
    }

    xc[ip-1] = (ic-1)*xlngth/(2*nx-2);
    yc[ip-1] = (jc-1)*ylngth/(2*ny-2);
  }

  if (2 <= iwrite) {
    std::println(" ");
    std::println("    I      XC           YC");
    std::println(" ");
    for(int i = 0; i < np; ++i) {
      std::println(
        " {:>5}{:>12.5f}{:>12.5f}",
        i+1, xc[i], yc[i]
      );
    }
  }
}

//******************************************************
double ubdry (double y, double para) {
//? ubdry() sets the parabolic inflow in terms of the value of the parameter.

  return 4.0 * para * y * ( 3.0 - y ) / 9.0;
}

//subroutine uval ( g, indx, iquad, it, nelemn, neqn, nnodes, node, np, nquad, &
//  para, phi, un, unx, uny, yc )
//
//******************************************************
void uval(
    std::array<double, maxeqn> &g,
    std::array<std::array<int, 2>, np> &indx,
    int iquad,
    int it,
    int neqn,
    std::array<std::array<int, nnodes>, nelemn> &node,
    double para,
    std::array<
      std::array<
        std::array<std::array<double, 3>, nnodes>,
        nquad
      >, nelemn
    > &phi,
    std::array<double, 2> &un,    // ← tamaño fijo 2
    std::array<double, 2> &unx,   // ← tamaño fijo 2
    std::array<double, 2> &uny,   // ← tamaño fijo 2
    std::array<double, np> &yc
) {
//? uval() evaluates the velocities at a given point in a particular triangle.
//  integer nelemn
//  integer neqn
//  integer nnodes
//  integer np
//  integer nquad
//
//  real ( kind = rk8 ) bb
//  real ( kind = rk8 ) bx
//  real ( kind = rk8 ) by
//  real ( kind = rk8 ) g(neqn)
//  integer indx(np,2)
//  integer ip
//  integer iq
//  integer iquad
//  integer it
//  integer iuk
//  integer iun
//  integer node(nelemn,nnodes)
//  real ( kind = rk8 ) para
//  real ( kind = rk8 ) phi(nelemn,nquad,nnodes,3)
//  real ( kind = rk8 ) ubc
//  real ( kind = rk8 ) ubdry
//  real ( kind = rk8 ) un(2)
//  real ( kind = rk8 ) unx(2)
//  real ( kind = rk8 ) uny(2)
//  real ( kind = rk8 ) yc(np)

  int ip;
  int iun;
  double bb;
  double bx;
  double by;
  double ubc;

  un[0] = 0.0; un[1] = 0.0;
  uny[0] = 0.0; uny[1] = 0.0;
  unx[0] = 0.0; unx[1] = 0.0;

  //do iq = 1, nnodes
  for (int iq = 0; iq < nnodes; ++iq) {

    ip = node[it][iq];
    bb = phi[it][iquad][iq][1];
    bx = phi[it][iquad][iq][2];
    by = phi[it][iquad][iq][3];

    for (int iuk = 0; iuk < 2; ++iuk) {

      iun = indx[ip][iuk]-1;

      if ( 0 < iun ) {
        un[iuk] += bb * g[iun];
        unx[iuk] += bx * g[iun];
        uny[iuk] += by * g[iun];
      } else if ( iun < 0 ) {
        ip = node[it][iq]-1;
        ubc = ubdry(yc[ip],para);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }
}

//subroutine uv_plot3d (f,indx,insc,ivunit,long,mx,my, &
//  nelemn,neqn,nnodes,node,np,para,press,reynld,yc)
//
//******************************************************
//!
//!! uv_plot3d() creates a velocity file for use by PLOT3D.
//!
//!  Given the following set of nodes:
//!
//!    A  B  C
//!    D  E  F
//!    G  H  I
//!
//!  the file will have the form:
//!
//!    D, U(G), V(G), P
//!    D, U(H), V(H), P
//!    D, U(I), V(I), P
//!    D, U(D), V(D), P
//!    D, U(E), V(E), P
//!    D, U(F), V(F), P
//!    D, U(A), V(A), P
//!    D, U(B), V(B), P
//!    D, U(C), V(C), P
//!
//!  Here both D and P are set to 1 for now, representing dummy values
//!  of density and pressure.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//!    Input, real ( kind = rk8 ) PARA, the value of the parameter.
//!
//!    Workspace, real ( kind = rk8 ) PRESS(MX,MY), used to hold the
//!    computed pressures.
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer nelemn
//  integer neqn
//  integer nnodes
//  integer np
//
//  real ( kind = rk8 ) alpha
//  real ( kind = rk8 ) dval
//  real ( kind = rk8 ) f(neqn)
//  real ( kind = rk8 ) fsmach
//  integer i
//  integer ii
//  integer indx(np,2)
//  integer insc(np)
//  integer ip
//  integer, save :: iset = 0
//  integer ivunit
//  integer j
//  integer k
//  logical long
//  integer mx
//  integer my
//  integer node(nelemn,nnodes)
//  real ( kind = rk8 ) para
//  real ( kind = rk8 ) press(mx,my)
//  real ( kind = rk8 ) reynld
//  real ( kind = rk8 ) time
//  real ( kind = rk8 ) ubdry
//  real ( kind = rk8 ) uval
//  real ( kind = rk8 ) vval
//  real ( kind = rk8 ) yc(np)
//
//  iset = iset+1
//
//  dval = 1.0
//  fsmach = 1.0
//  alpha = 1.0
//  time = 1.0
//
//  call pval (f,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)
//!
//!  If NY < NX, then nodes with a constant Y value are numbered consecutively.
//!
//  if ( long ) then
//    write(ivunit,'(2I5)')mx,my
//    write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time
//    do ii = 1,4
//      do j = 1,my
//        do i = 1,mx
//          ip = (i-1)*my+j
//          if ( ii == 1 ) then
//            write(ivunit,'(G15.5)')dval
//          else if ( ii == 2 ) then
//            k = indx(ip,1)
//            if (k == 0) then
//              uval = 0.0
//            else if (k < 0) then
//              uval = ubdry(yc(ip),para)
//            else
//              uval = f(k)
//            end if
//            write(ivunit,'(G15.5)')uval
//          else if ( ii == 3 ) then
//            k = indx(ip,2)
//            if (k == 0) then
//              vval = 0.0
//            else
//              vval = f(k)
//            end if
//            write(ivunit,'(G15.5)')vval
//          else
//            write(ivunit,'(G15.5)')press(i,j)
//          end if
//        end do
//      end do
//    end do
//!
//!  If NX < NY, then nodes with a constant X value are numbered consecutively.
//!
//  else
//    write(ivunit,'(2I5)')mx,my
//    write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time
//    do ii = 1,4
//      do i = 1,mx
//        do j = 1,my
//          if ( ii == 1 ) then
//             write(ivunit,'(G15.5)')dval
//          else if ( ii == 2 ) then
//            ip = (i-1)*my+j
//            k = indx(ip,1)
//            if (k == 0) then
//              uval = 0.0
//            else if (k < 0) then
//              uval = ubdry(yc(i),para)
//            else
//              uval = f(k)
//            end if
//            write(ivunit,'(G15.5)')uval
//          else if ( ii == 3 ) then
//            k = indx(ip,2)
//            if (k == 0) then
//              vval = 0.0
//            else
//              vval = f(k)
//            end if
//            write(ivunit,'(G15.5)')vval
//          else
//            write(ivunit,'(G15.5)')press(i,j)
//          end if
//        end do
//      end do
//    end do
//  end if
//
//  std::println("UV_PLOT3D wrote data set ',iset,' to file.")
//
//  return
//end
//subroutine uv_table ( f, indx, ivunit, neqn, np, para, yc )
//
//******************************************************
//!
//!! uv_table() creates a velocity table file.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    28 February 2006
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer neqn
//  integer np
//
//  real ( kind = rk8 ) f(neqn)
//  integer indx(np,2)
//  integer ip
//  integer ivunit
//  integer k
//  real ( kind = rk8 ) para
//  real ( kind = rk8 ) ubdry
//  real ( kind = rk8 ) uval
//  real ( kind = rk8 ) vval
//  real ( kind = rk8 ) yc(np)
//
//  do ip = 1, np
//
//    k = indx(ip,1)
//    if ( k == 0 ) then
//      uval = 0.0
//    else if ( k < 0 ) then
//      uval = ubdry ( yc(ip), para )
//    else
//      uval = f(k)
//    end if
//
//    k = indx(ip,2)
//    if ( k == 0 ) then
//      vval = 0.0
//    else
//      vval = f(k)
//    end if
//
//    write ( ivunit, '(2x,g14.6,2x,g14.6)' ) uval, vval
//
//  end do
//
//  return
//end
//subroutine xy_plot3d ( ixunit, long, np, nx, ny, xc, yc )
//
//******************************************************
//!
//!! xy_plot3d() creates a grid file for use by PLOT3D.
//!
//!  Given the following set of nodes:
//!
//!    A  B  C
//!    D  E  F
//!    G  H  I
//!
//!  the file will have the form:
//!
//!    X(G), X(H), X(I), X(D), X(E), X(F), X(A), X(B), X(C),
//!    Y(G), Y(H), Y(I), Y(D), Y(E), Y(F), Y(A), Y(B), Y(C).
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    20 January 2007
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer np
//
//  integer i
//  integer ip
//  integer, save :: iset = 0
//  integer ixunit
//  integer j
//  logical long
//  integer mx
//  integer my
//  integer nx
//  integer ny
//  real ( kind = rk8 ) xc(np)
//  real ( kind = rk8 ) yc(np)
//
//  iset = iset+1
//
//  mx = 2*nx-1
//  my = 2*ny-1
//!
//!  If NY < NX, then nodes with a constant Y value are numbered consecutively.
//!
//  if ( long ) then
//    write(ixunit,'(2I15)')mx,my
//    do i = 1,my
//      do j = 1,mx
//        ip = (j-1)*my+i
//        write(ixunit,'(G15.5)')xc(ip)
//      end do
//    end do
//
//    do i = 1,my
//      do j = 1,mx
//        ip = (j-1)*my+i
//        write(ixunit,'(G15.5)')yc(ip)
//      end do
//    end do
//!
//!  If NX < NY, then nodes with a constant X value are numbered consecutively.
//!
//  else
//
//    write(ixunit,'(2I15)')my,mx
//
//    do j = 1,mx
//      do i = 1,my
//        ip = (j-1)*my+i
//        write(ixunit,'(G15.5)')xc(ip)
//      end do
//    end do
//
//    do j = 1,mx
//      do i = 1,my
//        ip = (j-1)*my+i
//        write(ixunit,'(G15.5)')yc(ip)
//      end do
//    end do
//
//  end if
//
//  std::println("XYDUMP wrote data set ',iset,' to file.")
//
//  return
//end
//subroutine xy_table ( ixunit, np, xc, yc )
//
//******************************************************
//!
//!! xy_table() creates an XY table file.
//!
//!  Licensing:
//!
//!    This code is distributed under the MIT license.
//!
//!  Modified:
//!
//!    28 February 2006
//!
//!  Author:
//!
//!    John Burkardt
//!
//!  Parameters:
//!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer np
//
//  integer ip
//  integer ixunit
//  real ( kind = rk8 ) xc(np)
//  real ( kind = rk8 ) yc(np)
//
//  do ip = 1, np
//    write ( ixunit, '(2x,g14.6,2x,g14.6)' ) xc(ip), yc(ip)
//  end do
//
//  return
//end