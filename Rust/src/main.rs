use chrono::Local;
use chrono::Datelike;

//******************************************************
//? Constantes
const NX: usize = 21;
const NY: usize = 7;

const MAXROW: usize = 27 * NY;
const NELEMN: usize = 2 * (NX - 1) * (NY - 1);
const MX: usize = 2 * NX - 1;
const MY: usize = 2 * NY - 1;
const NP: usize = MX * MY;
const MAXEQN: usize = 2 * MX * MY + NX * NY;
const NNODES: usize = 6;
const NQUAD: usize = 3;

fn fmt_e(val: f64, prec: usize) -> String {
  if val == 0.0 {
    return format!(" {:.*}E+00", prec, 0.0);
  }
  let abs = val.abs();
  let exp = (abs.log10().floor() as i32) + 1;
  let mant = abs / 10.0_f64.powi(exp);
  let sign = if val < 0.0 { "-" } else { " " };
  format!("{}{:.*}E{:+.02}", sign, prec, mant, exp)
}

////******************************************************
fn timestamp() {
  let now = Local::now();
  println!("Fecha actual: {:02}/{:02}/{:04}", now.day(), now.month(), now.year());
}

fn main() {
  let mut a: Vec<Vec<f64>> = vec![vec![0.0; MAXEQN]; MAXROW];
  let mut a2: f64;
  let mut abound: f64;
  let mut anew: f64;
  let mut aold: f64;
  let mut area: Vec<f64> = vec![0.0; NELEMN];
  let mut dcda: Vec<f64> = vec![0.0; MY];
  let mut f: Vec<f64> = vec![0.0; MAXEQN];

  let mut g = vec![0.0; MAXEQN];
  let mut gr = vec![vec![0.0; MY]; MY];

  let mut iline: Vec<i32> = vec![-1; MY];
  let mut indx: Vec<Vec<i32>> = vec![vec![-1; 2]; NP];
  let mut insc: Vec<i32> = vec![-1; NP];

  let mut ipivot: Vec<i32> = vec![0; MAXEQN];
  let mut isotri: Vec<i32> = vec![0; NELEMN];

  let mut _long: bool = false;

  let mut nband: usize = 0;
  let mut neqn: usize = 0;
  let mut nlband: usize = 0;
  let mut node = vec![vec![0; NNODES]; NELEMN];

  let mut nodex0: usize = 0;

  let mut nrow = 0;

  let mut para: f64;
  let mut phi: Vec<Vec<Vec<Vec<f64>>>> = vec![vec![vec![vec![0.0; 3]; NNODES]; NQUAD]; NELEMN];
  let mut psi = vec![vec![vec![0.0; NNODES]; NQUAD]; NELEMN];

  let mut r = vec![0.0; MY];
  let mut res = vec![0.0; MAXEQN];

  let mut rjpold: f64;
  let mut test: f64;

  let mut ui = vec![0.0; MY];
  let mut unew = vec![0.0; MY];
  let mut xc = vec![0.0; NP];

  let mut xm = vec![vec![0.0; NQUAD]; NELEMN];
  let mut yc = vec![0.0; NP];

  let mut ym = vec![vec![0.0; NQUAD]; NELEMN];

  let inicio: std::time::Instant = std::time::Instant::now();

  timestamp();
  println!(" ");
  println!("channel():");
  println!("  Fortran90 version");
  println!("  Channel flow control problem");
  println!(" ");
  println!("Last modified:");
  println!("  Sabra dios");
  println!("");
  println!("  Flow control problem:");
  println!("    Inflow controlled by one parameter.");
  println!("    Velocities measured along vertical line.");
  println!("    Try to match specified velocity profile.");
  //?  Set input data
  let fileg: &str = "display.txt";
  let fileu: &str = "uv.dat";
  let filex: &str = "xy.dat";

  //let iounit: usize = 2;
  //let ivunit: usize = 4;
  let iwrite: usize = 10;
  //let ixunit: usize = 3;
  let maxnew: usize = 10;
  let maxsec: usize = 8;
  //let npara: usize = 1;
  let mut numnew: usize = 0;
  let mut numsec: usize = 0;
  let reynld: f64 = 1.0;
  let mut rjpnew: f64;
  let tolnew: f64 = 1.0E-04;
  let tolsec: f64 = 1.0E-06;
  let xlngth: f64 = 10.0;
  let ylngth: f64 = 3.0;

  println!("");
  println!("NX = {}", NX);
  println!("NY = {}", NY);
  println!("Number of elements = {}", NELEMN);
  println!("Reynolds number = {}", reynld);
  println!("Secant tolerance = {}", tolsec);
  println!("Newton tolerance = {}", tolnew);
  println!("");
  //?  SETGRD constructs grid, numbers unknowns, calculates areas,
  //?  and points for midpoint quadrature rule.
  setgrd(
    &mut indx,
    &mut insc,
    &mut isotri,
    iwrite,
    &mut _long,
    &mut neqn,
    &mut node,
  );
  //?  Compute the bandwidth
  setban(
    &indx,
    &insc,
    &mut nband,
    &mut nlband,
    &node,
    &mut nrow,
  );
  //?  Record variable numbers along profile sampling line.
  setlin(
    &mut iline,
    &indx,
    iwrite,
    &_long,
    &mut nodex0,
    xlngth,
  );
  //?  Set the coordinates of grid points.
  setxy(
    iwrite,
    &_long,
    &mut xc,
    xlngth,
    &mut yc,
    ylngth,
  );
  ////?  Set quadrature points
  setqud(&mut area, &node, &xc, &mut xm, &yc, &mut ym);
  //?  Evaluate basis functions at quadrature points
  setbas(
    &node, &xc, &yc, &mut phi, &mut psi, &xm, &ym,
  );
  //?  NSTOKE now solves the Navier Stokes problem for an inflow
  ////?  parameter of 1.0.
  para = 1.0;
  println!(" ");
  println!("Solve Navier Stokes problem with parameter = {}", para);
  println!("for profile at x = {}", xc[nodex0]);
  for i in 0..neqn {
    g[i] = 1.0;
  }

  nstoke(
    &mut a,
    &area,
    &mut f,
    &mut g,
    &indx,
    &insc,
    &mut ipivot,
    iwrite,
    maxnew,
    neqn,
    nlband,
    &node,
    nrow,
    &mut numnew,
    para,
    &phi,
    &psi,
    reynld,
    tolnew,
    &yc,
  );
  //?  RESID computes the residual at the given solution
  if 1 <= iwrite {
    resid(
      &area, &f, &indx, &insc, iwrite, neqn, &node, para, &phi, &psi, &mut res, reynld,
      &yc,
    );
  }
  //?  GETG computes the internal velocity profile at X = XC(NODEX0), which will
  //?  be used to measure the goodness-of-fit of the later solutions.
  getg(&f, &iline, MY, neqn, &mut ui);

  if 1 <= iwrite {
    println!(" ");
    println!("U profile:");
    println!(" ");
    for chunk in ui.chunks(5) {
      for &v in chunk {
        std::print!("{:>14.6}", v);
      }
      println!("");
    }
  }
  //?  GRAM generates the Gram matrix GR and the vector
  //?  R = line integral of ui*phi
  gram(
    &mut gr, &iline, &indx, iwrite, &node, nodex0, para, &mut r, &ui, &xc, &yc,
  );

  if let Err(e) = xy_table(&xc, &yc, filex) {
    eprintln!("Error writing xy.txt: {}", e);
  }
  if let Err(e) = uv_table(&f, &indx, para, &yc, fileu) {
    eprintln!("Error writing uv.txt: {}", e);
  }
  if let Err(e) = xy_plot3d(_long, &xc, &yc, fileg) {
    eprintln!("Error writing xy_plot3d.dat: {}", e);
  }

  //?  Destroy information about true solution
  for i in 0..neqn {
    f[i] = 0.0;
    g[i] = 0.0;
  }
  //?  Secant iteration loop
  aold = 0.0;
  rjpold = 0.0;
  anew = 0.1;

  for iter in 1..=maxsec {
    numsec += 1;
    println!(" ");
    println!("Secant iteration {}", iter);
    //?  Solve for unew at new value of parameter anew
    println!(" ");
    println!("Solving Navier Stokes problem for parameter = {}", anew);
    //?  Use solution F at previous value of parameter for starting point.
    for i in 0..neqn {
      g[i] = f[i];
    }
    para = anew;

    nstoke(
      &mut a,
      &area,
      &mut f,
      &mut g,
      &indx,
      &insc,
      &mut ipivot,
      iwrite,
      maxnew,
      neqn,
      nlband,
      &node,
      nrow,
      &mut numnew,
      para,
      &phi,
      &psi,
      reynld,
      tolnew,
      &yc,
    );
    //?  Get velocity profile
    getg(&f, &iline, MY, neqn, &mut unew);

    if 1 <= iwrite {
      println!(" ");
      println!("Velocity profile:");
      println!(" ");
      for chunk in unew.chunks(5) {
        for &v in chunk {
          std::print!("{:>14.6}", v);
        }
        println!("");
      }
    }
    //?  Solve linear system for du/da
    para = anew;
    abound = 1.0;
    linsys(
      &mut a,
      &mut area,
      &mut g,
      &mut f,
      &mut indx,
      &mut insc,
      &mut ipivot,
      neqn,
      nlband,
      &mut node,
      nrow,
      para,
      abound,
      &mut phi,
      &mut psi,
      reynld,
      &mut yc,
    );
    //?  Output in DCDA
    getg(&g, &iline, MY, neqn, &mut dcda);

    if 2 <= iwrite {
      println!(" ");
      println!("Sensitivities:");
      println!(" ");
      for chunk in dcda.chunks(5) {
        for &v in chunk {
          std::print!("{:>14.6}", v);
        }
        println!("");
      }
    }
    //?  Evaluate J prime at current value of parameter where J is
    //?  functional to be minimized.
    //?  JPRIME = 2.0 * DCDA(I) * (GR(I,J)*UNEW(J)-R(I))
    rjpnew = 0.0;
    for i in 0..MY {
      let mut temp = -r[i];
      for j in 0..MY {
        temp += gr[i][j] * unew[j];
      }
      rjpnew += 2.0 * dcda[i] * temp;
    }

    println!(" ");
    println!("Parameter  = {} J prime = {}", anew, rjpnew);
    //?  Update the estimate of the parameter using the secant step
    if iter == 1 {
      a2 = 0.5;
    } else {
      a2 = aold - rjpold * (anew - aold) / (rjpnew - rjpold);
    }

    aold = anew;
    anew = a2;
    rjpold = rjpnew;
    test = (anew - aold).abs() / anew.abs();

    println!("New value of parameter = {}", anew);
    println!("Convergence test = {}", test);

    if (anew - aold).abs() < anew.abs() * tolsec {
      println!("Secant iteration converged.");
      break;
    }
  }

  if maxsec < 1 || (anew - aold).abs() >= anew.abs() * tolsec {
    println!("Secant iteration failed to converge.");
  }

  println!("Number of secant steps = {}", numsec);
  println!("Number of Newton steps = {}", numnew);
  //?  Terminate.
  println!(" ");
  println!("CHANNEL:");
  println!("  Normal end of execution.");
  println!(" ");
  timestamp();
  let duracion = inicio.elapsed();
  println!("Tiempo transcurrido: {:?}", duracion);
}
//******************************************************
fn bsp(
  xq: f64,
  yq: f64,
  it: usize,
  iq: usize,
  id: usize,
  node: &[Vec<usize>],
  xc: &[f64],
  yc: &[f64],
) -> f64 {
  //? bsp() evaluates the linear basis functions associated with pressure.
  let iq1 = iq;
  let iq2: usize = (iq + 1) % 3;
  let iq3: usize = (iq + 2) % 3;
  let i1: usize = node[it][iq1];
  let i2: usize = node[it][iq2];
  let i3: usize = node[it][iq3];

  let d: f64 = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);

  if id == 1 {
    return 1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
  } else if id == 2 {
    return (yc[i2] - yc[i3]) / d;
  } else if id == 3 {
    return (xc[i3] - xc[i2]) / d;
  } else {
    println!("BSP - fatal error!");
    println!("unknown value of id = {}", id);
    std::process::exit(1);
  }
}

//subroutine daxpy ( n, da, dx, incx, dy, incy )
//
//******************************************************
//?! daxpy() computes constant times a vector plus a vector.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    Jack Dongarra
//?  Parameters:
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) dx(*),dy(*),da
//  integer i,incx,incy,ix,iy,m,n
//  if ( n <= 0)return
//  if (da  ==  0.0 ) return
//  if ( incx == 1.and.incy == 1)go to 20
//?        code for unequal increments or equal increments
//?          not equal to 1
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
//?        code for both increments equal to 1
//?        clean-up loop
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
//?! dcopy() copies a vector, x, to a vector, y.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    Jack Dongarra
//?  Parameters:
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) dx(*),dy(*)
//  integer i,incx,incy,ix,iy,m,n
//
//  if ( n <= 0)return
//  if ( incx == 1.and.incy == 1)go to 20
//?        code for unequal increments or equal increments
//?          not equal to 1
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
//?        code for both increments equal to 1
//?        clean-up loop
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
//?! ddot() forms the dot product of two vectors.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    Jack Dongarra
//?  Parameters:
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
//?        code for unequal increments or equal increments
//?          not equal to 1
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
//?        code for both increments equal to 1
//?        clean-up loop
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
//?! delete() deletes a file.
//?  Licensing:
//?    This code is distributed under the MIT license.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    John Burkardt
//?  Parameters:
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

//******************************************************
fn dgbfa(
  abd: &mut [Vec<f64>],
  n: usize,
  ml: usize,
  mu: usize,
  ipvt: &mut [i32],
  info: &mut usize,
) {
  let m = ml + mu + 1;
  *info = 0;

  let j0 = mu + 2;
  let j1 = std::cmp::min(n, m) - 1;

  if j0 <= j1 {
    for jz in j0..=j1 {
      let i0 = if m + 1 > jz { m + 1 - jz } else { 1 };
      if i0 <= ml {
        for i in i0..=ml {
          abd[i - 1][jz - 1] = 0.0;
        }
      }
    }
  }

  let mut jz = j1;
  let mut ju = 0;

  for k in 1..=n - 1 {
    jz += 1;
    if jz <= n {
      for i in 0..ml {
        abd[i][jz - 1] = 0.0;
      }
    }

    let lm = std::cmp::min(ml, n - k);
    let mut l = m;
    let mut dmax = abd[m - 1][k - 1].abs();
    for i in 1..=lm {
      let abs_val = abd[m - 1 + i][k - 1].abs();
      if dmax < abs_val {
        dmax = abs_val;
        l = m + i;
      }
    }
    ipvt[k - 1] = (l + k - m) as i32;

    if abd[l - 1][k - 1] == 0.0 {
      *info = k;
    } else {
      if l != m {
        let t = abd[l - 1][k - 1];
        abd[l - 1][k - 1] = abd[m - 1][k - 1];
        abd[m - 1][k - 1] = t;
      }

      let t = -1.0 / abd[m - 1][k - 1];
      for i in 0..lm {
        abd[m + i][k - 1] *= t;
      }

      ju = std::cmp::max(ju, mu + ipvt[k - 1] as usize);
      if ju > n {
        ju = n;
      }
      let mut mm = m;
      let mut ll = l;
      for j in k + 1..=ju {
        ll -= 1;
        mm -= 1;
        let t = abd[ll - 1][j - 1];
        if ll != mm {
          abd[ll - 1][j - 1] = abd[mm - 1][j - 1];
          abd[mm - 1][j - 1] = t;
        }
        for i in 0..lm {
          abd[mm + i][j - 1] += t * abd[m + i][k - 1];
        }
      }
    }
  }

  ipvt[n - 1] = n as i32;
  if abd[m - 1][n - 1] == 0.0 {
    *info = n;
  }
}

//******************************************************
fn dgbsl(
  abd: &[Vec<f64>],
  _lda: usize,
  n: usize,
  ml: usize,
  mu: usize,
  ipvt: &[i32],
  b: &mut [f64],
  job: usize,
) {
  let m = mu + ml + 1;

  if job == 0 {
    // Solve a * x = b
    // First solve L * y = b
    if 0 < ml {
      for k in 0..n - 1 {
        let lm = std::cmp::min(ml, n - 1 - k);
        let l = (ipvt[k] - 1) as usize;
        let t = b[l];
        if l != k {
          b[l] = b[k];
          b[k] = t;
        }
        for i in 0..lm {
          b[k + 1 + i] += t * abd[m + i][k];
        }
      }
    }
    // Now solve U * x = y
    for k in (0..n).rev() {
      b[k] /= abd[m - 1][k];
      let lm = std::cmp::min(k + 1, m) - 1;
      if lm > 0 {
        let la = m - lm - 1;
        let lb = k - lm;
        let t = -b[k];
        for i in 0..lm {
          b[lb + i] += t * abd[la + i][k];
        }
      }
    }
  } else {
    // Solve trans(a) * x = b
    // First solve trans(U) * y = b
    for k in 0..n {
      let lm = std::cmp::min(k + 1, m) - 1;
      if lm > 0 {
        let la = m - lm - 1;
        let lb = k - lm;
        let mut t = 0.0;
        for i in 0..lm {
          t += abd[la + i][k] * b[lb + i];
        }
        b[k] = (b[k] - t) / abd[m - 1][k];
      } else {
        b[k] /= abd[m - 1][k];
      }
    }
    // Now solve trans(L) * x = y
    if 0 < ml {
      for k in (0..n - 1).rev() {
        let lm = std::cmp::min(ml, n - k - 1);
        if lm > 0 {
          let mut t = 0.0;
          for i in 0..lm {
            t += abd[m + i][k] * b[k + 1 + i];
          }
          b[k] += t;
        }
        let l = (ipvt[k] - 1) as usize;
        if l != k {
          let t = b[l];
          b[l] = b[k];
          b[k] = t;
        }
      }
    }
  }
}
//subroutine dscal ( n, da, dx, incx )
//
//******************************************************
//?! dscal() scales a vector by a constant.
//?     uses unrolled loops for increment equal to one.
//?     jack dongarra, linpack, 3/11/78.
//?     modified 3/93 to return if incx  <=  0.
//?     modified 12/3/93, array(1) declarations changed to array(*)
//?  Parameters:
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  real ( kind = rk8 ) da,dx(*)
//  integer i,incx,m,n,nincx
//  if (  n <= 0 .or. incx <= 0 )return
//  if ( incx == 1)go to 20
//?        code for increment not equal to 1
//  nincx = n*incx
//  do i = 1,nincx,incx
//    dx(i) = da*dx(i)
//  end do
//  return
//?        code for increment equal to 1
//?        clean-up loop
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
//?! gdump() writes information to a file.
//?  Discussion:
//?    The information can be used to create
//?    graphics images.  In order to keep things simple, exactly one
//?    value, real or integer, is written per record.
//?  Licensing:
//?    This code is distributed under the MIT license.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    John Burkardt
//?  Parameters:
//?    Input, integer NPARA, the number of parameters.  Fixed at 1
//?    for now.
//?    Input, real ( kind = rk8 ) PARA(MAXPAR), the parameters.
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
//?  Pressures
//  do i = 1, np
//    j = insc(i)
//    if (j <= 0) then
//      fval = 0.0
//    else
//      fval = f(j)
//    end if
//    write (iounit,*) fval
//  end do
//?  Horizontal velocities, U
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
//?  Vertical velocities, V
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
//  println!("GDUMP wrote data set ',iset,' to file.")
//
//  return
//end
//******************************************************
fn getg(f: &[f64], iline: &[i32], my: usize, _neqn: usize, u: &mut [f64]) {
  for j in 0..my {
    let k = iline[j];
    if k >= 0 {
      u[j] = f[k as usize];
    } else {
      u[j] = 0.0;
    }
  }
}
//******************************************************
fn gram(
  gr: &mut [Vec<f64>],
  iline: &[i32],
  indx: &[Vec<i32>],
  iwrite: usize,
  node: &[Vec<usize>],
  nodex0: usize,
  para: f64,
  r: &mut [f64],
  ui: &[f64],
  xc: &[f64],
  yc: &[f64],
) {
  let wt: [f64; 3] = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
  let yq: [f64; 3] = [-0.7745966692, 0.0, 0.7745966692];

  for i in 0..MY {
    r[i] = 0.0;
    for j in 0..MY {
      gr[i][j] = 0.0;
    }
  }

  let xzero = xc[nodex0];

  for it in 0..NELEMN {
    let k = node[it][0];
    let kk = node[it][1];

    if 1.0E-04 < (xc[k] - xzero).abs() {
      continue;
    }
    if 1.0E-04 < (xc[kk] - xzero).abs() {
      continue;
    }

    for iquad in 0..3 {
      let bma2 = (yc[kk] - yc[k]) / 2.0;
      let ar = bma2 * wt[iquad];
      let x = xzero;
      let y = yc[k] + bma2 * (yq[iquad] + 1.0);
      // Compute u internal at quadrature points
      let mut uiqdpt = 0.0;
      for iq in 0..6 {
        if iq > 3 {
          continue;
        }
        if iq == 2 {
          continue;
        }
        let mut bb: f64 = 0.0;
        let mut bx: f64 = 0.0;
        let mut by: f64 = 0.0;
        qbf(x, y, it, iq, &mut bb, &mut bx, &mut by, node, xc, yc);
        let ip = node[it][iq];
        let iun = indx[ip][0];
        if iun >= 0 {
          let ii = igetl(iun, iline, MY);
          uiqdpt += bb * ui[ii];
        } else if iun == -2 {
          let ubc = ubdry(yc[ip], para);
          uiqdpt += bb * ubc;
        }
      }
      // Only loop over nodes lying on line x = xzero
      for iq in 0..6 {
        if iq == 0 || iq == 1 || iq == 3 {
          let ip = node[it][iq];
          let mut bb: f64 = 0.0;
          let mut bx: f64 = 0.0;
          let mut by: f64 = 0.0;
          qbf(x, y, it, iq, &mut bb, &mut bx, &mut by, node, xc, yc);
          let i = indx[ip][0];
          if i < 0 {
            continue;
          }
          let ii = igetl(i, iline, MY);
          r[ii] += bb * uiqdpt * ar;

          for iqq in 0..6 {
            if iqq == 0 || iqq == 1 || iqq == 3 {
              let ipp = node[it][iqq];
              let mut bbb: f64 = 0.0;
              let mut bbx: f64 = 0.0;
              let mut bby: f64 = 0.0;
              qbf(x, y, it, iqq, &mut bbb, &mut bbx, &mut bby, node, xc, yc);
              let j = indx[ipp][0];
              if j >= 0 {
                let jj = igetl(j, iline, MY);
                gr[ii][jj] += bb * bbb * ar;
              }
            }
          }
        }
      }
    }
  }

  if 2 <= iwrite {
    println!(" ");
    println!("Gram matrix:");
    println!(" ");
    for i in 0..MY {
      for j in 0..MY {
        println!("{} {} {}", i + 1, j + 1, gr[i][j]);
      }
    }
    println!(" ");
    println!("R vector:");
    println!(" ");
    for i in 0..MY {
      println!("{} {}", i + 1, r[i]);
    }
  }
}

//******************************************************
fn idamax(n: usize, dx: &[f64], incx: usize) -> usize {
  //? idamax() finds the index of element having max. absolute alue.
  //? Si el tamaño es menor que 1 o el incremento no es positivo,
  //? retorna -1 (equivalente a 0 en Fortran para indicar "no encontrado")
  //if n < 1 || incx <= 0 { return -1 };
  let mut idamax_val: usize = 0;
  if n == 1 {
    return idamax_val;
  };

  if incx == 1 {
    // Código para incremento igual a 1
    let mut dmax = dx[0].abs();
    for i in 1..n {
      if dmax < dx[i].abs() {
        idamax_val = i;
        dmax = dx[i].abs();
      }
    }
  } else {
    // Código para incremento diferente de 1
    let mut ix: usize = 0; // Índice basado en 0
    let mut dmax = dx[ix].abs();
    ix += incx;

    for i in 1..n {
      if dx[ix].abs() > dmax {
        idamax_val = i;
        dmax = dx[ix].abs();
      }
      ix += incx;
    }
  }

  idamax_val
}

//******************************************************
fn igetl(k: i32, iline: &[i32], my: usize) -> usize {
  for j in 0..my {
    if iline[j] == k {
      return j;
    }
  }
  println!(" ");
  println!("IGETL - fatal error!");
  println!("  Unable to get local unknown number for ");
  println!("  Global variable number {}", k);
  std::process::exit(1);
}
//end
//******************************************************
fn linsys(
  a: &mut [Vec<f64>],
  area: &[f64],
  f: &mut [f64],
  g: &[f64],
  indx: &[Vec<i32>],
  insc: &[i32],
  ipivot: &mut [i32],
  neqn: usize,
  nlband: usize,
  node: &[Vec<usize>],
  nrow: usize,
  para1: f64,
  para2: f64,
  phi: &[Vec<Vec<Vec<f64>>>],
  psi: &[Vec<Vec<f64>>],
  reynld: f64,
  yc: &[f64],
) {
  //? linsys() sets up and solves the linear system.
  //?    The G array contains the previous solution.
  //?
  //?    The F array contains the right hand side initially and then the
  //?    current solution.

  let mut bb: f64;
  let mut bx: f64;
  let mut by: f64;
  let mut bbl: f64;
  let mut bbb: f64;
  let mut bbx: f64;
  let mut bby: f64;
  let mut bbbl: f64;
  let mut info: usize = 0;
  let mut ip;
  let mut ihor;
  let mut iver;
  let mut iprs;
  let mut ipp;
  let mut ju;
  let mut jv;
  let mut jp;
  let mut iuse;

  let mut un: Vec<f64> = vec![0.0; 2];
  let mut unx: Vec<f64> = vec![0.0; 2];
  let mut uny: Vec<f64> = vec![0.0; 2];

  let ioff: i32 = (nlband + nlband) as i32;
  let mut ar: f64;
  let mut uu: f64;
  let visc: f64 = 1.0 / reynld;
  for j in 0..neqn {
    f[j] = 0.0;
  }

  for k in 0..nrow {
    for j in 0..neqn {
      a[k][j] = 0.0;
    }
  }
  //?  For each element,
  for it in 0..NELEMN {
    ar = area[it] / 3.0;
    //?  and for each quadrature point in the element,
    for iquad in 0..NQUAD {
      //?  Evaluate velocities at quadrature point
      uval(
        &g, indx, iquad, it, node, para1, phi, &mut un, &mut unx, &mut uny, yc,
      );
      for iq in 0..NNODES {
        ip = node[it][iq];
        bb = phi[it][iquad][iq][0];
        bx = phi[it][iquad][iq][1];
        by = phi[it][iquad][iq][2];
        bbl = psi[it][iquad][iq];
        ihor = indx[ip][0];
        iver = indx[ip][1];
        iprs = insc[ip];

        if ihor >= 0 {
          f[ihor as usize] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
        }

        if iver >= 0 {
          f[iver as usize] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
        }

        for iqq in 0..NNODES {
          ipp = node[it][iqq];
          bbb = phi[it][iquad][iqq][0];
          bbx = phi[it][iquad][iqq][1];
          bby = phi[it][iquad][iqq][2];
          bbbl = psi[it][iquad][iqq];
          ju = indx[ipp][0];
          jv = indx[ipp][1];
          jp = insc[ipp];
          if ju >= 0 {
            let ju_idx = ju as usize;
            if ihor >= 0 {
              iuse = ihor - ju + ioff;
              a[iuse as usize][ju_idx] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }

            if iver >= 0 {
              iuse = iver - ju + ioff;
              a[iuse as usize][ju_idx] += ar * bb * bbb * unx[1];
            }

            if iprs >= 0 {
              iuse = iprs - ju + ioff;
              a[iuse as usize][ju_idx] += ar * bbx * bbl;
            }
          } else if ju == -2 {
            uu = ubdry(yc[ipp], para2);
            if ihor >= 0 {
              f[ihor as usize] -= ar
                * uu
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }

            if iver >= 0 {
              f[iver as usize] -= ar * uu * bb * bbb * unx[1];
            }

            if iprs >= 0 {
              f[iprs as usize] -= ar * uu * bbx * bbl;
            }
          }
          if jv >= 0 {
            let jv_idx = jv as usize;
            if ihor >= 0 {
              iuse = ihor - jv + ioff;
              a[iuse as usize][jv_idx] += ar * bb * bbb * uny[0];
            }

            if iver >= 0 {
              iuse = iver - jv + ioff;
              a[iuse as usize][jv_idx] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
            }

            if iprs >= 0 {
              iuse = iprs - jv + ioff;
              a[iuse as usize][jv_idx] += ar * bby * bbl;
            }
          }
          if jp >= 0 {
            let jp_idx = jp as usize;
            if ihor >= 0 {
              iuse = ihor - jp + ioff;
              a[iuse as usize][jp_idx] -= ar * bx * bbbl;
            }
            if iver >= 0 {
              iuse = iver - jp + ioff;
              a[iuse as usize][jp_idx] -= ar * by * bbbl;
            }
          }
        }
      }
    }
  }
  //?  To avoid singularity of the pressure system, the last pressure
  //?  is simply assigned a value of 0.
  f[neqn - 1] = 0.0;
  for j in (neqn - nlband - 1)..(neqn - 1) {
    let i = (neqn - 1) - j + (ioff as usize);
    a[i][j] = 0.0;
  }
  a[ioff as usize][neqn - 1] = 1.0;
  //?  Factor the matrix
  dgbfa(a, neqn, nlband, nlband, ipivot, &mut info);

  if info != 0 {
    println!(" ");
    println!("LINSYS - fatal error!");
    println!("DGBFA returns INFO = {}", info);
    std::process::exit(1);
  }
  //?  Solve the linear system
  let job: usize = 0;
  dgbsl(a, MAXROW, neqn, nlband, nlband, ipivot, f, job);
}

//******************************************************
fn nstoke(
  a: &mut [Vec<f64>],
  area: &[f64],
  f: &mut [f64],
  g: &mut [f64],
  indx: &[Vec<i32>],
  insc: &[i32],
  ipivot: &mut [i32],
  iwrite: usize,
  maxnew: usize,
  neqn: usize,
  nlband: usize,
  node: &[Vec<usize>],
  nrow: usize,
  numnew: &mut usize,
  para: f64,
  phi: &[Vec<Vec<Vec<f64>>>],
  psi: &[Vec<Vec<f64>>],
  reynld: f64,
  tolnew: f64,
  yc: &[f64],
) {
  //?! nstoke() solves the Navier Stokes equation using //Taylor-Hood elements.
  //?  The G array contains the previous iterate.
  //?  The F array contains the right hand side initially and //then the current iterate.
  for iter in 0..maxnew {
    *numnew += 1;

    linsys(
      a, area, f, g, indx, insc, ipivot, neqn, nlband, node, nrow, para, para, phi, psi,
      reynld, yc,
    );
    //? Check for convergence
    for i in 0..neqn {
      g[i] = g[i] - f[i];
    }
    let diff = g[idamax(neqn, g, 1)].abs();

    if 1 <= iwrite {
      println!("NSTOKE iteration {} Mnorm = {}", iter + 1, diff);
    }

    for i in 0..neqn {
      g[i] = f[i];
    }

    if diff <= tolnew {
      println!("Navier Stokes iteration converged in {} iterations.", iter + 1);
      return;
    }
  }
  println!("Navier Stokes solution did not converge!");
}
//subroutine pval (g,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)
//
//******************************************************
//?! pval() computes a table of pressures.
//?  Licensing:
//?    This code is distributed under the MIT license.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    John Burkardt
//?  Parameters:
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
//?  Read the pressures where they are computed.
//?  These are "(odd, odd)" points.
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
//?  Interpolate the pressures at points (even, odd) and (odd, even).
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
//?  Interpolate the pressures at points (even,even).
//  do j = 2,my-1,2
//    do i = 2,mx-1,2
//      press(i,j) = 0.5*(press(i-1,j-1)+press(i+1,j+1))
//    end do
//  end do
//
//  return
//end

//******************************************************
fn qbf(
  x: f64,
  y: f64,
  it: usize,
  _in: usize,
  bb: &mut f64,
  bx: &mut f64,
  by: &mut f64,
  node: &[Vec<usize>],
  xc: &[f64],
  yc: &[f64],
) {
  //? qbf() evaluates a quadratic basis function in a triangle.
  let in1;
  let in2;
  let in3;
  let i1;
  let i2;
  let i3;
  let inn;
  let j1;
  let j2;
  let j3;
  let d;
  let t;
  let c;
  let s;

  if _in < 3 {
    in1 = _in;
    in2 = (_in + 1) % 3;
    in3 = (_in + 2) % 3;
    i1 = node[it][in1];
    i2 = node[it][in2];
    i3 = node[it][in3];
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    *bb = t * (2.0 * t - 1.0);
    *bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    *by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    inn = _in - 3;
    in1 = inn;
    in2 = (inn + 1) % 3;
    in3 = (inn + 2) % 3;
    i1 = node[it][in1];
    i2 = node[it][in2];
    i3 = node[it][in3];
    j1 = i2;
    j2 = i3;
    j3 = i1;
    d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    s = 1.0 + ((yc[j2] - yc[j3]) * (x - xc[j1]) + (xc[j3] - xc[j2]) * (y - yc[j1])) / c;
    *bb = 4.0 * s * t;
    *bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    *by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
}

//******************************************************
fn resid(
  area: &[f64],
  g: &[f64],
  indx: &[Vec<i32>],
  insc: &[i32],
  iwrite: usize,
  neqn: usize,
  node: &[Vec<usize>],
  para: f64,
  phi: &[Vec<Vec<Vec<f64>>>],
  psi: &[Vec<Vec<f64>>],
  res: &mut [f64],
  reynld: f64,
  yc: &[f64],
) {
  let visc = 1.0 / reynld;

  for i in 0..neqn {
    res[i] = 0.0;
  }

  for it in 0..NELEMN {
    let ar = area[it] / 3.0;
    for iquad in 0..NQUAD {
      let mut un = vec![0.0; 2];
      let mut unx = vec![0.0; 2];
      let mut uny = vec![0.0; 2];
      uval(
        &g, indx, iquad, it, node, para, phi, &mut un, &mut unx, &mut uny, yc,
      );
      for iq in 0..NNODES {
        let ip = node[it][iq];
        let bb = phi[it][iquad][iq][0];
        let bx = phi[it][iquad][iq][1];
        let by = phi[it][iquad][iq][2];
        let bbl = psi[it][iquad][iq];
        let iprs = insc[ip];
        let ihor = indx[ip][0];
        let iver = indx[ip][1];

        if ihor >= 0 {
          res[ihor as usize] += (un[0] * unx[0] + un[1] * uny[0]) * bb * ar;
        }

        if iver >= 0 {
          res[iver as usize] += (un[0] * unx[1] + un[1] * uny[1]) * bb * ar;
        }

        for iqq in 0..NNODES {
          let ipp = node[it][iqq];
          let bbb = phi[it][iquad][iqq][0];
          let bbx = phi[it][iquad][iqq][1];
          let bby = phi[it][iquad][iqq][2];
          let bbbl = psi[it][iquad][iqq];
          let ju = indx[ipp][0];
          let jv = indx[ipp][1];
          let jp = insc[ipp];

          if ju >= 0 {
            let ju_idx = ju as usize;
            if ihor >= 0 {
              let aijuu =
                visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              res[ihor as usize] += aijuu * ar * g[ju_idx];
            }
            if iver >= 0 {
              let aijvu = bb * bbb * unx[1];
              res[iver as usize] += aijvu * ar * g[ju_idx];
            }
            if iprs >= 0 {
              let aijpu = bbx * bbl;
              res[iprs as usize] += aijpu * ar * g[ju_idx];
            }
          } else if ju == -2 {
            let uu = ubdry(yc[ipp], para);
            if ihor >= 0 {
              let aijuu =
                visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              res[ihor as usize] += ar * aijuu * uu;
            }
            if iver >= 0 {
              let aijvu = bb * bbb * unx[1];
              res[iver as usize] += ar * aijvu * uu;
            }
            if iprs >= 0 {
              let aijpu = bbx * bbl;
              res[iprs as usize] += ar * aijpu * uu;
            }
          }

          if jv >= 0 {
            let jv_idx = jv as usize;
            if ihor >= 0 {
              let aijuv = bb * bbb * uny[0];
              res[ihor as usize] += aijuv * ar * g[jv_idx];
            }
            if iver >= 0 {
              let aijvv =
                visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]);
              res[iver as usize] += aijvv * ar * g[jv_idx];
            }
            if iprs >= 0 {
              let aijpv = bby * bbl;
              res[iprs as usize] += aijpv * ar * g[jv_idx];
            }
          }

          if jp >= 0 {
            let jp_idx = jp as usize;
            if ihor >= 0 {
              let aijup = -bx * bbbl;
              res[ihor as usize] += aijup * ar * g[jp_idx];
            }
            if iver >= 0 {
              let aijvp = -by * bbbl;
              res[iver as usize] += aijvp * ar * g[jp_idx];
            }
          }
        }
      }
    }
  }

  res[neqn - 1] = g[neqn - 1];

  let mut rmax = 0.0;
  let mut imax = 0;
  let mut ibad = 0;

  for i in 0..neqn {
    let test = res[i].abs();
    if rmax < test {
      rmax = test;
      imax = i;
    }
    if 1.0E-03 < test {
      ibad += 1;
    }
  }

  if 1 <= iwrite {
    println!(" ");
    println!("RESIDUAL INFORMATION:");
    println!(" ");
    println!("Worst residual is number {:5}", imax + 1);
    println!("of magnitude {}", fmt_e(rmax, 15));
    println!(" ");
    println!("Number of bad residuals is {} out of {}", ibad, neqn);
    println!(" ");
  }

  if 2 <= iwrite {
    println!("Raw residuals:");
    println!(" ");
    let mut i = 0;
    for j in 0..NP {
      if indx[j][0] >= 0 {
        i += 1;
        let (mark, label) = if res[i - 1].abs() <= 1.0E-03 {
          (" ", "U")
        } else {
          ("*", "U")
        };
        println!("{}{} {:5} {:5} {}", mark, label, i, j + 1, fmt_e(res[i - 1], 6));
      }
      if indx[j][1] >= 0 {
        i += 1;
        let (mark, label) = if res[i - 1].abs() <= 1.0E-03 {
          (" ", "V")
        } else {
          ("*", "V")
        };
        println!("{}{} {:5} {:5} {}", mark, label, i, j + 1, fmt_e(res[i - 1], 6));
      }
      if insc[j] >= 0 {
        i += 1;
        let (mark, label) = if res[i - 1].abs() <= 1.0E-03 {
          (" ", "P")
        } else {
          ("*", "P")
        };
        println!("{}{} {:5} {:5} {}", mark, label, i, j + 1, fmt_e(res[i - 1], 6));
      }
    }
  }
}
//******************************************************
fn setban(
  indx: &[Vec<i32>],
  insc: &[i32],
  nband: &mut usize,
  nlband: &mut usize,
  node: &[Vec<usize>],
  nrow: &mut usize,
) {
  //? setban() computes the half band width.

  *nlband = 0;

  //do it = 1, nelemn
  for it in 0..NELEMN {
    //do iq = 1, nnodes
    for iq in 0..NNODES {
      let ip = node[it][iq];
      for iuk in 0..3 {
        let i: i32 = if iuk == 2 {
          insc[ip]
        } else {
          indx[ip][iuk]
        };

        if i >= 0 {
          for iqq in 0..NNODES {
            let ipp = node[it][iqq];
            for iukk in 0..3 {
              let j = if iukk == 2 {
                insc[ipp]
              } else {
                indx[ipp][iukk]
              };
              if j > i {
                let diff = (j - i) as usize;
                *nlband = (*nlband).max(diff);
              }
            }
          }
        }
      }
    }
  }

  *nband = *nlband + *nlband + 1;
  *nrow = *nlband + *nlband + *nlband + 1;

  println!("Lower bandwidth = {}", *nlband);
  println!("Total bandwidth = {}", *nband);
  println!("NROW  = {}", *nrow);
  if MAXROW < *nrow {
    println!("SETBAN - NROW is too large!");
    println!("The maximum allowed is {}", MAXROW);
    panic!("NROW exceeds maximum allowed");
  }
}

//******************************************************
fn setbas(
  node: &[Vec<usize>],
  xc: &[f64],
  yc: &[f64],
  phi: &mut [Vec<Vec<Vec<f64>>>],
  psi: &mut [Vec<Vec<f64>>],
  xm: &[Vec<f64>],
  ym: &[Vec<f64>],
) {
  //? setbas() computes the basis functions at each integration point.

  let mut bb: f64 = 0.0;
  let mut bx: f64 = 0.0;
  let mut by: f64 = 0.0;

  for it in 0..NELEMN {
    for j in 0..NQUAD {
      let x = xm[it][j];
      let y = ym[it][j];
      for iq in 0..6 {
        psi[it][j][iq] = bsp(x, y, it, iq, 1, node, xc, yc);
        qbf(x, y, it, iq, &mut bb, &mut bx, &mut by, node, xc, yc);
        phi[it][j][iq][0] = bb;
        phi[it][j][iq][1] = bx;
        phi[it][j][iq][2] = by;
      }
    }
  }
}

//******************************************************
fn setgrd(
  indx: &mut [Vec<i32>],
  insc: &mut [i32],
  isotri: &mut [i32],
  iwrite: usize,
  _long: &mut bool,
  neqn: &mut usize,
  node: &mut [Vec<usize>],
) {
  //? setgrd() sets up the grid for the problem..

  //?  Determine whether region is long or skinny.
  //?  This will determine how we number the nodes and elements.

  if NY < NX {
    *_long = true;
    println!("Using vertical ordering.");
  } else {
    *_long = false;
    println!("Using horizontal ordering.");
  }
  //?  Set parameters for Taylor Hood element
  println!(" ");
  println!("SETGRD: Taylor Hood element");
  //?  Construct grid coordinates, elements, and ordering of unknowns
  *neqn = 0;
  let mut ielemn = 0;

  for ip in 0..NP {
    let ic;
    let jc;

    if *_long {
      ic = (ip / MY) + 1;
      jc = (ip % MY) + 1;
    } else {
      ic = (ip % MX) + 1;
      jc = (ip / MX) + 1;
    }

    let icnt = ic % 2;
    let jcnt = jc % 2;

    //?  If both the row count and the column count are odd,
    //?  and we're not in the last row or top column,
    //?  then we can define two new triangular elements based at the node.
    //?  For horizontal ordering,
    //?  given the following arrangement of nodes, for instance:
    //?    21 22 23 24 25
    //?    16 17 18 19 20
    //?    11 12 13 14 15
    //?    06 07 08 09 10
    //?    01 02 03 04 05
    //?  when we arrive at node 13, we will define
    //?  element 7: (13, 23, 25, 18, 24, 19)
    //?  element 8: (13, 25, 15, 19, 20, 14)
    //?  For vertical ordering,
    //?  given the following arrangement of nodes, for instance:
    //?    05 10 15 20 25
    //?    04 09 14 19 24
    //?    03 08 13 18 23
    //?    02 07 12 17 22
    //?    01 06 11 16 21
    //?  when we arrive at node 13, we will define
    //?  element 7: (13, 25, 23, 19, 24, 18)
    //?  element 8: (13, 15, 25, 14, 20, 19)
    if icnt == 1 && jcnt == 1 && ic != MX && jc != MY {
      if *_long {
        let ip1 = ip + MY;
        let ip2 = ip + MY + MY;

        // Primer elemento
        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2 + 2;
        node[ielemn - 1][2] = ip2;
        node[ielemn - 1][3] = ip1 + 1;
        node[ielemn - 1][4] = ip2 + 1;
        node[ielemn - 1][5] = ip1;
        isotri[ielemn - 1] = 0;

        // Segundo elemento
        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip + 2;
        node[ielemn - 1][2] = ip2 + 2;
        node[ielemn - 1][3] = ip + 1;
        node[ielemn - 1][4] = ip1 + 2;
        node[ielemn - 1][5] = ip1 + 1;
        isotri[ielemn - 1] = 0;
      } else {
        // .NOT. long
        let ip1 = ip + MX;
        let ip2 = ip + MX + MX;

        // Primer elemento
        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2;
        node[ielemn - 1][2] = ip2 + 2;
        node[ielemn - 1][3] = ip1;
        node[ielemn - 1][4] = ip2 + 1;
        node[ielemn - 1][5] = ip1 + 1;
        isotri[ielemn - 1] = 0;

        // Segundo elemento
        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2 + 2;
        node[ielemn - 1][2] = ip + 2;
        node[ielemn - 1][3] = ip1 + 1;
        node[ielemn - 1][4] = ip1 + 2;
        node[ielemn - 1][5] = ip + 1;
        isotri[ielemn - 1] = 0;
      }
    }

    //?  Consider whether velocity unknowns should be associated with this node.
    //?  If we are in column 1, horizontal velocities are specified, and
    //?  vertical velocities are zero.
    if ic == 1 && jc > 1 && jc < MY {
      indx[ip][0] = -2;
      indx[ip][1] = -1;
    }
    //?  If we are in column MX, horizontal velocities are unknown, and
    //?  vertical velocities are zero.
    else if ic == MX && jc > 1 && jc < MY {
      *neqn += 1;
      indx[ip][0] = *neqn as i32 - 1;
      indx[ip][1] = -1;
    }
    //?  Otherwise, if we are in row 1 or row MY, both horizontal and
    //?  vertical velocities are zero.
    else if jc == 1 || jc == MY {
      indx[ip][0] = -1;
      indx[ip][1] = -1;
    }
    //?  Otherwise, we are at an interior node
    else {
      *neqn += 2;
      indx[ip][0] = *neqn as i32 - 2;
      indx[ip][1] = *neqn as i32 - 1;
    }
    //?  Consider whether a pressure unknown should be associated with this node.
    //?  The answer is yes if both nodes are odd.
    if jcnt == 1 && icnt == 1 {
      *neqn += 1;
      insc[ip] = *neqn as i32 - 1;
    } else {
      insc[ip] = -1;
    }
  }
  //?  If debugging is requested, print out data.
  if 2 <= iwrite {
    println!(" ");
    println!("    I      INDX 1 & 2, INSC");
    println!(" ");
    //do i = 1,np
    for i in 0..NP {
      //write (*,'(2xi6,2x,i6,2x,i6,2x,i6)') i,indx(i,1:2),insc(i)
      let dsp = |v: i32| if v >= 0 { v + 1 } else if v == -2 { 0 } else { v };
      println!(
        "{:6} {:6} {:6} {:6}",
        i + 1,
        dsp(indx[i][0]),
        dsp(indx[i][1]),
        dsp(insc[i])
      );
    }
    println!(" ");
    println!("    IT    NODE(IT,1:6)");
    println!(" ");
    for it in 0..NELEMN {
      println!(
        "{:6} {:6} {:6} {:6} {:6} {:6} {:6}",
        it + 1,
        node[it][0] + 1,
        node[it][1] + 1,
        node[it][2] + 1,
        node[it][3] + 1,
        node[it][4] + 1,
        node[it][5] + 1
      );
    }
  }

  println!("Number of unknowns = {}", neqn);
  if MAXEQN < *neqn {
    println!("SETGRD - Too many unknowns!");
    println!("The maximum allowed is MAXEQN = {}", MAXEQN);
    println!("This problem requires NEQN = {}", neqn);
    std::process::exit(1);
  }
}

//******************************************************
fn setlin(
  iline: &mut [i32],
  indx: &[Vec<i32>],
  iwrite: usize,
  _long: &bool,
  nodex0: &mut usize,
  xlngth: f64,
) {
  //? setlin() gets the unknown indices along the profile line.
  //?  Determine the number of a node on the profile line
  let itemp: usize = ((18.0_f64 * (NX as f64 - 1.0)) / xlngth).round() as usize;
  *nodex0 = if *_long {
    itemp * (2 * NY - 1)
  } else {
    itemp
  };

  for i in 0..MY {
    let ip = if *_long {
      *nodex0 + i
    } else {
      *nodex0 + MX * i
    };
    iline[i] = indx[ip][0];
  }

  if 1 <= iwrite {
    println!(" ");
    println!("SETLIN: unknown numbers along line:");
    println!(" ");
    for i in 0..MY {
      std::print!("{:>5}", if iline[i] >= 0 { iline[i] + 1 } else { 0 });
      if (i + 1) % 15 == 0 && (i + 1) < MY {
        std::println!("");
      }
    }
    println!(" ");
  }
}

//******************************************************
fn setqud(
  area: &mut [f64],
  node: &[Vec<usize>],
  xc: &[f64],
  xm: &mut [Vec<f64>],
  yc: &[f64],
  ym: &mut [Vec<f64>],
) {
  //? setqud() sets midpoint quadrature rule information.
  for it in 0..NELEMN {
    let ip1 = node[it][0];
    let ip2 = node[it][1];
    let ip3 = node[it][2];
    let x1 = xc[ip1];
    let x2 = xc[ip2];
    let x3 = xc[ip3];
    let y1 = yc[ip1];
    let y2 = yc[ip2];
    let y3 = yc[ip3];
    xm[it][0] = 0.5 * (x1 + x2);
    xm[it][1] = 0.5 * (x2 + x3);
    xm[it][2] = 0.5 * (x3 + x1);
    ym[it][0] = 0.5 * (y1 + y2);
    ym[it][1] = 0.5 * (y2 + y3);
    ym[it][2] = 0.5 * (y3 + y1);
    area[it] = 0.5 * ((y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)).abs();
  }
}

////******************************************************
fn setxy(
  iwrite: usize,
  _long: &bool,
  xc: &mut [f64],
  xlngth: f64,
  yc: &mut [f64],
  ylngth: f64,
) {
  //? setxy() sets the X, Y coordinates of grid points.
  //?  Construct grid coordinates
  let mut ic;
  let mut jc;

  for ip in 0..NP {
    if *_long {
      ic = ip / MY;
      jc = ip % MY;
    } else {
      ic = ip % MX;
      jc = ip / MX;
    }

    xc[ip] = ic as f64 * xlngth / ((2 * NX - 2) as f64);
    yc[ip] = jc as f64 * ylngth / ((2 * NY - 2) as f64);
  }

  if 2 <= iwrite {
    println!(" ");
    println!("    I      XC           YC");
    println!(" ");
    for i in 0..NP {
      println!("{:>5} {:>12.5} {:>12.5}", i + 1, xc[i], yc[i]);
    }
  }
}

//******************************************************
fn ubdry(y: f64, para: f64) -> f64 {
  //? ubdry() sets the parabolic inflow in terms of the value of the parameter.
  4.0 * para * y * (3.0 - y) / 9.0
}

//******************************************************
fn uval(
  g: &[f64],
  indx: &[Vec<i32>],
  iquad: usize,
  it: usize,
  node: &[Vec<usize>],
  para: f64,
  phi: &[Vec<Vec<Vec<f64>>>],
  un: &mut [f64],
  unx: &mut [f64],
  uny: &mut [f64],
  yc: &[f64],
) {
  //? uval() evaluates the velocities at a given point in a //particular triangle.

  un[0] = 0.0;
  un[1] = 0.0;
  uny[0] = 0.0;
  uny[1] = 0.0;
  unx[0] = 0.0;
  unx[1] = 0.0;

  for iq in 0..NNODES {
    let ip = node[it][iq];
    let bb = phi[it][iquad][iq][0];
    let bx = phi[it][iquad][iq][1];
    let by = phi[it][iquad][iq][2];

    for iuk in 0..2 {
      let iun = indx[ip][iuk];

      if iun >= 0 {
        un[iuk] += bb * g[iun as usize];
        unx[iuk] += bx * g[iun as usize];
        uny[iuk] += by * g[iun as usize];
      } else if iun == -2 {
        let ubc: f64 = ubdry(yc[ip], para);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }
}

fn uv_table(
  f: &[f64],
  indx: &[Vec<i32>],
  para: f64,
  yc: &[f64],
  filename: &str,
) -> std::io::Result<()> {
  use std::io::Write;
  let mut file = std::fs::File::create(filename)?;
  for ip in 0..NP {
    let k = indx[ip][0];
    let uval = if k == -1 {
      0.0
    } else if k == -2 {
      ubdry(yc[ip], para)
    } else {
      f[k as usize]
    };
    let k = indx[ip][1];
    let vval = if k == -1 { 0.0 } else { f[k as usize] };
    writeln!(file, "  {:>14.6}  {:>14.6}", uval, vval)?;
  }
  Ok(())
}

fn xy_table(xc: &[f64], yc: &[f64], filename: &str) -> std::io::Result<()> {
  use std::io::Write;
  let mut file = std::fs::File::create(filename)?;
  for ip in 0..NP {
    writeln!(file, "  {:>14.6}  {:>14.6}", xc[ip], yc[ip])?;
  }
  Ok(())
}

fn xy_plot3d(
  long: bool,
  xc: &[f64],
  yc: &[f64],
  filename: &str,
) -> std::io::Result<()> {
  use std::io::Write;
  let mut file = std::fs::File::create(filename)?;
  if long {
    writeln!(file, "{} {}", MX, MY)?;
    for i in 0..MY {
      for j in 0..MX {
        let ip = j * MY + i;
        writeln!(file, "{:>14.6}", xc[ip])?;
      }
    }
    for i in 0..MY {
      for j in 0..MX {
        let ip = j * MY + i;
        writeln!(file, "{:>14.6}", yc[ip])?;
      }
    }
  } else {
    writeln!(file, "{} {}", MY, MX)?;
    for j in 0..MX {
      for i in 0..MY {
        let ip = j * MY + i;
        writeln!(file, "{:>14.6}", xc[ip])?;
      }
    }
    for j in 0..MX {
      for i in 0..MY {
        let ip = j * MY + i;
        writeln!(file, "{:>14.6}", yc[ip])?;
      }
    }
  }
  Ok(())
}