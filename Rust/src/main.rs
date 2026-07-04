use chrono::Local;

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

////******************************************************
fn timestamp() {
    let now = Local::now();
    let formatted_date = now.format("%d-%m-%Y").to_string();
    println!("Current date: {}", formatted_date);
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

    let mut iline: Vec<i32> = vec![0; MY];
    let mut indx: Vec<Vec<i32>> = vec![vec![0; 2]; NP];
    let mut insc: Vec<i32> = vec![0; NP];

    let mut ipivot: Vec<i32> = vec![0; MAXEQN];
    let mut isotri: Vec<i32> = vec![0; NELEMN];
    let mut iter: i32 = 0;

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
    let mut temp: f64;
    let mut test: f64;

    let mut ui = vec![0.0; MY];
    let mut unew = vec![0.0; MY];
    let mut xc = vec![0.0; NP];

    let mut xm = vec![vec![0.0; NQUAD]; NELEMN];
    let mut yc = vec![0.0; NP];

    let mut ym = vec![vec![0.0; NQUAD]; NELEMN];

    timestamp();
    println!(" ");
    println!("channel():");
    println!("  C++ version");
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

    let iounit: usize = 2;
    let ivunit: usize = 4;
    let iwrite: usize = 10;
    let ixunit: usize = 3;
    let mut maxnew: usize = 10;
    let maxsec: usize = 8;
    let npara: usize = 1;
    let mut numnew: usize = 0;
    let mut numsec: usize = 0;
    let mut reynld: f64 = 1.0;
    let rjpnew: f64 = 0.0;
    let mut tolnew: f64 = 1.0E-04;
    let mut tolsec: f64 = 1.0E-06;
    let mut xlngth: f64 = 10.0;
    let mut ylngth: f64 = 3.0;

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
        &mut indx,
        &mut insc,
        &mut nband,
        &mut nlband,
        &mut node,
        &mut nrow,
    );
    //?  Record variable numbers along profile sampling line.
    setlin(
        &mut iline,
        &mut indx,
        iwrite,
        &mut _long,
        &mut nodex0,
        &mut xlngth,
    );
    //?  Set the coordinates of grid points.
    setxy(
        iwrite,
        &mut _long,
        &mut xc,
        &mut xlngth,
        &mut yc,
        &mut ylngth,
    );
    ////?  Set quadrature points
    setqud(&mut area, &mut node, &mut xc, &mut xm, &mut yc, &mut ym);
    //?  Evaluate basis functions at quadrature points
    setbas(
        &mut node, &mut xc, &mut yc, &mut phi, &mut psi, &mut xm, &mut ym,
    );
    //?  NSTOKE now solves the Navier Stokes problem for an inflow
    ////?  parameter of 1.0.
    para = 1.0;
    println!(" ");
    println!("Solve Navier Stokes problem with parameter = {}", para);
    println!("for profile at x = {}", xc[nodex0 - 1]);
    for i in 0..neqn {
        g[i] = 1.0;
    }

    nstoke(
        &mut a,
        &mut area,
        &mut f,
        &mut g,
        &mut indx,
        &mut insc,
        &mut ipivot,
        iwrite,
        &mut maxnew,
        neqn,
        nlband,
        &mut node,
        &mut nrow,
        &mut numnew,
        &mut para,
        &mut phi,
        &mut psi,
        &mut reynld,
        &mut tolnew,
        &mut yc,
    );
    //?  RESID computes the residual at the given solution
    if 1 <= iwrite {
        resid(
            &mut area, &f, &indx, &insc, iwrite, &mut node, &mut phi, &mut psi, &mut res, reynld,
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
        &mut gr, &iline, &indx, iwrite, MY, &node, nodex0, para, &mut r, &ui, &xc, &yc,
    );
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
            &mut area,
            &mut f,
            &mut g,
            &mut indx,
            &mut insc,
            &mut ipivot,
            iwrite,
            &mut maxnew,
            neqn,
            nlband,
            &mut node,
            &mut nrow,
            &mut numnew,
            &mut para,
            &mut phi,
            &mut psi,
            &mut reynld,
            &mut tolnew,
            &mut yc,
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
}

//******************************************************
fn bsp(
    xq: f64,
    yq: f64,
    it: usize,
    iq: usize,
    id: usize,
    node: &Vec<Vec<usize>>,
    xc: &Vec<f64>,
    yc: &Vec<f64>,
) -> f64 {
    //? bsp() evaluates the linear basis functions associated with pressure.
    let iq1 = iq;
    let iq2: usize = iq % 3;
    let iq3: usize = (iq + 1) % 3;
    let i1: usize = node[it][iq1] - 1;
    let i2: usize = node[it][iq2] - 1;
    let i3: usize = node[it][iq3] - 1;

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
    abd: &mut Vec<Vec<f64>>,
    n: usize,
    ml: usize,
    mu: usize,
    ipvt: &mut Vec<i32>,
    info: &mut usize,
) {
    let m = ml + mu + 1;
    *info = 0;

    let j0 = mu + 1;
    let j1 = if n < m { n } else { m };

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
            for i in 1..=ml {
                abd[i - 1][jz - 1] = 0.0;
            }
        }

        let lm = std::cmp::min(ml, n - k);
        let l = idamax(lm + 1, &abd[m - 1][k - 1..], 1) + m - 1;
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
    abd: &mut Vec<Vec<f64>>,
    _lda: usize,
    n: usize,
    ml: usize,
    mu: usize,
    ipvt: &Vec<i32>,
    b: &mut Vec<f64>,
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
fn getg(f: &Vec<f64>, iline: &Vec<i32>, my: usize, _neqn: usize, u: &mut Vec<f64>) {
    for j in 0..my {
        let k = iline[j];
        if 0 < k {
            u[j] = f[(k - 1) as usize];
        } else {
            u[j] = 0.0;
        }
    }
}
//******************************************************
fn gram(
    gr: &mut Vec<Vec<f64>>,
    iline: &Vec<i32>,
    indx: &Vec<Vec<i32>>,
    iwrite: usize,
    node: &Vec<Vec<usize>>,
    nodex0: usize,
    para: f64,
    r: &mut Vec<f64>,
    ui: &Vec<f64>,
    xc: &Vec<f64>,
    yc: &Vec<f64>,
) {
    let wt: [f64; 3] = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
    let yq: [f64; 3] = [-0.7745966692, 0.0, 0.7745966692];

    for i in 0..MY {
        r[i] = 0.0;
        for j in 0..MY {
            gr[i][j] = 0.0;
        }
    }

    let xzero = xc[nodex0 - 1];

    for it in 0..NELEMN {
        let k = node[it][0];
        let kk = node[it][1];

        if 1.0E-04 < (xc[k - 1] - xzero).abs() {
            continue;
        }
        if 1.0E-04 < (xc[kk - 1] - xzero).abs() {
            continue;
        }

        for iquad in 0..3 {
            let bma2 = (yc[kk - 1] - yc[k - 1]) / 2.0;
            let ar = bma2 * wt[iquad];
            let x = xzero;
            let y = yc[k - 1] + bma2 * (yq[iquad] + 1.0);
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
                let iun = indx[ip - 1][0];
                if 0 < iun {
                    let ii = igetl(iun, iline, MY);
                    uiqdpt += bb * ui[ii];
                } else if iun < 0 {
                    let ubc = ubdry(yc[ip - 1], para);
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
                    let i = indx[ip - 1][0];
                    if i <= 0 {
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
                            let j = indx[ipp - 1][0];
                            if j != 0 {
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
fn igetl(k: i32, iline: &Vec<i32>, my: usize) -> usize {
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
    a: &mut Vec<Vec<f64>>,
    area: &mut Vec<f64>,
    f: &mut Vec<f64>,
    g: &mut Vec<f64>,
    indx: &mut Vec<Vec<i32>>,
    insc: &mut Vec<i32>,
    ipivot: &mut Vec<i32>,
    neqn: usize,
    nlband: usize,
    node: &mut Vec<Vec<usize>>,
    nrow: usize,
    para1: f64,
    para2: f64,
    phi: &mut Vec<Vec<Vec<Vec<f64>>>>,
    psi: &mut Vec<Vec<Vec<f64>>>,
    reynld: f64,
    yc: &mut Vec<f64>,
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
                g, indx, iquad, it, node, para1, phi, &mut un, &mut unx, &mut uny, yc,
            );
            //?  For each basis function,
            for iq in 0..NNODES {
                ip = node[it][iq] - 1;
                bb = phi[it][iquad][iq][0];
                bx = phi[it][iquad][iq][1];
                by = phi[it][iquad][iq][2];
                bbl = psi[it][iquad][iq];
                ihor = indx[ip][0] - 1;
                iver = indx[ip][1] - 1;
                iprs = insc[ip] - 1;

                if 0 <= ihor {
                    f[ihor as usize] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
                }

                if 0 <= iver {
                    f[iver as usize] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
                }
                //?  For another basis function,
                //do iqq = 1, nnodes
                for iqq in 0..NNODES {
                    ipp = node[it][iqq] - 1;
                    bbb = phi[it][iquad][iqq][0];
                    bbx = phi[it][iquad][iqq][1];
                    bby = phi[it][iquad][iqq][2];
                    bbbl = psi[it][iquad][iqq];
                    ju = indx[ipp][0] - 1;
                    jv = indx[ipp][1] - 1;
                    jp = insc[ipp] - 1;
                    //?  Horizontal velocity variable
                    if 0 <= ju {
                        if 0 <= ihor {
                            iuse = ihor - ju + ioff;
                            a[iuse as usize][ju as usize] += ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
                        }

                        if 0 <= iver {
                            iuse = iver - ju + ioff;
                            a[iuse as usize][ju as usize] += ar * bb * bbb * unx[1];
                        }

                        if 0 <= iprs {
                            iuse = iprs - ju + ioff;
                            a[iuse as usize][ju as usize] += ar * bbx * bbl;
                        }
                    } else if ju <= 0 {
                        uu = ubdry(yc[ipp - 1], para2);
                        if 0 <= ihor {
                            f[ihor as usize] -= ar
                                * uu
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
                        }

                        if 0 <= iver {
                            f[iver as usize] -= ar * uu * bb * bbb * unx[1];
                        }

                        if 0 <= iprs {
                            f[iprs as usize] -= ar * uu * bbx * bbl;
                        }
                    }
                    //?  Vertical velocity variable
                    if 0 <= jv {
                        if 0 <= ihor {
                            iuse = ihor - jv + ioff;
                            a[iuse as usize][jv as usize] += ar * bb * bbb * uny[0];
                        }

                        if 0 <= iver {
                            iuse = iver - jv + ioff;
                            a[iuse as usize][jv as usize] += ar
                                * (visc * (by * bby + bx * bbx)
                                    + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
                        }

                        if 0 <= iprs {
                            iuse = iprs - jv + ioff;
                            a[iuse as usize][jv as usize] += ar * bby * bbl;
                        }
                    }
                    //?  Pressure variable
                    if 0 <= jp {
                        if 0 <= ihor {
                            iuse = ihor - jp + ioff;
                            a[iuse as usize][jp as usize] -= ar * bx * bbbl;
                        }
                        if 0 < iver {
                            iuse = iver - jp + ioff;
                            a[iuse as usize][jp as usize] -= ar * by * bbbl;
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
    a: &mut Vec<Vec<f64>>,
    area: &mut Vec<f64>,
    f: &mut Vec<f64>,
    g: &mut Vec<f64>,
    indx: &mut Vec<Vec<i32>>,
    insc: &mut Vec<i32>,
    ipivot: &mut Vec<i32>,
    iwrite: usize,
    maxnew: &mut usize,
    neqn: usize,
    nlband: usize,
    node: &mut Vec<Vec<usize>>,
    nrow: &mut usize,
    numnew: &mut usize,
    para: &mut f64,
    phi: &mut Vec<Vec<Vec<Vec<f64>>>>,
    psi: &mut Vec<Vec<Vec<f64>>>,
    reynld: &mut f64,
    tolnew: &mut f64,
    yc: &mut Vec<f64>,
) {
    //?! nstoke() solves the Navier Stokes equation using //Taylor-Hood elements.
    //?  The G array contains the previous iterate.
    //?  The F array contains the right hand side initially and //then the current iterate.
    for iter in 0..*maxnew {
        *numnew += 1;

        linsys(
            a, area, f, g, indx, insc, ipivot, neqn, nlband, node, *nrow, *para, *para, phi, psi,
            *reynld, yc,
        );
        //? Check for convergence
        for i in 0..neqn {
            g[i] = g[i] - f[i];
        }
        let diff = g[idamax(neqn, g, 1)].abs();

        if 1 <= iwrite {
            println!("NSTOKE iteration {} Mnorm = {}", iter, diff);
        }

        for i in 0..neqn {
            g[i] = f[i];
        }

        if diff <= *tolnew {
            println!("Navier Stokes iteration converged in {}iterations.", iter);
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
    node: &Vec<Vec<usize>>,
    xc: &Vec<f64>,
    yc: &Vec<f64>,
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
        in2 = _in % 3;
        in3 = (_in + 1) % 3;
        i1 = node[it][in1] - 1;
        i2 = node[it][in2] - 1;
        i3 = node[it][in3] - 1;
        d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
        t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
        *bb = t * (2.0 * t - 1.0);
        *bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
        *by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
    } else {
        inn = _in - 3;
        in1 = inn;
        in2 = inn % 3;
        in3 = (inn + 1) % 3;
        i1 = node[it][in1] - 1;
        i2 = node[it][in2] - 1;
        i3 = node[it][in3] - 1;
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

//subroutine resid (area,g,indx,insc,iwrite,nelemn,neqn,nnodes, &
//  node,np,nquad,para,phi,psi,res,reynld,yc)
//
//******************************************************
//?! resid() computes the residual.
//?  Discussion:
//?    The G array contains the current iterate.
//?    The RES array will contain the value of the residual.
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
//?  For each element,
//  do 90 it = 1, nelemn
//
//    ar = area(it) / 3.0
//?  and for each quadrature point in the element,
//    do 80 iquad = 1, nquad
//?  Evaluate velocities at quadrature point
//      call uval (g,indx,iquad,it,nelemn,neqn,nnodes,node, &
//        np,nquad,para,phi,un,unx,uny,yc)
//?  For each basis function,
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
//?  For another basis function,
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
//    println!(" ")
//    println!("RESIDUAL INFORMATION:")
//    println!(" ")
//    println!("Worst residual is number ',IMAX
//    println!("of magnitude ',RMAX
//    println!(" ")
//    println!("Number of "bad" residuals is ',IBAD,' out of ',NEQN
//    println!(" ")
//  end if
//
//  if ( 2 <= iwrite ) then
//    println!("Raw residuals:")
//    println!(" ")
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
fn setban(
    indx: &mut Vec<Vec<i32>>,
    insc: &mut Vec<i32>,
    nband: &mut usize,
    nlband: &mut usize,
    node: &mut Vec<Vec<usize>>,
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
                    insc[ip - 1]
                } else {
                    indx[ip - 1][iuk]
                };

                if i > 0 {
                    for iqq in 0..NNODES {
                        let ipp = node[it][iqq];
                        for iukk in 0..3 {
                            let j = if iukk == 2 {
                                insc[ipp - 1]
                            } else {
                                indx[ipp - 1][iukk]
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
    node: &mut Vec<Vec<usize>>,
    xc: &mut Vec<f64>,
    yc: &mut Vec<f64>,
    phi: &mut Vec<Vec<Vec<Vec<f64>>>>,
    psi: &mut Vec<Vec<Vec<f64>>>,
    xm: &mut Vec<Vec<f64>>,
    ym: &mut Vec<Vec<f64>>,
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
    indx: &mut Vec<Vec<i32>>,
    insc: &mut Vec<i32>,
    isotri: &mut Vec<i32>,
    iwrite: usize,
    _long: &mut bool,
    neqn: &mut usize,
    node: &mut Vec<Vec<usize>>,
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

    for ip in 1..=NP {
        let ic;
        let jc;

        if *_long {
            ic = ((ip - 1) / MY) + 1;
            jc = ((ip - 1) % MY) + 1;
        } else {
            ic = ((ip - 1) % MX) + 1;
            jc = ((ip - 1) / MX) + 1;
        }

        let icnt = ic % 2; // igual que mod(ic,2)
        let jcnt = jc % 2; // igual que mod(jc,2)
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
                let ip2 = ip + MY + MY; // o 2 * my

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
                let ip2 = ip + MX + MX; // o 2 * mx

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
            indx[ip - 1][0] = -1;
            indx[ip - 1][1] = 0;
        }
        //?  If we are in column MX, horizontal velocities are unknown, and
        //?  vertical velocities are zero.
        else if ic == MX && jc > 1 && jc < MY {
            *neqn += 1;
            indx[ip - 1][0] = *neqn as i32;
            indx[ip - 1][1] = 0;
        }
        //?  Otherwise, if we are in row 1 or row MY, both horizontal and
        //?  vertical velocities are zero.
        else if jc == 1 || jc == MY {
            indx[ip - 1][0] = 0;
            indx[ip - 1][1] = 0;
        }
        //?  Otherwise, we are at an interior node
        else {
            *neqn += 2;
            indx[ip - 1][0] = *neqn as i32 - 1;
            indx[ip - 1][1] = *neqn as i32;
        }
        //?  Consider whether a pressure unknown should be associated with this node.
        //?  The answer is yes if both nodes are odd.
        if jcnt == 1 && icnt == 1 {
            *neqn += 1;
            insc[ip - 1] = *neqn as i32;
        } else {
            insc[ip - 1] = 0;
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
            println!(
                "{:6} {:6} {:6} {:6}",
                i + 1,
                indx[i][0],
                indx[i][1],
                insc[i]
            );
        }
        println!(" ");
        println!("    IT    NODE(IT,1:6)");
        println!(" ");
        for it in 0..NELEMN {
            println!(
                "{:6} {:6} {:6} {:6} {:6} {:6} {:6}",
                it + 1,
                node[it][0],
                node[it][1],
                node[it][2],
                node[it][3],
                node[it][4],
                node[it][5]
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
    iline: &mut Vec<i32>,
    indx: &mut Vec<Vec<i32>>,
    iwrite: usize,
    _long: &mut bool,
    nodex0: &mut usize,
    xlngth: &mut f64,
) {
    //? setlin() gets the unknown indices along the profile line.
    //?  Determine the number of a node on the profile line
    let itemp: usize = ((18.0_f64 * (NX as f64 - 1.0)) / *xlngth).round() as usize;
    *nodex0 = if *_long {
        itemp * (2 * NY - 1) + 1
    } else {
        itemp + 1
    };

    for i in 0..MY {
        let ip = if *_long {
            *nodex0 + i
        } else {
            *nodex0 + MX * i
        };
        iline[i] = indx[ip - 1][0];
    }

    if 1 <= iwrite {
        println!(" ");
        println!("SETLIN: unknown numbers along line:");
        println!(" ");
        for i in 0..MY {
            std::print!("{:>5}", iline[i]);
            if (i + 1) % 15 == 0 && (i + 1) < MY {
                std::println!("");
            }
        }
        println!(" ");
    }
}

//******************************************************
fn setqud(
    area: &mut Vec<f64>,
    node: &mut Vec<Vec<usize>>,
    xc: &mut Vec<f64>,
    xm: &mut Vec<Vec<f64>>,
    yc: &mut Vec<f64>,
    ym: &mut Vec<Vec<f64>>,
) {
    //? setqud() sets midpoint quadrature rule information.
    for it in 0..NELEMN {
        let ip1 = node[it][0];
        let ip2 = node[it][2];
        let ip3 = node[it][3];
        let x1 = xc[ip1 - 1];
        let x2 = xc[ip2 - 1];
        let x3 = xc[ip3 - 1];
        let y1 = yc[ip1 - 1];
        let y2 = yc[ip2 - 1];
        let y3 = yc[ip3 - 1];
        xm[it][0] = 0.5 * (x1 + x2);
        xm[it][1] = 0.5 * (x2 + x3);
        xm[it][2] = 0.5 * (x3 + x1);
        ym[it][0] = 0.5 * (y1 + y2);
        ym[it][1] = 0.5 * (y2 + y3);
        ym[it][2] = 0.5 * (y3 + y1);
        area[it] =
            0.5 * ((y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)).abs();
    }
}

////******************************************************
fn setxy(
    iwrite: usize,
    _long: &mut bool,
    xc: &mut Vec<f64>,
    xlngth: &mut f64,
    yc: &mut Vec<f64>,
    ylngth: &mut f64,
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

        xc[ip] = ic as f64 * *xlngth / ((2 * NX - 2) as f64);
        yc[ip] = jc as f64 * *ylngth / ((2 * NY - 2) as f64);
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
    g: &mut Vec<f64>,
    indx: &Vec<Vec<i32>>,
    iquad: usize,
    it: usize,
    node: &Vec<Vec<usize>>,
    para: f64,
    phi: &Vec<Vec<Vec<Vec<f64>>>>,
    un: &mut Vec<f64>,
    unx: &mut Vec<f64>,
    uny: &mut Vec<f64>,
    yc: &Vec<f64>,
) {
    //? uval() evaluates the velocities at a given point in a //particular triangle.

    un[0] = 0.0;
    un[1] = 0.0;
    uny[0] = 0.0;
    uny[1] = 0.0;
    unx[0] = 0.0;
    unx[1] = 0.0;

    //do iq = 1, nnodes
    for iq in 0..NNODES {
        let mut ip = node[it][iq] - 1;
        let bb = phi[it][iquad][iq][1];
        let bx = phi[it][iquad][iq][2];
        let by = phi[it][iquad][iq][3];

        for iuk in 0..2 {
            let iun = indx[ip][iuk] - 1;

            if 0 <= iun {
                un[iuk] += bb * g[iun as usize];
                unx[iuk] += bx * g[iun as usize];
                uny[iuk] += by * g[iun as usize];
            } else if iun < 0 {
                ip = node[it][iq] - 1;
                let ubc: f64 = ubdry(yc[ip], para);
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
//?! uv_plot3d() creates a velocity file for use by PLOT3D.
//?  Given the following set of nodes:
//?    A  B  C
//?    D  E  F
//?    G  H  I
//?  the file will have the form:
//?    D, U(G), V(G), P
//?    D, U(H), V(H), P
//?    D, U(I), V(I), P
//?    D, U(D), V(D), P
//?    D, U(E), V(E), P
//?    D, U(F), V(F), P
//?    D, U(A), V(A), P
//?    D, U(B), V(B), P
//?    D, U(C), V(C), P
//?  Here both D and P are set to 1 for now, representing dummy values
//?  of density and pressure.
//?  Licensing:
//?    This code is distributed under the MIT license.
//?  Modified:
//?    20 January 2007
//?  Author:
//?    John Burkardt
//?  Parameters:
//?    Input, real ( kind = rk8 ) PARA, the value of the parameter.
//?    Workspace, real ( kind = rk8 ) PRESS(MX,MY), used to hold the
//?    computed pressures.
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
//?  If NY < NX, then nodes with a constant Y value are numbered consecutively.
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
//?  If NX < NY, then nodes with a constant X value are numbered consecutively.
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
//  println!("UV_PLOT3D wrote data set ',iset,' to file.")
//
//  return
//end
//subroutine uv_table ( f, indx, ivunit, neqn, np, para, yc )
//
//******************************************************
//?! uv_table() creates a velocity table file.
//?  Licensing:
//?    This code is distributed under the MIT license.
//?  Modified:
//?    28 February 2006
//?  Author:
//?    John Burkardt
//?  Parameters:
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
//?! xy_plot3d() creates a grid file for use by PLOT3D.
//?  Given the following set of nodes:
//?    A  B  C
//?    D  E  F
//?    G  H  I
//?  the file will have the form:
//?    X(G), X(H), X(I), X(D), X(E), X(F), X(A), X(B), X(C),
//?    Y(G), Y(H), Y(I), Y(D), Y(E), Y(F), Y(A), Y(B), Y(C).
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
//?  If NY < NX, then nodes with a constant Y value are numbered consecutively.
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
//?  If NX < NY, then nodes with a constant X value are numbered consecutively.
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
//  println!("XYDUMP wrote data set ',iset,' to file.")
//
//  return
//end
//subroutine xy_table ( ixunit, np, xc, yc )
//
//******************************************************
//?! xy_table() creates an XY table file.
//?  Licensing:
//?    This code is distributed under the MIT license.
//?  Modified:
//?    28 February 2006
//?  Author:
//?    John Burkardt
//?  Parameters:
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
