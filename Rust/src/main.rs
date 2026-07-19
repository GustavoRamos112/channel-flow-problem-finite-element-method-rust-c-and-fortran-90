use chrono::Datelike;
use chrono::Local;
use std::fs;
use std::fs::OpenOptions;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

mod program;
use program::flow_struct::Flow_struct;

fn timestamp() {
  let now = Local::now();
  println!(
    "Fecha actual: {:02}/{:02}/{:04}",
    now.day(),
    now.month(),
    now.year()
  );
}

fn main() {
  let mut flow: Flow_struct = Flow_struct::new(21, 7);
  let inicio: std::time::Instant = std::time::Instant::now();
  let dir_dat = Path::new(flow.data_dir);

  timestamp();
  println!(" ");
  println!("NX = {}", flow.nx);
  println!("NY = {}", flow.ny);
  println!("Number of elements = {}", flow.nelemn);
  println!("Reynolds number = {}", flow.reynld);
  println!("Secant tolerance = {}", flow.tolsec);
  println!("Newton tolerance = {}", flow.tolnew);
  println!(" ");
  //?  SETGRD constructs grid, numbers unknowns, calculates areas, and points for midpoint quadrature rule.
  setgrd(&mut flow);
  flow.f = vec![0.0; flow.neqn];
  flow.g = vec![1.0; flow.neqn];
  flow.res = vec![0.0; flow.neqn];
  //?  Compute the bandwidth
  setban(&mut flow);
  flow.a = vec![0.0; flow.nrow * flow.neqn];
  //?  Record variable numbers along profile sampling line.
  setlin(&mut flow);
  //?  Set the coordinates of grid points.
  setxy(&mut flow);
  //?  Set quadrature points
  setqud(&mut flow);
  //?  Evaluate basis functions at quadrature points
  setbas(&mut flow);
  //?  NSTOKE now solves the Navier Stokes problem for an inflow
  //?  parameter of 1.0.
  flow.para1 = 1.0;
  flow.para2 = 1.0;
  println!(" ");
  println!(
    "Solve Navier Stokes problem with parameter = {}",
    flow.para1
  );
  println!("for profile at x = {}", flow.xc[flow.nodex0]);
  //flow.g.fill(1.0); ver linea 29
  nstoke(&mut flow);
  //r!  RESID computes the residual at the given solution
  if 1 <= flow.iwrite {
    resid(&mut flow);
  }
  //r!  GETG computes the internal velocity profile at X = XC//(NODEX0), which will
  //r!  be used to measure the goodness-of-fit of the later solutions.
  getg(&flow.f, &flow.iline, flow.my, flow.neqn, &mut flow.ui);

  if 1 <= flow.iwrite {
    println!("\nU profile:\n");
    for i in 0..flow.my {
      print!("{}, ", flow.ui[i]);
    }
    println!(" ");
  }
  //r!  GRAM generates the Gram matrix GR && the vector
  //r!  R = line integral of ui*phi
  gram(&mut flow);
  if flow.save_data {
    println!("Writing graphics data to file {}", flow.fileg);
    if !dir_dat.exists() {
      std::fs::create_dir_all(dir_dat).unwrap();
    }
    flow.rjpnew = 0.0;
    //r!  GDUMP dumps information for graphics display by DISPLAY.
    gdump(&mut flow).expect("Error al ejecutar gdump");
    //r!  Write the XY data to a file.
    xy_plot3d(&flow).expect("Erro al ejecutar xy_plot3d");
    //r!  Write the velocity data to a file.
    uv_plot3d(&flow);
  } else {
    if !dir_dat.exists() {
      std::fs::create_dir_all(dir_dat).unwrap();
    }
    xy_table(&flow).expect("Error al ejecutar xy_table");
    uv_table(&flow).expect("Error al ejecutar uv_table");
  }
  //r!  Destroy information about true solution
  flow.f.fill(0.0);
  flow.g.fill(0.0);
  //r!  Secant iteration loop
  flow.aold = 0.0;
  flow.rjpold = 0.0;
  flow.anew = 0.1;
  let mut temp;
  //  do iter = 1, maxsec
  for iter in 1..=flow.maxsec {
    flow.numsec += 1;
    println!(" ");
    println!("Secant iteration {}", iter);
    //r!  Solve for unew at new value of parameter anew
    println!(" ");
    println!(
      "Solving Navier Stokes problem for parameter = {}",
      flow.anew
    );
    //r!  Use solution F at previous value of parameter for starting point.
    flow.g.copy_from_slice(&flow.f);
    flow.para1 = flow.anew;
    flow.para2 = flow.anew;

    nstoke(&mut flow);
    //r!  Get velocity profile
    getg(&flow.f, &flow.iline, flow.my, flow.neqn, &mut flow.unew);
    if 1 <= flow.iwrite {
      println!(" ");
      println!("Velocity profile:");
      println!(" ");
      for i in 0..flow.my {
        print!("{}, ", flow.unew[i]);
      }
    }
    //r!  Solve linear system for du/da
    flow.para1 = flow.anew;
    flow.para2 = 1.0;
    linsys(
      &mut flow.a,
      &flow.area,
      &mut flow.g,
      &flow.f,
      &flow.indx,
      &flow.insc,
      &mut flow.ipivot,
      flow.maxrow,
      flow.nelemn,
      flow.neqn,
      flow.nlband,
      flow.nnodes,
      &flow.node,
      flow.np,
      flow.nquad,
      flow.nrow,
      flow.para1,
      flow.para2,
      &flow.phi,
      &flow.psi,
      flow.reynld,
      &flow.yc,
    );
    //r!  Output in DCDA
    getg(&flow.g, &flow.iline, flow.my, flow.neqn, &mut flow.dcda);

    if 2 <= flow.iwrite {
      println!("\n ");
      println!("Sensitivities:");
      println!(" ");
      for i in 0..flow.my {
        print!("{}, ", flow.dcda[i]);
      }
    }
    //r!  Evaluate J prime at current value of parameter where J is
    //r!  functional to be minimized.
    //r!  JPRIME = 2.0 * DCDA(I) * (GR(I,J)*UNEW(J)-R(I))
    flow.rjpnew = 0.0;
    //do i = 1, my
    for i in 0..flow.my {
      temp = -flow.r[i];
      //do j = 1, my
      for j in 0..flow.my {
        temp += flow.gr[i + j * flow.my] * flow.unew[j];
      }
      flow.rjpnew += 2.0 * flow.dcda[i] * temp;
    }
    println!("\n ");
    println!("Parameter  = {}, J prime = {}", flow.anew, flow.rjpnew);
    //r!  Dump information for graphics
    if flow.save_data {
      flow.para1 = flow.anew;
      gdump(&mut flow).expect("Error al ejecutar gdump");
    }
    //r!  Update the estimate of the parameter using the secant step
    if iter == 1 {
      flow.a2 = 0.5;
    } else {
      flow.a2 = flow.aold - flow.rjpold * (flow.anew - flow.aold) / (flow.rjpnew - flow.rjpold);
    }

    flow.aold = flow.anew;
    flow.anew = flow.a2;
    flow.rjpold = flow.rjpnew;
    flow.test = (flow.anew - flow.aold).abs() / (flow.anew).abs();

    println!("New value of parameter = {}", flow.anew);
    println!("Convergence test = {}", flow.test);
    if (flow.anew - flow.aold).abs() < flow.anew.abs() * flow.tolsec {
      flow.converged = true;
      break;
    }
  }
  if flow.converged {
    println!("Secant iteration converged.");
  } else {
    println!("Secant iteration failed to converge.");
  }
  // 40 continue

  println!("Number of secant steps = {}", flow.numsec);
  println!("Number of Newton steps = {}", flow.numnew);
  //?  Terminate.
  println!(" ");
  println!("CHANNEL:");
  println!("  Normal of execution.");
  println!(" ");

  timestamp();
  let duracion = inicio.elapsed();
  println!("Tiempo transcurrido: {:?}", duracion);
  if flow.save_times {
    let mut file = OpenOptions::new()
      .create(true)
      .append(true)
      .open("times.txt")
      .unwrap();
    write!(file, "{}\n", duracion.as_secs_f64()).unwrap();
    file.flush().unwrap();
  }
}

fn bsp(
  xq: f64,
  yq: f64,
  it: usize,
  iq: usize,
  id: usize,
  nelemn: usize,
  node: &Vec<usize>,
  xc: &Vec<f64>,
  yc: &Vec<f64>,
) -> f64 {
  //r!! bsp() evaluates the linear basis functions associated with pressure.
  let iq1 = iq;
  //* let iq2 = (iq % 3) + 1; -> iq estaba en 1-index
  //* pero ahora esta en 0-index, asi que es necesario sumarle 1
  let iq2 = (iq + 1) % 3;
  //* let iq3 = ((iq + 1) % 3) + 1;
  let iq3 = (iq + 2) % 3;
  let i1 = node[it + iq1 * nelemn];
  let i2 = node[it + iq2 * nelemn];
  let i3 = node[it + iq3 * nelemn];
  let d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);

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

/// daxpy para vectores 1D planos. Fiel a la definición original de LINPACK.
/// Usa `isize` para los incrementos para soportar valores negativos.
pub fn daxpy_1d(n: usize, da: f64, dx: &[f64], incx: isize, dy: &mut [f64], incy: isize) {
  if n == 0 || da == 0.0 {
    return;
  }

  // Calcular índices iniciales si los incrementos son negativos
  let mut ix: isize = if incx < 0 {
    ((n - 1) as isize) * (-incx)
  } else {
    0
  };
  let mut iy: isize = if incy < 0 {
    ((n - 1) as isize) * (-incy)
  } else {
    0
  };

  for _ in 0..n {
    dy[iy as usize] += da * dx[ix as usize];
    ix += incx;
    iy += incy;
  }
}

/// Computa dy = da * dx + dy para columnas de una matriz Vec<f64> (flat).
fn daxpy_same_matrix(
  n: usize,
  da: f64,
  abd: &mut [f64],
  start_row_x: usize,
  col_x: usize,
  incx: usize,
  start_row_y: usize,
  col_y: usize,
  incy: usize,
  nrow: usize,
) {
  // Optimización de Fortran: si el escalar es 0, no hacemos nada.
  if n == 0 || da == 0.0 {
    return;
  }

  let mut rx = start_row_x;
  let mut ry = start_row_y;

  for _ in 0..n {
    // 1. Leemos el valor de la columna X (préstamo inmutable)
    let val_x = abd[rx + col_x * nrow];

    // 2. Escribimos en la columna Y (préstamo mutable)
    // Como val_x ya es una copia (f64), el préstamo inmutable de abd ya terminó.
    abd[ry + col_y * nrow] += da * val_x;

    rx += incx;
    ry += incy;
  }
}

/// daxpy para columnas de DOS MATRICES DIFERENTES Vec<f64> (flat).
pub fn daxpy_diff_matrices(
  n: usize,
  da: f64,
  mat_x: &[f64],
  start_row_x: usize,
  col_x: usize,
  incx: usize,
  mat_y: &mut [f64],
  start_row_y: usize,
  col_y: usize,
  incy: usize,
  nrow: usize,
) {
  if n == 0 || da == 0.0 {
    return;
  }

  let mut rx = start_row_x;
  let mut ry = start_row_y;

  for _ in 0..n {
    mat_y[ry + col_y * nrow] += da * mat_x[rx + col_x * nrow];
    rx += incx;
    ry += incy;
  }
}

pub fn daxpy_vec_to_matrix(
  n: usize,
  da: f64,
  vec_x: &[f64],
  start_idx_x: usize,
  incx: usize,
  mat_y: &mut [f64],
  start_row_y: usize,
  col_y: usize,
  incy: usize,
  nrow: usize,
) {
  if n == 0 || da == 0.0 {
    return;
  }

  let mut ix = start_idx_x;
  let mut ry = start_row_y;

  for _ in 0..n {
    mat_y[ry + col_y * nrow] += da * vec_x[ix];
    ix += incx;
    ry += incy;
  }
}

pub fn daxpy_matrix_to_vec(
  n: usize,
  da: f64,
  mat_x: &[f64],
  start_row_x: usize,
  col_x: usize,
  incx: usize,
  vec_y: &mut [f64],
  start_idx_y: usize,
  incy: usize,
  nrow: usize,
) {
  if n == 0 || da == 0.0 {
    return;
  }
  let mut rx = start_row_x;
  let mut iy = start_idx_y;
  for _ in 0..n {
    vec_y[iy] += da * mat_x[rx + col_x * nrow];
    rx += incx;
    iy += incy;
  }
}

pub fn ddot_matrix_to_vec(
  n: usize,
  mat_x: &[f64],
  start_row_x: usize,
  col_x: usize,
  incx: usize,
  vec_y: &[f64],
  start_idx_y: usize,
  incy: usize,
  nrow: usize,
) -> f64 {
  if n == 0 {
    return 0.0;
  }

  let mut rx = start_row_x;
  let mut iy = start_idx_y;
  let mut dtemp = 0.0;

  for _ in 0..n {
    dtemp += mat_x[rx + col_x * nrow] * vec_y[iy];
    rx += incx;
    iy += incy;
  }

  dtemp
}

fn dgbfa(
  abd: &mut [f64],
  _lda: usize,
  n: usize,
  ml: usize,
  mu: usize,
  ipvt: &mut Vec<usize>,
  info: &mut usize,
) {
  //? dgbfa() factors a band matrix by elimination.
  //?   dgbfa is usually called by dgbco, but it can be called
  //?   directly with a saving in time if  rcond  is not needed.
  //?      abd     real ( kind = rk8 )(lda, n)
  //?              contains the matrix in band storage.  the //columns
  //?              of the matrix are stored in the columns of  //abd  and
  //?              the diagonals of the matrix are stored in rows
  //?              ml+1 through 2*ml+mu+1 of  abd .
  //?              see the comments below for details.
  //?
  //?      lda     integer
  //?              the leading dimension of the array  abd .
  //?              2*ml + mu + 1 <= LDA.
  //?
  //?      n       integer
  //?              the order of the original matrix.
  //?
  //?      ml      integer
  //?              number of diagonals below the main diagonal.
  //?              0 <= ml < n .
  //?
  //?      mu      integer
  //?              number of diagonals above the main diagonal.
  //?              0 <= mu < n .
  //?              more efficient if  ml <= mu .
  //?   on return
  //?
  //?      abd     an upper triangular matrix in band storage and
  //?              the multipliers which were used to obtain it.
  //?              the factorization can be written  a = l*u  //where
  //?              l  is a product of permutation && unit lower
  //?              triangular matrices &&  u  is upper //triangular.
  //?
  //?      ipvt    integer(n)
  //?              an integer vector of pivot indices.
  //?
  //?      info    integer
  //?              = 0  normal value.
  //?              = k  if  u(k,k) == 0.0 .  this is not an //error
  //?                   condition for this subroutine, but it does
  //?                   indicate that dgbsl will divide by zero if
  //?                   called.  use  rcond  in dgbco for a //reliable
  //?                   indication of singularity.
  //?
  //?   band storage
  //?
  //?         if  a  is a band matrix, the following program //segment
  //?         will set up the input.
  //?
  //?                 ml = (band width below the diagonal)
  //?                 mu = (band width above the diagonal)
  //?                 m = ml + mu + 1
  //?                 do j = 1, n
  //?                    i1 = max ( 1, j-mu )
  //?                    i2 = min ( n, j+ml )
  //?                    do i = i1, i2
  //?                       k = i - j + m
  //?                       abd(k,j) = a(i,j)
  //?                    }
  //?                 }
  //?
  //?         this uses rows  ml+1  through  2*ml+mu+1  of  abd .
  //?         in addition, the first  ml  rows in  abd  are used //for
  //?         elements generated during the triangularization.
  //?         the total number of rows needed in  abd  is  2*ml+mu//+1 .
  //?         the  ml+mu by ml+mu  upper left triangle && the
  //?         ml by ml  lower right triangle are not referenced.
  //?
  let mut t;
  let m = ml + mu;
  let nrow = _lda;
  *info = 0;
  //r!  Zero initial fill-in columns.
  let j0 = mu + 1;
  let j1 = n.min(m + 1) - 1;

  //do jz = j0, j1
  for jz in j0..j1 {
    //let i0 = m + 1 - jz;
    let i0 = m - jz;
    //do i = i0, ml
    for i in i0..ml {
      abd[i + jz * nrow] = 0.0;
    }
  }

  let mut jz = j1;
  let mut ju = 0;
  //r!  Gaussian elimination with partial pivoting.
  //do k = 1, n-1
  for k in 0..n - 1 {
    //?  Zero next fill-in column.
    jz += 1;
    if jz < n {
      for i in 0..ml {
        abd[i + jz * nrow] = 0.0;
      }
    }
    //  Find L = pivot index.
    //* Convertimos n de 1-index a 0-index restandole 1
    let lm = ((ml as i32).min((n - 1 - k) as i32)) as usize;
    //* m ya esta en 0-index
    let mut l = idamax_matrix(lm + 1, abd, m, k, 1, nrow) + m;
    ipvt[k] = l + k - m;
    //?  Zero pivot implies this column already triangularized.
    if abd[l + k * nrow] == 0.0 {
      *info = k;
    //?  Interchange if necessary.
    } else {
      if l != m {
        t = abd[l + k * nrow];
        abd[l + k * nrow] = abd[m + k * nrow];
        abd[m + k * nrow] = t;
      }
      //r!  Compute multipliers.
      t = -1.0 / abd[m + k * nrow];
      dscal_matrix(lm, t, abd, m + 1, k, 1, nrow);
      //r!  Row elimination with column indexing.
      //* La comparacion se hace en 1-index
      ju = ju.max(mu + ipvt[k] as usize + 1).min(n);
      let mut mm = m;

      //do j = k+1, ju
      for j in (k + 1)..ju {
        l -= 1;
        mm -= 1;
        let t = abd[l + j * nrow];
        if l != mm {
          abd[l + j * nrow] = abd[mm + j * nrow];
          abd[mm + j * nrow] = t;
        }
        daxpy_same_matrix(lm, t, abd, m + 1, k, 1, mm + 1, j, 1, nrow);
      }
    }
  }

  //* n esta en 1-index
  ipvt[n - 1] = n - 1;
  //* Hacemos una comprobacion para evitar un desvordameinto
  if abd.len() > m + (n - 1) * nrow {
    if abd[m + (n - 1) * nrow] == 0.0 {
      *info = n;
    }
  }
}

fn dgbsl(
  abd: &mut [f64],
  _lda: usize,
  n: usize,
  ml: usize,
  mu: usize,
  ipvt: &mut [usize],
  b: &mut [f64],
  job: usize,
) {
  //r!! dgbsl() solves a banded system factored by DGBFA.
  //r!
  //r!  Discussion:
  //r!
  //r!    SGBSL can solve either a * x = b  or  trans(a) * x = b.
  //r!
  //r!  Parameters:
  //r!
  //r!     on entry
  //r!
  //r!        abd     real ( kind = rk8 )(lda, n)
  //r!                the output from dgbco or dgbfa.
  //r!
  //r!        lda     integer
  //r!                the leading dimension of the array  abd .
  //r!
  //r!        n       integer
  //r!                the order of the original matrix.
  //r!
  //r!        ml      integer
  //r!                number of diagonals below the main diagonal.
  //r!
  //r!        mu      integer
  //r!                number of diagonals above the main diagonal.
  //r!
  //r!        ipvt    integer(n)
  //r!                the pivot vector from dgbco or dgbfa.
  //r!
  //r!        b       real ( kind = rk8 )(n)
  //r!                the right hand side vector.
  //r!
  //r!        job     integer
  //r!                = 0         to solve  a*x = b ,
  //r!                = nonzero   to solve  trans(a)*x = b , where
  //r!                            trans(a)  is the transpose.
  //r!
  //r!     on return
  //r!
  //r!        b       the solution vector  x .
  //r!
  //r!     error condition
  //r!
  //r!        a division by zero will occur if the input factor //contains a
  //r!        zero on the diagonal.  technically this indicates //singularity
  //r!        but it is often caused by improper arguments or //improper
  //r!        setting of lda .  it will not occur if the subroutines //are
  //r!        called correctly && if dgbco has set 0.0 < RCOND
  //r!        or dgbfa has set info == 0 .
  //r!
  //r!     to compute  inverse(a) * c  where  c  is a matrix
  //r!     with  p  columns
  //r!           call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
  //r!           ifrcond is too small) go to ...
  //r!           do j = 1, p
  //r!              call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), //0 )
  //r!           }
  //r!
  //r!     linpack. this version dated 08/14/78 .
  //r!     cleve moler, university of new mexico, argonne national //lab.
  //r!
  let m = mu + ml;
  let nrow = _lda;
  let mut t;
  //r!  JOB = 0, Solve  a * x = b.
  //r!  First solve l*y = b.
  if job == 0 {
    if 0 < ml {
      //do k = 1, n-1
      for k in 0..n - 1 {
        let lm = ml.min(n - k - 1);
        let l = ipvt[k];
        let t = b[l];
        if l != k {
          b[l] = b[k];
          b[k] = t;
        }
        daxpy_matrix_to_vec(lm, t, abd, m + 1, k, 1, b, k + 1, 1, nrow);
      }
    }
    //r!  Now solve u*x = y.
    //do k = n, 1, -1
    for k in (0..n).rev() {
      b[k] /= abd[m + k * nrow];
      let lm = k.min(m);
      let la = m - lm;
      let lb = k - lm;
      t = -b[k];
      daxpy_matrix_to_vec(lm, t, abd, la, k, 1, b, lb, 1, nrow);
    }
    //r!  JOB nonzero, solve  trans(a) * x = b.
    //r!  First solve  trans(u)*y = b.
  } else {
    //do k = 1, n
    for k in 0..n {
      let lm = k.min(m);
      let la = m - lm;
      let lb = k - lm;
      //let t = ddot( lm, abd(la,k), 1, b(lb), 1 );
      let t = ddot_matrix_to_vec(lm, &abd, la, k, 1, &b, lb, 1, nrow);
      b[k] = (b[k] - t) / abd[m + k * nrow];
    }
    //r!  Now solve trans(l)*x = y
    if 0 < ml {
      //do k = n-1, 1, -1
      for k in (0..n).rev() {
        let lm = ml.min(n - 1 - k);
        //b[k] += ddot(lm, abd(m+1,k), 1, b(k+1), 1);
        b[k] += ddot_matrix_to_vec(lm, &abd, m + 1, k, 1, &b, k + 1, 1, nrow);
        let l = ipvt[k];
        if l != k {
          t = b[l];
          b[l] = b[k];
          b[k] = t;
        }
      }
    }
  }
}

/// Escala un vector 1D por un escalar (dx = da * dx).
pub fn dscal(n: usize, da: f64, dx: &mut [f64], incx: usize) {
  if n == 0 || incx == 0 {
    return;
  }

  if incx == 1 {
    // El compilador de Rust optimizará esto automáticamente (SIMD/Unrolling)
    for i in 0..n {
      dx[i] *= da;
    }
  } else {
    let mut idx = 0;
    for _ in 0..n {
      dx[idx] *= da;
      idx += incx;
    }
  }
}

/// Escala una columna de una matriz Vec<f64> (flat) por un escalar.
pub fn dscal_matrix(
  n: usize,
  da: f64,
  abd: &mut [f64],
  start_row: usize,
  col: usize,
  incx: usize,
  nrow: usize,
) {
  if n == 0 || incx == 0 {
    return;
  }

  let mut current_row = start_row;
  for _ in 0..n {
    abd[current_row + col * nrow] *= da;
    current_row += incx;
  }
}

fn gdump(flow: &mut Flow_struct) -> std::io::Result<()> {
  //? gdump() writes information to a file.
  //r!  Discussion:
  //r!
  //r!    The information can be used to create
  //r!    graphics images.  In order to keep things simple, exactly one
  //r!    value, real or integer, is written per record.
  if flow.json {
    gdump_json(flow)?;
    Ok(())
  } else {
    let mut j;
    let mut fval: f64;

    let mut archivo = fs::File::create(format!("{}/{}.dat", flow.data_dir, flow.fileg)).unwrap();

    writeln!(archivo, "long: {}", flow._long)?;
    writeln!(archivo, "nelemn: {}", flow.nelemn)?;
    writeln!(archivo, "np: {}", flow.np)?;
    writeln!(archivo, "npara: {}", flow.npara)?;
    writeln!(archivo, "nx: {}", flow.nx)?;
    writeln!(archivo, "ny: {}", flow.ny)?;
    //r!
    //r!  Pressures
    //r!
    //do i = 1, np
    for i in 0..flow.np {
      j = flow.insc[i];
      if j <= 0 {
        fval = 0.0
      } else {
        fval = flow.f[(j - 1) as usize];
      }
      writeln!(archivo, "{}", fval)?;
    }
    //r!
    //r!  Horizontal velocities, U
    //r!
    //do i = 1, np
    for i in 0..flow.np {
      j = flow.indx[i];
      if j == 0 {
        fval = 0.0
      } else if j < 0 {
        fval = ubdry(flow.yc[i], flow.para1);
      } else {
        fval = flow.f[(j - 1) as usize];
      }
      writeln!(archivo, "{}", fval)?;
    }
    //r!
    //r!  Vertical velocities, V
    //r!
    for i in 0..flow.np {
      j = flow.indx[i + flow.np];
      if j <= 0 {
        fval = 0.0
      } else {
        fval = flow.f[(j - 1) as usize];
      }
      writeln!(archivo, "{}", fval)?;
    }

    for i in 0..flow.np {
      writeln!(archivo, "{}", flow.indx[i])?;
      writeln!(archivo, "{}", flow.indx[i + flow.np])?;
    }

    for i in 0..flow.np {
      writeln!(archivo, "{}", flow.insc[i])?;
    }

    for i in 0..flow.nelemn {
      writeln!(archivo, "{}", flow.isotri[i])?;
    }

    for i in 0..flow.nelemn {
      writeln!(
        archivo,
        "{}, {}, {}, {}, {}, {}",
        flow.node[i + 0 * flow.nelemn],
        flow.node[i + 1 * flow.nelemn],
        flow.node[i + 2 * flow.nelemn],
        flow.node[i + 3 * flow.nelemn],
        flow.node[i + 4 * flow.nelemn],
        flow.node[i + 5 * flow.nelemn]
      )?;
    }

    writeln!(archivo, "parametro = {}", flow.para1)?;
    writeln!(archivo, "reynold = {}", flow.reynld)?;
    writeln!(archivo, "rjpnew = {}", flow.rjpnew)?;

    for i in 0..flow.np {
      writeln!(archivo, "{}", flow.xc[i])?;
    }

    for i in 0..flow.np {
      writeln!(archivo, "{}", flow.yc[i])?;
    }

    println!(
      "GDUMP wrote data set to file {}/{}",
      flow.data_dir, flow.fileg
    );
    Ok(())
  }
}

fn gdump_json(flow: &mut Flow_struct) -> std::io::Result<()> {
  let mut j;
  let mut fval: f64;

  let archivo = fs::File::create(format!("{}/{}.json", flow.data_dir, flow.fileg)).unwrap();
  let mut buffer = BufWriter::new(archivo);

  writeln!(buffer, "{{")?;

  writeln!(buffer, "  \"long\": {},", flow._long)?;
  writeln!(buffer, "  \"nelemn\": {},", flow.nelemn)?;
  writeln!(buffer, "  \"np\": {},", flow.np)?;
  writeln!(buffer, "  \"npara\": {},", flow.npara)?;
  writeln!(buffer, "  \"nx\": {},", flow.nx)?;
  writeln!(buffer, "  \"ny\": {},", flow.ny)?;
  //r!
  //r!  Pressures
  //r!
  //do i = 1, np
  write!(buffer, "  \"p\": [")?;
  for i in 0..flow.np {
    j = flow.insc[i];
    if j <= 0 {
      fval = 0.0
    } else {
      fval = flow.f[(j - 1) as usize];
    }
    write!(buffer, "{}", fval)?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;
  //r!
  //r!  Horizontal velocities, U
  //r!
  write!(buffer, "  \"h_v_u\": [")?;
  for i in 0..flow.np {
    j = flow.indx[i];
    if j == 0 {
      fval = 0.0
    } else if j < 0 {
      fval = ubdry(flow.yc[i], flow.para1);
    } else {
      fval = flow.f[(j - 1) as usize];
    }
    write!(buffer, "{}", fval)?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;
  //r!
  //r!  Vertical velocities, V
  //r!
  write!(buffer, "  \"v_v_v\": [")?;
  for i in 0..flow.np {
    j = flow.indx[i + flow.np];
    if j <= 0 {
      fval = 0.0
    } else {
      fval = flow.f[(j - 1) as usize];
    }
    write!(buffer, "{}", fval)?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  write!(buffer, "  \"elements\": [")?;
  for i in 0..flow.np {
    write!(buffer, "[{}, {}]", flow.indx[i], flow.indx[i + flow.np])?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  write!(buffer, "  \"insc\": [")?;
  for i in 0..flow.np {
    write!(buffer, "{}", flow.insc[i])?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  write!(buffer, "  \"isotri\": [")?;
  for i in 0..flow.nelemn {
    write!(buffer, "{}", flow.isotri[i])?;
    if i < flow.nelemn - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  write!(buffer, "  \"node\": [")?;
  for i in 0..flow.nelemn {
    write!(
      buffer,
      "[{}, {}, {}, {}, {}, {}]",
      flow.node[i + 0 * flow.nelemn],
      flow.node[i + 1 * flow.nelemn],
      flow.node[i + 2 * flow.nelemn],
      flow.node[i + 3 * flow.nelemn],
      flow.node[i + 4 * flow.nelemn],
      flow.node[i + 5 * flow.nelemn]
    )?;
    if i < flow.nelemn - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  writeln!(buffer, "  \"parametro\": {},", flow.para1)?;
  writeln!(buffer, "  \"reynold\": {},", flow.reynld)?;
  writeln!(buffer, "  \"rjpnew\": {},", flow.rjpnew)?;

  write!(buffer, "  \"xc\": [")?;
  for i in 0..flow.np {
    write!(buffer, "{}", flow.xc[i])?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  write!(buffer, "  \"xm\": [")?;
  for i in 0..flow.np {
    write!(buffer, "{}", flow.yc[i])?;
    if i < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "]\n")?;

  writeln!(buffer, "}}")?;
  buffer.flush()?;

  println!(
    "GDUMP wrote data set to file {}/{}.json",
    flow.data_dir, flow.fileg
  );
  Ok(())
}

fn getg(f: &[f64], iline: &[i32], my: usize, _neqn: usize, u: &mut [f64]) {
  //? getg() outputs field values along the profile line X = XZERO.
  //do j = 1, my
  for j in 0..my {
    let k = iline[j];
    if 0 < k {
      u[j] = f[(k - 1) as usize];
    } else {
      u[j] = 0.0;
    }
  }
}

fn gram(flow: &mut Flow_struct) {
  //? gram() computes the Gram matrix, GR(I,J) = INTEGRAL PHI(I)//*PHI(J).
  //?
  //?  && the vector R(I) = INTEGRAL UI*PHI(I).
  //?
  //?  The integrals are computed along the line where the profile //is
  //?  specified.  The three point Gauss quadrature rule is used //for the
  //?  line integral.

  //?
  //?  Input values for 3 point Gauss quadrature
  //?
  let mut wt: [f64; 3] = [0.0; 3];
  let mut yq: [f64; 3] = [0.0; 3];
  wt[0] = 5.0 / 9.0;
  wt[1] = 8.0 / 9.0;
  wt[2] = wt[0];
  yq[0] = -0.7745966692;
  yq[1] = 0.0;
  yq[2] = -yq[0];

  let valid_nodes = [0, 1, 3];

  let mut bb: f64 = 0.0;
  let mut bx: f64 = 0.0;
  let mut by: f64 = 0.0;

  let mut bbb: f64 = 0.0;
  let mut bbx: f64 = 0.0;
  let mut bby: f64 = 0.0;
  //r!
  //r!  zero arrays
  //r!
  for i in 0..flow.my {
    flow.r[i] = 0.0;
    for j in 0..flow.my {
      flow.gr[i + j * flow.my] = 0.0;
    }
  }
  //r!
  //r!  Compute line integral by looping over intervals along line
  //r!  using three point Gauss quadrature
  //r!
  let xzero = flow.xc[flow.nodex0];
  //do 70 it = 1, nelemn
  for it in 0..flow.nelemn {
    //?
    //?  Check to see if we are in a triangle with a side along line
    //?  x = xzero.  If not, skip out
    //?
    let k = flow.node[it + 0 * flow.nelemn];
    let kk = flow.node[it + 1 * flow.nelemn];

    if 1.0e-04 < (flow.xc[k] - xzero).abs() || 1.0e-04 < (flow.xc[kk] - xzero).abs() {
      continue;
    }

    //do 60 iquad = 1, 3
    for iquad in 0..3 {
      let bma2 = (flow.yc[kk] - flow.yc[k]) / 2.0;
      let ar = bma2 * wt[iquad];
      let x = xzero;
      let y = flow.yc[k] + bma2 * (yq[iquad] + 1.0);
      //r!
      //r!  Compute u internal at quadrature points
      //r!
      let mut uiqdpt: f64 = 0.0;
      //do 30 iq = 1, nnodes
      for &iq in &valid_nodes {
        qbf(
          x,
          y,
          it,
          iq,
          flow.nelemn,
          &mut bb,
          &mut bx,
          &mut by,
          &flow.node,
          &flow.xc,
          &flow.yc,
        );
        let ip = flow.node[it + iq * flow.nelemn];
        let iun = flow.indx[ip];
        if 0 < iun {
          let ii = igetl(iun, &flow.iline, flow.my);
          uiqdpt = uiqdpt + bb * flow.ui[ii];
        } else {
          let ubc = ubdry(flow.yc[ip], flow.para1);
          uiqdpt += bb * ubc;
        }
      }
      //r!
      //r!  Only loop over flow.nodes lying on line x = xzero
      //r!
      //do 50 iq = 1, nnodes
      for &iq in &valid_nodes {
        let ip = flow.node[it + iq * flow.nelemn];
        qbf(
          x,
          y,
          it,
          iq,
          flow.nelemn,
          &mut bb,
          &mut bx,
          &mut by,
          &flow.node,
          &flow.xc,
          &flow.yc,
        );
        let i = flow.indx[ip];
        if i <= 0 {
          continue;
        };
        let ii = igetl(i, &flow.iline, flow.my);
        flow.r[ii] += bb * uiqdpt * ar;

        //do iqq = 1, nnodes
        for &iqq in &valid_nodes {
          let ipp = flow.node[it + iqq * flow.nelemn];
          qbf(
            x,
            y,
            it,
            iqq,
            flow.nelemn,
            &mut bbb,
            &mut bbx,
            &mut bby,
            &flow.node,
            &flow.xc,
            &flow.yc,
          );
          let j = flow.indx[ipp];
          if j != 0 {
            let jj = igetl(j, &flow.iline, flow.my);
            flow.gr[ii + jj * flow.my] += bb * bbb * ar;
          }
        }
      }
    }
  }

  if 2 <= flow.iwrite {
    println!("\nGram matrix:\n");
    //do i = 1,my
    for i in 0..flow.my {
      //do j = 1,my
      for j in 0..flow.my {
        println!("{}, {}, {}", i + 1, j + 1, flow.gr[i + j * flow.my]);
      }
    }
    println!(" ");
    println!("R vector:");
    println!(" ");
    //do i = 1,my
    for i in 0..flow.my {
      println!("{}, {}", i, flow.r[i]);
    }
  }
}

pub fn idamax(n: usize, dx: &[f64], incx: usize) -> usize {
  assert!(n >= 1, "idamax: n debe ser >= 1");
  assert!(incx > 0, "idamax: incx debe ser > 0");
  debug_assert!(
    dx.len() >= (n - 1) * incx + 1,
    "idamax: el slice es demasiado corto"
  );

  let mut max_abs = dx[0].abs();
  let mut max_idx = 0;

  if incx == 1 {
    for i in 1..n {
      let val = dx[i].abs();
      if val > max_abs {
        max_abs = val;
        max_idx = i;
      }
    }
  } else {
    let mut idx = incx; // primer índice a evaluar (0-based) después del primero
    for i in 1..n {
      let val = dx[idx].abs();
      if val > max_abs {
        max_abs = val;
        max_idx = i;
      }
      idx += incx;
    }
  }
  max_idx
}

/// Busca el índice del máximo valor absoluto en una columna de una matriz Vec<f64> (flat).
/// Asume que abd está indexada como abd[fila + columna * nrow].
pub fn idamax_matrix(
  n: usize,
  abd: &[f64],
  start_row: usize,
  col: usize,
  incx: usize,
  nrow: usize,
) -> usize {
  if n == 0 || incx == 0 {
    return 0;
  }

  // Empezamos con el primer elemento, igual que en Fortran
  let mut max_abs = abd[start_row + col * nrow].abs();
  let mut max_idx = 0;

  for i in 1..n {
    let current_row = start_row + i * incx;
    let val = abd[current_row + col * nrow].abs();

    if val > max_abs {
      max_abs = val;
      max_idx = i;
    }
  }
  max_idx // Devuelve el índice relativo (0-based) dentro del rango escaneado
}

fn igetl(k: i32, iline: &Vec<i32>, my: usize) -> usize {
  //r!! igetl() gets the local unknown number along the profile line.
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

fn linsys(
  a: &mut Vec<f64>,
  area: &Vec<f64>,
  f: &mut Vec<f64>,
  g: &Vec<f64>,
  indx: &Vec<i32>,
  insc: &Vec<i32>,
  ipivot: &mut Vec<usize>,
  maxrow: usize,
  nelemn: usize,
  neqn: usize,
  nlband: usize,
  nnodes: usize,
  node: &Vec<usize>,
  np: usize,
  nquad: usize,
  nrow: usize,
  para1: f64,
  para2: f64,
  phi: &Vec<f64>,
  psi: &Vec<f64>,
  reynld: f64,
  yc: &Vec<f64>,
) {
  //? linsys() sets up && solves the linear system.
  let mut info: usize = 0;
  let ioff: i32 = (nlband + nlband) as i32;
  let visc = 1.0 / reynld;
  let mut un: Vec<f64> = vec![0.0; 2];
  let mut unx: Vec<f64> = vec![0.0; 2];
  let mut uny: Vec<f64> = vec![0.0; 2];

  f.fill(0.0);
  a.fill(0.0);
  //?  For each element,
  for it in 0..nelemn {
    let ar = area[it] / 3.0;
    //?  and for each quadrature point in the element,
    for iquad in 0..nquad {
      //?  Evaluate velocities at quadrature point
      uval(
        &g, &indx, iquad, it, nelemn, neqn, nnodes, &node, np, nquad, para1, &phi, &mut un,
        &mut unx, &mut uny, &yc,
      );
      //?  For each basis function,
      //do iq = 1, nnodes
      for iq in 0..nnodes {
        let ip = node[it + iq * nelemn];
        let bb = phi[it + iquad * nelemn + iq * nelemn * nquad];
        let bx = phi[it + iquad * nelemn + iq * nelemn * nquad + nelemn * nquad * nnodes];
        let by = phi[it + iquad * nelemn + iq * nelemn * nquad + 2 * nelemn * nquad * nnodes];
        let bbl = psi[it + iquad * nelemn + iq * nelemn * nquad];
        let ihor = indx[ip] - 1;
        let iver = indx[ip + np] - 1;
        let iprs = insc[ip] - 1;

        if 0 <= ihor {
          f[ihor as usize] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
        }
        if 0 <= iver {
          f[iver as usize] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
        }
        //?  For another basis function,
        //do iqq = 1, nnodes
        for iqq in 0..nnodes {
          let ipp = node[it + iqq * nelemn];
          let bbb = phi[it + iquad * nelemn + iqq * nelemn * nquad];
          let bbx = phi[it + iquad * nelemn + iqq * nelemn * nquad + nelemn * nquad * nnodes];
          let bby = phi[it + iquad * nelemn + iqq * nelemn * nquad + 2 * nelemn * nquad * nnodes];
          let bbbl = psi[it + iquad * nelemn + iqq * nelemn * nquad];
          let ju = indx[ipp] - 1;
          let jv = indx[ipp + np] - 1;
          let jp = insc[ipp] - 1;
          //?  Horizontal velocity variable
          if 0 <= ju {
            if 0 <= ihor {
              let iuse = ihor - ju + ioff;
              //print!("{}-{}+{}={}, ",ihor, ju, ioff, iuse);
              a[iuse as usize + ju as usize * nrow] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]))
            }
            if 0 <= iver {
              let iuse = iver - ju + ioff;
              a[iuse as usize + ju as usize * nrow] += ar * bb * bbb * unx[1];
            }
            if 0 <= iprs {
              let iuse = iprs - ju + ioff;
              a[iuse as usize + ju as usize * nrow] += ar * bbx * bbl;
            }
          } else if ju == -2 {
            let uu = ubdry(yc[ipp], para2);
            if 0 <= ihor {
              f[ihor as usize] -= ar
                * uu
                * (visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]))
            }
            if 0 <= iver {
              f[iver as usize] -= ar * uu * bb * bbb * unx[1]
            }
            if 0 <= iprs {
              f[iprs as usize] -= ar * uu * bbx * bbl
            }
          }
          //?  Vertical velocity variable
          if 0 <= jv {
            if 0 <= ihor {
              let iuse = ihor - jv + ioff;
              a[iuse as usize + jv as usize * nrow] += ar * bb * bbb * uny[0]
            }
            if 0 <= iver {
              let iuse = iver - jv + ioff;
              a[iuse as usize + jv as usize * nrow] += ar
                * (visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]))
            }
            if 0 <= iprs {
              let iuse = iprs - jv + ioff;
              a[iuse as usize + jv as usize * nrow] += ar * bby * bbl;
            }
          }
          //?  Pressure variable
          if 0 <= jp {
            if 0 <= ihor {
              let iuse = ihor - jp + ioff;
              a[iuse as usize + jp as usize * nrow] -= ar * bx * bbbl;
            }
            if 0 <= iver {
              let iuse = iver - jp + ioff;
              a[iuse as usize + jp as usize * nrow] -= ar * by * bbbl;
            }
          }
        }
      }
    }
  }
  //?  To avoid singularity of the pressure system, the last pressure is simply assigned a value of 0.
  f[neqn - 1] = 0.0;
  //do j = neqn-nlband, neqn-1
  for j in (neqn - nlband - 1)..(neqn - 1) {
    let i = neqn - 1 - j + ioff as usize;
    a[i + j * nrow] = 0.0;
  }
  a[ioff as usize + (neqn - 1) * nrow] = 1.0;
  //?  Factor the matrix
  dgbfa(a, nrow, neqn, nlband, nlband, ipivot, &mut info);

  if info != 0 {
    println!(" ");
    println!("LINSYS - fatal error!");
    println!("DGBFA returns INFO = {}", info);
    std::process::exit(1);
  }
  //r!  Solve the linear system
  let job: usize = 0;
  dgbsl(a, nrow, neqn, nlband, nlband, ipivot, f, job);
}

fn nstoke(flow: &mut Flow_struct) {
  //? nstoke() solves the Navier Stokes equation using Taylor-Hood elements.
  //?  The G array contains the previous iterate.
  //?  The F array contains the right hand side initially && then the current iterate.
  for iter in 1..=flow.maxnew {
    flow.numnew += 1;
    linsys(
      &mut flow.a,
      &flow.area,
      &mut flow.f,
      &flow.g,
      &flow.indx,
      &flow.insc,
      &mut flow.ipivot,
      flow.maxrow,
      flow.nelemn,
      flow.neqn,
      flow.nlband,
      flow.nnodes,
      &flow.node,
      flow.np,
      flow.nquad,
      flow.nrow,
      flow.para1,
      flow.para2,
      &flow.phi,
      &flow.psi,
      flow.reynld,
      &flow.yc,
    );
    //?  Check for convergence
    flow
      .g
      .iter_mut()
      .zip(flow.f.iter())
      .for_each(|(a, b)| *a -= *b);
    let diff = flow.g[idamax(flow.neqn, &flow.g, 1)].abs();

    if 1 <= flow.iwrite {
      println!("NSTOKE iteration {} Mnorm = {}", iter, diff);
    }
    flow.g.copy_from_slice(&flow.f);
    if diff <= flow.tolnew {
      println!("Navier Stokes iteration converged in {} iteration.", iter);
      return;
    }
  }
  println!("Navier Stokes solution did not converge!");
}

//subroutine pval (g,insc,long,mx,my,nelemn,neqn,nnodes,node,np,//press)
//
//r!***************************************************************//**************80
//r!
//r!! pval() computes a table of pressures.
//r!
//r!  Licensing:
//r!
//r!    This code is distributed under the MIT license.
//r!
//r!  Modified:
//r!
//r!    20 January 2007
//r!
//r!  Author:
//r!
//r!    John Burkardt
//r!
//r!  Parameters:
//r!
//  implicit none
//
//  integer, parameter :: rk8 = kind ( 1.0 )
//
//  integer nelemn
//  integer flow.neqn
//  integer nnodes
//  integer np
//
//  real ( kind = rk8 ) g(neqn)
//  integer i
//  integerflow.insc(np)
//  integer ip
//  integer iq
//  integer it
//  integer ivar
//  integer j
//  logical long
//  integer mx
//  integer my
//  integer flow.node(nelemn,nnodes)
//  real ( kind = rk8 ) press(mx,my)
//
//  press(1:mx,1:my) = 0.0
//r!
//r!  Read the pressures where they are computed.
//r!  These are "(odd, odd)" points.
//r!
//  do it = 1, nelemn
//    do iq = 1, 3
//
//      ip = flow.node(it,iq)
//      ivar =flow.insc(ip)
//
//      if flow._long {
//        i = ((ip-1)/my)+1
//        j = mod(ip-1,my)+1
//      } else {
//        i = mod(ip-1,mx)+1
//        j = ((ip-1)/mx)+1
//      }
//
//      if 0 < ivar {
//        press(i,j) = g(ivar)
//      } else {
//        press(i,j) = 0.0
//      }
//
//    }
//  }
//r!
//r!  Interpolate the pressures at points (even, odd) && (odd, //even).
//r!
//  do i = 2,mx-1,2
//    do j = 1,my,2
//      press(i,j) = 0.5*(press(i-1,j)+press(i+1,j))
//    }
//  }
//
//  do j = 2,my-1,2
//    do i = 1,mx,2
//      press(i,j) = 0.5*(press(i,j-1)+press(i,j+1))
//    }
//  }
//r!
//r!  Interpolate the pressures at points (even,even).
//r!
//  do j = 2,my-1,2
//    do i = 2,mx-1,2
//      press(i,j) = 0.5*(press(i-1,j-1)+press(i+1,j+1))
//    }
//  }
//
//  return
//}

fn qbf(
  x: f64,
  y: f64,
  it: usize,
  _in: usize,
  nelemn: usize,
  bb: &mut f64,
  bx: &mut f64,
  by: &mut f64,
  node: &Vec<usize>,
  xc: &Vec<f64>,
  yc: &Vec<f64>,
) {
  //? qbf() evaluates a quadratic basis function in a triangle.
  if _in <= 2 {
    let in1 = _in;
    //* let in2 = (_in % 3) + 1; -> esta operacion espera
    //* que _in este en 1-index, asiq ue es necesario convertir
    let in2 = (_in + 1) % 3;
    //* let in3 = (_in + 1) % 3;
    let in3 = (_in + 2) % 3;
    let i1 = node[it + in1 * nelemn];
    let i2 = node[it + in2 * nelemn];
    let i3 = node[it + in3 * nelemn];
    let d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    let t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    *bb = t * (2.0 * t - 1.0);
    *bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    *by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    let inn = _in - 3;
    let in1 = inn;
    //* Estas operaciones esperan inn en 1-index
    //* asi que convertimos
    let in2 = (inn + 1) % 3;
    let in3 = (inn + 2) % 3;
    let i1 = node[it + in1 * nelemn];
    let i2 = node[it + in2 * nelemn];
    let i3 = node[it + in3 * nelemn];
    let j1 = i2;
    let j2 = i3;
    let j3 = i1;
    let d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    let c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    let t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    let s = 1.0 + ((yc[j2] - yc[j3]) * (x - xc[j1]) + (xc[j3] - xc[j2]) * (y - yc[j1])) / c;
    *bb = 4.0 * s * t;
    *bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    *by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
}

fn resid(flow: &mut Flow_struct) {
  //? resid() computes the residual.
  //?
  //?  Discussion:
  //?
  //?    The G array contains the current iterate.
  //?
  //?    The RES array will contain the value of the residual.
  let visc: f64 = 1.0 / flow.reynld;
  flow.res.fill(0.0);
  let mut un: Vec<f64> = vec![0.0; 2];
  let mut unx: Vec<f64> = vec![0.0; 2];
  let mut uny: Vec<f64> = vec![0.0; 2];
  //?  For each element,
  //do 90 it = 1, nelemn
  for it in 0..flow.nelemn {
    let ar = flow.area[it] / 3.0;
    //?  && for each quadrature point in the element,
    //do 80 iquad = 1, nquad
    for iquad in 0..flow.nquad {
      //?  Evaluate velocities at quadrature point
      uval(
        &flow.g,
        &flow.indx,
        iquad,
        it,
        flow.nelemn,
        flow.neqn,
        flow.nnodes,
        &flow.node,
        flow.np,
        flow.nquad,
        flow.para1,
        &flow.phi,
        &mut un,
        &mut unx,
        &mut uny,
        &flow.yc,
      );
      //?  For each basis function,
      //do 70 iq = 1, nnodes
      for iq in 0..flow.nnodes {
        let ip = flow.node[it + iq * flow.nelemn];
        let bb = flow.phi[it + iquad * flow.nelemn + iq * flow.nelemn * flow.nquad];
        let bx = flow.phi[it
          + iquad * flow.nelemn
          + iq * flow.nelemn * flow.nquad
          + flow.nelemn * flow.nquad * flow.nnodes];
        let by = flow.phi[it
          + iquad * flow.nelemn
          + iq * flow.nelemn * flow.nquad
          + 2 * flow.nelemn * flow.nquad * flow.nnodes];
        let bbl = flow.psi[it + iquad * flow.nelemn + iq * flow.nelemn * flow.nquad];
        let iprs = flow.insc[ip] - 1;
        let ihor = flow.indx[ip] - 1;
        let iver = flow.indx[ip + flow.np] - 1;

        if 0 <= ihor {
          flow.res[ihor as usize] += (un[0] * unx[0] + un[1] * uny[0]) * bb * ar;
        }

        if 0 <= iver {
          flow.res[iver as usize] += (un[0] * unx[1] + un[1] * uny[1]) * bb * ar;
        }
        //?  For another basis function,
        //do iqq = 1, nnodes
        for iqq in 0..flow.nnodes {
          let ipp = flow.node[it + iqq * flow.nelemn];
          let bbb = flow.phi[it + iquad * flow.nelemn + iqq * flow.nelemn * flow.nquad];
          let bbx = flow.phi[it
            + iquad * flow.nelemn
            + iqq * flow.nelemn * flow.nquad
            + flow.nelemn * flow.nquad * flow.nnodes];
          let bby = flow.phi[it
            + iquad * flow.nelemn
            + iqq * flow.nelemn * flow.nquad
            + 2 * flow.nelemn * flow.nquad * flow.nnodes];
          let bbbl = flow.psi[it + iquad * flow.nelemn + iqq * flow.nelemn * flow.nquad];
          let ju = flow.indx[ipp] - 1;
          let jv = flow.indx[ipp + flow.np] - 1;
          let jp = flow.insc[ipp] - 1;

          if 0 <= ju {
            if 0 <= ihor {
              let aijuu =
                visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              flow.res[ihor as usize] += aijuu * ar * flow.g[ju as usize];
            }
            if 0 <= iver {
              let aijvu = bb * bbb * unx[1];
              flow.res[iver as usize] += aijvu * ar * flow.g[ju as usize];
            }
            if 0 <= iprs {
              let aijpu = bbx * bbl;
              flow.res[iprs as usize] += aijpu * ar * flow.g[ju as usize];
            }
          } else if ju == -2 {
            let uu = ubdry(flow.yc[ipp], flow.para1);
            if 0 <= ihor {
              let aijuu =
                visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              flow.res[ihor as usize] += ar * aijuu * uu;
            }
            if 0 <= iver {
              let aijvu = bb * bbb * unx[1];
              flow.res[iver as usize] += ar * aijvu * uu;
            }
            if 0 <= iprs {
              let aijpu = bbx * bbl;
              flow.res[iprs as usize] += ar * aijpu * uu;
            }
          }

          if 0 <= jv {
            if 0 <= ihor {
              let aijuv = bb * bbb * uny[0];
              flow.res[ihor as usize] += aijuv * ar * flow.g[jv as usize];
            }
            if 0 <= iver {
              let aijvv =
                visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]);
              flow.res[iver as usize] += aijvv * ar * flow.g[jv as usize];
            }
            if 0 <= iprs {
              let aijpv = bby * bbl;
              flow.res[iprs as usize] += aijpv * ar * flow.g[jv as usize];
            }
          }

          if 0 <= jp {
            if 0 <= ihor {
              let aijup = -bx * bbbl;
              flow.res[ihor as usize] += aijup * ar * flow.g[jp as usize];
            }
            if 0 <= iver {
              let aijvp = -by * bbbl;
              flow.res[iver as usize] += aijvp * ar * flow.g[jp as usize];
            }
          }
        }
      }
    }
  }
  flow.res[flow.neqn - 1] = flow.g[flow.neqn - 1];

  let mut rmax = 0.0;
  let mut imax = 0;
  let mut ibad = 0;

  //do i = 1,neqn
  for i in 0..flow.neqn {
    let test = flow.res[i].abs();

    if rmax < test {
      rmax = test;
      imax = i;
    }

    if 1.0E-03 < test {
      ibad = ibad + 1;
    }
  }

  // ======== PRIMER BLOQUE: if (1 <= iwrite) ========
  if flow.iwrite >= 1 {
    println!();
    println!("RESIDUAL INFORMATION:");
    println!();
    println!("Worst residual is number {}", imax);
    println!("of magnitude {}", rmax);
    println!();
    println!(
      "Number of \"bad\" residuals is {} out of {}",
      ibad, flow.neqn
    );
    println!();
  }

  // ======== SEGUNDO BLOQUE: if (2 <= iwrite) ========
  if flow.iwrite >= 2 {
    println!("Raw residuals:");
    println!();

    let mut i = 0; // índice 0‑based para res (en Fortran es 1‑based)

    // Asumimos que indx1, indx2, insc tienen misma longitud (np)
    for j in 0..flow.np {
      let j_display = j + 1; // para mostrar en la salida (1‑based)

      // --- indx(j,1) ---
      if flow.indx[j] > 0 {
        let val = flow.res[i];
        if val.abs() <= 1e-3 {
          // formato: espacio, 'U', i, j, val (con g14.6)
          // Usamos {:>14.6} para un ancho de 14 con 6 decimales (similar a g14.6)
          println!(" U{:5}, {:5}, {:e}", i, j_display, val);
        } else {
          println!("*U{:5}, {:5}, {:e}", i, j_display, val);
        }
        i += 1;
      }

      // --- indx(j,2) ---
      if flow.indx[j + flow.np] > 0 {
        let val = flow.res[i];
        if val.abs() <= 1e-3 {
          println!(" V{:5}, {:5}, {:e}", i, j_display, val);
        } else {
          println!("*V{:5}, {:5}, {:e}", i, j_display, val);
        }
        i += 1;
      }

      // --- insc(j) ---
      if flow.insc[j] > 0 {
        let val = flow.res[i];
        if val.abs() <= 1e-3 {
          println!(" P{:5}, {:5}, {:e}", i, j_display, val);
        } else {
          println!("*P{:5}, {:5}, {:e}", i, j_display, val);
        }
        i += 1;
      }
    }
  }
}

//*******************************************************
fn setban(flow: &mut Flow_struct) {
  //? setban() computes the half band width.
  flow.nlband = 0;

  //do it = 1, nelemn
  for it in 0..flow.nelemn {
    //do iq = 1, nnodes
    for iq in 0..flow.nnodes {
      let ip: usize = flow.node[it + iq * flow.nelemn];
      //do iuk = 1, 3
      for iuk in 0..3 {
        let i: i32 = if iuk == 2 {
          flow.insc[ip]
        } else {
          flow.indx[ip + iuk * flow.np]
        };

        if 0 < i {
          //do iqq = 1, nnodes
          for iqq in 0..flow.nnodes {
            let ipp = flow.node[it + iqq * flow.nelemn];
            //do iukk = 1, 3
            for iukk in 0..3 {
              let j: i32 = if iukk == 2 {
                flow.insc[ipp]
              } else {
                flow.indx[ipp + iukk * flow.np]
              };
              flow.nlband = ((flow.nlband as i32).max((j - i) as i32)) as usize;
            }
          }
        }
      }
    }
  }

  flow.nband = flow.nlband + flow.nlband + 1;
  flow.nrow = flow.nlband + flow.nlband + flow.nlband + 1;

  println!("Lower bandwidth = {}", flow.nlband);
  println!("Total bandwidth = {}", flow.nband);
  println!("NROW  = {}", flow.nrow);
  if flow.maxrow < flow.nrow {
    println!("SETBAN - NROW is too large!");
    println!("The maximum allowed is {}", flow.maxrow);
    std::process::exit(1);
  }

  return;
}

fn setbas(flow: &mut Flow_struct) {
  //? computes the basis functions at each integration point.
  let mut bb: f64 = 0.0;
  let mut bx: f64 = 0.0;
  let mut by: f64 = 0.0;

  for it in 0..flow.nelemn {
    //do j = 1,nquad
    for j in 0..flow.nquad {
      let x = flow.xm[it + j * flow.nelemn];
      let y = flow.ym[it + j * flow.nelemn];
      for iq in 0..flow.nnodes {
        flow.psi[it + j * flow.nelemn + iq * flow.nelemn * flow.nquad] =
          bsp(x, y, it, iq, 1, flow.nelemn, &flow.node, &flow.xc, &flow.yc);
        qbf(
          x,
          y,
          it,
          iq,
          flow.nelemn,
          &mut bb,
          &mut bx,
          &mut by,
          &flow.node,
          &flow.xc,
          &flow.yc,
        );
        flow.phi[it + j * flow.nelemn + iq * flow.nelemn * flow.nquad] = bb;
        flow.phi[it
          + j * flow.nelemn
          + iq * flow.nelemn * flow.nquad
          + flow.nelemn * flow.nquad * flow.nnodes] = bx;
        flow.phi[it
          + j * flow.nelemn
          + iq * flow.nelemn * flow.nquad
          + 2 * flow.nelemn * flow.nquad * flow.nnodes] = by;
      }
    }
  }
}

//*********************************************************
fn setgrd(flow: &mut Flow_struct) {
  //? setgrd() sets up the grid for the problem..
  //?  Determine whether region is long or skinny.  This will  determine
  //?  how we number the flow.nodes && elements.
  if flow.ny < flow.nx {
    flow._long = true;
    println!("Using vertical ordering.");
  } else {
    flow._long = false;
    println!("Using horizontal ordering.");
  }
  //?  Set parameters for Taylor Hood element
  println!(" ");
  println!("SETGRD: Taylor Hood element");
  //?  Construct grid coordinates, elements, && ordering of //unknowns
  flow.neqn = 0;
  let mut ielemn = 0;
  let mut ic;
  let mut jc;

  for ip in 0..flow.np {
    if flow._long {
      ic = ((ip) / flow.my) + 1;
      jc = ((ip) % flow.my) + 1;
    } else {
      ic = ((ip) % flow.mx) + 1;
      jc = ((ip) / flow.mx) + 1;
    }

    let icnt = ic % 2;
    let jcnt = jc % 2;
    //?  If both the row count && the column count are odd,
    //?  && we're not in the last row or top column,
    //?  then we can define two new triangular elements based at the //node.
    //?
    //?  For horizontal ordering,
    //?  given the following arrangement of flow.nodes, for instance:
    //?
    //?    21 22 23 24 25
    //?    16 17 18 19 20
    //?    11 12 13 14 15
    //?    06 07 08 09 10
    //?    01 02 03 04 05
    //?
    //?  when we arrive at flow.node 13, we will define
    //?
    //?  element 7: (13, 23, 25, 18, 24, 19)
    //?  element 8: (13, 25, 15, 19, 20, 14)
    //?
    //?
    //?  For vertical ordering,
    //?  given the following arrangement of flow.nodes, for instance:
    //?
    //?    05 10 15 20 25
    //?    04 09 14 19 24
    //?    03 08 13 18 23
    //?    02 07 12 17 22
    //?    01 06 11 16 21
    //?
    //?  when we arrive at flow.node 13, we will define
    //?
    //?  element 7: (13, 25, 23, 19, 24, 18)
    //?  element 8: (13, 15, 25, 14, 20, 19)
    //?
    if icnt == 1 && jcnt == 1 && ic != flow.mx && jc != flow.my {
      if flow._long {
        let ip1 = ip + flow.my;
        let ip2 = ip + flow.my + flow.my;
        flow.node[ielemn + 0 * flow.nelemn] = ip;
        flow.node[ielemn + 1 * flow.nelemn] = ip2 + 2;
        flow.node[ielemn + 2 * flow.nelemn] = ip2;
        flow.node[ielemn + 3 * flow.nelemn] = ip1 + 1;
        flow.node[ielemn + 4 * flow.nelemn] = ip2 + 1;
        flow.node[ielemn + 5 * flow.nelemn] = ip1;
        flow.isotri[ielemn] = 0;
        ielemn = ielemn + 1;
        flow.node[ielemn + 0 * flow.nelemn] = ip;
        flow.node[ielemn + 1 * flow.nelemn] = ip + 2;
        flow.node[ielemn + 2 * flow.nelemn] = ip2 + 2;
        flow.node[ielemn + 3 * flow.nelemn] = ip + 1;
        flow.node[ielemn + 4 * flow.nelemn] = ip1 + 2;
        flow.node[ielemn + 5 * flow.nelemn] = ip1 + 1;
        flow.isotri[ielemn] = 0;
        ielemn = ielemn + 1;
      } else {
        let ip1 = ip + flow.mx;
        let ip2 = ip + flow.mx + flow.mx;
        ielemn += 1;
        flow.node[ielemn + 1 * flow.nelemn] = ip - 1;
        flow.node[ielemn + 2 * flow.nelemn] = ip2 - 1;
        flow.node[ielemn + 3 * flow.nelemn] = ip2 + 1;
        flow.node[ielemn + 4 * flow.nelemn] = ip1 - 1;
        flow.node[ielemn + 5 * flow.nelemn] = ip2;
        flow.node[ielemn + 6 * flow.nelemn] = ip1;
        flow.isotri[ielemn] = 0;
        ielemn += 1;
        flow.node[ielemn + 1 * flow.nelemn] = ip - 1;
        flow.node[ielemn + 2 * flow.nelemn] = ip2 + 1;
        flow.node[ielemn + 3 * flow.nelemn] = ip + 1;
        flow.node[ielemn + 4 * flow.nelemn] = ip1;
        flow.node[ielemn + 5 * flow.nelemn] = ip1 + 1;
        flow.node[ielemn + 6 * flow.nelemn] = ip + 1 - 1;
        flow.isotri[ielemn] = 0;
      };
    }
    if ic == 1 && 1 < jc && jc < flow.my {
      //?  Consider whether velocity unknowns should be associated with this flow.node.
      //? If we are in column 1, horizontal velocities are specified, and vertical velocities are zero.
      flow.indx[ip] = -1;
      flow.indx[ip + flow.np] = 0;
    } else if ic == flow.mx && 1 < jc && jc < flow.my {
      //?  If we are in column MX, horizontal velocities are unknown, and vertical velocities are zero.
      flow.neqn += 1;
      flow.indx[ip] = flow.neqn as i32;
      flow.indx[ip + flow.np] = 0;
    } else if jc == 1 || jc == flow.my {
      //?  Otherwise, if we are in row 1 or row MY, both horizontal and vertical velocities are zero.
      flow.indx[ip] = 0;
      flow.indx[ip + flow.np] = 0;
    } else {
      //?  Otherwise, we are at an interior flow.node
      flow.neqn += 2;
      flow.indx[ip] = flow.neqn as i32 - 1;
      flow.indx[ip + flow.np] = flow.neqn as i32;
    }
    //r!  Consider whether a pressure unknown should be associated with this flow.node.
    //r!  The answer is yes if both flow.nodes are odd.
    if jcnt == 1 && icnt == 1 {
      flow.neqn = flow.neqn + 1;
      flow.insc[ip] = flow.neqn as i32;
    } else {
      flow.insc[ip] = 0;
    }
  }
  //r!
  //r!  If debugging is requested, print out data.
  //r!
  if 2 <= flow.iwrite {
    println!(" ");
    println!("    I     flow.indx 1 & 2,flow.insc");
    println!(" ");
    for i in 0..flow.np {
      println!(
        "     {}, {}, {}, {}",
        i + 1,
        flow.indx[i],
        flow.indx[i + flow.np],
        flow.insc[i]
      )
    }
    println!(" ");
    println!("    IT    NODE(IT,1:6)");
    println!(" ");
    for it in 0..flow.nelemn {
      println!(
        "    {}, {}, {}, {}, {}, {}, {}",
        it + 1,
        flow.node[it + 0 * flow.nelemn] + 1,
        flow.node[it + 1 * flow.nelemn] + 1,
        flow.node[it + 2 * flow.nelemn] + 1,
        flow.node[it + 3 * flow.nelemn] + 1,
        flow.node[it + 4 * flow.nelemn] + 1,
        flow.node[it + 5 * flow.nelemn] + 1
      );
    }
  }

  println!("Number of unknowns = {}", flow.neqn);
  if flow.maxeqn < flow.neqn {
    println!("SETGRD - Too many unknowns!");
    println!("The maximum allowed is MAXEQN = {}", flow.maxeqn);
    println!("This problem requires NEQN = {}", flow.neqn);
    std::process::exit(1);
  }
}

fn setlin(flow: &mut Flow_struct) {
  //? setlin() gets the unknown indices along the profile line.
  //?  Determine the number of a flow.node on the profile line
  let itemp = ((2.0 * (flow.nx - 1) as f64 * 9.0) / flow.xlngth).round() as usize;

  flow.nodex0 = if flow._long {
    itemp * (2 * flow.ny - 1)
  } else {
    itemp
  };

  //do i = 1, my
  for i in 0..flow.my {
    let ip = if flow._long {
      flow.nodex0 + i
    } else {
      flow.nodex0 + flow.mx * i
    };
    flow.iline[i] = flow.indx[ip]
  }

  if 1 <= flow.iwrite {
    println!(" ");
    println!("SETLIN: unknown numbers along line:");
    println!(" ");
    for i in &flow.iline {
      print!(" {}", i);
    }
    println!(" ");
  }
}

fn setqud(flow: &mut Flow_struct) {
  //? setqud() sets midpoint quadrature rule information.
  for it in 0..flow.nelemn {
    let ip1 = flow.node[it + 0 * flow.nelemn];
    let ip2 = flow.node[it + 1 * flow.nelemn];
    let ip3 = flow.node[it + 2 * flow.nelemn];
    let x1 = flow.xc[ip1];
    let x2 = flow.xc[ip2];
    let x3 = flow.xc[ip3];
    let y1 = flow.yc[ip1];
    let y2 = flow.yc[ip2];
    let y3 = flow.yc[ip3];
    flow.xm[it + 0 * flow.nelemn] = 0.5 * (x1 + x2);
    flow.xm[it + 1 * flow.nelemn] = 0.5 * (x2 + x3);
    flow.xm[it + 2 * flow.nelemn] = 0.5 * (x3 + x1);
    flow.ym[it + 0 * flow.nelemn] = 0.5 * (y1 + y2);
    flow.ym[it + 1 * flow.nelemn] = 0.5 * (y2 + y3);
    flow.ym[it + 2 * flow.nelemn] = 0.5 * (y3 + y1);
    flow.area[it] =
      0.5 * ((y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3)).abs();
  }
  return;
}

fn setxy(flow: &mut Flow_struct) {
  //? setxy() sets the X, Y coordinates of grid points.
  //?  Construct grid coordinates
  //do ip = 1,np
  let mut ic;
  let mut jc;

  for ip in 0..flow.np {
    if flow._long {
      ic = (ip / flow.my) + 1;
      jc = (ip % flow.my) + 1;
    } else {
      ic = (ip % flow.mx) + 1;
      jc = (ip / flow.mx) + 1;
    }

    flow.xc[ip] = (ic - 1) as f64 * flow.xlngth / (2 * flow.nx - 2) as f64;
    flow.yc[ip] = (jc - 1) as f64 * flow.ylngth / (2 * flow.ny - 2) as f64;
  }

  if 2 <= flow.iwrite {
    println!(" ");
    println!("    I      XC           YC");
    println!(" ");
    for i in 0..flow.np {
      println!("{}, {}, {}", i + 1, flow.xc[i], flow.yc[i]);
    }
  }

  return;
}

fn ubdry(y: f64, para: f64) -> f64 {
  //? ubdry() sets the parabolic inflow in terms of the value of the parameter.
  4.0 * para * y * (3.0 - y) / 9.0
}

fn uval(
  g: &Vec<f64>,
  indx: &Vec<i32>,
  iquad: usize,
  it: usize,
  nelemn: usize,
  _neqn: usize,
  nnodes: usize,
  node: &Vec<usize>,
  np: usize,
  nquad: usize,
  para: f64,
  phi: &Vec<f64>,
  un: &mut Vec<f64>,
  unx: &mut Vec<f64>,
  uny: &mut Vec<f64>,
  yc: &Vec<f64>,
) {
  //? uval() evaluates the velocities at a given point in a //particular triangle.
  un[0] = 0.0;
  un[1] = 0.0;
  unx[0] = 0.0;
  unx[1] = 0.0;
  uny[0] = 0.0;
  uny[1] = 0.0;

  //do iq = 1, nnodes
  for iq in 0..nnodes {
    let ip = node[it + iq * nelemn];
    let bb = phi[it + iquad * nelemn + iq * nelemn * nquad];
    let bx = phi[it + iquad * nelemn + iq * nelemn * nquad + nelemn * nquad * nnodes];
    let by = phi[it + iquad * nelemn + iq * nelemn * nquad + 2 * nelemn * nquad * nnodes];

    //do iuk = 1, 2
    for iuk in 0..2 {
      let iun = indx[ip + iuk * np];

      if 0 < iun {
        un[iuk] += bb * g[(iun - 1) as usize];
        unx[iuk] += bx * g[(iun - 1) as usize];
        uny[iuk] += by * g[(iun - 1) as usize];
      } else if iun < 0 {
        let ip = node[it + iq * nelemn];
        let ubc = ubdry(yc[ip], para);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }
}

fn uv_plot3d(flow: &Flow_struct) {
  //r!! uv_plot3d() creates a velocity file for use by PLOT3D.
  //r!
  //r!  Given the following set of flow.nodes:
  //r!
  //r!    A  B  C
  //r!    D  E  F
  //r!    G  H  I
  //r!
  //r!  the file will have the form:
  //r!
  //r!    D, U(G), V(G), P
  //r!    D, U(H), V(H), P
  //r!    D, U(I), V(I), P
  //r!    D, U(D), V(D), P
  //r!    D, U(E), V(E), P
  //r!    D, U(F), V(F), P
  //r!    D, U(A), V(A), P
  //r!    D, U(B), V(B), P
  //r!    D, U(C), V(C), P
  //r!
  //r!  Here both D && P are set to 1 for now, representing dummy //values
  //r!  of density && pressure.

  let dval = 1.0;
  let fsmach = 1.0;
  let alpha = 1.0;
  let time = 1.0;

  //  call pval (f,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)
  //r!
  //r!  If NY < NX, then flow.nodes with a constant Y value are numbered //consecutively.
  //r!
  //  if flow._long {
  //    write(ivunit,'(2I5)')mx,my
  //    write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time
  //    do ii = 1,4
  //      do j = 1,my
  //        do i = 1,mx
  //          ip = (i-1)*my+j
  //          if ii == 1 {
  //            write(ivunit,'(G15.5)')dval
  //          } else if ii == 2 {
  //            k =flow.indx(ip,1)
  //            ifk == 0{
  //              uval = 0.0
  //            } else ifk < 0{
  //              uval = ubdry(yc(ip),para)
  //            } else {
  //              uval = f(k)
  //            }
  //            write(ivunit,'(G15.5)')uval
  //          } else if ii == 3 {
  //            k =flow.indx(ip,2)
  //            ifk == 0{
  //              vval = 0.0
  //            } else {
  //              vval = f(k)
  //            }
  //            write(ivunit,'(G15.5)')vval
  //          } else {
  //            write(ivunit,'(G15.5)')press(i,j)
  //          }
  //        }
  //      }
  //    }
  //r!
  //r!  If NX < NY, then flow.nodes with a constant X value are numbered //consecutively.
  //r!
  //  } else {
  //    write(ivunit,'(2I5)')mx,my
  //    write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time
  //    do ii = 1,4
  //      do i = 1,mx
  //        do j = 1,my
  //          if ii == 1 {
  //             write(ivunit,'(G15.5)')dval
  //          } else if ii == 2 {
  //            ip = (i-1)*my+j
  //            k =flow.indx(ip,1)
  //            ifk == 0{
  //              uval = 0.0
  //            } else ifk < 0{
  //              uval = ubdry(yc(i),para)
  //            } else {
  //              uval = f(k)
  //            }
  //            write(ivunit,'(G15.5)')uval
  //          } else if ii == 3 {
  //            k =flow.indx(ip,2)
  //            ifk == 0{
  //              vval = 0.0
  //            } else {
  //              vval = f(k)
  //            }
  //            write(ivunit,'(G15.5)')vval
  //          } else {
  //            write(ivunit,'(G15.5)')press(i,j)
  //          }
  //        }
  //      }
  //    }
  //  }
  //
  //  println!("UV_PLOT3D wrote data set {}", iset,' to file.");
}

fn uv_table(flow: &Flow_struct) -> std::io::Result<()> {
  //r!! uv_table() creates a velocity table file.
  if flow.json {
    uv_table_json(flow)?;
    Ok(())
  } else {
    let mut uval;
    let mut vval;
    let archivo = fs::File::create(format!("{}/{}.dat", flow.data_dir, flow.fileu))?;
    let mut buffer = BufWriter::new(archivo);
    //do ip = 1, np
    for ip in 0..flow.np {
      let mut k = flow.indx[ip];
      if k == 0 {
        uval = 0.0;
      } else if k < 0 {
        uval = ubdry(flow.yc[ip], flow.para1);
      } else {
        uval = flow.f[(k - 1) as usize];
      }

      k = flow.indx[ip + flow.np];
      if k == 0 {
        vval = 0.0;
      } else {
        vval = flow.f[(k - 1) as usize];
      }
      writeln!(buffer, "{}, {}", uval, vval)?;
    }
    buffer.flush()?;
    Ok(())
  }
}

fn uv_table_json(flow: &Flow_struct) -> Result<(), std::io::Error> {
  let mut uval;
  let mut vval;
  let archivo = fs::File::create(format!("{}/{}.json", flow.data_dir, flow.fileu))?;
  let mut buffer = BufWriter::new(archivo);
  //do ip = 1, np
  let mut k: i32;
  writeln!(buffer, "{{")?;
  write!(buffer, "  \"u_vals\": [")?;
  for ip in 0..flow.np {
    k = flow.indx[ip];
    if k == 0 {
      uval = 0.0;
    } else if k < 0 {
      uval = ubdry(flow.yc[ip], flow.para1);
    } else {
      uval = flow.f[(k - 1) as usize];
    }
    write!(buffer, "{}", uval)?;
    if ip < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;

  write!(buffer, "  \"v_vals\": [")?;
  for ip in 0..flow.np {
    k = flow.indx[ip + flow.np];
    if k == 0 {
      vval = 0.0;
    } else {
      vval = flow.f[(k - 1) as usize];
    }
    write!(buffer, "{}", vval)?;
    if ip < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "]\n")?;
  write!(buffer, "}}")?;
  buffer.flush()?;
  Ok(())
}

fn xy_plot3d(flow: &Flow_struct) -> std::io::Result<()> {
  //r!! xy_plot3d() creates a grid file for use by PLOT3D.
  //r!
  //r!  Given the following set of flow.nodes:
  //r!
  //r!    A  B  C
  //r!    D  E  F
  //r!    G  H  I
  //r!
  //r!  the file will have the form:
  //r!
  //r!    X(G), X(H), X(I), X(D), X(E), X(F), X(A), X(B), X(C),
  //r!    Y(G), Y(H), Y(I), Y(D), Y(E), Y(F), Y(A), Y(B), Y(C).
  //r!
  //r!  If NY < NX, then flow.nodes with a constant Y value are numbered //consecutively.
  //r!

  if flow.json {
    xy_plot3d_json(flow)?;
  } else {
    let archivo = fs::File::create(format!("{}/{}.dat", flow.data_dir, flow.filex))?;
    let mut buffer = BufWriter::new(archivo);
    let mut ip;
    if flow._long {
      writeln!(buffer, "{} {}", flow.mx, flow.my)?;
      //do i = 1,my
      for i in 0..flow.my {
        //do j = 1,mx
        for j in 0..flow.mx {
          ip = j * flow.my + i;
          writeln!(buffer, "{}", flow.xc[ip])?;
        }
      }

      for i in 0..flow.my {
        //do j = 1,mx
        for j in 0..flow.mx {
          ip = j * flow.my + i;
          writeln!(buffer, "{}", flow.yc[ip])?;
        }
      }
    //r!
    //r!  If NX < NY, then flow.nodes with a constant X value are numbered //consecutively.
    //r!
    } else {
      writeln!(buffer, "{} {}", flow.my, flow.mx)?;

      for i in 0..flow.mx {
        //do j = 1,mx
        for j in 0..flow.my {
          ip = j * flow.my + i;
          writeln!(buffer, "{}", flow.xc[ip])?;
        }
      }

      for i in 0..flow.mx {
        //do j = 1,mx
        for j in 0..flow.my {
          ip = j * flow.my + i;
          writeln!(buffer, "{}", flow.yc[ip])?;
        }
      }
    }
  }

  println!("XYDUMP wrote data set to file.");
  Ok(())
}

fn xy_plot3d_json(flow: &Flow_struct) -> std::io::Result<()> {
  let archivo = fs::File::create(format!("{}/{}.json", flow.data_dir, flow.filex))?;
  let mut buffer = BufWriter::new(archivo);
  let mut ip;
  let total_elements = flow.my * flow.mx;
  let mut current_index = 0;

  writeln!(buffer, "{{")?;
  writeln!(buffer, "  \"mx\": {},\n  \"my\": {},", flow.mx, flow.my)?;
  if flow._long {
    write!(buffer, "  \"xc\": [")?;
    //do i = 1,my
    for i in 0..flow.my {
      //do j = 1,mx
      for j in 0..flow.mx {
        ip = j * flow.my + i;
        write!(buffer, "{}", flow.xc[ip])?;

        if current_index < total_elements - 1 {
          write!(buffer, ", ")?;
        }
        current_index += 1;
      }
    }
    write!(buffer, "],\n")?;
    current_index = 0;
    write!(buffer, "  \"yc\": [")?;
    for i in 0..flow.my {
      //do j = 1,mx
      for j in 0..flow.mx {
        ip = j * flow.my + i;
        write!(buffer, "{}", flow.yc[ip])?;
        // Coma solo si NO es el último elemento
        if current_index < total_elements - 1 {
          write!(buffer, ", ")?;
        }
        current_index += 1;
      }
    }
    write!(buffer, "]\n")?;
  //r!
  //r!  If NX < NY, then flow.nodes with a constant X value are numbered //consecutively.
  //r!
  } else {
    write!(buffer, "  \"xc\": [")?;
    for i in 0..flow.mx {
      //do j = 1,mx
      for j in 0..flow.my {
        ip = j * flow.my + i;
        write!(buffer, "{}", flow.xc[ip])?;
        if current_index < total_elements - 1 {
          write!(buffer, ", ")?;
        }
        current_index += 1;
      }
    }
    write!(buffer, "],\n")?;
    current_index = 0;
    write!(buffer, "  \"yc\": [")?;
    for i in 0..flow.mx {
      //do j = 1,mx
      for j in 0..flow.my {
        ip = j * flow.my + i;
        write!(buffer, "{}", flow.yc[ip])?;
        if current_index < total_elements - 1 {
          write!(buffer, ", ")?;
        }
        current_index += 1;
      }
    }
    write!(buffer, "],\n")?;
  }
  writeln!(buffer, "}}")?;
  Ok(())
}

fn xy_table(flow: &Flow_struct) -> std::io::Result<()> {
  //r!! xy_table() creates an XY table file.
  if flow.json {
    xy_table_json(flow)?;
    Ok(())
  } else {
    let archivo = fs::File::create(format!("{}/{}.dat", flow.data_dir, flow.filex))?;
    let mut buffer = BufWriter::new(archivo);
    //do ip = 1, np
    for ip in 0..flow.np {
      writeln!(buffer, "{}, {}", flow.xc[ip], flow.yc[ip])?;
    }
    buffer.flush()?;
    Ok(())
  }
}

fn xy_table_json(flow: &Flow_struct) -> std::io::Result<()> {
  let archivo = fs::File::create(format!("{}/{}.json", flow.data_dir, flow.filex))?;
  let mut buffer = BufWriter::new(archivo);
  //do ip = 1, np

  writeln!(buffer, "{{")?;
  write!(buffer, "  \"xc\": [")?;
  for ip in 0..flow.np {
    write!(buffer, "{}", flow.xc[ip])?;
    if ip < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "],\n")?;
  write!(buffer, "  \"yc\": [")?;
  for ip in 0..flow.np {
    write!(buffer, "{}", flow.yc[ip])?;
    if ip < flow.np - 1 {
      write!(buffer, ", ")?;
    }
  }
  write!(buffer, "]\n")?;
  write!(buffer, "}}")?;

  buffer.flush()?;
  Ok(())
}
