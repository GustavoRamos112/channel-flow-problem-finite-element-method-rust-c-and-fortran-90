use std::vec;

#[allow(non_camel_case_types)]
pub struct Flow_struct<'c> {
  pub nx: usize,
  pub ny: usize,
  pub maxrow: usize,
  pub nelemn: usize,
  pub mx: usize,
  pub my: usize,
  pub np: usize,
  pub maxeqn: usize,
  pub nnodes: usize,
  pub nquad: usize,

  pub data_dir: &'c str,
  pub fileg: &'c str,
  pub filex: &'c str,
  pub fileu: &'c str,

  pub a2: f64,
  pub anew: f64,
  pub aold: f64,
  pub _long: bool,
  pub nband: usize,
  pub neqn: usize,
  pub nlband: usize,
  pub nodex0: usize,
  pub nrow: usize,

  pub para1: f64,
  pub para2: f64,
  pub rjpold: f64,
  pub test: f64,
  pub iwrite: usize,
  pub maxnew: usize,
  pub maxsec: usize,
  pub npara: usize,
  pub numnew: usize,
  pub numsec: usize,
  pub reynld: f64,
  pub rjpnew: f64,
  pub tolnew: f64,
  pub tolsec: f64,
  pub xlngth: f64,
  pub ylngth: f64,
  pub converged: bool,

  pub a: Vec<f64>,
  pub area: Vec<f64>,
  pub dcda: Vec<f64>,
  pub f: Vec<f64>,
  pub g: Vec<f64>,
  pub gr: Vec<f64>,
  pub iline: Vec<i32>,
  pub indx: Vec<i32>,
  pub insc: Vec<i32>,
  pub ipivot: Vec<usize>,
  pub isotri: Vec<i32>,
  pub node: Vec<usize>,
  pub phi: Vec<f64>,
  pub psi: Vec<f64>,
  pub r: Vec<f64>,
  pub res: Vec<f64>,
  pub ui: Vec<f64>,
  pub unew: Vec<f64>,
  pub xc: Vec<f64>,
  pub xm: Vec<f64>,
  pub yc: Vec<f64>,
  pub ym: Vec<f64>,

  pub save_times: bool,
  pub save_data: bool,
  pub json: bool,
}

impl Flow_struct<'_> {
  pub fn new(nx: usize, ny: usize) -> Self {
    let maxrow = 27 * ny;
    let nelemn = 2 * (nx - 1) * (ny - 1);
    let mx = 2 * nx - 1;
    let my = 2 * ny - 1;
    let np = mx * my;
    let maxeqn = 2 * mx * my + nx * ny;
    let nnodes = 6;
    let nquad = 3;
    Self {
      nx: nx,
      ny: ny,
      maxrow: maxrow,
      nelemn: nelemn,
      mx: mx,
      my: my,
      np: np,
      maxeqn: maxeqn,
      nnodes: nnodes,
      nquad: nquad,

      data_dir: "data",
      fileg: "display",
      filex: "xy",
      fileu: "uv",

      a2: 0.0,
      anew: 0.0,
      aold: 0.0,
      _long: false,
      nband: 0,
      neqn: 0,
      nlband: 0,
      nodex0: 0,
      nrow: 0,

      para1: 0.0,
      para2: 0.0,
      rjpold: 0.0,
      test: 0.0,
      iwrite: 10,
      maxnew: 10,
      maxsec: 8,
      npara: 1,
      numnew: 0,
      numsec: 0,
      reynld: 1.0,
      rjpnew: 0.0,
      tolnew: 1.0E-04,
      tolsec: 1.0E-06,
      xlngth: 10.0,
      ylngth: 3.0,
      converged: false,

      //* a: a(1:nrow,1:neqn) = 0.0D+00
      a: vec![0.0; 1 * 2],
      area: vec![0.0; nelemn],
      dcda: vec![0.0; my],
      //* f: f(1:neqn) = 0.0D+00
      f: vec![0.0],
      //* g: g(1:neqn) = 0.0D+00
      g: vec![0.0],
      gr: vec![0.0; my * my],
      iline: vec![0; my],
      indx: vec![0; np * 2],
      insc: vec![0; np],
      ipivot: vec![0; maxeqn],
      isotri: vec![0; nelemn],
      node: vec![0; nelemn * nnodes],
      phi: vec![0.0; nelemn * nquad * nnodes * 3],
      psi: vec![0.0; nelemn * nquad * nnodes],
      r: vec![0.0; my],
      //* res: res(1:neqn) = 0.0
      res: vec![0.0],
      ui: vec![0.0; my],
      unew: vec![0.0; my],
      xc: vec![0.0; np],
      xm: vec![0.0; nelemn * nquad],
      yc: vec![0.0; np],
      ym: vec![0.0; nelemn * nquad],

      save_times: true,
      save_data: false,
      json: false,
    }
  }
}
