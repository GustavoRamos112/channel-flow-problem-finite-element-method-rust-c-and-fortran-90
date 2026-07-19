#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <print>
#include <format>
#include <fstream>
#include <filesystem>

struct Flow_struct {
  size_t nx;
  size_t ny;
  size_t maxrow;
  size_t nelemn;
  size_t mx;
  size_t my;
  size_t np;
  size_t maxeqn;
  size_t nnodes;
  size_t nquad;

  std::string data_dir;
  std::string fileg;
  std::string filex;
  std::string fileu;

  double a2;
  double anew;
  double aold;
  bool _long;
  size_t nband;
  size_t neqn;
  size_t nlband;
  size_t nodex0;
  size_t nrow;

  double para1;
  double para2;
  double rjpold;
  double test;
  size_t iwrite;
  size_t maxnew;
  size_t maxsec;
  size_t npara;
  size_t numnew;
  size_t numsec;
  double reynld;
  double rjpnew;
  double tolnew;
  double tolsec;
  double xlngth;
  double ylngth;
  bool converged;

  std::vector<double> a;
  std::vector<double> area;
  std::vector<double> dcda;
  std::vector<double> f;
  std::vector<double> g;
  std::vector<double> gr;
  std::vector<int> iline;
  std::vector<int> indx;
  std::vector<int> insc;
  std::vector<size_t> ipivot;
  std::vector<int> isotri;
  std::vector<std::vector<size_t>> node;
  std::vector<std::vector<std::vector<std::vector<double>>>> phi;
  std::vector<std::vector<std::vector<double>>> psi;
  std::vector<double> r;
  std::vector<double> res;
  std::vector<double> ui;
  std::vector<double> unew;
  std::vector<double> xc;
  std::vector<double> xm;
  std::vector<double> yc;
  std::vector<double> ym;

  bool save_times;
  bool save_data;
  bool json;

  Flow_struct(size_t nx, size_t ny)
    : nx(nx), ny(ny) {
    maxrow = 27 * ny;
    nelemn = 2 * (nx - 1) * (ny - 1);
    mx = 2 * nx - 1;
    my = 2 * ny - 1;
    np = mx * my;
    maxeqn = 2 * mx * my + nx * ny;
    nnodes = 6;
    nquad = 3;

    data_dir = "data";
    fileg = "display";
    filex = "xy";
    fileu = "uv";

    a2 = 0.0;
    anew = 0.0;
    aold = 0.0;
    _long = false;
    nband = 0;
    neqn = 0;
    nlband = 0;
    nodex0 = 0;
    nrow = 0;
    para1 = 0.0;
    para2 = 0.0;
    rjpold = 0.0;
    test = 0.0;
    iwrite = 10;
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
    converged = false;

    a = std::vector<double>(1 * 2, 0.0);
    area = std::vector<double>(nelemn, 0.0);
    dcda = std::vector<double>(my, 0.0);
    f = std::vector<double>(1, 0.0);
    g = std::vector<double>(1, 0.0);
    gr = std::vector<double>(my * my, 0.0);
    iline = std::vector<int>(my, 0);
    indx = std::vector<int>(np * 2, 0);
    insc = std::vector<int>(np, 0);
    ipivot = std::vector<size_t>(maxeqn, 0);
    isotri = std::vector<int>(nelemn, 0);
    node = std::vector<std::vector<size_t>>(nelemn, std::vector<size_t>(nnodes, 0));
    phi = std::vector<std::vector<std::vector<std::vector<double>>>>(
      nelemn, std::vector<std::vector<std::vector<double>>>(
        nquad, std::vector<std::vector<double>>(
          nnodes, std::vector<double>(3, 0.0))));
    psi = std::vector<std::vector<std::vector<double>>>(
      nelemn, std::vector<std::vector<double>>(
        nquad, std::vector<double>(nnodes, 0.0)));
    r = std::vector<double>(my, 0.0);
    res = std::vector<double>(1, 0.0);
    ui = std::vector<double>(my, 0.0);
    unew = std::vector<double>(my, 0.0);
    xc = std::vector<double>(np, 0.0);
    xm = std::vector<double>(nelemn * nquad, 0.0);
    yc = std::vector<double>(np, 0.0);
    ym = std::vector<double>(nelemn * nquad, 0.0);

    save_times = true;
    save_data = false;
    json = true;
  }
};

double bsp(double xq, double yq, size_t it, size_t iq, size_t id,
  const std::vector<std::vector<size_t>>& node,
  const std::vector<double>& xc, const std::vector<double>& yc) {
  size_t iq1 = iq;
  size_t iq2 = (iq + 1) % 3;
  size_t iq3 = (iq + 2) % 3;
  size_t i1 = node[it][iq1];
  size_t i2 = node[it][iq2];
  size_t i3 = node[it][iq3];
  double d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
  if (id == 1) {
    return 1.0 + ((yc[i2] - yc[i3]) * (xq - xc[i1]) + (xc[i3] - xc[i2]) * (yq - yc[i1])) / d;
  } else if (id == 2) {
    return (yc[i2] - yc[i3]) / d;
  } else if (id == 3) {
    return (xc[i3] - xc[i2]) / d;
  } else {
    std::println("BSP - fatal error!");
    std::println("unknown value of id = {}", id);
    std::exit(1);
  }
  return 0.0;
}

void daxpy_1d(size_t n, double da, const double* dx, ptrdiff_t incx, double* dy, ptrdiff_t incy) {
  if (n == 0 || da == 0.0) { return; }
  ptrdiff_t ix = (incx < 0) ? ((static_cast<ptrdiff_t>(n) - 1) * (-incx)) : 0;
  ptrdiff_t iy = (incy < 0) ? ((static_cast<ptrdiff_t>(n) - 1) * (-incy)) : 0;
  for (size_t i = 0; i < n; i++) {
    dy[iy] += da * dx[ix];
    ix += incx;
    iy += incy;
  }
}

void daxpy_same_matrix(size_t n, double da, std::vector<double>& abd,
  size_t start_row_x, size_t col_x, size_t incx,
  size_t start_row_y, size_t col_y, size_t incy, size_t nrow) {
  if (n == 0 || da == 0.0) { return; }
  size_t rx = start_row_x;
  size_t ry = start_row_y;
  for (size_t i = 0; i < n; i++) {
    double val_x = abd[rx + col_x * nrow];
    abd[ry + col_y * nrow] += da * val_x;
    rx += incx;
    ry += incy;
  }
}

size_t idamax(size_t n, const std::vector<double>& dx, size_t incx);
size_t idamax_matrix(size_t n, const std::vector<double>& abd,
  size_t start_row, size_t col, size_t incx, size_t nrow);
void dscal_matrix(size_t n, double da, std::vector<double>& abd,
  size_t start_row, size_t col, size_t incx, size_t nrow);

void daxpy_diff_matrices(
  size_t n, double da,
  const std::vector<double>& mat_x, size_t start_row_x, size_t col_x, size_t incx,
  std::vector<double>& mat_y, size_t start_row_y, size_t col_y, size_t incy, size_t nrow)
{
  if (n == 0 || da == 0.0) { return; }
  size_t rx = start_row_x;
  size_t ry = start_row_y;
  for (size_t i = 0; i < n; i++) {
    mat_y[ry + col_y * nrow] += da * mat_x[rx + col_x * nrow];
    rx += incx;
    ry += incy;
  }
}

void daxpy_vec_to_matrix(
  size_t n, double da,
  const std::vector<double>& vec_x, size_t start_idx_x, size_t incx,
  std::vector<double>& mat_y, size_t start_row_y, size_t col_y, size_t incy, size_t nrow)
{
  if (n == 0 || da == 0.0) { return; }
  size_t ix = start_idx_x;
  size_t ry = start_row_y;
  for (size_t i = 0; i < n; i++) {
    mat_y[ry + col_y * nrow] += da * vec_x[ix];
    ix += incx;
    ry += incy;
  }
}

void daxpy_matrix_to_vec(
  size_t n, double da,
  const std::vector<double>& mat_x, size_t start_row_x, size_t col_x, size_t incx,
  std::vector<double>& vec_y, size_t start_idx_y, size_t incy, size_t nrow)
{
  if (n == 0 || da == 0.0) { return; }
  size_t rx = start_row_x;
  size_t iy = start_idx_y;
  for (size_t i = 0; i < n; i++) {
    vec_y[iy] += da * mat_x[rx + col_x * nrow];
    rx += incx;
    iy += incy;
  }
}

double ddot_matrix_to_vec(
  size_t n,
  const std::vector<double>& mat_x, size_t start_row_x, size_t col_x, size_t incx,
  const std::vector<double>& vec_y, size_t start_idx_y, size_t incy, size_t nrow)
{
  if (n == 0) { return 0.0; }
  size_t rx = start_row_x;
  size_t iy = start_idx_y;
  double dtemp = 0.0;
  for (size_t i = 0; i < n; i++) {
    dtemp += mat_x[rx + col_x * nrow] * vec_y[iy];
    rx += incx;
    iy += incy;
  }
  return dtemp;
}

void dgbfa(std::vector<double>& abd, size_t _lda, size_t n, size_t ml, size_t mu,
  std::vector<size_t>& ipvt, size_t& info)
{
  double t;
  size_t m = ml + mu;
  size_t nrow = _lda;
  info = 0;
  size_t j0 = mu + 1;
  size_t j1 = (n < m + 1 ? n : m + 1) - 1;
  for (size_t jz = j0; jz < j1; jz++) {
    size_t i0 = m - jz;
    for (size_t i = i0; i < ml; i++) {
      abd[i + jz * nrow] = 0.0;
    }
  }
  size_t jz = j1;
  size_t ju = 0;
  for (size_t k = 0; k < n - 1; k++) {
    jz += 1;
    if (jz < n) {
      for (size_t i = 0; i < ml; i++) {
        abd[i + jz * nrow] = 0.0;
      }
    }
    size_t lm = (ml < (n - 1 - k)) ? ml : (n - 1 - k);
    size_t l = idamax_matrix(lm + 1, abd, m, k, 1, nrow) + m;
    ipvt[k] = l + k - m;
    if (abd[l + k * nrow] == 0.0) {
      info = k;
    } else {
      if (l != m) {
        t = abd[l + k * nrow];
        abd[l + k * nrow] = abd[m + k * nrow];
        abd[m + k * nrow] = t;
      }
      t = -1.0 / abd[m + k * nrow];
      dscal_matrix(lm, t, abd, m + 1, k, 1, nrow);
      size_t ju_next = mu + (size_t)ipvt[k] + 1;
      if (ju < ju_next) { ju = ju_next; }
      if (ju > n) { ju = n; }
      size_t mm = m;
      for (size_t j = k + 1; j < ju; j++) {
        l -= 1;
        mm -= 1;
        double t = abd[l + j * nrow];
        if (l != mm) {
          abd[l + j * nrow] = abd[mm + j * nrow];
          abd[mm + j * nrow] = t;
        }
        daxpy_same_matrix(lm, t, abd, m + 1, k, 1, mm + 1, j, 1, nrow);
      }
    }
  }
  ipvt[n - 1] = n - 1;
  if (abd.size() > m + (n - 1) * nrow) {
    if (abd[m + (n - 1) * nrow] == 0.0) {
      info = n;
    }
  }
}

void dgbsl(std::vector<double>& abd, size_t _lda, size_t n, size_t ml, size_t mu,
  std::vector<size_t>& ipvt, std::vector<double>& b, size_t job)
{
  size_t m = mu + ml;
  size_t nrow = _lda;
  double t;
  if (job == 0) {
    if (0 < ml) {
      for (size_t k = 0; k < n - 1; k++) {
        size_t lm = ml;
        if (lm > n - k - 1) { lm = n - k - 1; }
        size_t l = ipvt[k];
        t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        daxpy_matrix_to_vec(lm, t, abd, m + 1, k, 1, b, k + 1, 1, nrow);
      }
    }
    for (size_t k = n; k > 0; k--) {
      b[k - 1] /= abd[m + (k - 1) * nrow];
      size_t lm = (k - 1 < m ? k - 1 : m);
      size_t la = m - lm;
      size_t lb = (k - 1) - lm;
      t = -b[k - 1];
      daxpy_matrix_to_vec(lm, t, abd, la, k - 1, 1, b, lb, 1, nrow);
    }
  } else {
    for (size_t k = 0; k < n; k++) {
      size_t lm = (k < m ? k : m);
      size_t la = m - lm;
      size_t lb = k - lm;
      double t = ddot_matrix_to_vec(lm, abd, la, k, 1, b, lb, 1, nrow);
      b[k] = (b[k] - t) / abd[m + k * nrow];
    }
    if (0 < ml) {
      for (size_t k = n; k > 0; k--) {
        size_t lm = ml;
        if (lm > n - 1 - (k - 1)) { lm = n - 1 - (k - 1); }
        b[k - 1] += ddot_matrix_to_vec(lm, abd, m + 1, k - 1, 1, b, (k - 1) + 1, 1, nrow);
        size_t l = ipvt[k - 1];
        if (l != k - 1) {
          t = b[l];
          b[l] = b[k - 1];
          b[k - 1] = t;
        }
      }
    }
  }
}

void dscal(size_t n, double da, std::vector<double>& dx, size_t incx) {
  if (n == 0 || incx == 0) { return; }
  if (incx == 1) {
    for (size_t i = 0; i < n; i++) { dx[i] *= da; }
  } else {
    size_t idx = 0;
    for (size_t i = 0; i < n; i++) { dx[idx] *= da; idx += incx; }
  }
}

void dscal_matrix(size_t n, double da, std::vector<double>& abd,
  size_t start_row, size_t col, size_t incx, size_t nrow)
{
  if (n == 0 || incx == 0) { return; }
  size_t current_row = start_row;
  for (size_t i = 0; i < n; i++) {
    abd[current_row + col * nrow] *= da;
    current_row += incx;
  }
}

size_t idamax(size_t n, const std::vector<double>& dx, size_t incx) {
  double max_abs = dx[0] >= 0.0 ? dx[0] : -dx[0];
  size_t max_idx = 0;
  if (incx == 1) {
    for (size_t i = 1; i < n; i++) {
      double val = dx[i] >= 0.0 ? dx[i] : -dx[i];
      if (val > max_abs) { max_abs = val; max_idx = i; }
    }
  } else {
    size_t idx = incx;
    for (size_t i = 1; i < n; i++) {
      double val = dx[idx] >= 0.0 ? dx[idx] : -dx[idx];
      if (val > max_abs) { max_abs = val; max_idx = i; }
      idx += incx;
    }
  }
  return max_idx;
}

size_t idamax_matrix(size_t n, const std::vector<double>& abd,
  size_t start_row, size_t col, size_t incx, size_t nrow)
{
  if (n == 0 || incx == 0) { return 0; }
  double max_abs = abd[start_row + col * nrow] >= 0.0 ? abd[start_row + col * nrow] : -abd[start_row + col * nrow];
  size_t max_idx = 0;
  for (size_t i = 1; i < n; i++) {
    size_t current_row = start_row + i * incx;
    double val = abd[current_row + col * nrow] >= 0.0 ? abd[current_row + col * nrow] : -abd[current_row + col * nrow];
    if (val > max_abs) { max_abs = val; max_idx = i; }
  }
  return max_idx;
}

void timestamp() {
  auto ahora = std::chrono::system_clock::now();
  std::println("{0:%Y-%m-%d %H:%M:%S}", ahora);
}

double ubdry(double y, double para) {
  return 4.0 * para * y * (3.0 - y) / 9.0;
}

size_t igetl(int k, const std::vector<int>& iline, size_t my) {
  for (size_t j = 0; j < my; j++) {
    if (iline[j] == k) return j;
  }
  std::println("IGETL - fatal error!");
  std::println("Unable to get local unknown number for");
  std::println("Global variable number {}", k);
  std::exit(1);
  return 0;
}

void gdump_json(Flow_struct& flow) {
  int j;
  double fval;
  std::ofstream archivo(flow.data_dir + "/" + flow.fileg + ".json");
  archivo << "{\n";
  archivo << "  \"long\": " << flow._long << ",\n";
  archivo << "  \"nelemn\": " << flow.nelemn << ",\n";
  archivo << "  \"np\": " << flow.np << ",\n";
  archivo << "  \"npara\": " << flow.npara << ",\n";
  archivo << "  \"nx\": " << flow.nx << ",\n";
  archivo << "  \"ny\": " << flow.ny << ",\n";
  archivo << "  \"p\": [";
  for (size_t i = 0; i < flow.np; i++) {
    j = flow.insc[i];
    if (j <= 0) { fval = 0.0; } else { fval = flow.f[j - 1]; }
    archivo << fval;
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"h_v_u\": [";
  for (size_t i = 0; i < flow.np; i++) {
    j = flow.indx[i];
    if (j == 0) { fval = 0.0; }
    else if (j < 0) { fval = ubdry(flow.yc[i], flow.para1); }
    else { fval = flow.f[j - 1]; }
    archivo << fval;
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"v_v_v\": [";
  for (size_t i = 0; i < flow.np; i++) {
    j = flow.indx[i + flow.np];
    if (j <= 0) { fval = 0.0; } else { fval = flow.f[j - 1]; }
    archivo << fval;
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"elements\": [";
  for (size_t i = 0; i < flow.np; i++) {
    archivo << "[" << flow.indx[i] << ", " << flow.indx[i + flow.np] << "]";
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"insc\": [";
  for (size_t i = 0; i < flow.np; i++) {
    archivo << flow.insc[i];
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"isotri\": [";
  for (size_t i = 0; i < flow.nelemn; i++) {
    archivo << flow.isotri[i];
    if (i < flow.nelemn - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"node\": [";
  for (size_t i = 0; i < flow.nelemn; i++) {
    archivo << "[" << flow.node[i][0] << ", " << flow.node[i][1] << ", " << flow.node[i][2] << ", "
            << flow.node[i][3] << ", " << flow.node[i][4] << ", " << flow.node[i][5] << "]";
    if (i < flow.nelemn - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"parametro\": " << flow.para1 << ",\n";
  archivo << "  \"reynold\": " << flow.reynld << ",\n";
  archivo << "  \"rjpnew\": " << flow.rjpnew << ",\n";
  archivo << "  \"xc\": [";
  for (size_t i = 0; i < flow.np; i++) {
    archivo << flow.xc[i];
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";
  archivo << "  \"xm\": [";
  for (size_t i = 0; i < flow.np; i++) {
    archivo << flow.yc[i];
    if (i < flow.np - 1) { archivo << ", "; }
  }
  archivo << "]\n";
  archivo << "}\n";
  archivo.flush();
  std::println("GDUMP wrote data set to file {}/{}", flow.data_dir, flow.fileg);
}

bool gdump(Flow_struct& flow) {
  if (flow.json) { gdump_json(flow); return true; }
  int j;
  double fval;
  std::ofstream archivo(flow.data_dir + "/" + flow.fileg + ".dat");
  if (!archivo.is_open()) { return false; }
  archivo << "long: " << flow._long << "\n";
  archivo << "nelemn: " << flow.nelemn << "\n";
  archivo << "np: " << flow.np << "\n";
  archivo << "npara: " << flow.npara << "\n";
  archivo << "nx: " << flow.nx << "\n";
  archivo << "ny: " << flow.ny << "\n";
  for (size_t i = 0; i < flow.np; i++) {
    j = flow.insc[i];
    if (j <= 0) { fval = 0.0; } else { fval = flow.f[j - 1]; }
    archivo << fval << "\n";
  }
  for (size_t i = 0; i < flow.np; i++) {
    j = flow.indx[i];
    if (j == 0) { fval = 0.0; }
    else if (j < 0) { fval = ubdry(flow.yc[i], flow.para1); }
    else { fval = flow.f[j - 1]; }
    archivo << fval << "\n";
  }
  for (size_t i = 0; i < flow.np; i++) {
    j = flow.indx[i + flow.np];
    if (j <= 0) { fval = 0.0; } else { fval = flow.f[j - 1]; }
    archivo << fval << "\n";
  }
  for (size_t i = 0; i < flow.np; i++) {
    archivo << flow.indx[i] << "\n";
    archivo << flow.indx[i + flow.np] << "\n";
  }
  for (size_t i = 0; i < flow.np; i++) { archivo << flow.insc[i] << "\n"; }
  for (size_t i = 0; i < flow.nelemn; i++) { archivo << flow.isotri[i] << "\n"; }
  for (size_t i = 0; i < flow.nelemn; i++) {
    archivo << flow.node[i][0] << ", " << flow.node[i][1] << ", " << flow.node[i][2] << ", "
            << flow.node[i][3] << ", " << flow.node[i][4] << ", " << flow.node[i][5] << "\n";
  }
  archivo << "parametro = " << flow.para1 << "\n";
  archivo << "reynold = " << flow.reynld << "\n";
  archivo << "rjpnew = " << flow.rjpnew << "\n";
  for (size_t i = 0; i < flow.np; i++) { archivo << flow.xc[i] << "\n"; }
  for (size_t i = 0; i < flow.np; i++) { archivo << flow.yc[i] << "\n"; }
  std::println("GDUMP wrote data set to file {}/{}", flow.data_dir, flow.fileg);
  return true;
}

void uv_plot3d(const Flow_struct& flow) {
}

void xy_plot3d_json(const Flow_struct& flow) {
  std::ofstream archivo(flow.data_dir + "/" + flow.filex + ".json");
  size_t ip;
  size_t total_elements = flow.my * flow.mx;
  size_t current_index = 0;
  archivo << "{\n";
  archivo << "  \"mx\": " << flow.mx << ",\n  \"my\": " << flow.my << ",\n";
  if (flow._long) {
    archivo << "  \"xc\": [";
    for (size_t i = 0; i < flow.my; i++) {
      for (size_t j = 0; j < flow.mx; j++) {
        ip = j * flow.my + i;
        archivo << flow.xc[ip];
        if (current_index < total_elements - 1) { archivo << ", "; }
        current_index++;
      }
    }
    archivo << "],\n";
    current_index = 0;
    archivo << "  \"yc\": [";
    for (size_t i = 0; i < flow.my; i++) {
      for (size_t j = 0; j < flow.mx; j++) {
        ip = j * flow.my + i;
        archivo << flow.yc[ip];
        if (current_index < total_elements - 1) { archivo << ", "; }
        current_index++;
      }
    }
    archivo << "]\n";
  } else {
    archivo << "  \"xc\": [";
    for (size_t i = 0; i < flow.mx; i++) {
      for (size_t j = 0; j < flow.my; j++) {
        ip = j * flow.my + i;
        archivo << flow.xc[ip];
        if (current_index < total_elements - 1) { archivo << ", "; }
        current_index++;
      }
    }
    archivo << "],\n";
    current_index = 0;
    archivo << "  \"yc\": [";
    for (size_t i = 0; i < flow.mx; i++) {
      for (size_t j = 0; j < flow.my; j++) {
        ip = j * flow.my + i;
        archivo << flow.yc[ip];
        if (current_index < total_elements - 1) { archivo << ", "; }
        current_index++;
      }
    }
    archivo << "],\n";
  }
  archivo << "}\n";
}

void xy_plot3d(const Flow_struct& flow) {
  if (flow.json) { xy_plot3d_json(flow); return; }
  std::ofstream archivo(flow.data_dir + "/" + flow.filex + ".dat");
  size_t ip;
  if (flow._long) {
    archivo << flow.mx << " " << flow.my << "\n";
    for (size_t i = 0; i < flow.my; i++) {
      for (size_t j = 0; j < flow.mx; j++) { ip = j * flow.my + i; archivo << flow.xc[ip] << "\n"; }
    }
    for (size_t i = 0; i < flow.my; i++) {
      for (size_t j = 0; j < flow.mx; j++) { ip = j * flow.my + i; archivo << flow.yc[ip] << "\n"; }
    }
  } else {
    archivo << flow.my << " " << flow.mx << "\n";
    for (size_t i = 0; i < flow.mx; i++) {
      for (size_t j = 0; j < flow.my; j++) { ip = j * flow.my + i; archivo << flow.xc[ip] << "\n"; }
    }
    for (size_t i = 0; i < flow.mx; i++) {
      for (size_t j = 0; j < flow.my; j++) { ip = j * flow.my + i; archivo << flow.yc[ip] << "\n"; }
    }
  }
  std::println("XYDUMP wrote data set to file.");
}

void getg(const std::vector<double>& f, const std::vector<int>& iline, size_t my, size_t /*neqn*/, std::vector<double>& u) {
  for (size_t j = 0; j < my; j++) {
    int k = iline[j];
    if (k > 0) {
      u[j] = f[(size_t)(k - 1)];
    } else {
      u[j] = 0.0;
    }
  }
}

void qbf(double x, double y, size_t it, size_t _in,
  double& bb, double& bx, double& by,
  const std::vector<std::vector<size_t>>& node,
  const std::vector<double>& xc, const std::vector<double>& yc) {
  if (_in <= 2) {
    size_t in1 = _in;
    size_t in2 = (_in + 1) % 3;
    size_t in3 = (_in + 2) % 3;
    size_t i1 = node[it][in1];
    size_t i2 = node[it][in2];
    size_t i3 = node[it][in3];
    double d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    double t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    bb = t * (2.0 * t - 1.0);
    bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    size_t inn = _in - 3;
    size_t in1 = inn;
    size_t in2 = (inn + 1) % 3;
    size_t in3 = (inn + 2) % 3;
    size_t i1 = node[it][in1];
    size_t i2 = node[it][in2];
    size_t i3 = node[it][in3];
    size_t j1 = i2;
    size_t j2 = i3;
    size_t j3 = i1;
    double d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    double c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    double t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    double s = 1.0 + ((yc[j2] - yc[j3]) * (x - xc[j1]) + (xc[j3] - xc[j2]) * (y - yc[j1])) / c;
    bb = 4.0 * s * t;
    bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
}

void setlin(Flow_struct& flow) {
  size_t itemp = (size_t)std::llround((2.0 * (double)(flow.nx - 1) * 9.0) / flow.xlngth);
  flow.nodex0 = flow._long ? itemp * (2 * flow.ny - 1) : itemp;
  for (size_t i = 0; i < flow.my; i++) {
    size_t ip = flow._long ? flow.nodex0 + i : flow.nodex0 + flow.mx * i;
    flow.iline[i] = flow.indx[ip];
  }
  if (1 <= flow.iwrite) {
    std::println(" \nSETLIN: unknown numbers along line:\n ");
    for (size_t i = 0; i < flow.my; i++) { std::print(" {}", flow.iline[i]); }
    std::println();
  }
}

void setqud(Flow_struct& flow) {
  for (size_t it = 0; it < flow.nelemn; it++) {
    size_t ip1 = flow.node[it][0];
    size_t ip2 = flow.node[it][1];
    size_t ip3 = flow.node[it][2];
    double x1 = flow.xc[ip1]; double x2 = flow.xc[ip2]; double x3 = flow.xc[ip3];
    double y1 = flow.yc[ip1]; double y2 = flow.yc[ip2]; double y3 = flow.yc[ip3];
    flow.xm[it + 0 * flow.nelemn] = 0.5 * (x1 + x2);
    flow.xm[it + 1 * flow.nelemn] = 0.5 * (x2 + x3);
    flow.xm[it + 2 * flow.nelemn] = 0.5 * (x3 + x1);
    flow.ym[it + 0 * flow.nelemn] = 0.5 * (y1 + y2);
    flow.ym[it + 1 * flow.nelemn] = 0.5 * (y2 + y3);
    flow.ym[it + 2 * flow.nelemn] = 0.5 * (y3 + y1);
    flow.area[it] = 0.5 * std::abs((y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3));
  }
}

void setxy(Flow_struct& flow) {
  int ic, jc;
  for (size_t ip = 0; ip < flow.np; ip++) {
    if (flow._long) {
      ic = (int)(ip / flow.my) + 1;
      jc = (int)(ip % flow.my) + 1;
    } else {
      ic = (int)(ip % flow.mx) + 1;
      jc = (int)(ip / flow.mx) + 1;
    }
    flow.xc[ip] = (double)(ic - 1) * flow.xlngth / (double)(2 * flow.nx - 2);
    flow.yc[ip] = (double)(jc - 1) * flow.ylngth / (double)(2 * flow.ny - 2);
  }
}

void setbas(Flow_struct& flow) {
  double bb = 0.0, bx = 0.0, by = 0.0;
  for (size_t it = 0; it < flow.nelemn; it++) {
    for (size_t j = 0; j < flow.nquad; j++) {
      double x = flow.xm[it + j * flow.nelemn];
      double y = flow.ym[it + j * flow.nelemn];
      for (size_t iq = 0; iq < flow.nnodes; iq++) {
        flow.psi[it][j][iq] = bsp(x, y, it, iq, 1, flow.node, flow.xc, flow.yc);
        qbf(x, y, it, iq, bb, bx, by, flow.node, flow.xc, flow.yc);
        flow.phi[it][j][iq][0] = bb;
        flow.phi[it][j][iq][1] = bx;
        flow.phi[it][j][iq][2] = by;
      }
    }
  }
}

void setban(Flow_struct& flow) {
  flow.nlband = 0;
  for (size_t it = 0; it < flow.nelemn; it++) {
    for (size_t iq = 0; iq < flow.nnodes; iq++) {
      size_t ip = flow.node[it][iq];
      for (size_t iuk = 0; iuk < 3; iuk++) {
        int i = (iuk == 2) ? flow.insc[ip] : flow.indx[ip + iuk * flow.np];
        if (i > 0) {
          for (size_t iqq = 0; iqq < flow.nnodes; iqq++) {
            size_t ipp = flow.node[it][iqq];
            for (size_t iukk = 0; iukk < 3; iukk++) {
              int j = (iukk == 2) ? flow.insc[ipp] : flow.indx[ipp + iukk * flow.np];
              int diff = j - i;
              if ((int)flow.nlband < diff) flow.nlband = (size_t)diff;
            }
          }
        }
      }
    }
  }
  flow.nband = flow.nlband + flow.nlband + 1;
  flow.nrow = flow.nlband + flow.nlband + flow.nlband + 1;
  std::println("Lower bandwidth = {}", flow.nlband);
  std::println("Total bandwidth = {}", flow.nband);
  std::println("NROW  = {}", flow.nrow);
  if (flow.maxrow < flow.nrow) {
    std::println("SETBAN - NROW is too large!");
    std::println("The maximum allowed is {}", flow.maxrow);
    std::exit(1);
  }
}

void setgrd(Flow_struct& flow);
void nstoke(Flow_struct& flow);
void resid(Flow_struct& flow);
void gram(Flow_struct& flow);
bool gdump(Flow_struct& flow);
void xy_plot3d(const Flow_struct& flow);
void uv_plot3d(const Flow_struct& flow);
bool xy_table(const Flow_struct& flow);
bool uv_table(const Flow_struct& flow);
void linsys(std::vector<double>& a, const std::vector<double>& area,
  std::vector<double>& g, const std::vector<double>& f,
  const std::vector<int>& indx, const std::vector<int>& insc,
  std::vector<size_t>& ipivot,
  size_t maxrow, size_t nelemn, size_t neqn, size_t nlband, size_t nnodes,
  const std::vector<std::vector<size_t>>& node,
  size_t np, size_t nquad, size_t nrow, double para1, double para2,
  const std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
  const std::vector<std::vector<std::vector<double>>>& psi,
  double reynld, const std::vector<double>& yc);

int main() {
  std::filesystem::create_directories("data");
  Flow_struct flow(21, 7);
  auto inicio = std::chrono::high_resolution_clock::now();

  timestamp();
  std::println(" ");
  std::println("NX = {}", flow.nx);
  std::println("NY = {}", flow.ny);
  std::println("Number of elements = {}", flow.nelemn);
  std::println("Reynolds number = {:g}", flow.reynld);
  std::println("Secant tolerance = {:g}", flow.tolsec);
  std::println("Newton tolerance = {:g}", flow.tolnew);
  std::println(" ");
  setgrd(flow);
  flow.f.assign(flow.neqn, 0.0);
  flow.g.assign(flow.neqn, 1.0);
  flow.res.assign(flow.neqn, 0.0);
  setban(flow);
  flow.a.assign(flow.nrow * flow.neqn, 0.0);
  setlin(flow);
  setxy(flow);
  setqud(flow);
  setbas(flow);
  flow.para1 = 1.0;
  flow.para2 = 1.0;
  std::println(" ");
  std::println("Solve Navier Stokes problem with parameter = {:g}", flow.para1);
  std::println("for profile at x = {:g}", flow.xc[flow.nodex0]);
  nstoke(flow);
  if (1 <= flow.iwrite) {
    resid(flow);
  }
  getg(flow.f, flow.iline, flow.my, flow.neqn, flow.ui);
  if (1 <= flow.iwrite) {
    std::println("\nU profile:\n");
    for (size_t i = 0; i < flow.my; i++) { std::print("{:g}, ", flow.ui[i]); }
    std::println(" ");
  }
  gram(flow);
  if (flow.save_data) {
    std::println("Writing graphics data to file {}", flow.fileg);
    flow.rjpnew = 0.0;
    if (!gdump(flow)) { std::println("Error al ejecutar gdump"); return 1; }
    xy_plot3d(flow);
    uv_plot3d(flow);
  } else {
    xy_table(flow);
    uv_table(flow);
  }
  std::fill(flow.f.begin(), flow.f.end(), 0.0);
  std::fill(flow.g.begin(), flow.g.end(), 0.0);
  flow.aold = 0.0;
  flow.rjpold = 0.0;
  flow.anew = 0.1;
  double temp;
  for (size_t iter = 1; iter <= flow.maxsec; iter++) {
    flow.numsec += 1;
    std::println(" ");
    std::println("Secant iteration {}", iter);
    std::println(" ");
    std::println("Solving Navier Stokes problem for parameter = {:g}", flow.anew);
    std::copy(flow.f.begin(), flow.f.end(), flow.g.begin());
    flow.para1 = flow.anew;
    flow.para2 = flow.anew;
    nstoke(flow);
    getg(flow.f, flow.iline, flow.my, flow.neqn, flow.unew);
    if (1 <= flow.iwrite) {
      std::println(" \nVelocity profile:\n ");
      for (size_t i = 0; i < flow.my; i++) { std::print("{:g}, ", flow.unew[i]); }
    }
    flow.para1 = flow.anew;
    flow.para2 = 1.0;
    linsys(flow.a, flow.area, flow.g, flow.f,
      flow.indx, flow.insc, flow.ipivot,
      flow.maxrow, flow.nelemn, flow.neqn, flow.nlband, flow.nnodes, flow.node,
      flow.np, flow.nquad, flow.nrow, flow.para1, flow.para2,
      flow.phi, flow.psi, flow.reynld, flow.yc);
    getg(flow.g, flow.iline, flow.my, flow.neqn, flow.dcda);
    if (2 <= flow.iwrite) {
      std::println("\n \nSensitivities:\n ");
      for (size_t i = 0; i < flow.my; i++) { std::print("{:g}, ", flow.dcda[i]); }
    }
    flow.rjpnew = 0.0;
    for (size_t i = 0; i < flow.my; i++) {
      temp = -flow.r[i];
      for (size_t j = 0; j < flow.my; j++) { temp += flow.gr[i + j * flow.my] * flow.unew[j]; }
      flow.rjpnew += 2.0 * flow.dcda[i] * temp;
    }
    std::println("\n \nParameter  = {:g}, J prime = {:g}", flow.anew, flow.rjpnew);
    if (flow.save_data) {
      flow.para1 = flow.anew;
      gdump(flow);
    }
    if (iter == 1) { flow.a2 = 0.5; }
    else { flow.a2 = flow.aold - flow.rjpold * (flow.anew - flow.aold) / (flow.rjpnew - flow.rjpold); }
    flow.aold = flow.anew;
    flow.anew = flow.a2;
    flow.rjpold = flow.rjpnew;
    flow.test = std::abs(flow.anew - flow.aold) / std::abs(flow.anew);
    std::println("New value of parameter = {:g}", flow.anew);
    std::println("Convergence test = {:g}", flow.test);
    if (std::abs(flow.anew - flow.aold) < std::abs(flow.anew) * flow.tolsec) {
      flow.converged = true;
      break;
    }
  }
  if (flow.converged) { std::println("Secant iteration converged."); }
  else { std::println("Secant iteration failed to converge."); }
  std::println("Number of secant steps = {}", flow.numsec);
  std::println("Number of Newton steps = {}", flow.numnew);
  std::println(" \nCHANNEL:");
  std::println("  Normal of execution.");
  std::println(" ");
  timestamp();
  auto duracion = std::chrono::high_resolution_clock::now() - inicio;
  auto secs = std::chrono::duration<double>(duracion).count();
  std::println("Tiempo transcurrido: {:f} seconds", secs);
  if (flow.save_times) {
    std::ofstream times_file("times.txt", std::ios::app);
    if (times_file.is_open()) { times_file << std::format("{:f}\n", secs); }
  }
  return 0;
}

void uval(
  const std::vector<double>& g,
  const std::vector<int>& indx,
  size_t iquad, size_t it, size_t nelemn, size_t neqn,
  size_t nnodes, const std::vector<std::vector<size_t>>& node,
  size_t np, size_t nquad, double para,
  const std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
  std::vector<double>& un, std::vector<double>& unx, std::vector<double>& uny,
  const std::vector<double>& yc)
{
  un[0] = 0.0; un[1] = 0.0;
  unx[0] = 0.0; unx[1] = 0.0;
  uny[0] = 0.0; uny[1] = 0.0;

  for (size_t iq = 0; iq < nnodes; iq++) {
    size_t ip = node[it][iq];
    double bb = phi[it][iquad][iq][0];
    double bx = phi[it][iquad][iq][1];
    double by = phi[it][iquad][iq][2];

    for (int iuk = 0; iuk < 2; iuk++) {
      int iun = indx[ip + iuk * np];

      if (0 < iun) {
        un[iuk] += bb * g[(size_t)(iun - 1)];
        unx[iuk] += bx * g[(size_t)(iun - 1)];
        uny[iuk] += by * g[(size_t)(iun - 1)];
      } else if (iun < 0) {
        double ubc = ubdry(yc[ip], para);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }
}

void linsys(
  std::vector<double>& a,
  const std::vector<double>& area,
  std::vector<double>& g,
  const std::vector<double>& f,
  const std::vector<int>& indx,
  const std::vector<int>& insc,
  std::vector<size_t>& ipivot,
  size_t maxrow, size_t nelemn, size_t neqn, size_t nlband, size_t nnodes,
  const std::vector<std::vector<size_t>>& node,
  size_t np, size_t nquad, size_t nrow, double para1, double para2,
  const std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
  const std::vector<std::vector<std::vector<double>>>& psi,
  double reynld, const std::vector<double>& yc)
{
  size_t info = 0;
  int ioff = (int)(nlband + nlband);
  double visc = 1.0 / reynld;
  std::vector<double> un(2, 0.0);
  std::vector<double> unx(2, 0.0);
  std::vector<double> uny(2, 0.0);

  std::fill(g.begin(), g.end(), 0.0);
  std::fill(a.begin(), a.end(), 0.0);

  for (size_t it = 0; it < nelemn; it++) {
    double ar = area[it] / 3.0;
    for (size_t iquad = 0; iquad < nquad; iquad++) {
      uval(f, indx, iquad, it, nelemn, neqn, nnodes, node, np, nquad, para1, phi,
        un, unx, uny, yc);

      for (size_t iq = 0; iq < nnodes; iq++) {
        size_t ip = node[it][iq];
        double bb = phi[it][iquad][iq][0];
        double bx = phi[it][iquad][iq][1];
        double by = phi[it][iquad][iq][2];
        double bbl = psi[it][iquad][iq];
        int ihor = indx[ip] - 1;
        int iver = indx[ip + np] - 1;
        int iprs = insc[ip] - 1;

        if (0 <= ihor) { g[(size_t)ihor] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]); }
        if (0 <= iver) { g[(size_t)iver] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]); }

        for (size_t iqq = 0; iqq < nnodes; iqq++) {
          size_t ipp = node[it][iqq];
          double bbb = phi[it][iquad][iqq][0];
          double bbx = phi[it][iquad][iqq][1];
          double bby = phi[it][iquad][iqq][2];
          double bbbl = psi[it][iquad][iqq];
          int ju = indx[ipp] - 1;
          int jv = indx[ipp + np] - 1;
          int jp = insc[ipp] - 1;

          if (0 <= ju) {
            if (0 <= ihor) {
              int iuse = ihor - ju + ioff;
              a[(size_t)iuse + (size_t)ju * nrow] += ar * (visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }
            if (0 <= iver) {
              int iuse = iver - ju + ioff;
              a[(size_t)iuse + (size_t)ju * nrow] += ar * bb * bbb * unx[1];
            }
            if (0 <= iprs) {
              int iuse = iprs - ju + ioff;
              a[(size_t)iuse + (size_t)ju * nrow] += ar * bbx * bbl;
            }
          } else if (ju == -2) {
            double uu = ubdry(yc[ipp], para2);
            if (0 <= ihor) {
              g[(size_t)ihor] -= ar * uu * (visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }
            if (0 <= iver) { g[(size_t)iver] -= ar * uu * bb * bbb * unx[1]; }
            if (0 <= iprs) { g[(size_t)iprs] -= ar * uu * bbx * bbl; }
          }

          if (0 <= jv) {
            if (0 <= ihor) {
              int iuse = ihor - jv + ioff;
              a[(size_t)iuse + (size_t)jv * nrow] += ar * bb * bbb * uny[0];
            }
            if (0 <= iver) {
              int iuse = iver - jv + ioff;
              a[(size_t)iuse + (size_t)jv * nrow] += ar * (visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
            }
            if (0 <= iprs) {
              int iuse = iprs - jv + ioff;
              a[(size_t)iuse + (size_t)jv * nrow] += ar * bby * bbl;
            }
          }

          if (0 <= jp) {
            if (0 <= ihor) {
              int iuse = ihor - jp + ioff;
              a[(size_t)iuse + (size_t)jp * nrow] -= ar * bx * bbbl;
            }
            if (0 <= iver) {
              int iuse = iver - jp + ioff;
              a[(size_t)iuse + (size_t)jp * nrow] -= ar * by * bbbl;
            }
          }
        }
      }
    }
  }

  g[neqn - 1] = 0.0;
  for (size_t j = neqn - nlband - 1; j < neqn - 1; j++) {
    size_t i = (size_t)((int)(neqn - 1) - (int)j + ioff);
    a[i + j * nrow] = 0.0;
  }
  a[(size_t)ioff + (neqn - 1) * nrow] = 1.0;

  dgbfa(a, nrow, neqn, nlband, nlband, ipivot, info);
  if (info != 0) {
    std::println(" \nLINSYS - fatal error!\nDGBFA returns INFO = {}", info);
    std::exit(1);
  }
  size_t job = 0;
  dgbsl(a, nrow, neqn, nlband, nlband, ipivot, g, job);
}

void nstoke(Flow_struct& flow) {
  for (size_t iter = 1; iter <= flow.maxnew; iter++) {
    flow.numnew += 1;
    linsys(flow.a, flow.area, flow.f, flow.g,
      flow.indx, flow.insc, flow.ipivot,
      flow.maxrow, flow.nelemn, flow.neqn, flow.nlband, flow.nnodes, flow.node,
      flow.np, flow.nquad, flow.nrow, flow.para1, flow.para2,
      flow.phi, flow.psi, flow.reynld, flow.yc);

    for (size_t i = 0; i < flow.neqn; i++) {
      flow.g[i] -= flow.f[i];
    }

    double diff = std::abs(flow.g[idamax(flow.neqn, flow.g, 1)]);

    if (1 <= flow.iwrite) {
      std::println("NSTOKE iteration {} Mnorm = {:g}", iter, diff);
    }

    std::copy(flow.f.begin(), flow.f.end(), flow.g.begin());

    if (diff <= flow.tolnew) {
      std::println("Navier Stokes solution converged in {} iteration.", iter);
      return;
    }
  }
  std::println("Navier Stokes solution did not converge!");
}

void uv_table_json(const Flow_struct& flow) {
  std::string path = flow.data_dir + "/" + flow.fileu + ".json";
  std::ofstream archivo(path);
  int k;
  double uval_v, vval_v;

  archivo << "{\n";
  archivo << "  \"u_vals\": [";
  for (size_t ip = 0; ip < flow.np; ip++) {
    k = flow.indx[ip];
    if (k == 0) { uval_v = 0.0; }
    else if (k < 0) { uval_v = ubdry(flow.yc[ip], flow.para1); }
    else { uval_v = flow.f[(size_t)(k - 1)]; }
    archivo << uval_v;
    if (ip < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";

  archivo << "  \"v_vals\": [";
  for (size_t ip = 0; ip < flow.np; ip++) {
    k = flow.indx[ip + flow.np];
    if (k == 0) { vval_v = 0.0; }
    else { vval_v = flow.f[(size_t)(k - 1)]; }
    archivo << vval_v;
    if (ip < flow.np - 1) { archivo << ", "; }
  }
  archivo << "]\n";
  archivo << "}";
  archivo.flush();
}

bool uv_table(const Flow_struct& flow) {
  if (flow.json) {
    uv_table_json(flow);
    return true;
  }

  std::string path = flow.data_dir + "/" + flow.fileu + ".dat";
  std::ofstream archivo(path);
  int k;
  double uval_v, vval_v;

  for (size_t ip = 0; ip < flow.np; ip++) {
    k = flow.indx[ip];
    if (k == 0) { uval_v = 0.0; }
    else if (k < 0) { uval_v = ubdry(flow.yc[ip], flow.para1); }
    else { uval_v = flow.f[(size_t)(k - 1)]; }

    k = flow.indx[ip + flow.np];
    if (k == 0) { vval_v = 0.0; }
    else { vval_v = flow.f[(size_t)(k - 1)]; }

    archivo << uval_v << ", " << vval_v << "\n";
  }
  archivo.flush();
  return true;
}

void xy_table_json(const Flow_struct& flow) {
  std::string path = flow.data_dir + "/" + flow.filex + ".json";
  std::ofstream archivo(path);

  archivo << "{\n";
  archivo << "  \"xc\": [";
  for (size_t ip = 0; ip < flow.np; ip++) {
    archivo << flow.xc[ip];
    if (ip < flow.np - 1) { archivo << ", "; }
  }
  archivo << "],\n";

  archivo << "  \"yc\": [";
  for (size_t ip = 0; ip < flow.np; ip++) {
    archivo << flow.yc[ip];
    if (ip < flow.np - 1) { archivo << ", "; }
  }
  archivo << "]\n";
  archivo << "}";
  archivo.flush();
}

bool xy_table(const Flow_struct& flow) {
  if (flow.json) {
    xy_table_json(flow);
    return true;
  }

  std::string path = flow.data_dir + "/" + flow.filex + ".dat";
  std::ofstream archivo(path);

  for (size_t ip = 0; ip < flow.np; ip++) {
    archivo << flow.xc[ip] << ", " << flow.yc[ip] << "\n";
  }
  archivo.flush();
  return true;
}

void setgrd(Flow_struct& flow) {
  if (flow.ny < flow.nx) {
    flow._long = true;
    std::println("Using vertical ordering.");
  } else {
    flow._long = false;
    std::println("Using horizontal ordering.");
  }
  std::println(" \nSETGRD: Taylor Hood element");
  flow.neqn = 0;
  size_t ielemn = 0;
  size_t ic;
  size_t jc;
  for (size_t ip = 0; ip < flow.np; ip++) {
    if (flow._long) {
      ic = (ip / flow.my) + 1;
      jc = (ip % flow.my) + 1;
    } else {
      ic = (ip % flow.mx) + 1;
      jc = (ip / flow.mx) + 1;
    }
    size_t icnt = ic % 2;
    size_t jcnt = jc % 2;
    if (icnt == 1 && jcnt == 1 && ic != flow.mx && jc != flow.my) {
      if (flow._long) {
        size_t ip1 = ip + flow.my;
        size_t ip2 = ip + flow.my + flow.my;
        flow.node[ielemn][0] = ip;
        flow.node[ielemn][1] = ip2 + 2;
        flow.node[ielemn][2] = ip2;
        flow.node[ielemn][3] = ip1 + 1;
        flow.node[ielemn][4] = ip2 + 1;
        flow.node[ielemn][5] = ip1;
        flow.isotri[ielemn] = 0;
        ielemn = ielemn + 1;
        flow.node[ielemn][0] = ip;
        flow.node[ielemn][1] = ip + 2;
        flow.node[ielemn][2] = ip2 + 2;
        flow.node[ielemn][3] = ip + 1;
        flow.node[ielemn][4] = ip1 + 2;
        flow.node[ielemn][5] = ip1 + 1;
        flow.isotri[ielemn] = 0;
        ielemn = ielemn + 1;
      } else {
        size_t ip1 = ip + flow.mx;
        size_t ip2 = ip + flow.mx + flow.mx;
        ielemn += 1;
        flow.node[ielemn][1] = ip - 1;
        flow.node[ielemn][2] = ip2 - 1;
        flow.node[ielemn][3] = ip2 + 1;
        flow.node[ielemn][4] = ip1 - 1;
        flow.node[ielemn][5] = ip2;
        flow.node[ielemn][6] = ip1;
        flow.isotri[ielemn] = 0;
        ielemn += 1;
        flow.node[ielemn][1] = ip - 1;
        flow.node[ielemn][2] = ip2 + 1;
        flow.node[ielemn][3] = ip + 1;
        flow.node[ielemn][4] = ip1;
        flow.node[ielemn][5] = ip1 + 1;
        flow.node[ielemn][6] = ip + 1 - 1;
        flow.isotri[ielemn] = 0;
      };
    }
    if (ic == 1 && 1 < jc && jc < flow.my) {
      flow.indx[ip] = -1;
      flow.indx[ip + flow.np] = 0;
    } else if (ic == flow.mx && 1 < jc && jc < flow.my) {
      flow.neqn += 1;
      flow.indx[ip] = (int)flow.neqn;
      flow.indx[ip + flow.np] = 0;
    } else if (jc == 1 || jc == flow.my) {
      flow.indx[ip] = 0;
      flow.indx[ip + flow.np] = 0;
    } else {
      flow.neqn += 2;
      flow.indx[ip] = (int)flow.neqn - 1;
      flow.indx[ip + flow.np] = (int)flow.neqn;
    }
    if (jcnt == 1 && icnt == 1) {
      flow.neqn = flow.neqn + 1;
      flow.insc[ip] = (int)flow.neqn;
    } else {
      flow.insc[ip] = 0;
    }
  }
  if (2 <= flow.iwrite) {
    std::println(" \n    I     flow.indx 1 & 2,flow.insc\n ");
    for (size_t i = 0; i < flow.np; i++) {
      std::println("     {}, {}, {}, {}", i + 1, flow.indx[i], flow.indx[i + flow.np], flow.insc[i]);
    }
    std::println(" \n    IT    NODE(IT,1:6)\n ");
    for (size_t it = 0; it < flow.nelemn; it++) {
      std::println("    {}, {}, {}, {}, {}, {}, {}", it + 1,
        flow.node[it][0] + 1, flow.node[it][1] + 1, flow.node[it][2] + 1,
        flow.node[it][3] + 1, flow.node[it][4] + 1, flow.node[it][5] + 1);
    }
  }
  std::println("Number of unknowns = {}", flow.neqn);
  if (flow.maxeqn < flow.neqn) {
    std::println("SETGRD - Too many unknowns!");
    std::println("The maximum allowed is MAXEQN = {}", flow.maxeqn);
    std::println("This problem requires NEQN = {}", flow.neqn);
    std::exit(1);
  }
}

void gram(Flow_struct& flow) {
  double wt[3];
  double yq[3];
  wt[0] = 5.0 / 9.0; wt[1] = 8.0 / 9.0; wt[2] = wt[0];
  yq[0] = -0.7745966692; yq[1] = 0.0; yq[2] = -yq[0];
  int valid_nodes[3] = {0, 1, 3};
  double bb = 0.0; double bx = 0.0; double by = 0.0;
  double bbb = 0.0; double bbx = 0.0; double bby = 0.0;
  for (size_t i = 0; i < flow.my; i++) {
    flow.r[i] = 0.0;
    for (size_t j = 0; j < flow.my; j++) { flow.gr[i + j * flow.my] = 0.0; }
  }
  double xzero = flow.xc[flow.nodex0];
  for (size_t it = 0; it < flow.nelemn; it++) {
    size_t k = flow.node[it][0];
    size_t kk = flow.node[it][1];
    if (1.0e-04 < std::abs(flow.xc[k] - xzero) || 1.0e-04 < std::abs(flow.xc[kk] - xzero)) { continue; }
    for (size_t iquad = 0; iquad < 3; iquad++) {
      double bma2 = (flow.yc[kk] - flow.yc[k]) / 2.0;
      double ar = bma2 * wt[iquad];
      double x = xzero;
      double y = flow.yc[k] + bma2 * (yq[iquad] + 1.0);
      double uiqdpt = 0.0;
      for (int vi = 0; vi < 3; vi++) {
        int iq = valid_nodes[vi];
        qbf(x, y, it, (size_t)iq, bb, bx, by, flow.node, flow.xc, flow.yc);
        size_t ip = flow.node[it][(size_t)iq];
        int iun = flow.indx[ip];
        if (0 < iun) {
          size_t ii = igetl(iun, flow.iline, flow.my);
          uiqdpt = uiqdpt + bb * flow.ui[ii];
        } else {
          double ubc_val = ubdry(flow.yc[ip], flow.para1);
          uiqdpt += bb * ubc_val;
        }
      }
      for (int vi = 0; vi < 3; vi++) {
        int iq = valid_nodes[vi];
        size_t ip = flow.node[it][(size_t)iq];
        qbf(x, y, it, (size_t)iq, bb, bx, by, flow.node, flow.xc, flow.yc);
        int i = flow.indx[ip];
        if (i <= 0) { continue; }
        size_t ii = igetl(i, flow.iline, flow.my);
        flow.r[ii] += bb * uiqdpt * ar;
        for (int vj = 0; vj < 3; vj++) {
          int iqq = valid_nodes[vj];
          size_t ipp = flow.node[it][(size_t)iqq];
          qbf(x, y, it, (size_t)iqq, bbb, bbx, bby, flow.node, flow.xc, flow.yc);
          int j = flow.indx[ipp];
          if (j != 0) {
            size_t jj = igetl(j, flow.iline, flow.my);
            flow.gr[ii + jj * flow.my] += bb * bbb * ar;
          }
        }
      }
    }
  }
  if (2 <= flow.iwrite) {
    std::println("\nGram matrix:\n");
    for (size_t i = 0; i < flow.my; i++) {
      for (size_t j = 0; j < flow.my; j++) { std::println("{}, {}, {:g}", i + 1, j + 1, flow.gr[i + j * flow.my]); }
    }
    std::println(" \nR vector:\n ");
    for (size_t i = 0; i < flow.my; i++) { std::println("{}, {:g}", i, flow.r[i]); }
  }
}

void resid(Flow_struct& flow) {
  double visc = 1.0 / flow.reynld;
  std::fill(flow.res.begin(), flow.res.end(), 0.0);
  std::vector<double> un(2, 0.0);
  std::vector<double> unx(2, 0.0);
  std::vector<double> uny(2, 0.0);
  for (size_t it = 0; it < flow.nelemn; it++) {
    double ar = flow.area[it] / 3.0;
    for (size_t iquad = 0; iquad < flow.nquad; iquad++) {
      uval(flow.g, flow.indx, iquad, it, flow.nelemn, flow.neqn,
        flow.nnodes, flow.node, flow.np, flow.nquad, flow.para1, flow.phi,
        un, unx, uny, flow.yc);
      for (size_t iq = 0; iq < flow.nnodes; iq++) {
        size_t ip = flow.node[it][iq];
        double bb = flow.phi[it][iquad][iq][0];
        double bx = flow.phi[it][iquad][iq][1];
        double by = flow.phi[it][iquad][iq][2];
        double bbl = flow.psi[it][iquad][iq];
        int iprs = flow.insc[ip] - 1;
        int ihor = flow.indx[ip] - 1;
        int iver = flow.indx[ip + flow.np] - 1;
        if (0 <= ihor) { flow.res[ihor] += (un[0] * unx[0] + un[1] * uny[0]) * bb * ar; }
        if (0 <= iver) { flow.res[iver] += (un[0] * unx[1] + un[1] * uny[1]) * bb * ar; }
        for (size_t iqq = 0; iqq < flow.nnodes; iqq++) {
          size_t ipp = flow.node[it][iqq];
          double bbb = flow.phi[it][iquad][iqq][0];
          double bbx = flow.phi[it][iquad][iqq][1];
          double bby = flow.phi[it][iquad][iqq][2];
          double bbbl = flow.psi[it][iquad][iqq];
          int ju = flow.indx[ipp] - 1;
          int jv = flow.indx[ipp + flow.np] - 1;
          int jp = flow.insc[ipp] - 1;
          if (0 <= ju) {
            if (0 <= ihor) {
              double aijuu = visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              flow.res[ihor] += aijuu * ar * flow.g[ju];
            }
            if (0 <= iver) {
              double aijvu = bb * bbb * unx[1];
              flow.res[iver] += aijvu * ar * flow.g[ju];
            }
            if (0 <= iprs) {
              double aijpu = bbx * bbl;
              flow.res[iprs] += aijpu * ar * flow.g[ju];
            }
          } else if (ju == -2) {
            double uu = ubdry(flow.yc[ipp], flow.para1);
            if (0 <= ihor) {
              double aijuu = visc * (by * bby + bx * bbx) + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              flow.res[ihor] += ar * aijuu * uu;
            }
            if (0 <= iver) {
              double aijvu = bb * bbb * unx[1];
              flow.res[iver] += ar * aijvu * uu;
            }
            if (0 <= iprs) {
              double aijpu = bbx * bbl;
              flow.res[iprs] += ar * aijpu * uu;
            }
          }
          if (0 <= jv) {
            if (0 <= ihor) {
              double aijuv = bb * bbb * uny[0];
              flow.res[ihor] += aijuv * ar * flow.g[jv];
            }
            if (0 <= iver) {
              double aijvv = visc * (by * bby + bx * bbx) + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]);
              flow.res[iver] += aijvv * ar * flow.g[jv];
            }
            if (0 <= iprs) {
              double aijpv = bby * bbl;
              flow.res[iprs] += aijpv * ar * flow.g[jv];
            }
          }
          if (0 <= jp) {
            if (0 <= ihor) {
              double aijup = -bx * bbbl;
              flow.res[ihor] += aijup * ar * flow.g[jp];
            }
            if (0 <= iver) {
              double aijvp = -by * bbbl;
              flow.res[iver] += aijvp * ar * flow.g[jp];
            }
          }
        }
      }
    }
  }
  flow.res[flow.neqn - 1] = flow.g[flow.neqn - 1];
  double rmax = 0.0;
  size_t imax = 0;
  size_t ibad = 0;
  for (size_t i = 0; i < flow.neqn; i++) {
    double test = std::abs(flow.res[i]);
    if (rmax < test) { rmax = test; imax = i; }
    if (1.0E-03 < test) { ibad = ibad + 1; }
  }
  if (flow.iwrite >= 1) {
    std::println("\nRESIDUAL INFORMATION:\n");
    std::println("Worst residual is number {}", imax);
    std::println("of magnitude {:g}\n", rmax);
    std::println("Number of \"bad\" residuals is {} out of {}\n", ibad, flow.neqn);
  }
  if (flow.iwrite >= 2) {
    std::println("Raw residuals:\n");
    int i = 0;
    for (size_t j = 0; j < flow.np; j++) {
      int j_display = (int)(j + 1);
      if (flow.indx[j] > 0) {
        double val = flow.res[i];
        if (std::abs(val) <= 1e-3) { std::println(" U{:5d}, {:5d}, {:e}", i, j_display, val); }
        else { std::println("*U{:5d}, {:5d}, {:e}", i, j_display, val); }
        i += 1;
      }
      if (flow.indx[j + flow.np] > 0) {
        double val = flow.res[i];
        if (std::abs(val) <= 1e-3) { std::println(" V{:5d}, {:5d}, {:e}", i, j_display, val); }
        else { std::println("*V{:5d}, {:5d}, {:e}", i, j_display, val); }
        i += 1;
      }
      if (flow.insc[j] > 0) {
        double val = flow.res[i];
        if (std::abs(val) <= 1e-3) { std::println(" P{:5d}, {:5d}, {:e}", i, j_display, val); }
        else { std::println("*P{:5d}, {:5d}, {:e}", i, j_display, val); }
        i += 1;
      }
    }
  }
}
