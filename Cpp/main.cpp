#include <print>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

constexpr int NX = 21;
constexpr int NY = 7;

constexpr int MAXROW = 27 * NY;
constexpr int NELEMN = 2 * (NX - 1) * (NY - 1);
constexpr int MX = 2 * NX - 1;
constexpr int MY = 2 * NY - 1;
constexpr int NP = MX * MY;
constexpr int MAXEQN = 2 * MX * MY + NX * NY;
constexpr int NNODES = 6;
constexpr int NQUAD = 3;

std::string fmt_e(double val, int prec) {
  std::string sign = (val < 0.0) ? "-" : " ";
  if (val == 0.0) {
    return sign + "0." + std::string(prec, '0') + "E+00";
  }
  double abs_val = std::abs(val);
  int exp = static_cast<int>(std::floor(std::log10(abs_val))) + 1;
  double mant = abs_val / std::pow(10.0, exp);
  auto s = std::format("{}{:.{}f}", sign, mant, prec);
  auto exp_str = std::format("E{:+03d}", exp);
  return s + exp_str;
}

void timestamp() {
  auto now = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(now);
  char buf[64];
  ctime_s(buf, sizeof(buf), &t);
  buf[std::strlen(buf) - 1] = '\0';
  std::println("Current time: {}", buf);
}

double bsp(double xq, double yq, int it, int iq, int id,
    const std::vector<std::vector<int>>& node,
    const std::vector<double>& xc,
    const std::vector<double>& yc) {
  int iq1 = iq;
  int iq2 = (iq + 1) % 3;
  int iq3 = (iq + 2) % 3;
  int i1 = node[it][iq1];
  int i2 = node[it][iq2];
  int i3 = node[it][iq3];
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
}

void qbf(double x, double y, int it, int in,
    double& bb, double& bx, double& by,
    const std::vector<std::vector<int>>& node,
    const std::vector<double>& xc,
    const std::vector<double>& yc) {
  if (in < 3) {
    int in1 = in;
    int in2 = (in + 1) % 3;
    int in3 = (in + 2) % 3;
    int i1 = node[it][in1];
    int i2 = node[it][in2];
    int i3 = node[it][in3];
    double d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    double t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    bb = t * (2.0 * t - 1.0);
    bx = (yc[i2] - yc[i3]) * (4.0 * t - 1.0) / d;
    by = (xc[i3] - xc[i2]) * (4.0 * t - 1.0) / d;
  } else {
    int inn = in - 3;
    int in1 = inn;
    int in2 = (inn + 1) % 3;
    int in3 = (inn + 2) % 3;
    int i1 = node[it][in1];
    int i2 = node[it][in2];
    int i3 = node[it][in3];
    int j1 = i2;
    int j2 = i3;
    int j3 = i1;
    double d = (xc[i2] - xc[i1]) * (yc[i3] - yc[i1]) - (xc[i3] - xc[i1]) * (yc[i2] - yc[i1]);
    double c = (xc[j2] - xc[j1]) * (yc[j3] - yc[j1]) - (xc[j3] - xc[j1]) * (yc[j2] - yc[j1]);
    double t = 1.0 + ((yc[i2] - yc[i3]) * (x - xc[i1]) + (xc[i3] - xc[i2]) * (y - yc[i1])) / d;
    double s = 1.0 + ((yc[j2] - yc[j3]) * (x - xc[j1]) + (xc[j3] - xc[j2]) * (y - yc[j1])) / c;
    bb = 4.0 * s * t;
    bx = 4.0 * (t * (yc[j2] - yc[j3]) / c + s * (yc[i2] - yc[i3]) / d);
    by = 4.0 * (t * (xc[j3] - xc[j2]) / c + s * (xc[i3] - xc[i2]) / d);
  }
}

double ubdry(double y, double para) {
  return 4.0 * para * y * (3.0 - y) / 9.0;
}

int idamax(int n, const std::vector<double>& dx, int incx) {
  int idamax_val = 0;
  if (n == 1) return idamax_val;

  if (incx == 1) {
    double dmax = std::abs(dx[0]);
    for (int i = 1; i < n; ++i) {
      if (dmax < std::abs(dx[i])) {
        idamax_val = i;
        dmax = std::abs(dx[i]);
      }
    }
  } else {
    int ix = 0;
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

int igetl(int k, const std::vector<int>& iline, int my) {
  for (int j = 0; j < my; ++j) {
    if (iline[j] == k) return j;
  }
  std::println(" ");
  std::println("IGETL - fatal error!");
  std::println("  Unable to get local unknown number for ");
  std::println("  Global variable number {}", k);
  std::exit(1);
}

void dgbfa(std::vector<std::vector<double>>& abd, int n, int ml, int mu,
    std::vector<int>& ipvt, int& info) {
  int m = ml + mu + 1;
  info = 0;

  int j_end = std::min(n, m);
  if (mu + 2 < j_end) {
    for (int col = mu + 1; col < j_end - 1; ++col) {
      int i0 = m - col;
      for (int row = i0 - 1; row < ml; ++row) {
        abd[row][col] = 0.0;
      }
    }
  }

  int jz = std::max(j_end - 2, 0);
  int ju = 0;

  for (int k = 0; k < n - 1; ++k) {
    jz += 1;
    if (jz < n) {
      for (int row = 0; row < ml; ++row) {
        abd[row][jz] = 0.0;
      }
    }

    int lm = std::min(ml, n - 1 - k);
    int l = m;
    double dmax = std::abs(abd[m - 1][k]);
    for (int i = 1; i <= lm; ++i) {
      double abs_val = std::abs(abd[m - 1 + i][k]);
      if (dmax < abs_val) {
        dmax = abs_val;
        l = m + i;
      }
    }
    ipvt[k] = l + k + 1 - m;

    if (abd[l - 1][k] == 0.0) {
      info = k + 1;
    } else {
      if (l != m) {
        double t = abd[l - 1][k];
        abd[l - 1][k] = abd[m - 1][k];
        abd[m - 1][k] = t;
      }

      double t = -1.0 / abd[m - 1][k];
      for (int i = 0; i < lm; ++i) {
        abd[m + i][k] *= t;
      }

      ju = std::max(ju, mu + static_cast<int>(ipvt[k]));
      if (ju > n) ju = n;
      int mm = m;
      int ll = l;
      for (int j = k + 1; j < ju; ++j) {
        ll -= 1;
        mm -= 1;
        double t2 = abd[ll - 1][j];
        if (ll != mm) {
          abd[ll - 1][j] = abd[mm - 1][j];
          abd[mm - 1][j] = t2;
        }
        for (int i = 0; i < lm; ++i) {
          abd[mm + i][j] += t2 * abd[m + i][k];
        }
      }
    }
  }

  ipvt[n - 1] = n;
  if (abd[m - 1][n - 1] == 0.0) {
    info = n;
  }
}

void dgbsl(const std::vector<std::vector<double>>& abd, int n, int ml, int mu,
    const std::vector<int>& ipvt, std::vector<double>& b, int job) {
  int m = mu + ml + 1;

  if (job == 0) {
    if (0 < ml) {
      for (int k = 0; k < n - 1; ++k) {
        int lm = std::min(ml, n - 1 - k);
        int l = ipvt[k] - 1;
        double t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
        for (int i = 0; i < lm; ++i) {
          b[k + 1 + i] += t * abd[m + i][k];
        }
      }
    }

    for (int k = n - 1; k >= 0; --k) {
      b[k] /= abd[m - 1][k];
      int lm = std::min(k + 1, m) - 1;
      if (lm > 0) {
        int la = m - lm - 1;
        int lb = k - lm;
        double t = -b[k];
        for (int i = 0; i < lm; ++i) {
          b[lb + i] += t * abd[la + i][k];
        }
      }
    }
  } else {
    for (int k = 0; k < n; ++k) {
      int lm = std::min(k + 1, m) - 1;
      if (lm > 0) {
        int la = m - lm - 1;
        int lb = k - lm;
        double t = 0.0;
        for (int i = 0; i < lm; ++i) {
          t += abd[la + i][k] * b[lb + i];
        }
        b[k] = (b[k] - t) / abd[m - 1][k];
      } else {
        b[k] /= abd[m - 1][k];
      }
    }

    if (0 < ml) {
      for (int k = n - 2; k >= 0; --k) {
        int lm = std::min(ml, n - k - 1);
        if (lm > 0) {
          double t = 0.0;
          for (int i = 0; i < lm; ++i) {
            t += abd[m + i][k] * b[k + 1 + i];
          }
          b[k] += t;
        }
        int l = ipvt[k] - 1;
        if (l != k) {
          double t = b[l];
          b[l] = b[k];
          b[k] = t;
        }
      }
    }
  }
}

void uval(const std::vector<double>& g,
    const std::vector<std::vector<int>>& indx,
    int iquad, int it,
    const std::vector<std::vector<int>>& node,
    double para,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
    std::vector<double>& un, std::vector<double>& unx, std::vector<double>& uny,
    const std::vector<double>& yc) {
  un[0] = 0.0; un[1] = 0.0;
  unx[0] = 0.0; unx[1] = 0.0;
  uny[0] = 0.0; uny[1] = 0.0;

  for (int iq = 0; iq < NNODES; ++iq) {
    int ip = node[it][iq];
    double bb = phi[it][iquad][iq][0];
    double bx = phi[it][iquad][iq][1];
    double by = phi[it][iquad][iq][2];

    for (int iuk = 0; iuk < 2; ++iuk) {
      int iun = indx[ip][iuk];
      if (iun >= 0) {
        un[iuk] += bb * g[iun];
        unx[iuk] += bx * g[iun];
        uny[iuk] += by * g[iun];
      } else if (iun == -2) {
        double ubc = ubdry(yc[ip], para);
        un[iuk] += bb * ubc;
        unx[iuk] += bx * ubc;
        uny[iuk] += by * ubc;
      }
    }
  }
}

void getg(const std::vector<double>& f, const std::vector<int>& iline,
    int my, int /*neqn*/, std::vector<double>& u) {
  for (int j = 0; j < my; ++j) {
    int k = iline[j];
    if (k >= 0) {
      u[j] = f[k];
    } else {
      u[j] = 0.0;
    }
  }
}

void gram(std::vector<std::vector<double>>& gr,
    const std::vector<int>& iline,
    const std::vector<std::vector<int>>& indx,
    int iwrite,
    const std::vector<std::vector<int>>& node,
    int nodex0, double para,
    std::vector<double>& r,
    const std::vector<double>& ui,
    const std::vector<double>& xc,
    const std::vector<double>& yc) {
  std::vector<double> wt = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  std::vector<double> yq = {-0.7745966692, 0.0, 0.7745966692};

  for (int i = 0; i < MY; ++i) {
    r[i] = 0.0;
    for (int j = 0; j < MY; ++j) {
      gr[i][j] = 0.0;
    }
  }

  double xzero = xc[nodex0];

  for (int it = 0; it < NELEMN; ++it) {
    int k = node[it][0];
    int kk = node[it][1];

    if (1.0E-04 < std::abs(xc[k] - xzero)) continue;
    if (1.0E-04 < std::abs(xc[kk] - xzero)) continue;

    for (int iquad = 0; iquad < 3; ++iquad) {
      double bma2 = (yc[kk] - yc[k]) / 2.0;
      double ar = bma2 * wt[iquad];
      double x = xzero;
      double y = yc[k] + bma2 * (yq[iquad] + 1.0);

      double uiqdpt = 0.0;
      for (int iq = 0; iq < 6; ++iq) {
        if (iq > 3) continue;
        if (iq == 2) continue;
        double bb = 0.0, bx = 0.0, by = 0.0;
        qbf(x, y, it, iq, bb, bx, by, node, xc, yc);
        int ip = node[it][iq];
        int iun = indx[ip][0];
        if (iun >= 0) {
          int ii = igetl(iun, iline, MY);
          uiqdpt += bb * ui[ii];
        } else if (iun == -2) {
          double ubc = ubdry(yc[ip], para);
          uiqdpt += bb * ubc;
        }
      }

      for (int iq = 0; iq < 6; ++iq) {
        if (iq == 0 || iq == 1 || iq == 3) {
          int ip = node[it][iq];
          double bb = 0.0, bx = 0.0, by = 0.0;
          qbf(x, y, it, iq, bb, bx, by, node, xc, yc);
          int i = indx[ip][0];
          if (i < 0) continue;
          int ii = igetl(i, iline, MY);
          r[ii] += bb * uiqdpt * ar;

          for (int iqq = 0; iqq < 6; ++iqq) {
            if (iqq == 0 || iqq == 1 || iqq == 3) {
              int ipp = node[it][iqq];
              double bbb = 0.0, bbx = 0.0, bby = 0.0;
              qbf(x, y, it, iqq, bbb, bbx, bby, node, xc, yc);
              int j = indx[ipp][0];
              if (j >= 0) {
                int jj = igetl(j, iline, MY);
                gr[ii][jj] += bb * bbb * ar;
              }
            }
          }
        }
      }
    }
  }

  if (2 <= iwrite) {
    std::println(" ");
    std::println("Gram matrix:");
    std::println(" ");
    for (int i = 0; i < MY; ++i) {
      for (int j = 0; j < MY; ++j) {
        std::println("{} {} {}", i + 1, j + 1, gr[i][j]);
      }
    }
    std::println(" ");
    std::println("R vector:");
    std::println(" ");
    for (int i = 0; i < MY; ++i) {
      std::println("{} {}", i + 1, r[i]);
    }
  }
}

void setgrd(std::vector<std::vector<int>>& indx,
    std::vector<int>& insc,
    std::vector<int>& isotri,
    int iwrite, bool& _long, int& neqn,
    std::vector<std::vector<int>>& node) {
  if (NY < NX) {
    _long = true;
    std::println("Using vertical ordering.");
  } else {
    _long = false;
    std::println("Using horizontal ordering.");
  }
  std::println(" ");
  std::println("SETGRD: Taylor Hood element");

  neqn = 0;
  int ielemn = 0;

  for (int ip = 0; ip < NP; ++ip) {
    int ic, jc;
    if (_long) {
      ic = (ip / MY) + 1;
      jc = (ip % MY) + 1;
    } else {
      ic = (ip % MX) + 1;
      jc = (ip / MX) + 1;
    }

    int icnt = ic % 2;
    int jcnt = jc % 2;

    if (icnt == 1 && jcnt == 1 && ic != MX && jc != MY) {
      if (_long) {
        int ip1 = ip + MY;
        int ip2 = ip + MY + MY;

        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2 + 2;
        node[ielemn - 1][2] = ip2;
        node[ielemn - 1][3] = ip1 + 1;
        node[ielemn - 1][4] = ip2 + 1;
        node[ielemn - 1][5] = ip1;
        isotri[ielemn - 1] = 0;

        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip + 2;
        node[ielemn - 1][2] = ip2 + 2;
        node[ielemn - 1][3] = ip + 1;
        node[ielemn - 1][4] = ip1 + 2;
        node[ielemn - 1][5] = ip1 + 1;
        isotri[ielemn - 1] = 0;
      } else {
        int ip1 = ip + MX;
        int ip2 = ip + MX + MX;

        ielemn += 1;
        node[ielemn - 1][0] = ip;
        node[ielemn - 1][1] = ip2;
        node[ielemn - 1][2] = ip2 + 2;
        node[ielemn - 1][3] = ip1;
        node[ielemn - 1][4] = ip2 + 1;
        node[ielemn - 1][5] = ip1 + 1;
        isotri[ielemn - 1] = 0;

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

    if (ic == 1 && jc > 1 && jc < MY) {
      indx[ip][0] = -2;
      indx[ip][1] = -1;
    } else if (ic == MX && jc > 1 && jc < MY) {
      indx[ip][0] = neqn;
      neqn += 1;
      indx[ip][1] = -1;
    } else if (jc == 1 || jc == MY) {
      indx[ip][0] = -1;
      indx[ip][1] = -1;
    } else {
      indx[ip][0] = neqn;
      indx[ip][1] = neqn + 1;
      neqn += 2;
    }

    if (jcnt == 1 && icnt == 1) {
      insc[ip] = neqn;
      neqn += 1;
    } else {
      insc[ip] = -1;
    }
  }

  if (2 <= iwrite) {
    std::println(" ");
    std::println("    I      INDX 1 & 2, INSC");
    std::println(" ");
    for (int i = 0; i < NP; ++i) {
      auto dsp = [](int v) -> int {
        if (v >= 0) return v + 1;
        else if (v == -2) return -1;
        else return 0;
      };
      std::println("{:6} {:6} {:6} {:6}", i + 1, dsp(indx[i][0]), dsp(indx[i][1]), dsp(insc[i]));
    }
    std::println(" ");
    std::println("    IT    NODE(IT,1:6)");
    std::println(" ");
    for (int it = 0; it < NELEMN; ++it) {
      std::println("{:6} {:6} {:6} {:6} {:6} {:6} {:6}",
        it + 1,
        node[it][0] + 1, node[it][1] + 1, node[it][2] + 1,
        node[it][3] + 1, node[it][4] + 1, node[it][5] + 1);
    }
  }

  std::println("Number of unknowns = {}", neqn);
  if (MAXEQN < neqn) {
    std::println("SETGRD - Too many unknowns!");
    std::println("The maximum allowed is MAXEQN = {}", MAXEQN);
    std::println("This problem requires NEQN = {}", neqn);
    std::exit(1);
  }
}

void setban(const std::vector<std::vector<int>>& indx,
    const std::vector<int>& insc,
    int& nband, int& nlband,
    const std::vector<std::vector<int>>& node,
    int& nrow) {
  nlband = 0;

  for (int it = 0; it < NELEMN; ++it) {
    for (int iq = 0; iq < NNODES; ++iq) {
      int ip = node[it][iq];
      for (int iuk = 0; iuk < 3; ++iuk) {
        int i = (iuk == 2) ? insc[ip] : indx[ip][iuk];
        if (i >= 0) {
          for (int iqq = 0; iqq < NNODES; ++iqq) {
            int ipp = node[it][iqq];
            for (int iukk = 0; iukk < 3; ++iukk) {
              int j = (iukk == 2) ? insc[ipp] : indx[ipp][iukk];
              if (j > i) {
                int diff = j - i;
                nlband = std::max(nlband, diff);
              }
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
  if (MAXROW < nrow) {
    std::println("SETBAN - NROW is too large!");
    std::println("The maximum allowed is {}", MAXROW);
    std::println("Aborting.");
    std::exit(1);
  }
}

void setlin(std::vector<int>& iline,
    const std::vector<std::vector<int>>& indx,
    int iwrite, bool _long,
    int& nodex0, double xlngth) {
  int itemp = static_cast<int>(std::lround((18.0 * (NX - 1)) / xlngth));
  nodex0 = _long ? itemp * (2 * NY - 1) : itemp;

  for (int i = 0; i < MY; ++i) {
    int ip = _long ? nodex0 + i : nodex0 + MX * i;
    iline[i] = indx[ip][0];
  }

  if (1 <= iwrite) {
    std::println(" ");
    std::println("SETLIN: unknown numbers along line:");
    std::println(" ");
    for (int i = 0; i < MY; ++i) {
      int val = (iline[i] >= 0) ? iline[i] + 1 : ((iline[i] == -2) ? -1 : 0);
      std::print("{:>5}", val);
      if ((i + 1) % 15 == 0 && (i + 1) < MY) {
        std::println("");
      }
    }
    std::println(" ");
  }
}

void setxy(int iwrite, bool _long,
    std::vector<double>& xc, double xlngth,
    std::vector<double>& yc, double ylngth) {
  for (int ip = 0; ip < NP; ++ip) {
    int ic, jc;
    if (_long) {
      ic = ip / MY;
      jc = ip % MY;
    } else {
      ic = ip % MX;
      jc = ip / MX;
    }
    xc[ip] = static_cast<double>(ic) * xlngth / (2.0 * NX - 2.0);
    yc[ip] = static_cast<double>(jc) * ylngth / (2.0 * NY - 2.0);
  }

  if (2 <= iwrite) {
    std::println(" ");
    std::println("    I      XC           YC");
    std::println(" ");
    for (int i = 0; i < NP; ++i) {
      std::println("{:>5} {:>12.5f} {:>12.5f}", i + 1, xc[i], yc[i]);
    }
  }
}

void setqud(std::vector<double>& area,
    const std::vector<std::vector<int>>& node,
    const std::vector<double>& xc,
    std::vector<std::vector<double>>& xm,
    const std::vector<double>& yc,
    std::vector<std::vector<double>>& ym) {
  for (int it = 0; it < NELEMN; ++it) {
    int ip1 = node[it][0];
    int ip2 = node[it][1];
    int ip3 = node[it][2];
    double x1 = xc[ip1], x2 = xc[ip2], x3 = xc[ip3];
    double y1 = yc[ip1], y2 = yc[ip2], y3 = yc[ip3];
    xm[it][0] = 0.5 * (x1 + x2);
    xm[it][1] = 0.5 * (x2 + x3);
    xm[it][2] = 0.5 * (x3 + x1);
    ym[it][0] = 0.5 * (y1 + y2);
    ym[it][1] = 0.5 * (y2 + y3);
    ym[it][2] = 0.5 * (y3 + y1);
    area[it] = 0.5 * std::abs(
      (y1 + y2) * (x2 - x1) + (y2 + y3) * (x3 - x2) + (y3 + y1) * (x1 - x3));
  }
}

void setbas(const std::vector<std::vector<int>>& node,
    const std::vector<double>& xc, const std::vector<double>& yc,
    std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
    std::vector<std::vector<std::vector<double>>>& psi,
    const std::vector<std::vector<double>>& xm,
    const std::vector<std::vector<double>>& ym) {
  for (int it = 0; it < NELEMN; ++it) {
    for (int j = 0; j < NQUAD; ++j) {
      double x = xm[it][j];
      double y = ym[it][j];
      for (int iq = 0; iq < 6; ++iq) {
        psi[it][j][iq] = bsp(x, y, it, iq, 1, node, xc, yc);
        double bb = 0.0, bx = 0.0, by = 0.0;
        qbf(x, y, it, iq, bb, bx, by, node, xc, yc);
        phi[it][j][iq][0] = bb;
        phi[it][j][iq][1] = bx;
        phi[it][j][iq][2] = by;
      }
    }
  }
}

void linsys(std::vector<std::vector<double>>& a,
    std::vector<double>& area,
    std::vector<double>& f,
    std::vector<double>& g,
    std::vector<std::vector<int>>& indx,
    std::vector<int>& insc,
    std::vector<int>& ipivot,
    int neqn, int nlband,
    std::vector<std::vector<int>>& node,
    int nrow,
    double para1, double para2,
    std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
    std::vector<std::vector<std::vector<double>>>& psi,
    double reynld,
    std::vector<double>& yc) {
  int ioff = nlband + nlband;
  double visc = 1.0 / reynld;

  std::fill(f.begin(), f.begin() + neqn, 0.0);
  for (int k = 0; k < nrow; ++k) {
    std::fill(a[k].begin(), a[k].begin() + neqn, 0.0);
  }

  for (int it = 0; it < NELEMN; ++it) {
    double ar = area[it] / 3.0;
    for (int iquad = 0; iquad < NQUAD; ++iquad) {
      std::vector<double> un(2, 0.0), unx(2, 0.0), uny(2, 0.0);
      uval(g, indx, iquad, it, node, para1, phi, un, unx, uny, yc);

      for (int iq = 0; iq < NNODES; ++iq) {
        int ip = node[it][iq];
        double bb = phi[it][iquad][iq][0];
        double bx = phi[it][iquad][iq][1];
        double by = phi[it][iquad][iq][2];
        double bbl = psi[it][iquad][iq];
        int ihor = indx[ip][0];
        int iver = indx[ip][1];
        int iprs = insc[ip];

        if (ihor >= 0) {
          f[ihor] += ar * bb * (un[0] * unx[0] + un[1] * uny[0]);
        }
        if (iver >= 0) {
          f[iver] += ar * bb * (un[0] * unx[1] + un[1] * uny[1]);
        }

        for (int iqq = 0; iqq < NNODES; ++iqq) {
          int ipp = node[it][iqq];
          double bbb = phi[it][iquad][iqq][0];
          double bbx = phi[it][iquad][iqq][1];
          double bby = phi[it][iquad][iqq][2];
          double bbbl = psi[it][iquad][iqq];
          int ju = indx[ipp][0];
          int jv = indx[ipp][1];
          int jp = insc[ipp];

          if (ju >= 0) {
            if (ihor >= 0) {
              int iuse = ihor - ju + ioff;
              a[iuse][ju] += ar * (visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }
            if (iver >= 0) {
              int iuse = iver - ju + ioff;
              a[iuse][ju] += ar * bb * bbb * unx[1];
            }
            if (iprs >= 0) {
              int iuse = iprs - ju + ioff;
              a[iuse][ju] += ar * bbx * bbl;
            }
          } else if (ju == -2) {
            double uu = ubdry(yc[ipp], para2);
            if (ihor >= 0) {
              f[ihor] -= ar * uu * (visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]));
            }
            if (iver >= 0) {
              f[iver] -= ar * uu * bb * bbb * unx[1];
            }
            if (iprs >= 0) {
              f[iprs] -= ar * uu * bbx * bbl;
            }
          }

          if (jv >= 0) {
            if (ihor >= 0) {
              int iuse = ihor - jv + ioff;
              a[iuse][jv] += ar * bb * bbb * uny[0];
            }
            if (iver >= 0) {
              int iuse = iver - jv + ioff;
              a[iuse][jv] += ar * (visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]));
            }
            if (iprs >= 0) {
              int iuse = iprs - jv + ioff;
              a[iuse][jv] += ar * bby * bbl;
            }
          }

          if (jp >= 0) {
            if (ihor >= 0) {
              int iuse = ihor - jp + ioff;
              a[iuse][jp] -= ar * bx * bbbl;
            }
            if (iver >= 0) {
              int iuse = iver - jp + ioff;
              a[iuse][jp] -= ar * by * bbbl;
            }
          }
        }
      }
    }
  }

  f[neqn - 1] = 0.0;
  for (int j = neqn - nlband - 1; j < neqn - 1; ++j) {
    int i = (neqn - 1) - j + ioff;
    a[i][j] = 0.0;
  }
  a[ioff][neqn - 1] = 1.0;

  int info = 0;
  dgbfa(a, neqn, nlband, nlband, ipivot, info);

  if (info != 0) {
    std::println(" ");
    std::println("LINSYS - fatal error!");
    std::println("DGBFA returns INFO = {}", info);
    std::exit(1);
  }

  int job = 0;
  dgbsl(a, neqn, nlband, nlband, ipivot, f, job);
}

void nstoke(std::vector<std::vector<double>>& a,
    std::vector<double>& area,
    std::vector<double>& f,
    std::vector<double>& g,
    std::vector<std::vector<int>>& indx,
    std::vector<int>& insc,
    std::vector<int>& ipivot,
    int iwrite, int maxnew, int neqn, int nlband,
    std::vector<std::vector<int>>& node,
    int nrow,
    int& numnew, double para,
    std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
    std::vector<std::vector<std::vector<double>>>& psi,
    double reynld, double tolnew,
    std::vector<double>& yc) {
  for (int iter = 0; iter < maxnew; ++iter) {
    numnew += 1;

    linsys(a, area, f, g, indx, insc, ipivot, neqn, nlband, node, nrow,
      para, para, phi, psi, reynld, yc);

    for (int i = 0; i < neqn; ++i) {
      g[i] = g[i] - f[i];
    }
    double diff = std::abs(g[idamax(neqn, g, 1)]);

    if (1 <= iwrite) {
      std::println("NSTOKE iteration {} Mnorm = {}", iter + 1, diff);
    }

    for (int i = 0; i < neqn; ++i) {
      g[i] = f[i];
    }

    if (diff <= tolnew) {
      std::println("Navier Stokes iteration converged in {} iterations.", iter + 1);
      return;
    }
  }
  std::println("Navier Stokes solution did not converge!");
}

void resid(const std::vector<double>& area,
    const std::vector<double>& g,
    const std::vector<std::vector<int>>& indx,
    const std::vector<int>& insc,
    int iwrite, int neqn,
    const std::vector<std::vector<int>>& node,
    double para,
    const std::vector<std::vector<std::vector<std::vector<double>>>>& phi,
    const std::vector<std::vector<std::vector<double>>>& psi,
    std::vector<double>& res,
    double reynld,
    const std::vector<double>& yc) {
  double visc = 1.0 / reynld;

  std::fill(res.begin(), res.begin() + neqn, 0.0);

  for (int it = 0; it < NELEMN; ++it) {
    double ar = area[it] / 3.0;
    for (int iquad = 0; iquad < NQUAD; ++iquad) {
      std::vector<double> un(2, 0.0), unx(2, 0.0), uny(2, 0.0);
      uval(g, indx, iquad, it, node, para, phi, un, unx, uny, yc);

      for (int iq = 0; iq < NNODES; ++iq) {
        int ip = node[it][iq];
        double bb = phi[it][iquad][iq][0];
        double bx = phi[it][iquad][iq][1];
        double by = phi[it][iquad][iq][2];
        double bbl = psi[it][iquad][iq];
        int iprs = insc[ip];
        int ihor = indx[ip][0];
        int iver = indx[ip][1];

        if (ihor >= 0) {
          res[ihor] += (un[0] * unx[0] + un[1] * uny[0]) * bb * ar;
        }
        if (iver >= 0) {
          res[iver] += (un[0] * unx[1] + un[1] * uny[1]) * bb * ar;
        }

        for (int iqq = 0; iqq < NNODES; ++iqq) {
          int ipp = node[it][iqq];
          double bbb = phi[it][iquad][iqq][0];
          double bbx = phi[it][iquad][iqq][1];
          double bby = phi[it][iquad][iqq][2];
          double bbbl = psi[it][iquad][iqq];
          int ju = indx[ipp][0];
          int jv = indx[ipp][1];
          int jp = insc[ipp];

          if (ju >= 0) {
            if (ihor >= 0) {
              double aijuu = visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              res[ihor] += aijuu * ar * g[ju];
            }
            if (iver >= 0) {
              double aijvu = bb * bbb * unx[1];
              res[iver] += aijvu * ar * g[ju];
            }
            if (iprs >= 0) {
              double aijpu = bbx * bbl;
              res[iprs] += aijpu * ar * g[ju];
            }
          } else if (ju == -2) {
            double uu = ubdry(yc[ipp], para);
            if (ihor >= 0) {
              double aijuu = visc * (by * bby + bx * bbx)
                + bb * (bbb * unx[0] + bbx * un[0] + bby * un[1]);
              res[ihor] += ar * aijuu * uu;
            }
            if (iver >= 0) {
              double aijvu = bb * bbb * unx[1];
              res[iver] += ar * aijvu * uu;
            }
            if (iprs >= 0) {
              double aijpu = bbx * bbl;
              res[iprs] += ar * aijpu * uu;
            }
          }

          if (jv >= 0) {
            if (ihor >= 0) {
              double aijuv = bb * bbb * uny[0];
              res[ihor] += aijuv * ar * g[jv];
            }
            if (iver >= 0) {
              double aijvv = visc * (by * bby + bx * bbx)
                + bb * (bbb * uny[1] + bby * un[1] + bbx * un[0]);
              res[iver] += aijvv * ar * g[jv];
            }
            if (iprs >= 0) {
              double aijpv = bby * bbl;
              res[iprs] += aijpv * ar * g[jv];
            }
          }

          if (jp >= 0) {
            if (ihor >= 0) {
              double aijup = -bx * bbbl;
              res[ihor] += aijup * ar * g[jp];
            }
            if (iver >= 0) {
              double aijvp = -by * bbbl;
              res[iver] += aijvp * ar * g[jp];
            }
          }
        }
      }
    }
  }

  res[neqn - 1] = g[neqn - 1];

  double rmax = 0.0;
  int imax = 0;
  int ibad = 0;

  for (int i = 0; i < neqn; ++i) {
    double test = std::abs(res[i]);
    if (rmax < test) {
      rmax = test;
      imax = i;
    }
    if (1.0E-03 < test) {
      ibad += 1;
    }
  }

  if (1 <= iwrite) {
    std::println(" ");
    std::println("RESIDUAL INFORMATION:");
    std::println(" ");
    std::println("Worst residual is number {:5}", imax + 1);
    std::println("of magnitude {}", fmt_e(rmax, 15));
    std::println(" ");
    std::println("Number of bad residuals is {} out of {}", ibad, neqn);
    std::println(" ");
  }

  if (2 <= iwrite) {
    std::println("Raw residuals:");
    std::println(" ");
    int ii = 0;
    for (int j = 0; j < NP; ++j) {
      if (indx[j][0] >= 0) {
        std::string mark = (std::abs(res[ii]) <= 1.0E-03) ? " " : "*";
        std::println("{}{} {:5} {:5} {}", mark, "U", ii + 1, j + 1, fmt_e(res[ii], 6));
        ii += 1;
      }
      if (indx[j][1] >= 0) {
        std::string mark = (std::abs(res[ii]) <= 1.0E-03) ? " " : "*";
        std::println("{}{} {:5} {:5} {}", mark, "V", ii + 1, j + 1, fmt_e(res[ii], 6));
        ii += 1;
      }
      if (insc[j] >= 0) {
        std::string mark = (std::abs(res[ii]) <= 1.0E-03) ? " " : "*";
        std::println("{}{} {:5} {:5} {}", mark, "P", ii + 1, j + 1, fmt_e(res[ii], 6));
        ii += 1;
      }
    }
  }
}

void uv_table(const std::vector<double>& f,
    const std::vector<std::vector<int>>& indx,
    double para, const std::vector<double>& yc,
    const std::string& filename) {
  std::ofstream file(filename);
  for (int ip = 0; ip < NP; ++ip) {
    int k = indx[ip][0];
    double uval = (k == -1) ? 0.0 : ((k == -2) ? ubdry(yc[ip], para) : f[k]);
    k = indx[ip][1];
    double vval = (k == -1) ? 0.0 : f[k];
    file << std::format("  {:>14.6f}  {:>14.6f}\n", uval, vval);
  }
}

void xy_table(const std::vector<double>& xc,
    const std::vector<double>& yc,
    const std::string& filename) {
  std::ofstream file(filename);
  for (int ip = 0; ip < NP; ++ip) {
    file << std::format("  {:>14.6f}  {:>14.6f}\n", xc[ip], yc[ip]);
  }
}

void xy_plot3d(bool _long,
    const std::vector<double>& xc,
    const std::vector<double>& yc,
    const std::string& filename) {
  std::ofstream file(filename);
  if (_long) {
    file << std::format("{} {}\n", MX, MY);
    for (int i = 0; i < MY; ++i) {
      for (int j = 0; j < MX; ++j) {
        int ip = j * MY + i;
        file << std::format("{:>14.6f}\n", xc[ip]);
      }
    }
    for (int i = 0; i < MY; ++i) {
      for (int j = 0; j < MX; ++j) {
        int ip = j * MY + i;
        file << std::format("{:>14.6f}\n", yc[ip]);
      }
    }
  } else {
    file << std::format("{} {}\n", MY, MX);
    for (int j = 0; j < MX; ++j) {
      for (int i = 0; i < MY; ++i) {
        int ip = j * MY + i;
        file << std::format("{:>14.6f}\n", xc[ip]);
      }
    }
    for (int j = 0; j < MX; ++j) {
      for (int i = 0; i < MY; ++i) {
        int ip = j * MY + i;
        file << std::format("{:>14.6f}\n", yc[ip]);
      }
    }
  }
}

int main() {
  auto inicio = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<double>> a(MAXROW, std::vector<double>(MAXEQN, 0.0));
  double a2;
  double abound;
  double anew;
  double aold;
  std::vector<double> area(NELEMN, 0.0);
  std::vector<double> dcda(MY, 0.0);
  std::vector<double> f(MAXEQN, 0.0);

  std::vector<double> g(MAXEQN, 0.0);
  std::vector<std::vector<double>> gr(MY, std::vector<double>(MY, 0.0));

  std::vector<int> iline(MY, -1);
  std::vector<std::vector<int>> indx(NP, std::vector<int>(2, -1));
  std::vector<int> insc(NP, -1);

  std::vector<int> ipivot(MAXEQN, 0);
  std::vector<int> isotri(NELEMN, 0);

  bool _long = false;

  int nband = 0;
  int neqn = 0;
  int nlband = 0;
  std::vector<std::vector<int>> node(NELEMN, std::vector<int>(NNODES, 0));

  int nodex0 = 0;

  int nrow = 0;

  double para;
  std::vector<std::vector<std::vector<std::vector<double>>>> phi(
    NELEMN, std::vector<std::vector<std::vector<double>>>(
      NQUAD, std::vector<std::vector<double>>(
        NNODES, std::vector<double>(3, 0.0))));
  std::vector<std::vector<std::vector<double>>> psi(
    NELEMN, std::vector<std::vector<double>>(
      NQUAD, std::vector<double>(NNODES, 0.0)));

  std::vector<double> r(MY, 0.0);
  std::vector<double> res(MAXEQN, 0.0);

  double rjpold = 0.0;
  double test = 0.0;

  std::vector<double> ui(MY, 0.0);
  std::vector<double> unew(MY, 0.0);
  std::vector<double> xc(NP, 0.0);

  std::vector<std::vector<double>> xm(NELEMN, std::vector<double>(NQUAD, 0.0));
  std::vector<double> yc(NP, 0.0);

  std::vector<std::vector<double>> ym(NELEMN, std::vector<double>(NQUAD, 0.0));

  timestamp();
  std::println(" ");
  std::println("channel():");
  std::println("  C++23 version");
  std::println("  Channel flow control problem");
  std::println(" ");
  std::println("Last modified:");
  std::println("  2026");
  std::println("");
  std::println("  Flow control problem:");
  std::println("    Inflow controlled by one parameter.");
  std::println("    Velocities measured along vertical line.");
  std::println("    Try to match specified velocity profile.");

  std::string fileg = "display.txt";
  std::string fileu = "uv.dat";
  std::string filex = "xy.dat";

  int iwrite = 10;
  int maxnew = 10;
  int maxsec = 8;
  int numnew = 0;
  int numsec = 0;
  double reynld = 1.0;
  double rjpnew = 0.0;
  double tolnew = 1.0E-04;
  double tolsec = 1.0E-06;
  double xlngth = 10.0;
  double ylngth = 3.0;

  std::println("");
  std::println("NX = {}", NX);
  std::println("NY = {}", NY);
  std::println("Number of elements = {}", NELEMN);
  std::println("Reynolds number =  {}.", reynld);
  std::println("Secant tolerance = {}", tolsec);
  std::println("Newton tolerance = {}", tolnew);
  std::println("");

  setgrd(indx, insc, isotri, iwrite, _long, neqn, node);

  setban(indx, insc, nband, nlband, node, nrow);

  setlin(iline, indx, iwrite, _long, nodex0, xlngth);

  setxy(iwrite, _long, xc, xlngth, yc, ylngth);

  setqud(area, node, xc, xm, yc, ym);

  setbas(node, xc, yc, phi, psi, xm, ym);

  para = 1.0;
  std::println(" ");
  std::println("Solve Navier Stokes problem with parameter = {}", para);
  std::println("for profile at x = {}", xc[nodex0]);
  for (int i = 0; i < neqn; ++i) {
    g[i] = 1.0;
  }

  nstoke(a, area, f, g, indx, insc, ipivot, iwrite, maxnew, neqn, nlband,
    node, nrow, numnew, para, phi, psi, reynld, tolnew, yc);

  if (1 <= iwrite) {
    resid(area, f, indx, insc, iwrite, neqn, node, para, phi, psi, res, reynld, yc);
  }

  getg(f, iline, MY, neqn, ui);

  if (1 <= iwrite) {
    std::println(" ");
    std::println("U profile:");
    std::println(" ");
    for (int i = 0; i < MY; i += 5) {
      for (int j = i; j < std::min(i + 5, MY); ++j) {
        std::print("{:>14.6f}", ui[j]);
      }
      std::println("");
    }
  }

  gram(gr, iline, indx, iwrite, node, nodex0, para, r, ui, xc, yc);

  xy_table(xc, yc, filex);
  uv_table(f, indx, para, yc, fileu);
  xy_plot3d(_long, xc, yc, fileg);

  for (int i = 0; i < neqn; ++i) {
    f[i] = 0.0;
    g[i] = 0.0;
  }

  aold = 0.0;
  rjpold = 0.0;
  anew = 0.1;

  for (int iter = 1; iter <= maxsec; ++iter) {
    numsec += 1;
    std::println(" ");
    std::println("Secant iteration {}", iter);
    std::println(" ");
    std::println("Solving Navier Stokes problem for parameter = {}", anew);

    for (int i = 0; i < neqn; ++i) {
      g[i] = f[i];
    }
    para = anew;

    nstoke(a, area, f, g, indx, insc, ipivot, iwrite, maxnew, neqn, nlband,
      node, nrow, numnew, para, phi, psi, reynld, tolnew, yc);

    getg(f, iline, MY, neqn, unew);

    if (1 <= iwrite) {
      std::println(" ");
      std::println("Velocity profile:");
      std::println(" ");
      for (int i = 0; i < MY; i += 5) {
        for (int j = i; j < std::min(i + 5, MY); ++j) {
          std::print("{:>14.6f}", unew[j]);
        }
        std::println("");
      }
    }

    para = anew;
    abound = 1.0;
    linsys(a, area, g, f, indx, insc, ipivot, neqn, nlband, node, nrow,
      para, abound, phi, psi, reynld, yc);

    getg(g, iline, MY, neqn, dcda);

    if (2 <= iwrite) {
      std::println(" ");
      std::println("Sensitivities:");
      std::println(" ");
      for (int i = 0; i < MY; i += 5) {
        for (int j = i; j < std::min(i + 5, MY); ++j) {
          std::print("{:>14.6f}", dcda[j]);
        }
        std::println("");
      }
    }

    rjpnew = 0.0;
    for (int i = 0; i < MY; ++i) {
      double temp = -r[i];
      for (int j = 0; j < MY; ++j) {
        temp += gr[i][j] * unew[j];
      }
      rjpnew += 2.0 * dcda[i] * temp;
    }

    std::println(" ");
    std::println("Parameter  = {} J prime = {}", anew, rjpnew);

    if (iter == 1) {
      a2 = 0.5;
    } else {
      a2 = aold - rjpold * (anew - aold) / (rjpnew - rjpold);
    }

    aold = anew;
    anew = a2;
    rjpold = rjpnew;
    test = std::abs(anew - aold) / std::abs(anew);

    std::println("New value of parameter = {}", anew);
    std::println("Convergence test = {}", test);

    if (std::abs(anew - aold) < std::abs(anew) * tolsec) {
      std::println("Secant iteration converged.");
      break;
    }
  }

  if (maxsec < 1 || std::abs(anew - aold) >= std::abs(anew) * tolsec) {
    std::println("Secant iteration failed to converge.");
  }

  std::println("Number of secant steps = {}", numsec);
  std::println("Number of Newton steps = {}", numnew);

  std::println(" ");
  std::println("CHANNEL:");
  std::println("  Normal end of execution.");
  std::println(" ");
  timestamp();
  auto duracion = std::chrono::high_resolution_clock::now() - inicio;
  std::ofstream file("times.txt", std::ios::app);
  file << std::chrono::duration<double>(duracion).count() << "\n";
  file.flush();
  std::println(
    "Elapsed time: {:f} s",
    std::chrono::duration<double>(duracion).count()
  );

}
