import re
import time
import subprocess as sp
import os
import shutil
from tabulate import tabulate

REPETITIONS = 11
MEASURED = 10

CONFIGS = [
  (21, 7, "XS"),
  (31, 11, "S"),
  (41, 15, "M"),
]

OUTPUT_FILES = ["display.txt", "uv.txt", "uv.dat", "xy.txt", "xy.dat", "times.txt"]

LANGUAGES = {
  "Fortran": {
    "dir": "Fortran90",
    "src": "main.f90",
    "build": ["flang", "main.f90", "-O3", "-o", "build/main.exe"],
    "run": ["build/main.exe"],
    "re_nx": re.compile(r"(integer,\s*parameter\s*::\s*nx\s*=)\s*\d+", re.IGNORECASE),
    "re_ny": re.compile(r"(integer,\s*parameter\s*::\s*ny\s*=)\s*\d+", re.IGNORECASE),
  },
  "C++": {
    "dir": "Cpp",
    "src": "main.cpp",
    "build": ["clang++", "main.cpp", "-std=c++23", "-O3", "-o", "build/main.exe"],
    "run": ["build/main.exe"],
    "re_nx": re.compile(r"(constexpr\s+int\s+NX\s*=)\s*\d+"),
    "re_ny": re.compile(r"(constexpr\s+int\s+NY\s*=)\s*\d+"),
  },
  "Rust": {
    "dir": "Rust",
    "src": "src/main.rs",
    "build": ["cargo", "build", "--release"],
    "run": ["./target/release/Rust.exe"],
    "re_nx": re.compile(r"(const\s+NX\s*:\s*usize\s*=)\s*\d+"),
    "re_ny": re.compile(r"(const\s+NY\s*:\s*usize\s*=)\s*\d+"),
  },
}


def set_nx_ny(filepath: str, nx: int, ny: int, re_nx, re_ny):
  with open(filepath, "r", encoding="utf-8") as f:
    content = f.read()
  content = re_nx.sub(rf"\g<1> {nx}", content)
  content = re_ny.sub(rf"\g<1> {ny}", content)
  with open(filepath, "w", encoding="utf-8") as f:
    f.write(content)


def _extract_int(text: str) -> int:
  return int(re.search(r"\d+", text).group(0)) # type: ignore


def get_nx_ny(filepath: str, re_nx, re_ny) -> tuple[int, int]:
  with open(filepath, "r", encoding="utf-8") as f:
    content = f.read()
  return _extract_int(re_nx.search(content).group(0)), _extract_int(
    re_ny.search(content).group(0)
  )


def compile_language(lang_info: dict, dirpath: str) -> bool:
  result = sp.run(
    lang_info["build"],
    cwd=dirpath,
    stdout=sp.PIPE,
    stderr=sp.PIPE,
    timeout=300,
  )
  if result.returncode != 0:
    err = result.stderr.decode("utf-8", errors="replace")
    print(f"    Compilacion FALLO:\n{err[:500]}")
    return False
  return True


def cleanup_outputs(dirpath: str):
  for fname in OUTPUT_FILES:
    fpath = os.path.join(dirpath, fname)
    try:
      if os.path.exists(fpath):
        os.remove(fpath)
    except OSError:
      pass


def benchmark(
  lang_info: dict, dirpath: str, timeout_s: int = 600
) -> list[float] | None:
  original_dir = os.getcwd()
  os.chdir(dirpath)

  cleanup_outputs(".")

  tiempos: list[float] = []
  for i in range(REPETITIONS):
    try:
      start = time.perf_counter()
      sp.run(
        lang_info["run"],
        stdout=sp.DEVNULL,
        stderr=sp.DEVNULL,
        timeout=timeout_s,
      )
      end = time.perf_counter()
      elapsed = end - start
      tiempos.append(elapsed)
      if i == 0:
        print(f"      Warm-up: {elapsed:.3f}s")
      cleanup_outputs(".")
    except sp.TimeoutExpired:
      print(f"      Timeout en ejecucion {i + 1}")
      os.chdir(original_dir)
      return None

  os.chdir(original_dir)
  return tiempos[1:]


def main():
  print("check_time.py - Benchmark multi-configuracion")
  print(f"{'=' * 70}")
  print(
    f"Repeticiones por configuracion: {REPETITIONS} (1 warm-up + {MEASURED} medidas)"
  )
  print(
    f"Configuraciones: {', '.join(f'{lbl} (NX={nx}, NY={ny})' for nx, ny, lbl in CONFIGS)}"
  )
  print()

  base_dir = os.getcwd()
  resultados: dict[tuple[int, int], dict[str, list[float] | None]] = {}
  tiempos_compilacion: dict[tuple[int, int], dict[str, float]] = {}

  for lang_name, lang_info in LANGUAGES.items():
    src_path = os.path.join(base_dir, lang_info["dir"], lang_info["src"])
    bak_path = src_path + ".bak"
    if not os.path.exists(bak_path):
      shutil.copy2(src_path, bak_path)

  try:
    for nx, ny, lbl in CONFIGS:
      print(f"{'=' * 70}")
      print(f"Configuracion: {lbl} - NX={nx}, NY={ny}")
      print(f"{'=' * 70}")
      resultados[(nx, ny)] = {}
      tiempos_compilacion[(nx, ny)] = {}

      for lang_name, lang_info in LANGUAGES.items():
        dirpath = os.path.join(base_dir, lang_info["dir"])
        src_path = os.path.join(dirpath, lang_info["src"])

        orig_nx, orig_ny = get_nx_ny(src_path, lang_info["re_nx"], lang_info["re_ny"])
        print(
          f"\n  [{lang_name}] (original NX={orig_nx}, NY={orig_ny}) -> NX={nx}, NY={ny}"
        )

        set_nx_ny(src_path, nx, ny, lang_info["re_nx"], lang_info["re_ny"])

        build_dir = os.path.join(dirpath, "build")
        os.makedirs(build_dir, exist_ok=True)

        t0 = time.perf_counter()
        ok = compile_language(lang_info, dirpath)
        t1 = time.perf_counter()
        tiempos_compilacion[(nx, ny)][lang_name] = t1 - t0

        if not ok:
          resultados[(nx, ny)][lang_name] = None
          continue

        t_comp = t1 - t0
        print(f"    Compilacion: {t_comp:.2f}s")

        times = benchmark(lang_info, dirpath)
        if times is None:
          resultados[(nx, ny)][lang_name] = None
        else:
          resultados[(nx, ny)][lang_name] = times
          media = sum(times) / len(times)
          print(f"    Media ({MEASURED} ejecuciones): {media:.4f}s")
          if len(times) > 1:
            print(f"    Min: {min(times):.4f}s  Max: {max(times):.4f}s")

  finally:
    for lang_name, lang_info in LANGUAGES.items():
      src_path = os.path.join(base_dir, lang_info["dir"], lang_info["src"])
      bak_path = src_path + ".bak"
      if os.path.exists(bak_path):
        shutil.copy2(bak_path, src_path)
        os.remove(bak_path)
        print(f"  Restaurado: {lang_info['src']}")

  print(f"\n{'=' * 70}")
  print("RESUMEN DE RESULTADOS")
  print(f"{'=' * 70}\n")

  headers = [
    "Config",
    "NX",
    "NY",
    "Lenguaje",
    "Media (s)",
    "Min (s)",
    "Max (s)",
    "Compil. (s)",
  ]
  rows = []
  for nx, ny, lbl in CONFIGS:
    for lang_name in LANGUAGES:
      times = resultados[(nx, ny)].get(lang_name)
      t_comp = tiempos_compilacion[(nx, ny)].get(lang_name, 0)
      if times is not None:
        media = sum(times) / len(times)
        rows.append(
          [
            lbl,
            nx,
            ny,
            lang_name,
            f"{media:.4f}",
            f"{min(times):.4f}",
            f"{max(times):.4f}",
            f"{t_comp:.2f}",
          ]
        )
      else:
        rows.append([lbl, nx, ny, lang_name, "FAIL", "FAIL", "FAIL", f"{t_comp:.2f}"])

  print(tabulate(rows, headers=headers, tablefmt="grid"))
  print()

  print("COMPARATIVA POR CONFIGURACION (media en segundos):\n")
  cmp_headers = ["Config"] + list(LANGUAGES)
  cmp_rows = []
  for nx, ny, lbl in CONFIGS:
    row = [f"{lbl} ({nx}x{ny})"]
    for lang_name in LANGUAGES:
      times = resultados[(nx, ny)].get(lang_name)
      row.append(f"{sum(times) / len(times):.4f}" if times else "FAIL")
    cmp_rows.append(row)
  print(tabulate(cmp_rows, headers=cmp_headers, tablefmt="grid", stralign="center"))


if __name__ == "__main__":
  main()
