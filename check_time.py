
import time
import subprocess
import os
# pyrefly: ignore [untyped-import]
from tabulate import tabulate

def check_times_file (dir_times: str) -> None:
  if os.path.exists(dir_times):
    os.remove(dir_times)

  if not os.path.exists(dir_times):
    with open(dir_times, "w") as archivo:
      archivo.write("")

def medir_tiempo(directorio: str, comando: list, repeticiones: int) -> tuple[list[float], float]:
    # Cambiar de directorio de forma segura
    original_dir: str = os.getcwd()
    os.chdir(directorio)

    # 1. Warm-up (Lanzar una vez para cargar en caché)
    subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    tiempos: list[float] = []
    for _ in range(repeticiones):
        # 2. Usar perf_counter para alta precisión
        start: float = time.perf_counter()
        subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        end: float = time.perf_counter()
        tiempos.append(end - start)

    os.chdir(original_dir)

    # Retornar el promedio
    return (tiempos, sum(tiempos) / len(tiempos))

if __name__ == "__main__":
  _repetitions = 10

  #! fortran
  _dir = "fortran"
  _command: list[str] = ["build/bump.exe"]

  archivo_times: str = "times.txt"
  dir_times: str = f"{_dir}/{archivo_times}"

  check_times_file(dir_times)

  tiempos_f90, media_f90 = medir_tiempo(
    directorio=_dir, comando=_command,
    repeticiones=_repetitions
  )

  #! C++
  _dir = "c++"
  _command: list[str] = ["build/main.exe"]

  archivo_times: str = "times.txt"
  dir_times: str = f"{_dir}/{archivo_times}"

  check_times_file(dir_times)

  medir_tiempo(
    directorio=_dir, comando=_command,
    repeticiones=_repetitions
  )

  with open(dir_times, "r") as archivo:
    tiempos_cpp: list[float] = [float(l) for l in archivo.readlines()]

  tiempos_cpp.pop(0)

  media_cpp: float = sum(tiempos_cpp) / len(tiempos_cpp)

  #! rust
  _dir = "rust"
  _command: list[str] = ["./target/release/rust.exe"]

  archivo_times: str = "times.txt"
  dir_times: str = f"{_dir}/{archivo_times}"

  check_times_file(dir_times)

  medir_tiempo(
    directorio=_dir, comando=_command,
    repeticiones=_repetitions
  )

  with open(dir_times, "r") as archivo:
    tiempos_rust: list[float] = [float(l) for l in archivo.readlines()]

  tiempos_rust.pop(0)

  media_rust: float = sum(tiempos_rust) / len(tiempos_rust)

  print("Repeticiones: ", _repetitions)

  tiempos_generales: list[list[float]] = [
    [tiempos_f90[i], tiempos_cpp[i], tiempos_rust[i]]
    for i in range(_repetitions)
  ]

  print(tabulate(
    tiempos_generales,
    headers=["Tiempo F90", "Tiempo Cpp", "Tiempo Rust"],
    tablefmt="grid")
  )
  print(tabulate(
    [
      ["f90", media_f90],
      ["c++", media_cpp],
      ["rust", media_rust]
    ],
    headers=["Lenguaje", "Media"],
    tablefmt="grid",
    stralign="center",
  ))
