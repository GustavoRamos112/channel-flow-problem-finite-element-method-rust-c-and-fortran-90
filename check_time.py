
import time
import subprocess
import os
# pyrefly: ignore [untyped-import]
from tabulate import tabulate
from tqdm import tqdm

def check_times_file (dir_times: str) -> None:
  if os.path.exists(dir_times):
    os.remove(dir_times)

  if not os.path.exists(dir_times):
    with open(dir_times, "w") as archivo:
      archivo.write("")

def medir_tiempo(directorio: str, comando: list, repeticiones: int, fortran: bool = False) -> tuple[list[float], float]:
    #? Cambiar de directorio de forma segura
    original_dir: str = os.getcwd()
    os.chdir(directorio)

    #? 1. Warm-up (Lanzar una vez para cargar en caché)
    if fortran:
      try:
        subprocess.run("del uv.txt")
        subprocess.run("del xy.txt")
      except:
        pass
    subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    tiempos: list[float] = []
    for _ in tqdm(range(repeticiones)):
        if fortran:
          subprocess.run("del uv.txt")
          subprocess.run("del xy.txt")
        #? 2. Usar perf_counter para alta precisión
        start: float = time.perf_counter()
        subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        end: float = time.perf_counter()
        tiempos.append(end - start)

    os.chdir(original_dir)

    # Retornar el promedio
    return (tiempos, sum(tiempos) / len(tiempos))

if __name__ == "__main__":
  #! Bandera para compilar
  compilar = False

  if compilar:
    print("Compilando...")
    original_dir: str = os.getcwd()
    
    print("Cpp: ", end="")
    start: float = time.perf_counter()
    os.chdir("Cpp")
    if not os.path.exists("build"):
      os.mkdir("build")
    subprocess.run("clang++ main.cpp -std=c++23 -O3 -o build\\main.exe", shell=True, stderr=subprocess.DEVNULL)
    os.chdir(original_dir)
    end: float = time.perf_counter()
    print(end - start)

    print("Fortran: ", end="")
    start: float = time.perf_counter()
    os.chdir("Fortran90")
    if not os.path.exists("build"):
      os.mkdir("build")
    subprocess.run("flang main.f90 -O3 -o build\\main.exe", shell=True, stderr=subprocess.DEVNULL)
    os.chdir(original_dir)
    end: float = time.perf_counter()
    print(end - start)

    print("Rust: ", end="")
    start: float = time.perf_counter()
    os.chdir("Rust")
    subprocess.run("cargo build --release", shell=True, stderr=subprocess.DEVNULL)
    os.chdir(original_dir)
    end: float = time.perf_counter()
    print(end - start)

  _repetitions = 10

  #! fortran
  print("fortran")
  _dir = "fortran90"
  _command: list[str] = ["build/main.exe"]

  tiempos_f90, media_f90 = medir_tiempo(
    directorio=_dir, comando=_command,
    repeticiones=_repetitions
  )

  #! C++
  print("C++")
  _dir = "cpp"
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
  print("rust")
  _dir = "Rust"
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
  tabla = [
    ["f90", media_f90],
    ["c++", media_cpp],
    ["rust", media_rust]
  ]

  tabla.sort(key=lambda x: x[1])
  print(tabla)
  for i in range(len(tabla)):
    tabla[i].append(f"X{tabla[i][1]/tabla[0][1]:.3}")
  print(tabulate(
    tabla,
    headers=["Lenguaje", "Media", "Comparativa"],
    tablefmt="grid",
    stralign="center",
  ))
  