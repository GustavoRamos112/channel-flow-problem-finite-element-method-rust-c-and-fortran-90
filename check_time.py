
import time
import subprocess
import os
# pyrefly: ignore [untyped-import]
from tabulate import tabulate
from tqdm import tqdm

class benchmark:
  def __init__(
    self, 
    lenguajes: dict[str, tuple[str, str]], 
    repeticiones: int = 10, 
    archivo_times: str = "times.txt",
    compilar: bool = True) -> None:
    self.lenguajes: dict[str, tuple[str, str]] = lenguajes
    self.compilar: bool = compilar
    self.repeticiones: int = repeticiones
    self.archivo_times: str = archivo_times

  def compilar_programa(self, dir: str, comando: str, rust = False) -> None:
    print(f"{dir}: ", end="")
    original_dir: str = os.getcwd()
    start: float = time.perf_counter()
    os.chdir(dir)
    if not os.path.exists("build") and not rust:
      os.mkdir("build")
    subprocess.run(comando, shell=True, stderr=subprocess.DEVNULL)
    os.chdir(original_dir)
    end: float = time.perf_counter()
    print(f"{(end - start):.2f}")

  def check_times_file(self, dir_times: str) -> None:
    if os.path.exists(dir_times):
      os.remove(dir_times)

    if not os.path.exists(dir_times):
      with open(dir_times, "w") as archivo:
        archivo.write("")

  def medir_tiempo(self, directorio: str, comando: list, repeticiones: int, fortran: bool = False) -> tuple[list[float], float]:
      #? Cambiar de directorio de forma segura
      original_dir: str = os.getcwd()
      os.chdir(directorio)

      #? 1. Warm-up (Lanzar una vez para cargar en caché)
      if fortran:
        try:
          os.remove("uv.txt")
          os.remove("xy.txt")
        except:
          pass
      subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

      tiempos: list[float] = []
      for _ in tqdm(range(repeticiones)):
          if fortran:
            try:
              os.remove("uv.txt")
              os.remove("xy.txt")
            except:
              pass
          #? 2. Usar perf_counter para alta precisión
          start: float = time.perf_counter()
          subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
          end: float = time.perf_counter()
          tiempos.append(end - start)

      os.chdir(original_dir)

      # Retornar el promedio
      return (tiempos, sum(tiempos) / len(tiempos))

  def realizar_benchmark(self) -> None:
    tiempos_generales: list[list[float]] = [[] for _ in range(self.repeticiones)]
    medias = {}
    i: int = 0

    if self.compilar:
      print("Compilando...")
      for key, value in self.lenguajes.items():
        nombre_l = value[2] == "rust"
        self.compilar_programa(
          key, value[0], 
          rust = nombre_l
        )
    
    for key, value in self.lenguajes.items():
      print(f"{key}:")
      directorio_programa: str = key
      comando: list[str] = [value[1]]

      if value[2] != "fortran":
        dir_times: str = os.path.join(directorio_programa, self.archivo_times)
        self.check_times_file(dir_times)
      
      es_fortran = value[2] == "fortran"
      tiempos_todos, media_python = self.medir_tiempo(
        directorio=directorio_programa, comando=comando,
        repeticiones=self.repeticiones,
        fortran = es_fortran
      )

      if value[2] == "fortran":
        i = 0
        for tiempo in tiempos_todos:
          tiempos_generales[i].append(tiempo)
          i += 1
        medias[key] = media_python
        continue

      with open(dir_times, "r") as archivo:
        tiempos_archivo: list[float] = [float(l) for l in archivo.readlines()]
      
      tiempos_archivo.pop(0)
      media_individual: float = sum(tiempos_archivo) / len(tiempos_archivo)
      
      i = 0
      for tiempo in tiempos_archivo:
        tiempos_generales[i].append(tiempo)
        i += 1
      
      medias[key] = media_individual

    print("Repeticiones: ", self.repeticiones)

    print(tabulate(
      tiempos_generales,
      headers=["Tiempo F90", "Tiempo Cpp", "Tiempo Rust"],
      tablefmt="grid")
    )

    tabla = [
      [key, medias[key]]
      for key in medias.keys()
    ]

    tabla.sort(key=lambda x: x[1])
    for i in range(len(tabla)):
      tabla[i].append(f"X{tabla[i][1]/tabla[0][1]:.3}")

    print(tabulate(
      tabla,
      headers=["Lenguaje", "Media", "Comparativa"],
      tablefmt="grid",
      stralign="center",
    ))


if __name__ == "__main__":
  #? Formato:
  #? "Directorio": (
  #?    "Comando de compilación",
  #?    "Ubicacion del ejecutable",
  #?    "Nombre del lenguaje"
  #?  )
  lenguajes: dict[str, tuple[str, str]] = {
    "Fortran90": (
      "flang main.f90 -O3 -o build/main.exe",
      "build/main.exe",
      "fortran"
    ),
    "Cpp": (
      "clang++ main.cpp -std=c++23 -O3 -o build/main.exe",
      "build/main.exe",
      "cpp"
    ),
    "Rust": (
      "cargo build --release",
      "./target/release/rust.exe",
      "rust"
    )
  }

  check = benchmark(
    lenguajes = lenguajes, 
    repeticiones = 10, 
    compilar = False
  )
  check.realizar_benchmark()