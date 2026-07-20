# Comparativa C++/Rust/Fortran90 a través de un programa que resuelve el flujo de un fluido a través de un canal usando elementos finitos.

<span style="color:red">
  Importante: Algunas de las explicaciones se toman de las originales descritas en el propio código original de Fortran90 en Fortran90/main.f90, créditos a: John Burkardt
</span>

## Problema:

Flujo independiente del tiempo de un fluido viscoso incompresible en un canal bidimensional usando elementos finitos.

El problema de flujo de fluidos se formula en términos de las variables primitivas $u$, $v$ y $p$.

Este código intenta ajustar un perfil corriente abajo modificando un parámetro: el valor del parámetro de flujo entrante en la altura media del canal.

Las funciones lineales a trozos sobre triángulos aproximan la presión, y las funciones cuadráticas sobre triángulos aproximan la velocidad.

Esta es la base de elementos finitos de "Taylor-Hood". La formulación en variables primitivas de las ecuaciones de Navier-Stokes involucra la velocidad horizontal $U$, la velocidad vertical $V$ y la presión $P$. Las ecuaciones son:
```math
U \frac{\partial U}{\partial x} + 
  V \frac{\partial U}{\partial y} + 
  \frac{\partial P}{\partial x} - 
  \mu\left(\frac{\partial^2 U}{\partial x^2} + 
  \frac{\partial^2 U}{\partial y^2}\right) = F_1 \\
```
```math
U \frac{\partial V}{\partial x} + 
  V \frac{\partial V}{\partial y} + 
  \frac{\partial P}{\partial y} - 
  \mu\left(\frac{\partial^2 V}{\partial x^2} + 
  \frac{\partial^2 V}{\partial y^2}\right) = F_2 \\
```
```math
\frac{\partial U}{\partial x} + 
  \frac{\partial V}{\partial y} = 0
```

Al reformularse en la forma de elementos finitos, donde $\phi_i$ es la función de base común $i$-ésima para $U$ y $V$, y $\psi_i$ es la función de base $i$-ésima para $P$, estas ecuaciones se convierten en:
```math
\int \left( 
  U \frac{dU}{dx} \phi_i + 
  V \frac{dU}{dy} \phi_i - 
  P \frac{d\phi_i}{dx} + 
  \mu \left( 
    \frac{dU}{dx} \frac{d\phi_i}{dx} + 
    \frac{dU}{dy} \frac{d\phi_i}{dy} 
  \right) 
\right) = \int F_1 \phi_i  \\
```
```math
\int \left( 
  U \frac{dV}{dx} \phi_i + 
  V \frac{dV}{dy} \phi_i - 
  P \frac{d\phi_i}{dy} + 
  \mu \left( 
    \frac{dV}{dx} \frac{d\phi_i}{dx} + 
    \frac{dV}{dy} \frac{d\phi_i}{dy} 
  \right) 
\right) = \int F_2 \phi_i  \\
```
```math
\int \left( 
  \frac{dU}{dx} \psi_i + 
  \frac{dV}{dx} \psi_i 
\right) = 0
```

## Código y diferencias con el original

### Tipos de datos
El código original hace uso de los tipos `integer` y `real`; sin embargo, usando otro enfoque, usaremos los tipos `int` (`i32` para Rust y `int` en C++), `unsigned` (`usize` en Rust y `size_t` en C++) y `double` (`f64` en Rust y `double` en C++).

### Estructura de la memoria
En Fortran, las matrices y tensores funcionan como bloques de memoria contiguos (eficiente en memoria); sin embargo, en C++ y Rust, cada matriz y tensor está separado por renglones (`std::vector<std::vector<...>>` en C++ y `Vec<Vec<...>>` en Rust), lo cual es más ineficiente en memoria. Para comparar en las mismas condiciones, es necesario usar memoria continua.

Para conseguir esto, en lugar de usar vectores anidados, se usan vectores planos. Así, una matriz `a` de tamaño `n*m` que se declararía, por ejemplo, en C++ como:
```c++
  std::vector<
    std::vector<double>
  > a(n, std::vector<double>(m, 0.0));
  // Para acceder a un elemento i,j arbitrario
  std::println("{}", a[i][j]);
```
Se declara como:
```c++
  std::vector<double> a(n*m, 0.0);
  // Para acceder a un elemento i,j arbitrario
  std::println("{}", a[(i*m) + j]);
```

### Datos: el uso de un Struct
El código en su versión original estructura dentro de `main` todas las variables y constantes que se usarán a lo largo del programa, usando todas estas como variables locales que se pasan como referencias de escritura y lectura en cada función.

En esta versión se opta por el uso de un tipo `struct` (`Flow_struct` en C++ y Rust) que contiene todas las variables y constantes que se usarán a lo largo del programa. Así pasamos de funciones como:
```fortran
  call nstoke (a,area,f,g,indx,insc,ipivot,iwrite, &
    maxnew,maxrow,nelemn,neqn,nlband,nnodes,node, &
    np,nquad,nrow,numnew,para,phi,psi,reynld,tolnew,yc)
```
a funciones como:
```rust
  // Rust
  nstoke(&mut flow);
```
```c++
  // C++
  nstoke(flow);
```

Esto permite tener una mayor legibilidad en el código.

### Variables principales dentro del struct

`double_matrix a`, la matriz en bandas utilizada en `nstoke` y `solvlin`. Se almacena de acuerdo con el modo de almacenamiento de banda general de LINPACK/LAPACK.

`double_matrix area`, contiene el área de cada elemento.

`double_vector f`. Después de la llamada a `nstoke`, `f` contiene la solución del problema de Navier-Stokes actual.

`string data_dir`, el directorio donde se guardarán los archivos de datos.*

`string fileg`, el nombre del archivo de datos gráficos que se volcará para el programa DISPLAY.

`string fileu`, el nombre del archivo de datos gráficos UV que se volcará para el programa PLOT3D.

`string filex`, el nombre del archivo de datos gráficos XY que se volcará para el programa PLOT3D.

`double_vec g`. Después de la llamada a `solvlin`, `g` contiene las sensibilidades.

`int_vec indx`, registra, para cada nodo, si hay incógnitas de velocidad asociadas al nodo. `indx[i][0]` registra información para las velocidades horizontales. Si es 0, no hay incógnita asociada al nodo y se asume una velocidad horizontal cero. Si es -1, no hay incógnita asociada al nodo, pero se especifica una velocidad mediante la rutina UBDRY. Si es positivo, el valor es el índice en los arreglos `f` y `g` del coeficiente desconocido. `indx[i][1]` registra información para las velocidades verticales. Si es 0, no hay incógnita asociada al nodo y se asume una velocidad vertical cero. Si es positivo, el valor es el índice en los arreglos `f` y `g` del coeficiente desconocido.

`int_vec insc`, registra, para cada nodo, si hay una incógnita de presión asociada al nodo. Si `insc[i]` es cero, no hay incógnita asociada al nodo y se asume una presión de 0. Si `insc[i]` es positivo, el valor es el índice en los arreglos `f` y `g` del coeficiente desconocido.

`unsigned_vec ipivot`, información de pivotes utilizada por el solucionador lineal.

`int_vec isotri`, registra si un elemento dado es isométrico o no. `isotri[i]` es 0 si el elemento no es isométrico, y 1 si es isométrico.

`unsigned iwrite`, controla la cantidad de salida generada por el programa. 0, salida mínima. 1, salida normal, más archivos de datos gráficos creados por GDUMP, UVDUMP y XYDUMP. 2, salida abundante.

`bool long`, `.TRUE.` si la región es "larga y delgada", y `.FALSE.` si la región es "alta y estrecha". Esto determina cómo se numeran los nodos, elementos y variables.

`unsigned maxeqn`, el número máximo de ecuaciones permitido.

`unsigned maxnew`, el número máximo de pasos de Newton por iteración.

`unsigned maxrow`, el número máximo de filas reales de la matriz de coeficientes que se permite.

`unsigned maxsec`, el número máximo de pasos de la secante permitidos.

`unsigned mx`, `mx = 2*nx - 1`, el número total de puntos de la malla en el lado horizontal de la región.

`unsigned my`, `my = 2*ny - 1`, el número total de puntos de la malla en el lado vertical de la región.

`unsigned nelemn`, el número de elementos utilizados.

`unsigned neqn`, el número de ecuaciones o funciones para el sistema completo.

`unsigned nlband`, el número de diagonales por debajo de la diagonal principal de la matriz `a` que son distintas de cero.

`unsigned nnodes`, el número de nodos por elemento, 6.

`unsigned_matrix node`, registra los números globales de los 6 nodos que componen cada elemento. `node[i][j]` (con `j` de 0 a 5) contiene el número del nodo local `j` dentro del elemento `i`.

`unsigned nodex0`, el nodo con el número más bajo en la columna de nodos donde se mide el perfil.

`unsigned np`, el número de nodos.

`unsigned nquad`, el número de puntos de cuadratura, actualmente establecido en 3.

`unsigned nrow`, la dimensión de filas utilizada de la matriz de coeficientes.

`unsigned numnew`, número total de iteraciones de Newton realizadas en la rutina NSTOKE durante toda la ejecución.

`unsigned numsec`, el número de pasos de la secante realizados.

`unsigned nx`, el número de puntos "principales" de la malla en el lado horizontal de la región.

`unsigned ny`, el número de puntos "principales" de la malla en el lado vertical de la región.

`double_tensor phi`, Cada entrada de `phi` contiene el valor de una función de base cuadrática o su derivada, evaluada en un punto de cuadratura. En particular, `phi[i][j][k][0]` es el valor de la función de base cuadrática asociada con el nodo local `k` en el elemento `i`, evaluada en el punto de cuadratura `j`. `phi[i][j][k][1]` es la derivada X de esa misma función de base, y `phi[i][j][k][2]` es la derivada Y.

`double_tensor psi`, Cada entrada de `psi` contiene el valor de una función de base lineal evaluada en un punto de cuadratura. `psi[i][j][k]` es el valor de la función de base lineal asociada con el nodo local `k` en el elemento `i`, evaluada en el punto de cuadratura `j`.

`double_vec res`, contiene los residuos.

`double reynld`, el valor del número de Reynolds. En el sistema de unidades del programa, la viscosidad = 1 / `reynld`.

`double rjpnew`, la derivada con respecto al parámetro A del funcional J.

`double tolnew`, la tolerancia de convergencia para la iteración de Newton en NSTOKE.

`double tolsec`, la tolerancia de convergencia para la iteración de la secante en el programa principal.

`double_vec xc`, `xc[i]` es la coordenada X del nodo `i`.

`double xlngth`, la longitud de la región.

`double_matrix xm`, `xm[it][i]` es la coordenada X del `i`-ésimo punto de cuadratura en el elemento `it`.

`double_vec yc`, `yc[i]` es la coordenada Y del nodo `i`.

`double ylngth`, la altura de la región.

`double_matrix ym`, `ym[it][i]` es la coordenada Y del `i`-ésimo punto de cuadratura en el elemento `it`.

`bool save_times`, si es verdadero, guarda los tiempos de ejecución.*

`bool save_data`, si es verdadero, guarda los archivos de datos.*

`bool json`, si es verdadero, guarda los archivos de datos en formato JSON.*

\* Constantes nuevas.

## Requisitos

### Fortran
Flang o, en su defecto, cualquier compilador de Fortran.

### C++
Clang con soporte para C++23 o, en su defecto, cualquier compilador con soporte para C++23.

### Rust
Cargo y Rust.

## Compilación y comparativa

Primero es necesario clonar el repositorio:
``` bash
git clone https://github.com/GustavoRamos112/channel-flow-problem-finite-element-method-rust-c-and-fortran-90.git
```

### Compilación
En caso de C++ y Fortran se hará uso de la bandera de compilación `O3`. Esto se debe a que, para Rust, se hará uso de `cargo build --release`, ya que en su modo de depuración (debug) Rust activa muchas banderas que ayudan a una mejor depuración, pero hacen más lento el código en ejecución.

#### C++
```bash
cd C++
mkdir build
clang++ main.cpp -std=c++23 -O3 -o build/main.exe
```

#### Fortran
```bash
cd Fortran90
mkdir build
flang main.f90 -O3 -o build/main.exe
```

#### Rust
```bash
cd Rust
cargo build --release
```

### Ejecución

#### C++
```bash
cd C++
./build/main.exe
```

#### Fortran
```bash
cd Fortran90
./build/main.exe
```

#### Rust
```bash
cd Rust
cargo run --release
```

### Comparación

El script de Python `check_time.py` se encarga de ejecutar los programas una cantidad de repeticiones y, luego de comparar los tiempos de ejecución de los diferentes programas, obtiene el promedio mínimo y, en base a ese, compara el resto.

Para que esto funcione, `save_times` debe ser `true` en el código fuente, ya que Python lee estos valores para sacar el promedio (a diferencia de Fortran, donde a través de Python se mide el tiempo de ejecución y se compara).

El script ya compila cada programa. En caso de compilar todo manualmente, es necesario establecer `w_compilar = False` dentro del script.

De igual forma, si es necesario cambiar de compilador (usar gfortran o g++), es necesario modificar el comando de compilación dentro del script:
```python
  if w_compilar:
    print("Compilando...")
    compilar("Cpp", "clang++ main.cpp -std=c++23 -O3 -o build\\main.exe")
    compilar("Fortran90", "flang main.f90 -O3 -o build\\main.exe")
    compilar("Rust", "cargo build --release", rust = True)
```

Para ejecutar el script:
```bash
python check_time.py
```

En mi caso, este es el output de Python:

```bash
Compilando...
Cpp: 41.29041140002664
Fortran90: 4.927102700108662
Rust: 10.300957300001755
fortran
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 10/10 [00:02<00:00,  4.08it/s]
C++
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 10/10 [00:02<00:00,  4.41it/s]
rust
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 10/10 [00:01<00:00,  6.81it/s]
Repeticiones:  10
+--------------+--------------+---------------+
|   Tiempo F90 |   Tiempo Cpp |   Tiempo Rust |
+==============+==============+===============+
|     0.238903 |     0.213527 |      0.134147 |
+--------------+--------------+---------------+
|     0.301028 |     0.407633 |      0.126831 |
+--------------+--------------+---------------+
|     0.282615 |     0.177025 |      0.120838 |
+--------------+--------------+---------------+
|     0.209329 |     0.16664  |      0.125682 |
+--------------+--------------+---------------+
|     0.253211 |     0.176433 |      0.127628 |
+--------------+--------------+---------------+
|     0.211235 |     0.166562 |      0.138872 |
+--------------+--------------+---------------+
|     0.194851 |     0.176813 |      0.139741 |
+--------------+--------------+---------------+
|     0.228553 |     0.204286 |      0.134469 |
+--------------+--------------+---------------+
|     0.234353 |     0.165549 |      0.120816 |
+--------------+--------------+---------------+
|     0.284882 |     0.180622 |      0.118879 |
+--------------+--------------+---------------+

+------------+----------+---------------+
|  Lenguaje  |    Media |  Comparativa  |
+============+==========+===============+
|    Rust    | 0.12879  |     X1.0      |
+------------+----------+---------------+
|    C++     | 0.203509 |     X1.58     |
+------------+----------+---------------+
|    F90     | 0.243896 |     X1.89     |
+------------+----------+---------------+
```

## Licencia

<div align="center">

[![GNU GPLv3 Image](https://www.licen.cc/images/license-logos/licen.cc-mit.png)](https://mit-license.org/)

</div>

<div align="left">

Este proyecto está bajo la Licencia MIT. Permite el uso, copia, modificación, distribución, sublicenciamiento y venta del software de forma gratuita y sin restricciones. El único requisito es incluir el aviso de derechos de autor original y este texto de licencia en todas las copias o partes sustanciales del software. El código se proporciona "tal cual", sin garantías de ningún tipo.

</div>

## Creditos:

Codigo Original den Fortran 90: John Burkardt

Mas detalles sobre funciones especificas o heredadas, revisar `Fortran90/main.f90`.