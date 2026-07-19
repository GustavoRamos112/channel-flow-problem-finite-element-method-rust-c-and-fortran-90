# Notas de la traduccion Fortran -> Rust

- Todas las variables que se tomaban como globales o se pasaban de una funcion a otra, ahora se guardan en un struct `flow_struct.rs`

- Muchas variables que originalmente se declaraban en el main ahora se declaran localmente cada que se usan es por ello que no aparecen directamente

- El vector `node` en el struct ya se ha modificado para que sea 0-index

- La matriz `index` se mantiene en 1-index ya que es mas dolor de abeza ambiarlo que dejarlo asi

- El vector `iline` esta en 1-index ya que se llena con la matriz `index`

- `fileg` -> gram

- `filex` -> `xy_plot3d` y `xy_table`

- `fileu` -> `uv_plot3d` y `uv_table`