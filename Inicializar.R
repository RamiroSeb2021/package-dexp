usethis::proj_activate("G:/Mi unidad/RSTUDIO_/Diseño_Experimentos/Proyecto_APP/Proyecto_DEXP/dexp")

devtools::document()

devtools::install()

library(dexp)

# Probar con ayuda
?calcular_r_MT

# Ejecutar ejemplo real

calcular_r_MT(
  T_ = 6,
  D = 20,
  ro = 6,
  S1 = sqrt(141.6),
  df1 = 40
)

# Verificar que todo este limpio
devtools::check()


usethis::use_readme_rmd()

devtools::build_readme()
usethis::use_pkgdown()
pkgdown::build_site()

install.packages("pkgbuild")  # si no lo tienes aún
pkgbuild::find_rtools()
pkgbuild::rtools_path()

# Agregar un readme con ejemplos

usethis::use_readme_rmd()
devtools::build_readme()

# Agrega documentación web con pkgdown

usethis::use_pkgdown()
pkgdown::build_site()

devtools::install_github("RamiroSeb2021/package-dexp/")

