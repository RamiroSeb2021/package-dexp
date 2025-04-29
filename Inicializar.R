usethis::proj_activate("G:/Mi unidad/RSTUDIO_/Diseño_Experimentos/Proyecto_APP/Proyecto_DEXP/dexp")


remove.packages("dexp")
devtools::clean_dll()

devtools::document()   # genera documentación
devtools::build()      # empaqueta
devtools::install()    # instala


library(help = dexp)

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

?encontrar_r_minimo

resultado <- encontrar_r_minimo(t = 5, rho = 0.4, potencia_objetivo = 0.8)
resultado$r_optimo
resultado$potencia
resultado$grafico
resultado$tabla

# Verificar que todo este limpio
devtools::check()
usethis::use_package("ggplot2")



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

devtools::install_github("RamiroSeb2021/package-dexp")

