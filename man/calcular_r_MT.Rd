% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metodo_tukey.R
\name{calcular_r_MT}
\alias{calcular_r_MT}
\title{Cálculo del número óptimo de réplicas en diseños con múltiples tratamientos}
\usage{
calcular_r_MT(T_, D, ro, S1, df1, alfa = 0.05, Beta = 0.1)
}
\arguments{
\item{T_}{Número de tratamientos.}

\item{D}{Diferencia mínima detectable entre medias (objetivo de longitud del intervalo: \code{D = 2A}).}

\item{ro}{Estimación inicial del número de réplicas.}

\item{S1}{Estimación de la desviación estándar experimental.}

\item{df1}{Grados de libertad del estimador de \code{S1}.}

\item{alfa}{Nivel de significancia (por defecto \code{0.05}).}

\item{Beta}{Error tipo II deseado (por defecto \code{0.1}).}
}
\value{
Una lista con los siguientes elementos:
\describe{
\item{r}{Número estimado de réplicas que cumple la potencia deseada.}
\item{A}{Valor de \code{A} asociado a \code{r}.}
\item{r_i}{Valor de \code{r} más cercano que logra que \code{2A ≈ D}.}
\item{A_i}{Valor de \code{A} correspondiente a \code{r_i}.}
\item{lista_r}{Vector de valores de \code{r} evaluados en la revisión del intervalo.}
\item{valores_A}{Vector de valores de \code{A} asociados a \code{lista_r}.}
\item{posicion_escogida}{Índice de \code{r_i} dentro de \code{lista_r}.}
}
}
\description{
Esta función estima el número óptimo de réplicas necesarias \code{r} en un diseño experimental
con múltiples tratamientos, de manera que la longitud del intervalo de confianza sea cercana
a un valor deseado (\code{D = 2A}). Utiliza un procedimiento iterativo basado en distribuciones
F y de Tukey para ajustar el número de réplicas y grados de libertad del error, y realiza una
revisión posterior para seleccionar el valor de \code{r} que produce un \code{A} más cercano a \code{D/2}.
}
\examples{
# Ejemplo básico
calcular_r_MT(
  T_ = 6, D = 20, ro = 5,
  S1 = sqrt(141.6), df1 = 40,
  alfa = 0.05, Beta = 0.1
)

# Prueba ejemplo 5.12
calcular_r_MT(
  T_ = 6,
  D = 20,
  ro = 6,
  S1 = sqrt(141.6),
  df1 = 40
)

# Prueba ejemplo 5.13
calcular_r_MT(
  T_ = 8,
  D = 500,
  ro = 4,
  S1 = sqrt(90000),
  df1 = 200,
  alfa = 0.1,
  Beta = 0.25
)

}
