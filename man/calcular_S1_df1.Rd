% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metodo_tukey.R
\name{calcular_S1_df1}
\alias{calcular_S1_df1}
\title{Estimación de la desviación estándar corregida y grados de libertad}
\usage{
calcular_S1_df1(
  desviacion_estandar,
  Si,
  Ss,
  max_error = 0.01,
  confianza = 0.9,
  maximum_df = 1000
)
}
\arguments{
\item{desviacion_estandar}{Desviación estándar inicial (en unidades de la variable estudiada).}

\item{Si}{Porcentaje relativo inferior aceptado (por ejemplo, 0.07 corresponde a un 7\%).}

\item{Ss}{Porcentaje relativo superior aceptado (por ejemplo, 0.12 corresponde a un 12\%).}

\item{max_error}{Error relativo máximo permitido entre el cociente observado y el cociente esperado (por defecto \code{0.01}).}

\item{confianza}{Nivel de confianza del intervalo considerado (por defecto \code{0.9}).}

\item{maximum_df}{Máximo número de grados de libertad que se evaluarán (por defecto \code{1000}).}
}
\value{
Una lista con:
\describe{
\item{S1}{Estimación corregida de la desviación estándar.}
\item{grados_libertad}{Número de grados de libertad estimado.}
\item{valor_x}{Cociente entre cuantiles de chi-cuadrado.}
\item{error_relativo}{Error relativo del cociente respecto al objetivo.}
}
}
\description{
Esta función estima una desviación estándar corregida (\code{S1}) y el número de grados de libertad
(\code{df1}) que reproduce, mediante cocientes de chi-cuadrado, la relación entre un límite superior
y un límite inferior dados. Se usa un criterio de error relativo máximo permitido para aceptar el resultado.
}
\details{
Se busca el número de grados de libertad tal que el cociente de los valores críticos de la distribución chi-cuadrado,
evaluados a \code{p = confianza} y \code{p = 1 - confianza}, sea cercano al cociente observado entre \code{Ss} y \code{Si}.
El proceso finaliza tan pronto como el error relativo esté por debajo del máximo permitido (\code{max_error}).

El parámetro \code{confianza} actualmente se fija en los cálculos internos.
}
\examples{
set.seed(123)
calcular_S1_df1(desviacion_estandar = 30, Si = 0.07, Ss = 0.12)

}
