% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Potencia.R
\name{encontrar_r_minimo_Potencia}
\alias{encontrar_r_minimo_Potencia}
\title{Determinar el número mínimo de réplicas para alcanzar una potencia deseada}
\usage{
encontrar_r_minimo_Potencia(
  t,
  rho,
  dif,
  potencia_objetivo = 0.8,
  sigma2_,
  alpha = 0.05,
  nsim = 1000,
  r_max = 50
)
}
\arguments{
\item{t}{Número de tratamientos \eqn{t}.}

\item{rho}{Razón de varianzas \eqn{\rho = \sigma_\tau^2 / \sigma^2}.}

\item{dif}{Diferencia mínima estandarizada de interés \eqn{d}.}

\item{potencia_objetivo}{Potencia deseada (1 - \eqn{\beta}). (por defecto: 0.8)}

\item{sigma2_}{Varianza del error \eqn{\sigma^2}. (por defecto: 1)}

\item{alpha}{Nivel de significancia. (por defecto: 0.05)}

\item{nsim}{Número de simulaciones por valor de \eqn{r}. (por defecto: 1000)}

\item{r_max}{Límite superior para la búsqueda del número de réplicas. (por defecto: 50)}
}
\value{
Una lista con:
\describe{
\item{r_optimo}{Primer valor de \eqn{r} que alcanza la potencia deseada.}
\item{potencia}{Potencia estimada alcanzada para \eqn{r_optimo}.}
\item{grafico}{Gráfico de la curva de potencia simulada.}
\item{tabla}{Data frame con los valores de \eqn{r} y sus respectivas potencias.}
}
Si no se alcanza la potencia deseada, se devuelve \code{NULL} y un \code{warning}.
}
\description{
Ejecuta una búsqueda secuencial del número de réplicas \eqn{r} requerido para alcanzar una potencia deseada
en un diseño con efectos aleatorios. Utiliza simulaciones de Monte Carlo y grafica la curva de potencia.
}
\examples{
# Encontrar r mínimo para potencia deseada 0.8
resultado <- encontrar_r_minimo_Potencia(t = 4, rho = 0.5, dif = 3, potencia_objetivo = 0.8, sigma2_ = 10.35)
resultado$r_optimo
resultado$grafico

}
