% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Potencia_Efectos_Aleatorios.R
\name{simular_potencia}
\alias{simular_potencia}
\title{Simular potencia para ANOVA de efectos aleatorios balanceado}
\usage{
simular_potencia(t, r, sigma2 = 1, rho = 0.5, alpha = 0.05, nsim = 1000)
}
\arguments{
\item{t}{Número de tratamientos \eqn{t}.}

\item{r}{Número de réplicas por tratamiento \eqn{r}.}

\item{sigma2}{Varianza del error \eqn{\sigma^2}. (default: 1)}

\item{rho}{Razón de varianzas \eqn{\rho = \sigma_\tau^2 / \sigma^2}. (default: 0.5)}

\item{alpha}{Nivel de significancia. (default: 0.05)}

\item{nsim}{Número de simulaciones de Monte Carlo. (default: 1000)}
}
\value{
Potencia empírica (valor entre 0 y 1).
}
\description{
Ejecuta una simulación de Monte Carlo para estimar la potencia de la prueba
\eqn{F} en un diseño a una vía con efectos aleatorios y \emph{r} réplicas por tratamiento.
Cada observación se genera como
\deqn{Y_{ij} = \tau_i + \varepsilon_{ij},}
donde \eqn{\tau_i \sim N(0,\sigma_\tau^2)} y \eqn{\varepsilon_{ij} \sim N(0,\sigma^2)}.
}
\details{
La función calcula la proporción de simulaciones donde el estadístico F observado
supera el valor crítico bajo la hipótesis nula. Es útil para verificar la potencia teórica.
}
\examples{
# Simular potencia con 5 tratamientos, 10 réplicas y rho = 0.4
simular_potencia(t = 5, r = 10, rho = 0.4, nsim = 1000)

}
