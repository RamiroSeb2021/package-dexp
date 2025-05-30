% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metodo_Harris_Hurvitz_Mood.R
\name{k_tabla_mc}
\alias{k_tabla_mc}
\title{Calcular el valor k = a/s1 para tablas de potencia (simulación Monte Carlo)}
\usage{
k_tabla_mc(n, m, beta = 0.8, alpha = 0.05, B = 2e+06, seed = 1)
}
\arguments{
\item{n}{Grados de libertad de la desviación estándar muestral nueva \eqn{s2}.}

\item{m}{Grados de libertad de la desviación estándar preliminar \eqn{s1}.}

\item{beta}{Potencia deseada de la prueba (\eqn{1 - \text{error tipo II}}), por defecto 0.80.}

\item{alpha}{Nivel de significancia de la prueba (cola superior), por defecto 0.05.}

\item{B}{Número de réplicas de Monte Carlo para estimar potencia, por defecto 2e6.}

\item{seed}{Semilla para reproducibilidad de la simulación, por defecto 1.}

\item{tol}{Tolerancia de la búsqueda de raíces con \code{uniroot()}, por defecto 0.001.}
}
\value{
Valor numérico \eqn{k = a/s1} que garantiza la potencia \code{beta} al nivel \code{alpha}.
}
\description{
Esta función estima el valor \eqn{k = a/s1} necesario para que una prueba de Student
unitaria (una sola cola), con \eqn{\alpha} nivel de significancia y \eqn{\beta} potencia,
rechace la hipótesis nula \eqn{\mu = 0} cuando la media verdadera es \eqn{\mu = a}.
Utiliza simulación Monte Carlo de las distribuciones muestrales de \eqn{s1} y \eqn{s2}.
}
\details{
\describe{
\item{1.}{Se simulan \eqn{B} valores de \eqn{s1 \sim \sqrt{\chi^2_m / m}} y \eqn{s2 \sim \sqrt{\chi^2_n / n}}.}
\item{2.}{Para cada \eqn{k} candidato, se fija \eqn{a = k \cdot s1} y se simula \eqn{\bar X \sim N(a, 1/(n+1))}.}
\item{3.}{Se computa \eqn{t = \sqrt{n+1} \cdot \bar X / s2}.}
\item{4.}{La potencia empírica es la proporción de réplicas donde \eqn{t > t_{1-\alpha}(n)}.}
\item{5.}{\code{uniroot()} busca el \eqn{k} tal que la potencia simulada sea igual a \code{beta}.}
}
}
\examples{
# Generar la tabla completa de k para beta = 0.95, alpha = 0.05
n_vals <- 1:10
m_vals <- c(1:6, 8, 12, 16, 24, 32)

# Inicializar matriz vacía con nombres de fila/columna
k_matrix <- matrix(NA,
                   nrow = length(n_vals),
                   ncol = length(m_vals),
                   dimnames = list(paste0("n=", n_vals),
                                   paste0("m=", m_vals)))

# Llenar matriz con k calculados
for (i in seq_along(n_vals)) {
  for (j in seq_along(m_vals)) {
    n <- n_vals[i]
    m <- m_vals[j]
    k_matrix[i, j] <- k_tabla_mc(n, m,
                                  beta = 0.95,
                                  alpha = 0.05,
                                  B = 5e5   # ajustar B si se desea más precisión
                                  )
  }
}

# Convertir a data.frame para mostrar
k_df <- as.data.frame(k_matrix, row.names = NULL)
k_df <- cbind(n = n_vals, k_df)
print(k_df, digits = 4)
View(round(k_df[-1], 3))

}
