#' @import stats
#' @import ggplot2
#' @importFrom stats qf qtukey ptukey qchisq qt
#' @importFrom utils install.packages

NULL

#' Simular potencia para ANOVA de efectos aleatorios balanceado
#'
#' Ejecuta una simulación de Monte Carlo para estimar la potencia de la prueba
#' \eqn{F} en un diseño a una vía con efectos aleatorios y \emph{r} réplicas por tratamiento.
#' Cada observación se genera como
#' \deqn{Y_{ij} = \tau_i + \varepsilon_{ij},}
#' donde \eqn{\tau_i \sim N(0,\sigma_\tau^2)} y \eqn{\varepsilon_{ij} \sim N(0,\sigma^2)}.
#'
#' @param t Número de tratamientos \eqn{t}.
#' @param r Número de réplicas por tratamiento \eqn{r}.
#' @param sigma2 Varianza del error \eqn{\sigma^2}. (default: 1)
#' @param rho Razón de varianzas \eqn{\rho = \sigma_\tau^2 / \sigma^2}. (default: 0.5)
#' @param alpha Nivel de significancia. (default: 0.05)
#' @param nsim Número de simulaciones de Monte Carlo. (default: 1000)
#'
#' @return Potencia empírica (valor entre 0 y 1).
#'
#' @details
#' La función calcula la proporción de simulaciones donde el estadístico F observado
#' supera el valor crítico bajo la hipótesis nula. Es útil para verificar la potencia teórica.
#'
#' @examples
#' # Simular potencia con 5 tratamientos, 10 réplicas y rho = 0.4
#' simular_potencia(t = 5, r = 10, rho = 0.4, nsim = 1000)
#'
#' @export
simular_potencia <- function(t, r, sigma2 = 1, rho = 0.5, alpha = 0.05, nsim = 1000) {
  sigma2_tau <- rho * sigma2
  f_crit <- qf(1 - alpha, df1 = t - 1, df2 = t * (r - 1))
  rechazos <- numeric(nsim)

  for (i in 1:nsim) {
    tau <- rnorm(t, mean = 0, sd = sqrt(sigma2_tau))
    datos <- vector()
    for (j in 1:t) {
      errores <- rnorm(r, mean = 0, sd = sqrt(sigma2))
      y <- tau[j] + errores
      datos <- c(datos, y)
    }
    grupo <- factor(rep(1:t, each = r))
    modelo <- aov(datos ~ grupo)
    F_obs <- summary(modelo)[[1]]["grupo", "F value"]
    rechazos[i] <- as.numeric(F_obs > f_crit)
  }

  mean(rechazos)
}

#' Encontrar el tamaño mínimo de réplica que logra la potencia objetivo
#'
#' Incrementa el número de réplicas \eqn{r} desde 2 hasta \code{r_max},
#' llamando a \code{simular_potencia()} en cada paso, hasta alcanzar la potencia
#' deseada. Además, construye el gráfico de la curva de potencia y devuelve
#' la tabla de resultados.
#'
#' @param t Número de tratamientos \eqn{t}.
#' @param rho Razón de varianzas \eqn{\rho = \sigma_\tau^2 / \sigma^2}.
#' @param potencia_objetivo Potencia deseada (1 - β). (default: 0.8)
#' @param sigma2 Varianza del error \eqn{\sigma^2}. (default: 1)
#' @param alpha Nivel de significancia. (default: 0.05)
#' @param nsim Número de simulaciones de Monte Carlo por valor de \eqn{r}. (default: 1000)
#' @param r_max Límite superior para la búsqueda de \eqn{r}. (default: 50)
#'
#' @return Lista con:
#' \describe{
#'   \item{r_optimo}{Primer \eqn{r} con potencia ≥ \code{potencia_objetivo}.}
#'   \item{potencia}{Potencia empírica alcanzada para \eqn{r_optimo}.}
#'   \item{grafico}{Gráfico de la potencia simulada vs \eqn{r}.}
#'   \item{tabla}{Data frame con los valores de \eqn{r} y sus respectivas potencias \eqn{1-\beta}.}
#' }
#' Devuelve \code{NULL} (y lanza \code{warning}) si no se alcanza la potencia deseada.
#'
#' @examples
#' # Encontrar el tamaño mínimo de réplicas para rho = 0.4
#' resultado <- encontrar_r_minimo(t = 5, rho = 0.4, potencia_objetivo = 0.8)
#' resultado$r_optimo
#' resultado$potencia
#' resultado$grafico
#' resultado$tabla
#'
#' @export
encontrar_r_minimo <- function(t, rho, potencia_objetivo = 0.8,
                               sigma2 = 1, alpha = 0.05, nsim = 1000,
                               r_max = 50) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)

  resultados <- data.frame(r = integer(0), potencia = numeric(0))

  for (r in 2:r_max) {
    potencia <- simular_potencia(t = t, r = r, sigma2 = sigma2,
                                 rho = rho, alpha = alpha, nsim = nsim)
    # cat("Potencia estimada para r =", r, ":", round(potencia, 4), "\n")
    resultados <- rbind(resultados, data.frame(r = r, potencia = potencia))

    if (potencia >= potencia_objetivo) {
      grafico <- ggplot2::ggplot(resultados, ggplot2::aes(x = r, y = potencia)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_text(ggplot2::aes(label = r), vjust = -1, size = 3) +
        ggplot2::geom_hline(yintercept = potencia_objetivo, linetype = "dashed", color = "red") +
        ggplot2::labs(
          title = "Curva de potencia estimada",
          x = "Número de réplicas (r)",
          y = expression("Potencia (1 - " * beta * ")")
        ) +
        ggplot2::theme_minimal()


      # cat("Se alcanza potencia >= ", potencia_objetivo, " con r =", r, "\n")
      return(list(
        r_optimo = r,
        potencia = potencia,
        grafico = grafico,
        tabla = resultados
      ))
    }
  }

  warning("No se alcanzó la potencia deseada con r_max = ", r_max)
  NULL
}
