#' @import stats
#' @import ggplot2
#' @importFrom stats qf qtukey ptukey qchisq qt
#' @importFrom utils install.packages

NULL

#' Simulación de potencia para ANOVA con efectos aleatorios (F no central)
#'
#' Ejecuta una simulación de Monte Carlo para estimar la potencia de la prueba
#' \eqn{F} en un diseño de una vía con efectos aleatorios balanceado y \emph{r} réplicas por tratamiento.
#' Cada observación se genera como:
#' \deqn{Y_{ij} = \tau_i + \varepsilon_{ij},}
#' donde \eqn{\tau_i \sim \mathcal{N}(0,\sigma_\tau^2)} representa el efecto aleatorio del tratamiento,
#' y \eqn{\varepsilon_{ij} \sim \mathcal{N}(0,\sigma^2)} es el error experimental.
#'
#' @param t Número de tratamientos \eqn{t}.
#' @param r Número de réplicas por tratamiento \eqn{r}.
#' @param d Diferencia mínima estandarizada de interés \eqn{d}, usada para el cálculo del parámetro de no centralidad.
#' @param sigma2 Varianza del error \eqn{\sigma^2}. (por defecto: 1)
#' @param rho Razón de varianzas \eqn{\rho = \sigma_\tau^2 / \sigma^2}. (por defecto: 0.5)
#' @param alpha Nivel de significancia. (por defecto: 0.05)
#' @param nsim Número de simulaciones. (por defecto: 1000)
#'
#' @return Valor numérico entre 0 y 1 que representa la potencia empírica.
#'
#' @details
#' El parámetro de no centralidad se calcula como:
#' \deqn{\phi^2 = \frac{r d^2}{2 \sigma^2 t}}.
#' Se compara el estadístico F observado con el cuantil correspondiente de la distribución
#' \eqn{F_{t-1,t(r-1)}(\phi)}. La función devuelve la proporción de rechazos de la hipótesis nula.
#'
#' @examples
#' # Simular potencia con 5 tratamientos, 10 réplicas y rho = 0.4
#' simular_potencia_FNC(t = 5, r = 10, d = 1, sigma2 = 1, rho = 0.4, nsim = 1000)
#'
#' @export

simular_potencia_FNC <- function(t, 
                                 r,
                                 d,
                                 sigma2,
                                 rho = 0.5,
                                 alpha = 0.05,
                                 nsim = 1000) {
  sigma2_tau <- rho * sigma2
  phi <- sqrt((r*d^2)/(2*t*sigma2))
  f_crit <- qf(1 - alpha, df1 = t - 1, df2 = t * (r - 1), ncp = phi)
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


#' Determinar el número mínimo de réplicas para alcanzar una potencia deseada
#'
#' Ejecuta una búsqueda secuencial del número de réplicas \eqn{r} requerido para alcanzar una potencia deseada
#' en un diseño con efectos aleatorios. Utiliza simulaciones de Monte Carlo y grafica la curva de potencia.
#'
#' @param t Número de tratamientos \eqn{t}.
#' @param rho Razón de varianzas \eqn{\rho = \sigma_\tau^2 / \sigma^2}.
#' @param dif Diferencia mínima estandarizada de interés \eqn{d}.
#' @param potencia_objetivo Potencia deseada (1 - \eqn{\beta}). (por defecto: 0.8)
#' @param sigma2_ Varianza del error \eqn{\sigma^2}. (por defecto: 1)
#' @param alpha Nivel de significancia. (por defecto: 0.05)
#' @param nsim Número de simulaciones por valor de \eqn{r}. (por defecto: 1000)
#' @param r_max Límite superior para la búsqueda del número de réplicas. (por defecto: 50)
#'
#' @return Una lista con:
#' \describe{
#'   \item{r_optimo}{Primer valor de \eqn{r} que alcanza la potencia deseada.}
#'   \item{potencia}{Potencia estimada alcanzada para \eqn{r_optimo}.}
#'   \item{grafico}{Gráfico de la curva de potencia simulada.}
#'   \item{tabla}{Data frame con los valores de \eqn{r} y sus respectivas potencias.}
#' }
#' Si no se alcanza la potencia deseada, se devuelve \code{NULL} y un \code{warning}.
#'
#' @examples
#' # Encontrar r mínimo para potencia deseada 0.8
#' resultado <- encontrar_r_minimo_Potencia(t = 4, rho = 0.5, dif = 3, potencia_objetivo = 0.8, sigma2_ = 10.35)
#' resultado$r_optimo
#' resultado$grafico
#'
#' @export

encontrar_r_minimo_Potencia <- function(t, 
                                rho, 
                                dif,
                                potencia_objetivo = 0.8,
                                sigma2_, 
                                alpha = 0.05, 
                                nsim = 1000,
                                r_max = 50) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  
  resultados <- data.frame(r = integer(0), potencia = numeric(0))
  
  for (r in 2:r_max) {
    potencia <- simular_potencia_FNC(t = t, r = r, d = dif, sigma2 = sigma2_,
                                 rho = rho, alpha = alpha, nsim = nsim)
    # cat("Potencia estimada para r =", r, ":", round(potencia, 4), "\n")
    resultados <- rbind(resultados, data.frame(r = r, potencia = round(potencia, 2)))
    
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



