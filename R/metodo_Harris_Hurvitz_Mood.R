#' Calcular el valor \(k = a/s_1\) para tablas de potencia (simulación Monte Carlo)
#'
#' Esta función estima el valor \(k = a/s_1\) necesario para que una prueba de Student
#' unitaria (una sola cola), con \(\alpha\) nivel de significancia y \(\beta\) potencia,
#' rechace la hipótesis nula \(\mu=0\) cuando la media verdadera es \(\mu=a\).
#' Utiliza simulación Monte Carlo de las distribuciones muestrales de \(s_1\) y \(s_2\).
#'
#' @param n  Integer. Grados de libertad de la desviación estándar muestral nueva \(s_2\).
#' @param m  Integer. Grados de libertad de la desviación estándar preliminar \(s_1\).
#' @param beta Numeric. Potencia deseada de la prueba (\(1-\) error tipo II), por defecto 0.80.
#' @param alpha Numeric. Nivel de significancia de la prueba (cola superior), por defecto 0.05.
#' @param B Integer. Número de réplicas de Monte Carlo para estimar potencias, por defecto 2e6.
#' @param seed Integer. Semilla para reproducibilidad de la simulación, por defecto 1.
#' @param tol  Numeric. Tolerancia de la búsqueda de raíces con \code{uniroot()}, por defecto 0.001.
#'
#' @return Numeric. Estimación de \(k = a/s_1\) que logra potencia \code{beta} al nivel \code{alpha}.
#'
#' @details
#' 1. Se simulan \code{B} valores de
#'    \(\,s_1 \sim \sqrt{\chi^2_m / m}\) y \(\,s_2 \sim \sqrt{\chi^2_n / n}\).
#' 2. Para cada \(k\) candidato, se fija \(a = k\,s_1\), se simulan \(\bar X \sim N(a,\,1/(n+1))\)
#'    y se computa la estadística \(t = \sqrt{n+1}\,\bar X / s_2\).
#' 3. La potencia empírica es la proporción de réplicas donde \(t > t_{1-\alpha}(n)\).
#' 4. \code{uniroot()} busca el \(k\) tal que la potencia simulada sea igual a \code{beta}.
#'
#' @examples
#' # Generar la tabla completa de k para beta = 0.95, alpha = 0.05
#' n_vals <- 1:10
#' m_vals <- c(1:6, 8, 12, 16, 24, 32)
#'
#' # Inicializar matriz vacía con nombres de fila/columna
#' k_matrix <- matrix(NA,
#'                    nrow = length(n_vals),
#'                    ncol = length(m_vals),
#'                    dimnames = list(paste0("n=", n_vals),
#'                                    paste0("m=", m_vals)))
#'
#' # Llenar matriz con k calculados
#' for (i in seq_along(n_vals)) {
#'   for (j in seq_along(m_vals)) {
#'     n <- n_vals[i]
#'     m <- m_vals[j]
#'     k_matrix[i, j] <- k_tabla_mc(n, m,
#'                                   beta = 0.95,
#'                                   alpha = 0.05,
#'                                   B = 5e5,   # ajustar B si se desea más precisión
#'                                   tol = 1e-3)
#'   }
#' }
#'
#' # Convertir a data.frame para mostrar
#' k_df <- as.data.frame(k_matrix, row.names = NULL)
#' k_df <- cbind(n = n_vals, k_df)
#' print(k_df, digits = 4)
#' View(round(k_df[-1], 3))
#'
#' @export

k_tabla_mc <- function(n, m, beta = 0.80, alpha = 0.05,
                       B = 2e6, seed = 1) {
  set.seed(seed)
  tcrit <- qt(1 - alpha, df = n)
  s1 <- sqrt(rchisq(B, m) / m)
  s2 <- sqrt(rchisq(B, n) / n)
  
  potencia <- function(k) {
    a <- k * s1
    xbar <- rnorm(B, mean = a, sd = 1 / sqrt(n + 1))
    mean((sqrt(n + 1) * xbar / s2) > tcrit)
  }
  uniroot(function(k) potencia(k) - beta,
          c(0, 100), tol = 0.001)$root
}


#' Cálculo iterativo del número mínimo de réplicas usando el método de Harris–Hurvitz–Mood (HHM)
#'
#' Esta función estima de forma iterativa el número mínimo de réplicas necesarias \code{r}
#' en un diseño experimental con múltiples tratamientos, empleando el método de Harris–Hurvitz–Mood (HHM).
#' La estimación se basa en el ajuste simultáneo de \code{df2}, el valor de \code{K'} proveniente de la tabla
#' de potencia (calculado mediante \code{k_tabla_mc}), y la fórmula:
#'
#' \deqn{
#'   r = 2(df_2 + 1)\left(\frac{K' \cdot S_1}{d}\right)^2
#' }
#'
#' El proceso se repite hasta que el cambio relativo entre iteraciones sucesivas sea menor a un umbral
#' especificado por \code{max_error}, o hasta alcanzar el número máximo de iteraciones.
#'
#' @param t Número de tratamientos.
#' @param d Diferencia mínima detectable entre medias.
#' @param r_inicial Valor inicial para el número de réplicas.
#' @param S1_sq Estimación de la varianza experimental (\code{S1^2}).
#' @param df1 Grados de libertad asociados al estimador de varianza experimental.
#' @param alpha_ Nivel de significancia (por defecto \code{0.05}).
#' @param beta_ Error tipo II (\code{1 - potencia deseada}, por defecto \code{0.8}).
#' @param max_error Error relativo máximo permitido para la convergencia (por defecto \code{0.01}).
#' @param max_iter Número máximo de iteraciones permitidas (por defecto \code{1000}).
#'
#' @return Una lista con los siguientes elementos:
#' \describe{
#'   \item{r}{Número mínimo de réplicas requeridas (redondeado hacia arriba).}
#'   \item{S1}{Desviación estándar estimada a partir de \code{S1_sq}.}
#'   \item{df1}{Grados de libertad del estimador de varianza experimental.}
#'   \item{df2}{Grados de libertad del error al finalizar la iteración.}
#'   \item{dif}{Diferencia mínima detectable utilizada.}
#'   \item{alfa}{Nivel de significancia.}
#'   \item{iteraciones}{Número de iteraciones realizadas.}
#'   \item{convergencia}{\code{TRUE} si se alcanzó convergencia; \code{FALSE} si se alcanzó el límite de iteraciones.}
#'   \item{K_}{Valor de \(k = a / s_1\) que garantiza la potencia deseada, obtenido por simulación Monte Carlo mediante \code{k_tabla_mc()}.}
#' }
#'
#' @details
#' Esta función es útil en contextos donde el valor de \code{K'} depende del número de grados de libertad
#' \code{df2}, que a su vez depende de \code{r}. El proceso iterativo ajusta ambos hasta que el número
#' de réplicas requerido estabiliza su valor.
#'
#' @seealso \code{\link{k_tabla_mc}}, \code{\link{obtener_K_A9}}, \code{\link{calcular_df2}}
#'
#' @examples
#' resultados <- calcular_r_minimo_HHM(
#'   t = 6,
#'   d = 20,
#'   r_inicial = 3,
#'   S1_sq = 141.6,
#'   df1 = 40,
#'   alpha_ = 0.05,
#'   beta_ = 0.8
#' )
#'
#' resultados$r      # número estimado de réplicas
#' resultados$K_     # valor de K = a / s1 obtenido por simulación
#' resultados$df2    # grados de libertad del error final
#'
#' @export


calcular_r_minimo_HHM <- function(t, 
                                  d, 
                                  r_inicial, 
                                  S1_sq,
                                  df1, 
                                  alpha_ = 0.05, 
                                  beta_ = 0.8, 
                                  max_error = 0.01, 
                                  max_iter = 1000) {
  # Inicializar
  r_est <- r_inicial
  df2 <- t * (r_est - 1)
  error <- 1
  iter <- 0
  # Iteración hasta convergencia o límite de iteraciones
  while (error > max_error && iter < max_iter) {
    iter <- iter + 1
    # Tamaño de efecto estandarizado
    K <- k_tabla_mc(n = df2, m = df1, beta = beta_, alpha = alpha_)
    # Cálculo de r usando la fórmula 5.24 corregida (sin redondear)
    r_new <- 2*(df2 + 1)*((K*sqrt(S1_sq))/d)^2
    # Error relativo
    error <- abs(r_new - r_est) / r_est
    # Actualizar
    r_est <- r_new
    df2 <- t * (r_est - 1)
  }
  return(list(
    r = ceiling(r_est),
    S1 = sqrt(S1_sq),
    df1 = df1,
    df2 = floor(df2),
    dif = d,
    alfa = alpha_,
    iteraciones = iter,
    convergencia = error <= max_error,
    K_ = K 
  ))
}


