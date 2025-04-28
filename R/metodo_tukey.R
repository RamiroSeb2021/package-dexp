#' @import stats
#' @importFrom stats qf qtukey ptukey qchisq qt
NULL

#' Estimación de la desviación estándar corregida y grados de libertad
#'
#' Esta función estima una desviación estándar corregida (\code{S1}) y el número de grados de libertad
#' (\code{df1}) que reproduce, mediante cocientes de chi-cuadrado, la relación entre un límite superior
#' y un límite inferior dados. Se usa un criterio de error relativo máximo permitido para la aceptación del resultado.
#'
#' @param desviacion_estandar Desviación estándar inicial (en unidades de la variable estudiada).
#' @param Si Porcentaje relativo inferior aceptado (por ejemplo, 0.07 corresponde a un 7\%).
#' @param Ss Porcentaje relativo superior aceptado (por ejemplo, 0.12 corresponde a un 12\%).
#' @param max_error Error relativo máximo permitido entre el cociente observado y el cociente esperado (por defecto \code{0.01}).
#' @param confianza Nivel de confianza del intervalo considerado (por defecto \code{0.9}). *Actualmente fijo en los cálculos internos.*
#' @param maximum_df Máximo número de grados de libertad que se evaluarán (por defecto \code{1000}).
#'
#' @return
#' Una lista con:
#' \describe{
#'   \item{S1}{Estimación corregida de la desviación estándar.}
#'   \item{grados_libertad}{Número de grados de libertad estimado.}
#'   \item{valor_x}{Cociente entre cuantiles de chi-cuadrado.}
#'   \item{error_relativo}{Error relativo del cociente respecto al objetivo.}
#' }

#'
#' @details
#' Se busca el número de grados de libertad tal que el cociente de los valores críticos de la distribución chi-cuadrado,
#' evaluados a \code{p = confianza} y \code{p = 1 - confianza}, sea cercano al cociente observado entre \code{Ss} y \code{Si}.
#' El proceso finaliza tan pronto como el error relativo esté por debajo del máximo permitido (\code{max_error}).
#'
#' @examples
#' set.seed(123)
#' calcular_S1_df1(desviacion_estandar = 30, Si = 0.07, Ss = 0.12)
#'
#' @export
calcular_S1_df1 <- function(desviacion_estandar,
                            Si,
                            Ss,
                            max_error = 0.01,
                            confianza = 0.9,
                            maximum_df = 1000) {

  SI <- desviacion_estandar * Si
  SS <- desviacion_estandar * Ss

  S1 <- (SI + SS) / 2

  Cociente <- SS / SI

  for (i in 1:maximum_df) {
    x_i <- sqrt(qchisq(p = confianza, df = i) / qchisq(p = 1 - confianza, df = i))
    error_relativo <- abs(x_i - Cociente) / Cociente

    if (error_relativo < max_error) {
      return(list(
        S1 = S1,
        grados_libertad = i,
        valor_x = x_i,
        error_relativo = error_relativo
      ))
    }
  }

  warning("No se encontró un número de grados de libertad que cumpla el criterio.")
  return(NULL)
}



# Calculo de A ------------------------------------------------------------

#' Cálculo del valor A en intervalos de confianza para comparación de medias
#'
#' Esta función calcula el valor \code{A}, que representa una estimación ajustada
#' de la mitad de la longitud esperada del intervalo de confianza para comparar dos medias
#' en un diseño experimental. Es útil para evaluar si el número de réplicas \code{r}
#' genera un intervalo suficientemente estrecho, es decir, si \code{2A ≈ D}.
#'
#' @param alfa Nivel de significancia deseado.
#' @param r Número de réplicas por tratamiento.
#' @param sigma Desviación estándar estimada de la variable de interés.
#'
#' @return Un valor numérico que representa \code{A}, la mitad de la longitud esperada
#' del intervalo de confianza.
#'
#' @examples
#' calcular_A(alfa = 0.05, r = 6, sigma = sqrt(141.6))
#'
#' @export

calcular_A <- function(alfa, r, sigma){
  # alfa: Nivel de significancia
  # t: Numero de tratamientos
  # sigma: Desviacion estandar

  numerador <- qt(p = 1 - alfa, df = r - 1) * 2 * sigma * gamma(r/2)
  denominador <- sqrt(r)*sqrt((r-1)) * gamma((r-1)/2)

  A <- numerador/denominador
  return(A)
}

# Calcular_r_basico -------------------------------------------------------

#' Cálculo del número de réplicas según fórmula básica corregida
#'
#' Esta función estima el número requerido de réplicas \code{r} en un diseño experimental
#' con múltiples tratamientos, usando una versión corregida de la fórmula 5.24. Se basa en
#' cuantiles de la distribución F (para potencia estadística) y de la distribución de Tukey
#' (para comparaciones múltiples), así como en la varianza experimental estimada.
#'
#' @param S1_sq Estimación de la varianza experimental (\code{S1^2}).
#' @param df1 Grados de libertad del estimador de varianza experimental.
#' @param df2 Grados de libertad del error.
#' @param beta Nivel de error tipo II (es decir, \code{1 - potencia deseada}).
#' @param alpha Nivel de significancia para las comparaciones múltiples.
#' @param t Número de tratamientos.
#' @param d Diferencia mínima detectable entre medias.
#'
#' @return Un valor numérico que representa el número estimado de réplicas \code{r} (no redondeado).
#'
#' @examples
#' # Ejemplo directo con valores definidos
#' calcular_r_B(
#'   S1_sq = 141.6,
#'   df1 = 40,
#'   df2 = 48,
#'   beta = 0.1,
#'   alpha = 0.1,
#'   t = 6,
#'   d = 20
#' )
#'
#' # Ejemplo extendido (como prueba)
#' S1_sq_ <- 141.6
#' Df1 <- 40
#' D <- 20
#' Df2 <- 48
#' Beta <- 0.1
#' Alpha <- 0.1
#' T_ <- 6
#'
#' qf(1 - 0.10, df1 = 48, df2 = 40)
#' qtukey(p = 0.90, nmeans = 5, df = 48)
#' ptukey(q = 4.2, nmeans = 5, df = 48)
#'
#' calcular_r_B(
#'   S1_sq = S1_sq_,
#'   df1 = Df1,
#'   d = D,
#'   df2 = Df2,
#'   beta = Beta,
#'   alpha = Alpha,
#'   t = T_
#' )
#'
#' @export


calcular_r_B <- function(S1_sq, df1, df2, beta, alpha, t, d){
  # Cuantil F para potencia (1 - beta)
  F_crit <- qf(1 - beta, df2, df1)
  # Cuantil de Tukey para comparaciones múltiples (nivel de confianza 1 - alpha)
  q_crit <- qtukey(1 - alpha, nmeans = t - 1, df = df2)
  # Cálculo de r usando la fórmula 5.24 corregida (sin redondear)
  r <- (F_crit * S1_sq * q_crit^2) / d^2
  return(r)
}

# Calcular df2 ------------------------------------------------------------

#' Calcular número de réplicas y grados de libertad del error
#'
#' Esta función estima de forma iterativa el número de réplicas \code{r} y
#' los grados de libertad del error \code{df2} para un diseño con comparaciones
#' múltiples entre tratamientos. El algoritmo se basa en cuantiles de la
#' distribución F y la distribución de Tukey, y se detiene cuando el error
#' relativo entre iteraciones es menor que un umbral especificado.
#'
#' @param t Número de tratamientos.
#' @param d Diferencia mínima detectable entre medias.
#' @param r_inicial Estimación inicial para el número de réplicas.
#' @param df1 Grados de libertad asociados al estimador de la varianza experimental.
#' @param alpha Nivel de significancia (error tipo I).
#' @param beta Nivel de error tipo II (1 - potencia deseada).
#' @param S1_sq Estimación de la varianza experimental (S1 al cuadrado).
#' @param max_error Error relativo máximo permitido para la convergencia (por defecto 0.01).
#' @param max_iter Número máximo de iteraciones permitidas (por defecto 1000).
#'
#' @return Una lista con los siguientes elementos:
#' \describe{
#'   \item{r}{Número estimado de réplicas (redondeado hacia arriba).}
#'   \item{S1}{Varianza experimental utilizada (S1^2).}
#'   \item{df1}{Grados de libertad del estimador de la varianza.}
#'   \item{df2}{Grados de libertad del error (redondeado hacia abajo).}
#'   \item{dif}{Diferencia mínima detectable (d).}
#'   \item{alfa}{Nivel de significancia.}
#'   \item{iteraciones}{Número de iteraciones realizadas.}
#'   \item{convergencia}{\code{TRUE} si el proceso convergió; \code{FALSE} en caso contrario.}
#' }
#'
#' @examples
#'
#' calcular_df2(
#'   t = 6, d = 20, r_inicial = 5, df1 = 40,
#'   alpha = 0.05, beta = 0.1, S1_sq = 141.6
#' )
#'
#' @export

calcular_df2 <- function(t, d, r_inicial, df1, alpha, beta, S1_sq, max_error = 0.01, max_iter = 1000) {
  # Inicializar
  r_est <- r_inicial
  df2 <- t * (r_est - 1)
  error <- 1
  iter <- 0
  # Iteración hasta convergencia o límite de iteraciones
  while (error > max_error && iter < max_iter) {
    iter <- iter + 1
    # Cuantil F para potencia (1 - beta)
    F_crit <- qf(1 - beta, df2, df1)
    # Cuantil de Tukey para comparaciones múltiples (nivel de confianza 1 - alpha)
    q_crit <- qtukey(1 - alpha, nmeans = t, df = df2)
    # Cálculo de r usando la fórmula 5.24 corregida (sin redondear)
    r_new <- (F_crit * S1_sq * q_crit^2) / d^2
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
    alfa = alpha,
    iteraciones = iter,
    convergencia = error <= max_error
  ))
}

# Revisión intervalo ------------------------------------------------------


#' Revisión de intervalo para número óptimo de réplicas
#'
#' Esta función evalúa un rango de valores alrededor del número estimado de
#' réplicas \code{r}, calculando el valor correspondiente de \code{A} para cada uno
#' y determinando cuál valor de \code{r} hace que la longitud del intervalo
#' de confianza (\code{2A}) sea más cercana a la diferencia mínima detectable \code{d}.
#'
#' @param info_r Lista resultante de la función \code{calcular_df2}, que debe contener
#' los elementos \code{r}, \code{dif}, \code{alfa}, y \code{S1}.
#' @param des Desviación entera para explorar alrededor del valor de \code{r}. Por defecto es 5.
#'
#' @return Una lista con los siguientes elementos:
#' \describe{
#'   \item{r_}{Valor de \code{r} que minimiza la diferencia entre \code{2A} y \code{d}.}
#'   \item{A}{Valor de \code{A} correspondiente a \code{r_}.}
#'   \item{r_l}{Vector de valores de \code{r} evaluados.}
#'   \item{A_l}{Vector de valores de \code{A} correspondientes a cada \code{r_l}.}
#'   \item{position}{Índice de \code{r_} dentro del vector evaluado.}
#' }
#'
#' @examples
#'
#' # Primero es necesario ejecutar calcular_r_MT
#' info <- calcular_r_MT(
#'   T_ = 6, D = 20, ro = 5,
#'   S1 = sqrt(141.6), df1 = 40,
#'   alfa = 0.05, Beta = 0.1
#' )
#' revision_intervalo(info_r = info, des = 3)

#'
#' @export

revision_intervalo <- function(info_r, des = 5){
  x <- seq(from = info_r$r - des, to = info_r$r + des)
  d_ <- rep(info_r$dif, length(x))
  x_ <- calcular_A(info_r$alfa, x, info_r$S1)
  diff <- abs(2*x_ - d_)

  i = which(diff == min(diff))
  return(list(r_ = x[i], A = x_[i], r_l = x, A_l = x_, position = i))
}

# Calcular el número de replicas ------------------------------------------

#' Cálculo del número óptimo de réplicas en diseños con múltiples tratamientos
#'
#' Esta función estima el número óptimo de réplicas necesarias \code{r} en un diseño experimental
#' con múltiples tratamientos, de manera que la longitud del intervalo de confianza sea cercana
#' a un valor deseado (\code{D = 2A}). Utiliza un procedimiento iterativo basado en distribuciones
#' F y de Tukey para ajustar el número de réplicas y grados de libertad del error, y realiza una
#' revisión posterior para seleccionar el valor de \code{r} que produce un \code{A} más cercano a \code{D/2}.
#'
#' @param T_ Número de tratamientos.
#' @param D Diferencia mínima detectable entre medias (objetivo de longitud del intervalo: \code{D = 2A}).
#' @param ro Estimación inicial del número de réplicas.
#' @param S1 Estimación de la desviación estándar experimental.
#' @param df1 Grados de libertad del estimador de \code{S1}.
#' @param alfa Nivel de significancia (por defecto \code{0.05}).
#' @param Beta Error tipo II deseado (por defecto \code{0.1}).
#'
#' @return Una lista con los siguientes elementos:
#' \describe{
#'   \item{r}{Número estimado de réplicas que cumple la potencia deseada.}
#'   \item{A}{Valor de \code{A} asociado a \code{r}.}
#'   \item{r_i}{Valor de \code{r} más cercano que logra que \code{2A ≈ D}.}
#'   \item{A_i}{Valor de \code{A} correspondiente a \code{r_i}.}
#'   \item{lista_r}{Vector de valores de \code{r} evaluados en la revisión del intervalo.}
#'   \item{valores_A}{Vector de valores de \code{A} asociados a \code{lista_r}.}
#'   \item{posicion_escogida}{Índice de \code{r_i} dentro de \code{lista_r}.}
#' }
#'
#' @examples
#' # Ejemplo básico
#' calcular_r_MT(
#'   T_ = 6, D = 20, ro = 5,
#'   S1 = sqrt(141.6), df1 = 40,
#'   alfa = 0.05, Beta = 0.1
#' )
#'
#' # Prueba ejemplo 5.12
#' calcular_r_MT(
#'   T_ = 6,
#'   D = 20,
#'   ro = 6,
#'   S1 = sqrt(141.6),
#'   df1 = 40
#' )
#'
#' # Prueba ejemplo 5.13
#' calcular_r_MT(
#'   T_ = 8,
#'   D = 500,
#'   ro = 4,
#'   S1 = sqrt(90000),
#'   df1 = 200,
#'   alfa = 0.1,
#'   Beta = 0.25
#' )
#'
#' @export


calcular_r_MT <- function(T_, D, ro, S1, df1, alfa = 0.05, Beta = 0.1){

  info <-
    calcular_df2(
      t = T_,
      d = D,
      r_inicial = ro,
      df1 = df1,
      alpha = alfa,
      beta = Beta,
      S1_sq = S1^2
    )
  # resultados de la iteracion
  info2 <- revision_intervalo(info)
  return(list(r = info$r,
              A = calcular_A(alfa, info$r, S1),
              r_i = info2$r_,
              A_i = info2$A,
              lista_r = info2$r_l,
              valores_A = info2$A_l,
              posicion_escogida =  info2$position))
}







