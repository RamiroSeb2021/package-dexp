#' Asignación proporcional de réplicas sin considerar costos ni tamaño de muestra
#'
#' Esta función calcula un esquema de asignación proporcional de réplicas por tratamiento,
#' basado únicamente en las desviaciones estándar estimadas de cada tratamiento.
#' El objetivo es redistribuir un número total fijo de observaciones iniciales (\code{n = r0 \times a}),
#' sin considerar costos individuales ni un tamaño de muestra objetivo distinto al inicial.
#'
#' @param a Número de tratamientos.
#' @param r0 Número inicial de réplicas por tratamiento antes de ajustar proporcionalmente.
#' @param sigmas Vector numérico de desviaciones estándar estimadas para cada tratamiento. Su longitud debe ser igual a \code{a}.
#'
#' @return
#' Un vector numérico de longitud \code{a} indicando el número de réplicas ajustadas proporcionalmente para cada tratamiento,
#' redondeadas al número entero más cercano.
#'
#' @details
#' El procedimiento de asignación es:
#' \deqn{n = r_0 \times a}
#' \deqn{r_i = \mathrm{round}\left(\frac{n \times \sigma_i}{\sum_{i=1}^a \sigma_i}\right)}
#' donde \eqn{\sigma_i} es la desviación estándar del tratamiento \eqn{i}.
#'
#' Si la longitud de \code{sigmas} no coincide con \code{a}, la función devuelve un error.
#'
#' @examples
#' sigmas <- c(6.27, 9.57, 12, 3.32)
#' proporcionalidad_sin_costo_ni_tamaño_de_muestra(a = 4, r0 = 5, sigmas = sigmas)
#'
#' @export
proporcionalidad_sin_costo_ni_tamaño_de_muestra <- function(a, r0, sigmas) {
  if (length(sigmas) != a) {
    stop("La cantidad de tratamientos no coincide con la longitud del vector de desviaciones")
  }
  n <- r0 * a
  r_prop <- round((n * sigmas) / sum(sigmas))
  return(r_prop)
}


#' Asignación proporcional de réplicas considerando costos y presupuesto total
#'
#' Esta función calcula un esquema de asignación proporcional de réplicas por tratamiento,
#' considerando tanto las desviaciones estándar como los costos unitarios de observación
#' de cada tratamiento. La asignación se ajusta para respetar un presupuesto total disponible (\code{costo_total}).
#'
#' @param a Número de tratamientos.
#' @param sigmas Vector numérico de desviaciones estándar estimadas para cada tratamiento (longitud igual a \code{a}).
#' @param costos Vector numérico de costos unitarios por tratamiento (longitud igual a \code{a}).
#' @param costo_total Presupuesto total disponible para realizar las observaciones.
#'
#' @return
#' Un vector numérico de longitud \code{a} indicando el número de réplicas ajustadas proporcionalmente para cada tratamiento,
#' redondeadas al número entero más cercano.
#'
#' @details
#' El procedimiento de asignación sigue los siguientes pasos:
#' \itemize{
#'   \item Se define \eqn{\lambda = 1/a}.
#'   \item Se calcula el factor \eqn{\phi}:
#'     \deqn{\phi = \left( \sum_{i=1}^a \lambda \sigma_i \sqrt{c_i} \right)^2 / (costo\_total)^2}
#'   \item Finalmente, el número de réplicas para cada tratamiento se estima como:
#'     \deqn{r_i = \mathrm{round}\left(\frac{\lambda \sigma_i}{\sqrt{\phi c_i}}\right)}
#' }
#'
#' Si la longitud de \code{sigmas} o \code{costos} no coincide con \code{a}, se detiene la ejecución con un mensaje de error.
#'
#' @examples
#' sigmas <- c(6.27, 9.57, 12, 3.32)
#' costos <- c(1000, 200, 700, 1100)
#' costo_total <- 50000
#' proporcionalidad_con_costo_ni_tamaño_de_muestra(a = 4, sigmas = sigmas, costos = costos, costo_total = costo_total)
#'
#' @export
proporcionalidad_con_costo_ni_tamaño_de_muestra <- function(a, sigmas, costos, costo_total) {
  if (length(sigmas) != a) {
    stop("La cantidad de tratamientos no coincide con la longitud del vector de desviaciones")
  }
  if (length(costos) != a) {
    stop("La cantidad de tratamientos no coincide con la longitud del vector de costos")
  }

  lambda <- 1 / a
  phi <- (sum(lambda * sigmas * sqrt(costos))^2) / (costo_total^2)
  r_prop <- round(lambda * sigmas / sqrt(phi * costos))

  return(r_prop)
}

#' Cálculo del número de tratamientos y réplicas bajo un modelo de efectos aleatorios
#'
#' Esta función estima el número óptimo de tratamientos y el número de réplicas por tratamiento
#' en un diseño experimental con modelo de efectos aleatorios. Se consideran los costos de tratamientos y unidades experimentales,
#' la varianza de error, la proporción de varianza atribuible a los efectos aleatorios (\code{rho}),
#' y un valor máximo permitido para la varianza relativa de las medias (\code{v_max}).
#'
#' @param costo_tratamiento Costo unitario de cada tratamiento (\eqn{C_1}).
#' @param costo_ue Costo unitario por cada unidad experimental (\eqn{C_2}).
#' @param sigma_cuadrado Varianza de los errores (\eqn{\sigma^2}).
#' @param rho Proporción de la varianza total atribuida a los efectos aleatorios (\eqn{\rho = \sigma_\tau^2 / \sigma^2}).
#' @param v_max Máximo valor permitido para la varianza relativa (\eqn{v_{\max}}).
#'
#' @return
#' Una lista con dos elementos:
#' \describe{
#'   \item{\code{num_de_tratamientos}}{Número estimado de tratamientos (entero redondeado).}
#'   \item{\code{num_de_replicas}}{Número estimado de réplicas por tratamiento (entero redondeado).}
#' }
#'
#' @details
#' El procedimiento utiliza las siguientes fórmulas:
#' \itemize{
#'   \item Se estima la varianza de efectos aleatorios:
#'     \deqn{\sigma^2_A = \rho \sigma^2}
#'   \item Número de tratamientos:
#'     \deqn{a = \left(\frac{1}{v_{\max}}\right) \left(\sigma^2_A + \sqrt{\frac{\sigma^2_A \sigma^2 C_2}{C_1}}\right)}
#'   \item Número de réplicas por tratamiento:
#'     \deqn{r = \sqrt{\frac{\sigma^2 C_1}{\sigma^2_A C_2}}}
#' }
#'
#' @examples
#' resultado <- numero_de_tratamientos_y_replicas_con_efectos_aleatorios(
#'   costo_tratamiento = 150000,
#'   costo_ue = 50000,
#'   sigma_cuadrado = 416.21,
#'   rho = 0.3796,
#'   v_max = 43.49
#' )
#' resultado
#'
#' @export
numero_de_tratamientos_y_replicas_con_efectos_aleatorios <- function(costo_tratamiento, costo_ue, sigma_cuadrado, rho, v_max) {
  sigma_A2 <- rho * sigma_cuadrado

  num_de_tratamientos <- round((1 / v_max) * (sigma_A2 + sqrt((sigma_A2 * sigma_cuadrado * costo_ue) / costo_tratamiento)))
  num_de_replicas <- round(sqrt((sigma_cuadrado * costo_tratamiento) / (sigma_A2 * costo_ue)))

  return(list(num_de_tratamientos = num_de_tratamientos, num_de_replicas = num_de_replicas))
}
