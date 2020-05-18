#' @keywords internal
"_PACKAGE"

## usethis namespace: start

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom pracma coth
#' @importFrom Rdpack reprompt
#' @import gsl
## usethis namespace: end
NULL





aquifer_predict_interface <- function(aquifer_lst, time_lst) {

   m3_to_bbl <- 6.289814
   m_to_ft <- 3.28084
   kPa_to_psi <- 0.145038
   month_to_day <- 30.41667
   year_to_day <- 365
   input_unit <- aquifer_lst$input_unit
   output_unit <- aquifer_lst$output_unit
   model <-  class(aquifer_lst)[1]
   mu_water <- aquifer_lst$mu_water
   c_water <- aquifer_lst$c_water
   c_rock <- aquifer_lst$c_rock
   pressure <- aquifer_lst$pressure
   time <- time_lst$t
   time_ref_date <- time_lst$reference_date
   time_unit <- class(time_lst)[1]
   if (time_unit == "month") {
      time_c <- time * month_to_day
   } else if (time_unit == "year") {
      time_c <- time * year_to_day
   } else {
      time_c <- time
   }
   if (length(time) != length(pressure)) {
      stop("pressure and time must have the same length.")
   }
   Date <- NULL



   if (model == "veh_rad_edge") {

      phi <- aquifer_lst$phi
      perm_h <- aquifer_lst$perm_h
      h_a <- aquifer_lst$h_a
      r_a <- aquifer_lst$r_a
      r_R <- aquifer_lst$r_R
      tetha <- aquifer_lst$tetha
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         r_a <- r_a * m_to_ft
         r_R <- r_R * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- veh_uss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time_c, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }




   if (model == "ykh_rad_bottom") {

      phi <- aquifer_lst$phi
      perm_h <- aquifer_lst$perm_h
      perm_v <- aquifer_lst$perm_v
      h_a <- aquifer_lst$h_a
      r_a <- aquifer_lst$r_a
      r_R <- aquifer_lst$r_R
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         r_a <- r_a * m_to_ft
         r_R <- r_R * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- yk_uss_rad_bottom(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time_c, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }



   if (model == "fetk_rad_edge") {

      phi <- aquifer_lst$phi
      perm_h <- aquifer_lst$perm_h
      perm_v <- aquifer_lst$perm_v
      h_a <- aquifer_lst$h_a
      r_a <- aquifer_lst$r_a
      r_R <- aquifer_lst$r_R
      tetha <- aquifer_lst$tetha
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         r_a <- r_a * m_to_ft
         r_R <- r_R * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- fetkovich_pss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time_c, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }



   if (model == "nb_lin_edge") {

      phi <- aquifer_lst$phi
      perm_h <- aquifer_lst$perm_h
      h_a <- aquifer_lst$h_a
      w_a <- aquifer_lst$w_a
      l_a <- aquifer_lst$l_a
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         w_a <- w_a * m_to_ft
         l_a <- l_a * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- nb_uss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }



   if (model == "nb_lin_bottom") {

      phi <- aquifer_lst$phi
      perm_v <- aquifer_lst$perm_v
      h_a <- aquifer_lst$h_a
      w_a <- aquifer_lst$w_a
      l_a <- aquifer_lst$l_a
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         w_a <- w_a * m_to_ft
         l_a <- l_a * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- nb_uss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }




   if (model == "fetk_lin_edge") {

      phi <- aquifer_lst$phi
      perm_h <- aquifer_lst$perm_h
      h_a <- aquifer_lst$h_a
      w_a <- aquifer_lst$w_a
      l_a <- aquifer_lst$l_a
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         w_a <- w_a * m_to_ft
         l_a <- l_a * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- fetk_pss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }




   if (model == "fetk_lin_bottom") {

      phi <- aquifer_lst$phi
      perm_v <- aquifer_lst$perm_v
      h_a <- aquifer_lst$h_a
      w_a <- aquifer_lst$w_a
      l_a <- aquifer_lst$l_a
      if (input_unit == "SI") {
         h_a <- h_a * m_to_ft
         w_a <- w_a * m_to_ft
         l_a <- l_a * m_to_ft
         c_water <- c_water / kPa_to_psi
         c_rock <- c_rock / kPa_to_psi
         pressure <- pressure * kPa_to_psi
      }
      We_table <- fetk_pss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure)
      if (time_unit == "day") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (time_unit == "month") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (months)", "We")
      }
      if (time_unit == "year") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (years)", "We")
      }
      if (time_unit == "date") {
         We_table$Date <- time_ref_date + (We_table$Time)
         We_table$Time <- time
         We_table <- We_table %>% dplyr::select(Date, dplyr::everything())
         colnames(We_table) <- c("Date", "Time (days)", "We")
      }
      if (output_unit == "SI") {
         We_table$We <- We_table$We * 1e6 / m3_to_bbl
         colnames(We_table)[3] <- "We (m3)"
      } else {
         colnames(We_table)[3] <- "We (MMbbl)"
      }
      return(We_table)
   }



}




#****************************************************************************************************
# Stehfest Algorithm

stehfest <- function(n) {
   # Coefficients for the Gaver-Stehfest algorithm
   v <- vector(length = n)
   m <- 0.5 * n
   for (i in 1:n) {
      sum <- 0
      j <- floor(0.5 * (i + 1))
      l <- min(i, m)
      for (k in j:l) {
         sum <- sum + (k ^ m * factorial(2 * k)) / (factorial(m - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i))
      }
      v[i] = (-1) ^ (i + m) * sum
   }
   return(v)
}


#****************************************************************************************************
# Van Everdingen and Hurst Unsteady State Model - Radial Flow - Edge Water - Laplace Transform

QDs_inf <- function(s) {

   a <- bessel_K0(sqrt(s)) / bessel_K1(sqrt(s)) / (s ^ 1.5)
   sol <- (1 / a) / (s * s * s)
   return(sol)
}

QD_inf <- function(tD, stehfest_size) {

   v <- stehfest(stehfest_size)
   sum <- 0
   for (i in 1:stehfest_size) {
      s = i * log(2) / tD
      sum <- sum + v[i] * QDs_inf(s)
   }
   sol <- sum * log(2) / tD
   return(sol)
}

QDs <- function(s, reD) {

   if (sqrt(s) > 742) {
      ratio <- 1.0
      f <- ratio / s / sqrt(s)
      g <- 1 / f / (s * s * s)
      return(g)
   } else {
      if ((reD * sqrt(s)) > 713) {
         if (sqrt(s) > 742) {
            ratio <- 1.0
            f <- ratio / s / sqrt(s)
            g <- 1 / f / (s * s * s)
            return(g)
         } else {
            ratio <- bessel_K0(sqrt(s)) / bessel_K1(sqrt(s))
            f <- ratio / s / sqrt(s)
            g <- 1 / f / (s * s * s)
            return(g)
         }
      } else {
         a <- bessel_K0(sqrt(s)) / bessel_I0(sqrt(s))
         b <- bessel_K1(sqrt(s)) / bessel_I0(sqrt(s))
         c <- bessel_I1(sqrt(s)) / bessel_I0(sqrt(s))
         d <- bessel_K1(reD * sqrt(s)) / bessel_I1(reD * sqrt(s))
         ratio <- (a + d) / (b - d * c)
         f <- ratio / s / sqrt(s)
         g <- 1 / f / (s * s * s)
         return(g)
      }
   }
}

QD <- function(tD, reD, stehfest_size) {

   v <- stehfest(stehfest_size)
   sum <- 0
   for (i in 1:stehfest_size) {
      s = i * log(2) / tD
      sum <- sum + v[i] * QDs(s, reD)
   }
   sol <- sum * log(2) / tD
   return(sol)
}

veh_uss_rad_edge_WeD <- function(tD, reD) {

   stehfest_size <- 8
   if (reD > 200) {
      vec <- QD_inf(tD, stehfest_size)
   } else {
      vec <- QD(tD, reD, stehfest_size)
   }
   return(vec)
}

veh_uss_rad_edge <- function(phi, perm, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      f <- tetha / 360
      B <- 1.119 * phi * c_total * r_R * r_R * h_a * f
      tD <- 6.328e-03 * perm * time / phi / mu_water / c_total / r_R / r_R
      raD <- r_a / r_R
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- veh_uss_rad_edge_WeD(tD_temp[j], raD)
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}


#*******************************************************************************

# Yildiz and Khosravi Unsteady State Model - Radial Flow - Bottom Water

QDs_yk <- function(s, raD, haD) {

   sai_1 <- coth(sqrt(s) * haD) / s ^ 1.5
   k <- 100
   zeros <- bessel_zero_J1(s = c(1:k))
   beta <- zeros / raD
   sai_2 <- vector(length = k)
   sum_sai_2 <- 0
   for (i in 1:k) {
      sai_2[i] <- (bessel_J1(beta[i]) ^ 2 *
                      coth(sqrt(beta[i] * beta[i] + s) * haD)) /
         (beta[i] * beta[i] * sqrt(beta[i] * beta[i] + s) *
             bessel_J0(beta[i] * raD) ^ 2)
      sum_sai_2 <- sum_sai_2 + sai_2[i]
   }
   pDs <- sai_1 / raD / raD + 4 * sum_sai_2 / raD / raD / s
   QDs <- 1 / pDs / s / s / s / 2 / haD
   return(QDs)
}

QD_yk <- function(tD, raD, haD, stehfest_size) {

   v <- stehfest(stehfest_size)
   sum <- 0
   for (i in 1:stehfest_size) {
      s <- i * log(2) / tD
      sum <- sum + v[i] * QDs_yk(s, raD, haD)
   }
   sol <- sum * log(2) / tD
   return(sol)
}

yk_uss_rad_bottom_WeD <- function(tD, raD, haD) {

   stehfest_size <- 8
   vec <- QD_yk(tD, raD, haD, stehfest_size)
   return(vec)
}


# Yildiz and Khosravi Unsteady State Model - Radial Flow - Bottom Water

yk_uss_rad_bottom <- function(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure) {

   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      WeD <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      Fk <- perm_v / perm_h
      B <- 1.119 * phi * c_total * r_R * r_R * h_a
      tD <- 6.328e-03 * perm_h * time / phi / mu_water / c_total / r_R / r_R
      raD <- r_a / r_R
      zD <- h_a / r_R / sqrt(Fk)
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- yk_uss_rad_bottom_WeD(tD_temp[j], raD, zD)
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}


#*******************************************************************************

# Fetkovich Pseudo-Steady State Model - Radial Flow - Edge Water

fetkovich_pss_rad_edge <- function(phi, perm, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      dt <- vector(length = len_t)
      p_r <- vector(length = len_t)
      p_a <- vector(length = len_t)
      dp_ar <- vector(length = len_t)
      dWe <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      f <- tetha / 360
      reD <- r_a / r_R
      Wi <- pi * (r_a * r_a - r_R * r_R) * h_a * phi / 5.615
      Wei <- c_total * Wi * pressure[1] * f
      for (i in seq_along(time)) {
         if (i == 1) {
            dt[i] <- 0
            p_r[i] <- pressure[1]
            p_a[i] <- pressure[1]
            dp_ar[i] <- 0
            dWe[i] <- 0
            We[i] <- 0
         } else {
            J <- 0.00708 * perm * h_a * f / mu_water / (log(reD) - 0.75)
            dt[i] <- time[i] - time[i - 1]
            p_r[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dp_ar[i] <- p_a[i - 1] - p_r[i]
            dWe[i] <- Wei * dp_ar[i] * (1 - exp(-J * pressure[1] * dt[i] / Wei)) / pressure[1]
            We[i] <- We[i - 1] + dWe[i]
            p_a[i] <- pressure[1] * (1 - We[i] / Wei)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}




#*******************************************************************************

# Nabor and Barham Unsteady State Model - Linear Flow - Edge Water - Exact Solution

nb_uss_lin_edge_WeD <- function(tD) {

   if (tD >= 10) {
      WeD <- 1
      return(WeD)
   } else {
      sum <- 0
      for (i in seq(1, 101, by = 2)) {
         sum <- sum + (1 / i / i) * exp(-1 * i * i * pi * pi * tD / 4)
      }
      WeD <- 1 - 8 * sum / pi / pi
      return(WeD)
   }
}


# Nabor and Barham Unsteady State Model - Linear Flow - Edge Water

nb_uss_lin_edge <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      WeD <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      V_a <- h_a * L_a * w_a * phi
      B <- V_a * c_total / 5.615
      tD <- 6.328e-03 * perm * time / phi / mu_water / c_total / L_a / L_a
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- nb_uss_lin_edge_WeD(tD_temp[j])
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}




#*******************************************************************************
# Nabor and Barham Unsteady State Model - Linear Flow - Bottom Water - Exact Solution

nb_uss_lin_bottom_WeD <- function(tD) {

   if (tD >= 10) {
      WeD <- 1
      return(WeD)
   } else {
      sum <- 0
      for (i in seq(1, 101, by = 2)) {
         sum <- sum + (1 / i / i) * exp(-1 * i * i * pi * pi * tD / 4)
      }
      WeD <- 1 - 8 * sum / pi / pi
      return(WeD)
   }
}


# Nabor and Barham Unsteady State Model - Linear Flow - Bottom Water

nb_uss_lin_bottom <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      WeD <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      V_a <- h_a * L_a * w_a * phi
      B <- V_a * c_total / 5.615
      tD <- 6.328e-03 * perm * time / phi / mu_water / c_total / h_a / h_a
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- nb_uss_lin_bottom_WeD(tD_temp[j])
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}




#*******************************************************************************

# Fetkovich Pseudo-Steady State Model - Linear Flow - Edge Water

fetk_pss_lin_edge <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      dt <- vector(length = len_t)
      p_r <- vector(length = len_t)
      p_a <- vector(length = len_t)
      dp_ar <- vector(length = len_t)
      dWe <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      Wi <- h_a * L_a * w_a * phi / 5.615
      Wei <- c_total * Wi * pressure[1]
      for (i in seq_along(time)) {
         if (i == 1) {
            dt[i] <- 0
            p_r[i] <- pressure[1]
            p_a[i] <- pressure[1]
            dp_ar[i] <- 0
            dWe[i] <- 0
            We[i] <- 0
         } else {
            J <- 0.003381 * perm * h_a * w_a / mu_water / L_a
            dt[i] <- time[i] - time[i - 1]
            p_r[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dp_ar[i] <- p_a[i - 1] - p_r[i]
            dWe[i] <- Wei * dp_ar[i] * (1 - exp(-J * pressure[1] * dt[i] / Wei)) / pressure[1]
            We[i] <- We[i - 1] + dWe[i]
            p_a[i] <- pressure[1] * (1 - We[i] / Wei)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}




#*******************************************************************************

# Fetkovich Pseudo-Steady State Model - Linear Flow - Bottom Water

fetk_pss_lin_bottom <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      dt <- vector(length = len_t)
      p_r <- vector(length = len_t)
      p_a <- vector(length = len_t)
      dp_ar <- vector(length = len_t)
      dWe <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      Wi <- h_a * L_a * w_a * phi / 5.615
      Wei <- c_total * Wi * pressure[1]
      for (i in seq_along(time)) {
         if (i == 1) {
            dt[i] <- 0
            p_r[i] <- pressure[1]
            p_a[i] <- pressure[1]
            dp_ar[i] <- 0
            dWe[i] <- 0
            We[i] <- 0
         } else {
            J <- 0.003381 * perm * w_a * L_a / mu_water / h_a
            dt[i] <- time[i] - time[i - 1]
            p_r[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dp_ar[i] <- p_a[i - 1] - p_r[i]
            dWe[i] <- Wei * dp_ar[i] * (1 - exp(-J * pressure[1] * dt[i] / Wei)) / pressure[1]
            We[i] <- We[i - 1] + dWe[i]
            p_a[i] <- pressure[1] * (1 - We[i] / Wei)
         }
      }
      We <- We / 1e6
      df <- data.frame(Time = time, We = We)
      return(df)
   }
}