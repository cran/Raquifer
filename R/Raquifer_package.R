
#' @keywords internal
"_PACKAGE"

## usethis namespace: start

#' @import magrittr
#' @importFrom dplyr select
#' @importFrom pracma coth
#' @importFrom Rdpack reprompt
#' @import gsl
## usethis namespace: end
NULL



#' A list object for aquifer prameters
#'
#' Create an object of class 'aquifer'
#'
#' @param input_unit a unit system for parameters, a character string either 'SI' or 'Field'
#' @param output_unit a unit system for properties, a character string either 'SI' or 'Field'
#' @param model state of flow in the aquifer, a character string either 'uss' for the un-steady state flow or 'pss' for the pseudo-steady state flow
#' @param flow_type a character string either 'radial' or 'linear'
#' @param water_drive a character string either 'edge' or 'bottom'
#' @param phi aquifer porosity, a numeric fraction
#' @param perm_h aquifer horizontal permeability in 'md' in both 'SI' and 'Field' input unit systems. A NULL value must be used for the combination of 'uss', 'linear', and 'bottom' flow
#' @param perm_v aquifer vertical permeability in 'md' in both 'SI' and 'Field' input unit systems. A NULL value must be used for the combination of 'uss', 'linear', 'edge' flow. A NULL value must be used for the combination of 'uss', 'radial', 'edge' flow. A NULL value must be used for the combination of 'pss', 'radial', 'edge' flow.
#' @param h_a aquifer height in 'm' or 'ft' in 'SI' and 'Field' input unit systems, respectively.
#' @param r_a aquifer radius in 'm' or 'ft' in 'SI' and 'Field' input unit systems, respectively. A NULL value must be used for the combination of 'uss', 'linear', 'edge' flow. A NULL value must be used for the combination of 'uss', 'linear', 'bottom' flow.
#' @param r_R reservoir radius in 'm' or 'ft' in 'SI' and 'Field' input unit systems, respectively. A NULL value must be used for the combination of 'uss', 'linear', 'edge' flow. A NULL value must be used for the combination of 'uss', 'linear', 'bottom' flow.
#' @param w_a aquifer width in 'm' or 'ft' in 'SI' and 'Field' input unit systems, respectively. A NULL value must be used for the combination of 'uss', 'radial', 'edge' flow. A NULL value must be used for the combination of 'uss', 'radial', 'bottom' flow. A NULL value must be used for the combination of 'pss', 'radial', 'edge' flow.
#' @param l_a aquifer length in 'm' or 'ft' in 'SI' and 'Field' input unit systems, respectively. A NULL value must be used for the combination of 'uss', 'radial', 'edge' flow. A NULL value must be used for the combination of 'uss', 'radial', 'bottom' flow. A NULL value must be used for the combination of 'pss', 'radial', 'edge' flow.
#' @param tetha fraction of reservoir encircled by the aquifer, reported in "degrees" in both 'SI' and 'Field' input unit systems. A NULL value must be used for the combination of 'uss', 'radial', 'bottom' flow. A NULL value must be used for the combination of 'uss', 'linear', 'edge' flow. A NULL value must be used for the combination of 'uss', 'linear', 'bottom' flow.
#' @param mu_water water viscosity in 'mPa.s' or 'cp' in 'SI' and 'Field' input unit systems, respectively
#' @param c_water water compressibility in '1/kPa' or '1/psi' in 'SI' and 'Field' input unit systems, respectively
#' @param c_rock rock compressibility in '1/kPa' or '1/psi' in 'SI' and 'Field' input unit systems, respectively
#' @param pressure a numeric vector of pressure data at the boundary of reservoir/aquifer. Must have the same length as the 'aquifer_time()' object
#'
#' @return a list of class 'aquifer' with all the required parameters for the aquifer_predict() S3 methods
#' @export
#'
#' @examples
#'
#' aquifer_param_01 <- aquifer_param(input_unit = "Field", output_unit = "Field",
#' model = "uss", flow_type = "radial", water_drive = "edge", phi = 0.2, perm_h = 100,
#' h_a = 47, r_a = 2e4, r_R = 2e3, tetha = 360, mu_water = 0.34, c_water = 4e-6,
#' c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312))
#'
#' aquifer_param_01
#'
#' aquifer_param_02 <- aquifer_param(input_unit = "SI", output_unit = "SI",
#' model = "uss", flow_type = "radial", water_drive = "bottom", phi = 0.2, perm_h = 100,
#' perm_v = 25, h_a = 25, r_a = 6000, r_R = 600, mu_water = 0.34, c_water = 6e-7,
#' c_rock = 4.5e-7, pressure = c(3456, 3425, 3387, 3350, 3312) * 6.895)
#'
#' aquifer_param_02
#'
#' aquifer_param_03 <- aquifer_param(input_unit = "Field", output_unit = "Field",
#' model = "pss", flow_type = "radial", water_drive = "edge", phi = 0.2, perm_h = 100,
#' h_a = 47, r_a = 2e4, r_R = 2e3, tetha = 360, mu_water = 0.34, c_water = 4e-6,
#' c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312))
#'
#' aquifer_param_03
#'
#' aquifer_param_04 <- aquifer_param(input_unit = "Field", output_unit = "Field",
#' model = "uss", flow_type = "linear", water_drive = "edge", phi = 0.2, perm_h = 100,
#' h_a = 47, w_a = 30000, l_a = 10000, mu_water = 0.34, c_water = 4e-6,
#' c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312))
#'
#' aquifer_param_04
#'
#' aquifer_param_05 <- aquifer_param(input_unit = "Field", output_unit = "Field",
#' model = "uss", flow_type = "linear", water_drive = "bottom", phi = 0.2, perm_v = 10,
#' h_a = 47, w_a = 4000, l_a = 4000, mu_water = 0.34, c_water = 4e-6,
#' c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312))
#'
#' aquifer_param_05

aquifer_param <- function(input_unit = NULL, output_unit = NULL, model = NULL, flow_type = NULL, water_drive = NULL, phi = NULL, perm_h = NULL, perm_v = NULL, h_a = NULL, r_a = NULL, r_R = NULL, w_a = NULL, l_a = NULL, tetha = NULL, mu_water = NULL, c_water = NULL, c_rock = NULL, pressure = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be a character string.")
   if (!is.character(output_unit)) stop("'output_unit' must be a character string.")
   if (!is.character(model)) stop("'model' must be a character string.")
   if (!(input_unit %in% c("SI","Field"))) stop("'input_unit' must be either 'SI' or 'Field'.")
   if (!(output_unit %in% c("SI","Field"))) stop("'output_unit' must be either 'SI' or 'Field'.")
   if (!(model %in% c("uss","pss"))) stop("The aquifer model must be either 'uss' or 'pss'.")
   if (!(flow_type %in% c("radial","linear"))) stop("'flow_type' must be either 'radial' or 'linear'.")
   if (!(water_drive %in% c("edge","bottom"))) stop("'water_drive' must be either 'edge' or 'bottom'.")
   if (length(model) > 1) stop("Only one model is acceptable for the aquifer modeling")
   if (model == "pss") {
      if (flow_type == "radial") {
         if (water_drive != "edge") stop("'water_drive' must be 'edge' for the 'radial' ,'pss' model.")
         }
   }
   if (!is.numeric(pressure)) stop("'pressure' must be a numeric vector.")

   if (model == "uss") {
      if (!is.numeric(phi)) stop("'phi' must be a numeric value less than one.")
      if (flow_type == "radial") {
         if (water_drive == "edge") {
            if (!is.numeric(perm_h)) stop("'perm_h' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
            if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
            if (!is.numeric(tetha)) stop("'tetha' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "veh_rad_edge", phi = phi,
                                perm_h = perm_h, h_a = h_a, r_a = r_a, r_R = r_R, tetha = tetha, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("veh_rad_edge", "aquifer")
            }
         if (water_drive == "bottom") {
            if (!is.numeric(perm_h)) stop("'perm_h' must be a numeric.")
            if (!is.numeric(perm_v)) stop("'perm_v' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
            if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "ykh_rad_bottom", phi = phi,
                                perm_h = perm_h, perm_v = perm_v, h_a = h_a, r_a = r_a, r_R = r_R, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("ykh_rad_bottom", "aquifer")
         }
      }
      if (flow_type == "linear") {
         if (water_drive == "edge") {
            if (!is.numeric(perm_h)) stop("'perm_h' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
            if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "nb_lin_edge", phi = phi,
                                perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("nb_lin_edge", "aquifer")
         }
         if (water_drive == "bottom") {
            if (!is.numeric(perm_v)) stop("'perm_v' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
            if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "nb_lin_bottom", phi = phi,
                                perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("nb_lin_bottom", "aquifer")
         }
      }
   }
   if (model == "pss") {
      if (!is.numeric(phi)) stop("'phi' must be a numeric value less than one.")
      if (flow_type == "radial") {
         if (water_drive == "edge") {
            if (!is.numeric(perm_h)) stop("'perm_h' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
            if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
            if (!is.numeric(tetha)) stop("'tetha' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "fetk_rad_edge", phi = phi,
                                perm_h = perm_h, h_a = h_a, r_a = r_a, r_R = r_R, tetha = tetha, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("fetk_rad_edge", "aquifer")
         }
      }
      if (flow_type == "linear") {
         if (water_drive == "edge") {
            if (!is.numeric(perm_h)) stop("'perm_h' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
            if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "fetk_lin_edge", phi = phi,
                                perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("fetk_lin_edge", "aquifer")
         }
         if (water_drive == "bottom") {
            if (!is.numeric(perm_v)) stop("'perm_v' must be a numeric.")
            if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
            if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
            if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
            if (!is.numeric(mu_water)) stop("'mu_water' must be a numeric.")
            if (!is.numeric(c_water)) stop("'c_water' must be a numeric.")
            if (!is.numeric(c_rock)) stop("'c_rock' must be a numeric.")
            aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "fetk_lin_bottom", phi = phi,
                                perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = mu_water,
                                c_water = c_water, c_rock = c_rock, pressure = pressure)
            class(aquifer_lst) <- c("fetk_lin_bottom", "aquifer")
         }
      }
   }
   return(aquifer_lst)
}



#' A list object of class 'time' for aquifer models
#'
#' Create an object of class 'time'
#'
#' @param x a vector of times or a daily sequence of dates
#' @param unit time/date unit of vector x
#'
#' @return a list of class 'time' with all the required parameters for the aquifer_predict() S3 methods
#' @export
#'
#' @examples
#' aquifer_time_1 <- aquifer_time(c(0:4) * 365, unit = "day")
#'
#' aquifer_time_1
#'
#' aquifer_time_2 <- aquifer_time(c(0:4), unit = "month")
#'
#' aquifer_time_2
#'
#' aquifer_time_3 <- aquifer_time(c(0:4), unit = "year")
#'
#' aquifer_time_3
#'
#' aquifer_time_4 <- aquifer_time(seq(as.Date("2020/1/1"), by = "year",
#' length.out = 5), unit = "date")
#'
#' aquifer_time_4


aquifer_time <- function(x, unit = "day") {

   if (!is.character(unit)) stop("'unit' must be a character string of type 'day', 'month', 'year', or 'date'.")
   if (!(unit %in% c("day", "month", "year", "date"))) stop("'unit' must be a character string of type 'day', 'month', 'year', or 'date'.")
   if (unit == "date") {
      if (class(x) != "Date") stop("'x' must be a sequence of type 'Date'.")
   }
   if (unit %in% c("day", "month", "year")) {
      if (!is.vector(x)) stop("'x' must be a vector of days, months, or years.")
      if (!is.numeric(x)) stop("'x' must be a numeric vector of days, months, or years.")
      if (any(duplicated(x))) stop("'x' must be a non-duplicated numeric vector of days, months, or years.")
      if (is.unsorted(x)) stop("'x' must be sorted in an ascending order.")
      time_lst <- list(t = x, unit = unit, reference_date = Sys.Date())
      class(time_lst) <- c(unit, "time")
   }
   if (unit == "date") {
      x <- as.Date(x)
      if (any(duplicated(x))) stop("'x' must be a non-duplicated sequence of dates.")
      if (is.unsorted(x)) stop("'x' must be sorted in an ascending order.")
      time_day <- as.numeric(x - x[1])
      time_lst <- list(t = time_day, unit = unit, reference_date = x[1])
      class(time_lst) <- c("day", "time")
   }
   return(time_lst)
}





#' Generic function for cumulative water influx predictions
#'
#' Generate a data frame of cumulative water influx estimates according to the class of 'aquifer_lst' and 'time_lst' objects
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates according to the class of 'aquifer_lst' and 'time_lst' objects
#'
#' @export
#'
#' @examples
#' aquifer_time_1 <- aquifer_time(c(0:4) * 365, unit = "day")
#' aquifer_param_01 <- aquifer_param(input_unit = "Field", output_unit = "Field",
#' model = "uss", flow_type = "radial", water_drive = "edge", phi = 0.2, perm_h = 100,
#' h_a = 47, r_a = 2e4, r_R = 2e3, tetha = 360, mu_water = 0.34, c_water = 4e-6,
#' c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312))
#' results_01 <- aquifer_predict(aquifer_param_01, aquifer_time_1)
#'
#' results_01
#'
#' @references
#' \insertRef{Yildiz2007}{Raquifer}
#'
#' \insertRef{Nabor1964}{Raquifer}
#'
#' \insertRef{Fetkovich1971}{Raquifer}
#'
#' \insertRef{VanEverdingen1949}{Raquifer}

aquifer_predict <- function(aquifer_lst, time_lst) {
   if ((inherits(aquifer_lst, "aquifer") == TRUE) & (inherits(time_lst, "time"))) {
      UseMethod("aquifer_predict")
   } else {
      if (!inherits(aquifer_lst, "aquifer")) {
         stop("A class of 'aquifer' must be assigned to the 'aquifer_lst' parameter of the aquifer_predict() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the aquifer_predict() function.")
      }
   }
}




#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Van Everdingen-Hurst un-steady state radial flow model, edge-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Van Everdingen-Hurst un-steady state radial flow model, edge-water-drive
#' @export


aquifer_predict.veh_rad_edge <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}


#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Yildiz-Khosravi un-steady state radial flow model, bottom-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Yildiz-Khosravi un-steady state radial flow model, bottom-water-drive
#' @export


aquifer_predict.ykh_rad_bottom <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}


#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Fetkovich pseudo-steady state radial flow model, edge-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Fetkovich pseudo-steady state radial flow model, edge-water-drive
#' @export


aquifer_predict.fetk_rad_edge <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}


#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Nabor-Barham un-steady state linear flow model, edge-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Nabor-Barham un-steady state linear flow model, edge-water-drive
#' @export


aquifer_predict.nb_lin_edge <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}


#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Nabor-Barham un-steady state linear flow model, bottom-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Nabor-Barham un-steady state linear flow model, bottom-water-drive
#' @export


aquifer_predict.nb_lin_bottom <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}



#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Fetkovich pseudo-steady state linear flow model, edge-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Fetkovich pseudo-steady state linear flow model, edge-water-drive
#' @export


aquifer_predict.fetk_lin_edge <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}



#' S3 method for class 'aquifer_predict'
#'
#' Return a data frame of estimated cumulative water influx for the Fetkovich pseudo-steady state linear flow model, bottom-water-drive
#'
#' @param aquifer_lst a list object of class 'aquifer'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of cumulative water influx estimates using the Fetkovich pseudo-steady state linear flow model, bottom-water-drive
#' @export


aquifer_predict.fetk_lin_bottom <- function(aquifer_lst, time_lst) {

   results <- aquifer_predict_interface(aquifer_lst, time_lst)
   return(results)
}

