## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE, warning = FALSE-------------------------------------------
#  install.packages("Raquifer")

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "uss", 
                            flow_type = "radial", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, r_a = 5 * 14892, r_R = 14892, tetha = 180,
                            mu_water = 0.485, c_water = 3.88e-6, c_rock = 2e-6, 
                            pressure = c(1640,1600,1400,1200,1000,800,600,400))

aqu_time

parameters

pred_veh <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_veh)

pred_veh %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "uss", 
                            flow_type = "radial", water_drive = "bottom", phi = 0.27, perm_h = 64.2, 
                            perm_v = 64.2, h_a = 20, r_a = 5 * 14892, r_R = 14892, 
                            mu_water = 0.485, c_water = 3.88e-6, c_rock = 2e-6, 
                            pressure = c(1640,1600,1400,1200,1000,800,600,400))

pred_ykh <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_ykh)

pred_ykh %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "pss", 
                            flow_type = "radial", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, r_a = 5 * 14892, r_R = 14892, tetha = 180,
                            mu_water = 0.485, c_water = 3.88e-6, c_rock = 2e-6, 
                            pressure = c(1640,1600,1400,1200,1000,800,600,400))

pred_fetk <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_fetk)

pred_fetk %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "uss", 
                            flow_type = "linear", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, w_a = 29784, l_a = 161145, mu_water = 0.485, c_water = 3.88e-6, 
                            c_rock = 2e-6, pressure = c(1640,1600,1400,1200,1000,800,600,400))

pred_nb_01 <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_nb_01)

pred_nb_01 %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "uss", 
                            flow_type = "linear", water_drive = "bottom", phi = 0.27, perm_v = 64.2, 
                            h_a = 20, w_a = 29784, l_a = 161145, mu_water = 0.485, c_water = 3.88e-6, 
                            c_rock = 2e-6, pressure = c(1640,1600,1400,1200,1000,800,600,400))

pred_nb_02 <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_nb_02)

pred_nb_02 %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "pss", 
                            flow_type = "linear", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, w_a = 29784, l_a = 161145, mu_water = 0.485, c_water = 3.88e-6, 
                            c_rock = 2e-6, pressure = c(1640,1600,1400,1200,1000,800,600,400))

parameters

pred_fetk_02 <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_fetk_02)

pred_fetk_02 %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = c(0,0.368,2.439,4.957,7.732,11.926,18.126,30.044) * 365, unit = "day")

parameters <- aquifer_param(input_unit = "Field", output_unit = "Field", model = "pss", 
                            flow_type = "linear", water_drive = "bottom", phi = 0.27, perm_v = 64.2, 
                            h_a = 20, w_a = 29784, l_a = 161145, mu_water = 0.485, c_water = 3.88e-6, 
                            c_rock = 2e-6, pressure = c(1640,1600,1400,1200,1000,800,600,400))

pred_fetk_03 <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_fetk_03)

pred_fetk_03 %>% ggplot(aes(x = `Time (days)`, y = `We (MMbbl)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = seq(as.Date("2020/1/1"), by = "year", length.out = 8), unit = "date")

parameters <- aquifer_param(input_unit = "Field", output_unit = "SI", model = "uss", 
                            flow_type = "radial", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, r_a = 5 * 14892, r_R = 14892, tetha = 180,
                            mu_water = 0.485, c_water = 3.88e-6, c_rock = 2e-6, 
                            pressure = c(1640,1600,1400,1200,1000,800,600,400))

aqu_time

parameters

pred_veh <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_veh)

pred_veh %>% ggplot(aes(x = `Time (days)`, y = `We (m3)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = TRUE------
library(Raquifer)
library(ggplot2)
library(magrittr)

aqu_time <- aquifer_time(x = 1:8, unit = "month")

parameters <- aquifer_param(input_unit = "Field", output_unit = "SI", model = "uss", 
                            flow_type = "radial", water_drive = "edge", phi = 0.27, perm_h = 64.2, 
                            h_a = 20, r_a = 5 * 14892, r_R = 14892, tetha = 180,
                            mu_water = 0.485, c_water = 3.88e-6, c_rock = 2e-6, 
                            pressure = c(1640,1600,1400,1200,1000,800,600,400))

aqu_time

parameters

pred_veh <- aquifer_predict(aquifer_lst = parameters, time_lst = aqu_time)

head(pred_veh)

pred_veh %>% ggplot(aes(x = `Time (months)`, y = `We (m3)`)) +
  geom_point(size = 3, color = "blue") +
  theme_bw()


