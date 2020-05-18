context("test-check")

test_that("error if aquifer_time unit is not 'day'", {

   expect_error(aquifer_time(c(1:7300), unit = "days")   )

})

test_that("error if aquifer_time unit is not 'month'", {

   expect_error(aquifer_time(c(1:360), unit = "months")   )

})


test_that("error if aquifer_time unit is not 'year'", {

   expect_error(aquifer_time(c(1:30), unit = "years")   )

})



test_that("error if aquifer_param phi value is 'NULL'", {

   expect_error(aquifer_param(input_unit = "Field", output_unit = "Field",
                              model = "uss", flow_type = "radial", water_drive = "edge", perm_h = 100,
                              h_a = 47, r_a = 2e4, r_R = 2e3, tetha = 360, mu_water = 0.34, c_water = 4e-6,
                              c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312)))

})

test_that("error if aquifer_param model is missing", {

   expect_error(aquifer_param(input_unit = "Field", output_unit = "Field",
                              flow_type = "radial", water_drive = "edge", phi = 0.27, perm_h = 100,
                              h_a = 47, r_a = 2e4, r_R = 2e3, tetha = 360, mu_water = 0.34, c_water = 4e-6,
                              c_rock = 3e-6, pressure = c(3456, 3425, 3387, 3350, 3312)))

})
