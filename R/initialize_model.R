library(SimInf)
source("R/events.R")
source("R/maintenance_fun.R")
source("R/ldata.R")


## Define the MRSA model.

## Cumulative mortality rates based on Winpig 2019:
## From birth to weaning: 17.7
## From weaning to finishing: 2.0
## From finishing to slaughter 1.7
## ** Annual average removal rate for sows based on Engblom et al. 2007: 49.5 %. However, this is mostly taken into account with
## the removal of animals due to not getting pregnant/at time of weaning, therefore the rate needs to be adjusted.
## Total of 14.9 % of the removed sows were either euthanized or found dead which is ~7.4 % of the sows removed annually
## However this isn't implemented currently as it would be problematic to end up having piglets without a sow.

## Times spent in units:
## Birth to weaning  currently timer for farrowing is 35 d because piglets are born "too" early
## Weaning to finishing 56 d
## Finishing to slaughter 92 d
## Daily mortality rates are calculated by using the exponential decay function N1 = Ne^(-xt)
## which transforms into x = -log(1-cumulative mortality probability)/time

## Daily mortality probabilities based on above:
## Birth to weaning: 0.006 (calculated based on 35 d, approximately the same if using 33 d)
## Weaning to finishing: 0.0004
## Finishing to slaughter: 0.0002

## half-life for LA-MRSA is dust is 5 days based on Feld at al. 2018
## based on exponential decay function Nt = N(1/2)^(t/t1/2)
## -> 1/2 ^(1/5) -> 0.8705506

## duration of LA-MRSA carriage 17.4d (Broens 2012 "Quantification of transmission...")

## indirect transmission rates here are the mean transmission rates from the medium parameter set obtained through approximate Bayesian computation.
## For manuscript, model was run with imported fitted values.

## starting point gdata, transmission based on Broens' TP values
gdata = c(Nt=0.8705506,
          beta_mature = 0.00,
          beta_piglets = 0.00,
          beta_growers = 0.00,
          beta_finish = 0.00,
          duration = 17.4,
          mortality_grower = 0.0004,
          mortality_finish = 0.0002,
          between_pen = 0.0,
          beta_phi_mature = 0.0018,
          beta_phi_piglets = 0.0051,
          beta_phi_growers = 0.0035,
          beta_phi_finish = 0.0028)

## transmission only through indirect route
transitions <- c("Ssows -> Ssows * phi * beta_phi_mature-> Isows",
                 "Sgilts -> Sgilts * phi * beta_phi_mature -> Igilts",
                 "Spiglets -> Spiglets * phi * beta_phi_piglets-> Ipiglets",
                 "Sgrowers -> Sgrowers * phi * beta_phi_growers-> Igrowers",
                 "Sfinish -> Sfinish * phi * beta_phi_finish-> Ifinish",
                 "Isows -> Isows/duration -> Ssows",
                 "Igilts -> Igilts/duration -> Sgilts",
                 "Ipiglets -> Ipiglets/duration -> Spiglets",
                 "Igrowers -> Igrowers/duration -> Sgrowers",
                 "Ifinish -> Ifinish/duration -> Sfinish",
                 "Sgrowers -> Sgrowers * mortality_grower -> @",
                 "Igrowers -> Igrowers * mortality_grower -> @",
                 "Sfinish -> Sfinish * mortality_finish-> @",
                 "Ifinish -> Ifinish * mortality_finish -> @")


compartments <- c("Ssows",
                  "Isows",
                  "Sgilts",
                  "Igilts",
                  "Spiglets",
                  "Ipiglets",
                  "Sgrowers",
                  "Igrowers",
                  "Sfinish",
                  "Ifinish",
                  "pentype",
                  "section",
                  "countdown")

##1 Move a sow from breeding to gestation/the opposite direction: select = 1 shift = 0
##2 Move a sow from gestation to farrowing: select = 1 shift = 0
##3 Move a sow from farrowing to breeding: select = 1 shift = 0
##4 Move a gilt grower to breeding: select = 2 shift = 0
##5 Move a gilt from breeding to gestation/the opposite direction: select = 2 shift = 0
##6 Move a gilt from gestation to farrowing: select = 2 shift = 1
##7 Piglets born: select = 3 shift = 0
##8 Piglets weaned to grower: select = 4 shift = 2
##9 Growers moved to finisher: select = 5 shift = 3
##10 Finishers to slaughter: select= 6 shift = 0
##11 piglets move to gilt grower (-> gilts): select = 4, shift = 4
##12 piglet mortality: select = 4, shift = 0
##13 grower mortality: select = 5, shift = 0
##14 sow mortality: select = 1, shift = 0
##15 gilt mortality: select = 2, shift = 0
##16 Purchase gilt: select = 7, shift = 0
##17 Gilt growers to finishing (intTrans): select = 5, shift = 3 ##NOT IN USE
##18 Growers moved to growing buffer: select = 5, shift = 0
##19 Finishers moved to finishing buffer: select = 6, shift = 0
##20 Purchase infected gilt: select = 8, shift = 0


E <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 1, 0,
              0, 1, 0, 0, 0, 0, 0, 1,
              0, 0, 1, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0),
            byrow = TRUE,
            ncol = 8,
            dimnames = list(compartments, c("Sows",
                                            "Gilts",
                                            "Susceptible piglets",
                                            "Piglets",
                                            "Growers",
                                            "Finishers",
                                            "Susceptible gilts",
                                            "Infected gilts")))

N <- matrix(c(0, 0, 0,  0,
              0, 0, 0,  0,
             -2, 0, 0,  0,
             -2, 0, 0,  0,
              0, 2, 0, -2,
              0, 2, 0, -2,
              0, 0, 2,  0,
              0, 0, 2,  0,
              0, 0, 0,  0,
              0, 0, 0,  0,
              0, 0, 0,  0,
              0, 0, 0,  0,
              0, 0, 0,  0),
            dimnames = list(compartments, c("gilts_to_farrowing",
                                            "weaning_of_piglets",
                                            "grow_to_finish",
                                            "piglets_to_gilts")),
            ncol = 4,
            byrow = TRUE)

u0 <- initialize_u0()

## feeding some growing gilts to the model
u0$Sgilts[u0$pentype == "Gilt breeding"][1] <- 26
u0$countdown[u0$pentype == "Gilt breeding"][1] <- 33
u0$Sgilts[u0$pentype == "Gilt growing"][1] <- 10
u0$countdown[u0$pentype == "Gilt growing"][1] <- 165
u0$Sgilts[u0$pentype == "Gilt growing"][2] <- 10
u0$countdown[u0$pentype == "Gilt growing"][2] <- 165
u0$Sgilts[u0$pentype == "Gilt growing"][3] <- 10
u0$countdown[u0$pentype == "Gilt growing"][3] <- 165


nPens <- nrow(u0)

v0 <- data.frame(phi = rep(0, nPens))

## environmental accumulation and decay
pts_fun <- c("v_new[0] = gdata[0] * v[0] + u[1] + u[3] + u[5] + u[7] + u[9];",
             "if (ldata[0] > 0)",
             "    v_new[0] += gdata[8] * (v[-1] - v[0]);",
             "if (ldata[1] > 0)",
             "    v_new[0] += gdata[8] * (v[1] - v[0]);",
             "return 1;")

model <- mparse(transitions = transitions,
                compartments = compartments,
                ldata = ldata,
                gdata = gdata,
                u0 = u0,
                v0 = v0,
                E = E,
                N = N,
                tspan = 1:2000,
                pts_fun = pts_fun)

## Install the MRSA model as an R package.
path <- tempdir()
unlink(path, recursive = TRUE)
package_skeleton(model = model, name = "MRSA", path = path)
pkg <- file.path(path, "MRSA")
install.packages(pkg, repos = NULL, type = "source")

## Load the new library:
try(detach("package:MRSA", unload = TRUE), silent = TRUE)
library(MRSA)

model <- MRSA(ldata = ldata,
              gdata = gdata,
              u0 = u0,
              v0 = v0,
              tspan = 1:2000,
              events = NULL)

## Iterate over all time-points in tspan.

final_result <- NULL
result <- u0
result <- cbind(data.frame(node = seq(nrow(u0)),
                           time = 0,
                           phi = 0), result)
result$pentype <- as.numeric(result$pentype)
class(result) <- c(class(result), "MRSA_single_step_trajectory")
result <- clean_trajectory(result)


