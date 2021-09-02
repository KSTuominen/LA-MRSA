########
library("SimInf")
source("R/init_model_d730.R")
source("R/maintenance_fun.R")
source("R/prevalence_percentile.R")

## run a trajectory

events <- read.csv("events/events-d731-3000-2002.csv", header = TRUE)

model@tspan <- as.double(731:3000)
model@events <- SimInf_events(model@events@E, model@events@N, events)

farr_pens <- 99:254
growing_pens <- 255:514 ## buffer pens excluded
gilt_pens <- 1131:1155 ## gilt growing pens
finishing_pens <- 541:1080 # buffer pens excluded

##### collect event times and nodes for animals entering into farrowing, growing and finishing

## collect events for sows (select 1)/gilts (select 2) entering to farrowing units
## we need to define the destination nodes as well because select isn't specific for these movements
farrowing<- events[events$event == "extTrans" & (events$select == 1 | events$select == 2) & events$dest %in% farr_pens & events$time > 2600,]

## weaned piglets
weaning <- events[events$event == "extTrans" & events$dest %in% growing_pens & events$time > 2600,]

##  growers to finishing
finishing <- events[events$event == "extTrans" & events$dest %in% finishing_pens & events$time > 2600,]

##### use the events to collect sampling times and destination nodes
## m1 sow entering farrowing unit
## m2 three days after entering
## m3 three weeks after
## m4 7 days after weaning (pigs 42 d old)
## m5 35 after weaning (70 d old)
## m6 84 days after moving to finishing

m1 <- farrowing[,c(2,4)]
m2 <- m1
m3 <- m1
m4 <- weaning[,c(2,4)]
m5 <- m4
m6 <- finishing[,c(2,4)]
m2$time <- m2$time+3
m3$time <- m3$time+21
m4$time <- m4$time+7
m5$time <- m5$time+35
m6$time <- m6$time+84
m1 <- m1[m2$time <= 3000,]
m2 <- m2[m2$time <= 3000,]
m3 <- m3[m3$time <= 3000,]
m4 <- m4[m4$time <= 3000,]
m5 <- m5[m5$time <= 3000,]
m6 <- m6[m6$time <= 3000,]

## expected values, 10%, 50%, 90% percentiles from prevalences in Broens et al 2012
low <- prev_percentiles$`10%`
med <- prev_percentiles$`50%`
high <- prev_percentiles$`90%`

ex_m1 <- high[1]
ex_m2_sow <- high[2]
ex_m2_piglet <- high[3]
ex_m3_sow <- high[4]
ex_m3_piglet <- high[5]
ex_m4 <- high[6]
ex_m5 <- high[7]
ex_m6 <- high[8]

expected <- c("ex_m1" = ex_m1,
              "ex_m2_sow" = ex_m2_sow,
              "ex_m2_piglet" = ex_m2_piglet,
              "ex_m3_sow" = ex_m3_sow,
              "ex_m3_piglet" = ex_m3_piglet,
              "ex_m4" = ex_m4,
              "ex_m5" = ex_m5,
              "ex_m6" = ex_m6)

ngen = 3
ptol = 0.9
npart = 10

result <- run(model)
model@tspan <- sort(as.double(c(731, unique(c(m1$time, m2$time, m3$time, m4$time, m5$time, m6$time)))))
m1 <- match(unique(m1$time), model@tspan)
m2 <- match(unique(m2$time), model@tspan)
m3 <- match(unique(m3$time), model@tspan)
m4 <- match(unique(m4$time), model@tspan)
m5 <- match(unique(m5$time), model@tspan)
m6 <- match(unique(m6$time), model@tspan)

acceptFun <- function(result, generation, tol, expected,
                      farr_pens, m1, m2, m3, growing_pens, m4, m5,
                      finishing_pens, m6){

    sows_farr <- prevalence(result, Isows~Ssows+Isows, index = farr_pens)
    sows_m1 <- sows_farr[m1, "prevalence"]
    sows_m2 <- sows_farr[m2, "prevalence"]
    sows_m3 <- sows_farr[m3, "prevalence"]

    piglet_farr <- prevalence(result, Ipiglets~Spiglets+Ipiglets, index = farr_pens)
    piglet_m2 <- piglet_farr[m2, "prevalence"]
    piglet_m3 <- piglet_farr[m3, "prevalence"]

    wean <- prevalence(result, Igrowers~Sgrowers+Igrowers, index = growing_pens)
    wean_m4 <- wean[m4, "prevalence"]
    wean_m5 <- wean[m5, "prevalence"]

    fin <- prevalence(result, Ifinish~Sfinish+Ifinish, index = finishing_pens)
    fin_m6 <- fin[m6, "prevalence"]

  ## calculate the distance of the observed values from expected
    dist <- sum((expected["ex_m1"] - sows_m1)^2,
                (expected["ex_m2_sow"] - sows_m2)^2,
                (expected["ex_m2_piglet"] - piglet_m2)^2,
                (expected["ex_m3_sow"] - sows_m3)^2,
                (expected["ex_m3_piglet"] - piglet_m3)^2,
                (expected["ex_m4"] - wean_m4)^2,
                (expected["ex_m5"] - wean_m5)^2,
                (expected["ex_m6"] - fin_m6)^2)

    tol <- tol[generation]
    abc_accept(dist < tol, tol)
}

tol <- c(10, 9.9, 9.8, 9.7, 9.6, 9.5, 9.4, 9.3, 9.2, 9.1, 9, 8.9, 8.8,
         8.7, 8.6, 8.5, 8.4, 8.3, 8.2, 8.1, 8, 7.9, 7.8, 7.7, 7.6,
         7.5, 7.4, 7.3, 7.2, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4,
         6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2,
         5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.9,
         3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3, 2.9, 2.8, 2.7,
         2.6, 2.5, 2.4, 2.3, 2.29, 2.28, 2.26, 2.24, 2.22, 2.2, 2.18,
         2.16, 2.14, 2.139, 2.138, 2.137, 2.136, 2.135, 2.134, 2.133,
         2.132, 2.131, 2.13, 2.129, 2.128, 2.127, 2.126, 2.125, 2.124,
         2.123, 2.122, 2.121, 2.12, 2.119, 2.118, 2.117, 2.116, 2.115,
         2.114, 2.113, 2.112, 2.111, 2.11, 2.109, 2.108, 2.107, 2.106,
         2.105, 2.104, 2.103, 2.102, 2.101, 2.1, 2.099, 2.098, 2.097,
         2.096, 2.095, 2.094, 2.093, 2.092, 2.091, 2.09, 2.089, 2.088,
         2.087, 2.086, 2.085, 2.084, 2.083, 2.082, 2.081, 2.08, 2.079,
         2.078, 2.077, 2.076, 2.075, 2.074, 2.073, 2.072, 2.071, 2.07,
         2.069, 2.068, 2.067, 2.066, 2.065, 2.064, 2.063, 2.062,
         2.061, 2.06, 2.059, 2.058, 2.057, 2.056, 2.055, 2.054, 2.053,
         2.052, 2.051, 2.05, 2.049, 2.048, 2.047, 2.046, 2.045, 2.044,
         2.043, 2.042, 2.041, 2.04, 2.039, 2.038, 2.037, 2.036, 2.035,
         2.034, 2.033, 2.032, 2.031, 2.03, 2.029, 2.028, 2.027, 2.026,
         2.025, 2.024, 2.023, 2.022, 2.021, 2.02, 2.019, 2.018, 2.017,
         2.016, 2.015, 2.014, 2.013, 2.012, 2.011, 2.01, 2.009, 2.008,
         2.007, 2.006, 2.005, 2.004, 2.003, 2.002, 2.001, 2)


fit <- abc(model,
           priors = c(beta_phi_mature~U(0.0001, 0.3),
                      beta_phi_piglets~U(0.0001, 0.3),
                      beta_phi_growers~U(0.0001, 0.3),
                      beta_phi_finish~U(0.0001, 0.3)),
           ngen = 1,
           npart = 200,
           fn = acceptFun,
           tol = tol,
           expected = expected,
           farr_pens = farr_pens,
           m1 = m1,
           m2 = m2,
           m3 = m3,
           growing_pens = growing_pens,
           m4 = m4,
           m5 = m5,
           finishing_pens = finishing_pens,
           m6 = m6,
           verbose = TRUE)
save(fit, file = "fit/high/fitted.Rda")
pdf("fit/high/fitted.pdf")
plot(fit, main = paste("Gen:", 1, " Tolerance:", tail(as.vector(fit@epsilon), 1)))
dev.off()

load("fit/high/fitted.Rda")
for(i in 2:200) {
    fit <- continue(fit,
                    tol = tol,
                    expected = expected,
                    farr_pens = farr_pens,
                    m1 = m1,
                    m2 = m2,
                    m3 = m3,
                    growing_pens = growing_pens,
                    m4 = m4,
                    m5 = m5,
                    finishing_pens = finishing_pens,
                    m6 = m6,
                    verbose = TRUE)
    save(fit, file = "fit/high/fitted.Rda")
    pdf("fit/high/fitted.pdf")
    plot(fit, main = paste("Gen:", i, " Tolerance:", tail(as.vector(fit@epsilon), 1)))
    dev.off()
}

