## run the model with premade events and all animals infected at d730
## NOTE! With the pre-generated events, the result -dataframe countdown column is unaltered (99999),
## even when there are animals in the pen

library("SimInf")
source("R/init_model_d730.R")
source("R/maintenance_fun.R")
source("R/prevalence_percentile.R")

## pre-generated events from day 731-3000. All nimals infected at day 730.

events <- read.csv("events/events-d731-3000-2002.csv", header = TRUE)

model@tspan <- as.double(731:3000)
model@events <- SimInf_events(model@events@E, model@events@N, events)

final_result <- run(model)

## Plot the whole herd prevalance over the study period in this trajectory
plot(prevalence(final_result, Isows + Igilts + Ipiglets + Igrowers +
                    Ifinish ~ Isows + Igilts + Ipiglets +
                    Igrowers + Ifinish + Ssows +
                    Sgilts + Spiglets + Sgrowers +
                    Sfinish), type = "l", ylim = c(0, 1))
