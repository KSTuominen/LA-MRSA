source("R/initialize_model.R")
source("R/daily.R")
events = NULL
result <- list(result=result, events = events)
final_result <- NULL
totals <- NULL
for(tspan in 1:3000) {
    result <- daily(result$result,
                    model,
                    tspan)
    final_result <- rbind(final_result,
                          result$result)
    events <- rbind(events,
                    result$events)

}
