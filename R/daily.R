daily <- function(result,
                  model,
                  tspan) {
    ## step the counter
    result$countdown <- as.integer(ifelse(result$countdown == 0 | result$countdown == 99999,
                                          result$countdown,
                                          result$countdown - 1))

    ## check that we don't have timer 99999 and animals in a pen at the same time
    null_timer <- result$countdown == 99999
    non_empty_pen <- result[, "npigs"] != 0
    stopifnot(!any(null_timer & non_empty_pen))

    ## Remove the events from the model object and set the tspan
    model@events <- SimInf_events(model@events@E, model@events@N, NULL)
    model@tspan <- as.numeric(tspan)

    ## Set the countdown to 99999 if there are no animals but the
    ## timer reached zero. This happens when the downtime in a pen
    ## expires.
    result[result$countdown == 0 & nanimals(result) == 0, "countdown"] <- 99999

    ## Get the rows that we need to work with 'today'
    residual <- result[result$countdown == 0,]

    ## Initialize the herd with gilts to breeding
    ob <- initialize_herd(result,
                          residual)

    ob <- abortion(ob$result,
                   ob$residual,
                   ob$events)

    ob <- piglet_mortality(ob$result,
                           ob$residual,
                           ob$events)

    ob <- cross_foster(ob$result,
                       ob$residual,
                       ob$events)

    ob <- finishing_mingle(ob$result,
                           ob$residual,
                           ob$events)

    ob <- finishing_slaughter(ob$result,
                              ob$residual,
                              ob$events)

    ## Start a series of events the while control flow just allows you
    ## to break out anytime you don't have anything left to do or you
    ## get to the end of the control flow
    while(TRUE) {
        if(nrow(residual) == 0) break

        ob <- weaning(ob$result,
                      ob$residual,
                      ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- growing_to_finishing(ob$result,
                                   ob$residual,
                                   ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- growing_buffer_to_finishing(ob$result,
                                          ob$residual,
                                          ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- farrowing(ob$result,
                        ob$residual,
                        ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- gilt_to_gestation(ob$result,
                                ob$residual,
                                ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- sow_to_gestation(ob$result,
                               ob$residual,
                               ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- gilts_to_slaughter(ob$result,
                                 ob$residual,
                                 ob$events)
        if(nrow(ob$residual) == 0) break

        ob <- buffer_culling(ob$result,
                             ob$residual,
                             ob$events)
        if(nrow(ob$residual) == 0) break

        ## If we still have rows left in the residual that we have not
        ## taken care of, we have a problem.
        stopifnot(nrow(ob$residual) == 0)
    }

    ## Finally run the trajectory and return the result
    model@u0 <- as_u_matrix(ob$result)
    model@v0 <- as_v_matrix(ob$result)

    model@events <- SimInf_events(model@events@E, model@events@N, ob$events)
    result <- trajectory(run(model))
    return(list(result = clean_trajectory(result), events = ob$events))
}

##' In this function we add gilts to the breeding at weekly intervals
##' in the beginning of the study to startup to herd.
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
initialize_herd <- function(result,
                            residual = NULL,
                            events = NULL) {

    ## calculate if there are animals in the gilt growing unit
    gilt.growing.ind <- result$node[result$pentype == "Gilt growing"]
    total.growing.gilts <- sum(result$Spiglets[gilt.growing.ind], result$Ipiglets[gilt.growing.ind])

    #  Generating the initial gilts, 1 batch to breeding each week until there is a flow of animals from elsewhere
    ## weeks for importing every other week
    weeks <-  1:21
    interval <- 7
    gilts.in <- (weeks*interval+1)

    if(result$time[1] %in% gilts.in){
        ## finding empty destination pen
        empty.gilt.breeding <- empty.pens(result, "Gilt breeding", events)[1]

        ## if there is none break out
        if(is.na(empty.gilt.breeding)) {
            return(list(result = result,
                        residual = residual,
                        events = events))
        }

        ## type 16 generates susceptible gilts and 20 infected gilts
        ## currently only feeding in susceptible ones
        event1 <- event(type = 16,
                        time = result$time[1],
                        node = empty.gilt.breeding,
                        dest = 0,
                        n = 22,
                        proportion = 0)

        # event2 <- event(type = 20,
        #                 time = result$time[1],
        #                 node = empty.gilt.breeding,
        #                 dest = 0,
        #                 n = 21,
        #                 proportion = 0)

        result$countdown[result$node %in% empty.gilt.breeding] <- timers("Breeding")

        events <- rbind(events, event1)
    }

    return(list(result = result,
                residual = residual,
                events = events))
}

#' Slaughtering finishing pigs in several batches
#'
#' @param result The result from the previous day or modified result
#'     from the current day.
#' @param residual The rows with timer == 0 that still need to be
#'     processed. The default is NULL if there is nothing to be done
#' @param events The events that have already been generated. The
#'     default is NULL in the case that no events have yet been
#'     generated for this day
#' @return A list of the (modified) result, the new residual and the
#'     new events appended to the events you fed to the function

finishing_slaughter <- function(result,
                                residual = NULL,
                                events = NULL) {
    ## gather pens that are "eligible for slaughter
    new_in <- result[(result$pentype == "Finishing") & (result$countdown %in% c(0, 7, 14)), ]

    ## residual only has lines with countdown 0, but these pens will need to be removed after creating events
    index <- residual$pentype == "Finishing" & residual$countdown == 0

    if (nrow(new_in) == 0) {
      return(list(result = result,
                  residual = residual,
                  events = events))
    }

    prop <- c(0.33, 0.5, 1.0)
    for (i in seq_len(nrow(new_in))){
        row <- new_in[i, ]
        if (row$countdown == 14){
            ## events to slaughter proportion from each pen
            event1 <- event(type = 10,
                            time = result$time[1],
                            node = row$node,
                            dest = 0,
                            n = 0,
                            proportion = prop[1])
            events <- rbind(events, event1)
        }

        else if (row$countdown == 7){
          event1 <- event(type = 10,
                          time = result$time[1],
                          node = row$node,
                          dest = 0,
                          n = 0,
                          proportion = prop[2])
          events <- rbind(events, event1)

        }

        else{
          event1 <- event(type = 10,
                          time = result$time[1],
                          node = row$node,
                          dest = 0,
                          n = 0,
                          proportion = prop[3])
          events <- rbind(events, event1)

        ## set timers when the pens are empty
        result$countdown[result$node == row$node] <- timers("Finishing dt")
        }
    }

    ## remove the rows with countdown 0 from residual

    residual <- residual[!index, ]
    return(list(result = result,
                residual = residual,
                events = events))
}

#' Slaughtering growing gilts when they haven't been needed to replace sows
#'
#' @param result The result from the previous day or modified result
#'     from the current day.
#' @param residual The rows with timer == 0 that still need to be
#'     processed. The default is NULL if there is nothing to be done
#' @param events The events that have already been generated. The
#'     default is NULL in the case that no events have yet been
#'     generated for this day
#' @return A list of the (modified) result, the new residual and the
#'     new events appended to the events you fed to the function

gilts_to_slaughter <- function(result,
                               residual = NULL,
                               events = NULL) {

    ## Index of growing pens that have pigs
    index <- residual$pentype == "Gilt growing" & (residual$Sgilts != 0 | residual$Igilts !=0)

    ## If there is nothing to do then just get out now
    if(sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }

    dep <- residual$node[index]
    event1 <- event(type = 15,
                    time = result$time[1],
                    node = dep,
                    dest = 0,
                    n = 0,
                    proportion = 1)

    events <- rbind(events, event1)
    ## timer for the pen to be empty before new animals.
    result$countdown[result$node %in% dep] <- timers("Gilt growing dt")
    residual <- residual[!index, ]

    return(list(result = result,
                residual = residual,
                events = events))
}

#' Piglet mortality outside siminf. Prob = 0.006 when piglets spend 33 d in the unit
#' (total 35 d - 1 day for weaning - 1 day to ignore the first day of life because of piglets mixing)
#' total mortality between birth and weaning according to 2019 statistics is 17.7
#'
#' @param result The result from the previous day or modified result
#'     from the current day.
#' @param residual The rows with timer == 0 that still need to be
#'     processed. The default is NULL if there is nothing to be done
#' @param events The events that have already been generated. The
#'     default is NULL in the case that no events have yet been
#'     generated for this day
#' @param prob probability of death for each piglet per day
#' @return A list of the (modified) result, the new residual and the
#'     new events appended to the events you fed to the function

piglet_mortality <- function(result,
                      residual = NULL,
                      events = NULL,
                      prob = 0.006) {

  ## excluding the first day of life because we don't want deaths happen the same day as mixing
  index <- (result$pentype == "Farrowing") & (result$npigs !=0) & (result$countdown > 0) & (result$countdown != 34)

  if (sum(index) == 0) {
    return(list(result = result,
                residual = residual,
                events = events))
  }

  ## Calculate the number of piglets to remove from each pen
  nodes <- result[index, "node"]
  n.pigs <- rowSums(result[index, c("Spiglets", "Ipiglets")])
  n.remove <- sapply(n.pigs, function(x) {
    rbinom(1, x, prob)
  })

  ## No reason to create events to remove 0 animals so drop those nodes
  nodes <- nodes[n.remove > 0]
  n.remove <- n.remove[n.remove > 0]
  if (length(nodes) == 0) {
    return(list(result = result,
                residual = residual,
                events = events))
  }

  ## now create events if needed
  event1 <- event(type = 12,
                  time = result$time[1],
                  node = nodes,
                  dest = 0,
                  n = n.remove,
                  proportion = 0)

  ## And concatenate these events to the old events and return the list
  events <- rbind(events, event1)
  return(list(result = result,
              residual = residual,
              events = events))
}

#' Mixing of piglets in the farrowing unit (= cross-fostering)
#'
#' @param result The result from the previous day or modified result
#'     from the current day.
#' @param residual The rows with timer == 0 that still need to be
#'     processed. The default is NULL if there is nothing to be done
#' @param events The events that have already been generated. The
#'     default is NULL in the case that no events have yet been
#'     generated for this day
#' @param prop proportion of piglets being mixed
#' @return A list of the (modified) result, the new residual and the
#'     new events appended to the events you fed to the function
cross_foster <- function(result,
                         residual = NULL,
                         events = NULL,
                         prop = 0.1){

  ## gather pens that have "newborn" piglets
  newborn <- result[(result$pentype == "Farrowing") & (result$Spiglets !=0 | result$Ipiglets !=0) & (result$countdown == 34),]

  if (nrow(newborn) == 0) {
    return(list(result = result,
                residual = residual,
                events = events))
  }

  sections <- unique(newborn$section)

  event1 <- do.call("rbind", lapply(sections, function(x) {

    piglets <- newborn[newborn$section == x, c("Spiglets", "Ipiglets")]
    pigs.per.pen <- rowSums(piglets)
    pens <- newborn$node[newborn$section == x]
    total.pigs <- sum(pigs.per.pen)
    total.mingle <- rpois(1, prop*total.pigs)

    ## don't try to move more pigs than the total
    if(total.mingle > total.pigs) total.mingle <- total.pigs

    ## if none then skip to next section
    if(total.mingle == 0) return(NULL)

    ## generate the source pens
    source_pens <- sample(unlist(mapply(rep,
                                        times = pigs.per.pen,
                                        x = names(pigs.per.pen))),
                          total.mingle,
                          replace = FALSE)

    ## Sample the dest pens one at a time to avoid a movement from a
    ## source pen to the same pen
    dest <- sapply(source_pens, function(x) {
      sample(pens[pens != x], 1)
    })

    ## Generate the events
    event(type = 21,
          time = result$time[1],
          node = source_pens,
          dest = dest,
          n = 1,
          proportion = 0)
  }))

  events <- rbind(events, event1)

  return(list(result = result,
              residual = residual,
              events = events))
}

#' Mixing of pigs in the finishing unit
#'
#' @param result The result from the previous day or modified result
#'     from the current day.
#' @param residual The rows with timer == 0 that still need to be
#'     processed. The default is NULL if there is nothing to be done
#' @param events The events that have already been generated. The
#'     default is NULL in the case that no events have yet been
#'     generated for this day
#' @param prop proportion of animals being mixed
#' @return A list of the (modified) result, the new residual and the
#'     new events appended to the events you fed to the function
finishing_mingle <- function(result,
                             residual = NULL,
                             events = NULL,
                             prop = 1.0) {

  ## gather pens that have "new" finishing pigs
  new_in <- result[(result$pentype == "Finishing") & (result$Sfinish != 0 | result$Ifinish != 0) & (result$countdown == 98), ]

  if (nrow(new_in) == 0) {
    return(list(result = result,
                residual = residual,
                events = events))
  }
  sections <- unique(new_in$section)

  ## make a vector of each animal in a pen and shuffle that?
  event1 <- do.call("rbind", lapply(sections, function(x) {

    pigs <- new_in[new_in$section == x, c("Sfinish", "Ifinish")]
    pigs.per.pen <- rowSums(pigs)
    total_finishers <- sum(pigs.per.pen)
    pens <- new_in$node[new_in$section == x]
    total_pens <- length(pens)
    events_section <- NULL

    ## because the capacity is smaller, add some extra pens
    ## how many new pens we need
    needed_pens <- total_finishers/pen_capacity("Finishing")
    balance <- needed_pens - total_pens

    if(balance > 0){
      ## pens that are still empty in the section
      empty <- result[result$pentype == "Finishing" & result$section == x & result$npigs == 0, "node"]
      ## do nothing if there are no free pens (means that we are overfilling pens)
      if(length(empty) > 0){
        ##  take only the empty pens we need, if not, take what is available
        if(balance <= length(empty)){
          empty <- empty[1:length(balance)]
        }
      }
      ## append the extra pens to pens-variable to get more destination pens
      dest_pens <- c(pens, empty)

    }

    ## if we don't add pens, use the departure pens also as destination pens
    else{
      dest_pens <- pens
    }

    ## the events in SimInf are processed in order of
    ## node then destination. Here we are taking care of the scaling of
    ## proportion of pigs to move out of each pen correctly along the
    ## destination pens but also adjusting for the fact that the pens
    ## are also recieving pigs. These recieved pigs should not be
    ## among those that are moved out when we come to that pen. For
    ## this to work we must assume a uniform population size in each
    ## node before the operation starts. This may not truely be the
    ## case but on average may be true since everything is based on
    ## litter size which is sampled from a Poisson. Below, the line
    ## marked with 'Note' will gradually reduce the proportion of pigs
    ## moved out of the pen as we move along the source pens
    ## vector. If there are 10 pens to move out of then the first pen
    ## will have prop_i == prop the next pen will have prop_i =
    ## prop/(1 + ((1/10) * 1)) the third prop/(1 + ((1/10) * 2)) and
    ## so on.

    ## generate the events
    for(i in seq_len(total_pens)) {

      prop_i <- prop / (1 + ((1 / length(dest_pens)) * (i - 1))) ## Note

      ## split the proportion over destination pens
      p <- prop_i / length(dest_pens)

      ## Scale the p to the order to the dest pen ID
      scaled_p <- scale_p(dest_pens, rep(p, length(dest_pens)))

      ## create the events:
      event_pen <- event(type = 22,
                         time = result$time[1],
                         node = pens[i],
                         dest = dest_pens,
                         n = 0,
                         proportion = scaled_p)

      events_section <- rbind(events_section, event_pen)
    }
    events_section

  }))
  events <- rbind(events, event1)

  ## set the timer for the newly occupied pens to same as with previously occupied pens
  result$countdown[result$node %in% event1$dest] <- timers("Finishing")-1

  return(list(result = result,
              residual = residual,
              events = events))
}

##' In this function we move the gilts to gestation from breeding and
##' deal with buffers in the gestation and breeding.
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
gilt_to_gestation <- function(result,
                              residual = NULL,
                              events = NULL) {

    index <- residual$pentype == "Gilt breeding"

    ## Check for anything to do, take only pens with animals in them

    if (sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }
    ## make this more general, so not row by row but for all gestation pens
    pigs.per.pen <- residual[index, "npigs"]
    total.to.move <-sum(pigs.per.pen)

    ## calculating how many animals in the pens are not pregnant and how many buffer pens they need
    empty <- pregnancy(total.to.move)

    ## save the departure pens
    dep <- residual[index, "node"]

    ## Find the number of pens that we need at the destination for buffer
    total.buffer.pens <- ceiling(empty / pen_capacity("Gilt breeding buffer"))

    ## Get the available gilt breeding buffer pens:
    empty.gilt.buffer.pens <- free.pens(result, "Gilt breeding buffer", residual)

    ## Check that there are enough pens in the buffer and throw error if not
    stopifnot(length(empty.gilt.buffer.pens) >= total.buffer.pens)
    empty.gilt.buffer.pens <- empty.gilt.buffer.pens[seq_len(total.buffer.pens)]

    ## dataframe with only pen ids of the pens to be sampled and
    ## sum of both susceptible and infected piglets
    sampling.gilts <- data.frame(node = dep,
                                 pigs = pigs.per.pen,
                                 stringsAsFactors = FALSE)

    ## Divide pigs going to buffer evenly among dest pens
    dest.pens.vec <- rep(empty.gilt.buffer.pens,
                         length.out = empty)

    ## samples the pigs from the source pens. Having only one pig in one pen needs to be handled
    ## separately, because otherwise R samples from 1:penID instead from than penId only
    if(nrow(sampling.gilts == 1) & sampling.gilts$pigs[1] == 1){
      source_pens <- sampling.gilts$node
    }

    else{
      source_pens <- sample(unlist(mapply(rep,
                                          times = sampling.gilts$pigs,
                                          x = sampling.gilts$node)),
                                          empty,
                                          replace = FALSE)
    }

    ## Keep track of those that have been removed
    removed <- table(source_pens)
    removed <- removed[match(sampling.gilts$node, names(removed))]
    removed[is.na(removed)] <- 0
    sampling.gilts$pigs <- sampling.gilts$pigs - removed

    if(empty > 0) {
      event1 <- event(type = 5,
                      time = result$time[1],
                      node = source_pens,
                      dest = dest.pens.vec,
                      n = 1,
                      proportion = 0)
      events <- rbind(events, event1)
    }

    ## Set the timers for the destination pens
    result$countdown[result$node %in% dest.pens.vec] <- timers("Breeding buffer")

    ## How many gestation pens do we need for the remainder that are
    ## going onto gestation?
    npens <- ceiling(sum(sampling.gilts$pigs) / pen_capacity("Gilt gestation"))

    ## Get some empty pens:
    empty.gilt.gest.pens <- sectioning(result, "Gilt gestation", npens, events)

    ## stop if there are fewer than needed
    stopifnot(length(empty.gilt.gest.pens) >= npens)

    ## keep just those needed
    dest.pens.vec <- empty.gilt.gest.pens[seq_len(npens)]

    ## determine where to put which animal:
    dest.pens.vec <- sort(rep(dest.pens.vec,
                         length.out = sum(sampling.gilts$pigs)))

    ## And which animals to move
    source.pens.vec <- unlist(mapply(rep,
                                     times = sampling.gilts$pigs,
                                     x = sampling.gilts$node))

    if(sum(sampling.gilts$pigs) > 0) {
        event1 <- event(type = 5,
                        time = result$time[1],
                        node = source.pens.vec,
                        dest = dest.pens.vec,
                        n = 1,
                        proportion = 0)
        events <- rbind(events, event1)
    }

    ## set timers
    result$countdown[result$node %in% dest.pens.vec] <- timers("Gestation")

    ## set the downtown in all the source pens at once
    result$countdown[result$node %in% sampling.gilts$node] <- timers("Breeding dt")

    ## modify the residual and return
    residual <- residual[!index, ]
    return(list(result = result,
                residual = residual,
                events = events))
}

##' In this function we move the sows to gestation from breeding and
##' deal with buffers in the gestation and breeding.
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
sow_to_gestation <- function(result,
                             residual = NULL,
                             events = NULL) {

    index <- residual$pentype == "Sow breeding"
    ## Check for anything to do, take only pens with animals in them

    if (sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }

    ## make this more general, so not row by row but for all gestation pens
    pigs.per.pen <- residual[index, "npigs"]
    total.to.move <-sum(pigs.per.pen)

    ## calculating how many animals in the pens are not pregnant and how many buffer pens they need
    empty <- pregnancy(total.to.move)

    ## save the departure pens
    dep <- residual[index, "node"]

    ## Find the number of pens that we need at the destination for buffer
    total.buffer.pens <- ceiling(empty / pen_capacity("Sow breeding buffer"))

    ## Get the availble sow breeding buffer pens:
    empty.sow.buffer.pens <- free.pens(result, "Sow breeding buffer", residual)

    ## Check that there are enough pens in the buffer and throw error if not
    stopifnot(length(empty.sow.buffer.pens) >= total.buffer.pens)
    empty.sow.buffer.pens <- empty.sow.buffer.pens[seq_len(total.buffer.pens)]

    ## dataframe with only pen ids of the pens to be sampled and
    ## sum of both suspectible and infected piglets
    sampling.sows <- data.frame(node = dep,
                                 pigs = pigs.per.pen,
                                 stringsAsFactors = FALSE)

    ## Divide pigs going to buffer evenly amoung dest pens
    dest.pens.vec <- rep(empty.sow.buffer.pens,
                         length.out = empty)
    if(nrow(sampling.sows == 1)& sampling.sows$pigs[1] == 1){
      source_pens <- sampling.sows$node
    }

    else{
    ## samples the pigs from the source pens.
    source_pens <- sample(unlist(mapply(rep,
                                        times = sampling.sows$pigs,
                                        x = sampling.sows$node)),
                          empty,
                          replace = FALSE)
    }

    ## Keep track of those that have been removed
    removed <- table(source_pens)
    removed <- removed[match(sampling.sows$node, names(removed))]
    removed[is.na(removed)] <- 0
    sampling.sows$pigs <- sampling.sows$pigs - removed

    if(empty > 0) {
        event1 <- event(type = 1,
                        time = result$time[1],
                        node = source_pens,
                        dest = dest.pens.vec,
                        n = 1,
                        proportion = 0)
        events <- rbind(events, event1)
    }

    ## Set the timers for the destination pens
    result$countdown[result$node %in% dest.pens.vec] <- timers("Breeding buffer")

    ## How many gestation pens do we need for the remainder that are
    ## going onto gestation?
    npens <- ceiling(sum(sampling.sows$pigs) / pen_capacity("Sow gestation"))

    ## Get some empty pens:
    empty.sow.gest.pens <- sectioning(result, "Sow gestation", npens, events)

    ## stop if there are fewer than needed
    stopifnot(length(empty.sow.gest.pens) >= npens)

    ## keep just those needed
    dest.pens.vec <- empty.sow.gest.pens[seq_len(npens)]

    ## determine where to put which animal:
    dest.pens.vec <- sort(rep(dest.pens.vec,
                         length.out = sum(sampling.sows$pigs)))

    ## And which animals to move
    source.pens.vec <- unlist(mapply(rep,
                                     times = sampling.sows$pigs,
                                     x = sampling.sows$node))

    if(sum(sampling.sows$pigs) > 0) {
        event1 <- event(type = 1,
                        time = result$time[1],
                        node = source.pens.vec,
                        dest = dest.pens.vec,
                        n = 1,
                        proportion = 0)
        events <- rbind(events, event1)
    }

    ## set timers
    result$countdown[result$node %in% dest.pens.vec] <- timers("Gestation")

    ## set the downtown in all the source pens at once
    result$countdown[result$node %in% sampling.sows$node] <- timers("Breeding dt")

    ## modify the residual and return
    residual <- residual[!index, ]
    return(list(result = result,
                residual = residual,
                events = events))
}

##' Farrowing: In this function we move the sows into the farrowing room and
##' generate piglet (birth events).
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
farrowing <- function(result,
                      residual = NULL,
                      events = NULL) {

    ## collecting all the gestation pens with countdown == 0 that have animals
    index <- (residual$pentype == "Sow gestation" | residual$pentype == "Gilt gestation")

    if (sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }

    pens.to.farrowing <- residual[index,]

    ## how many animals to farrowing in total
    total.to.farrow <- sum(pens.to.farrowing$npigs)

    gestation.pen.id <- unlist(mapply(rep,
                                      times = pens.to.farrowing$npigs,
                                      x = pens.to.farrowing$node))

    ## finding an empty section and return the node list
    empty.farrow.pens <- sectioning(result, "Farrowing", 0, events)

    ## stop if there are not enough pens in the farrowing room. There
    ## are several checks that determine what to do with extra sows
    ## and log them. Let's just throw an error for now.
    if(length(empty.farrow.pens) < total.to.farrow){
        excess <- total.to.farrow - length(empty.farrow.pens)
        culled <- paste("D", tspan, ":", excess, "animal(s) culled because we ran out of farrowing pens in a section")
        warning(culled)
        ## sample which animals we are going to cull
        pen.id <- pens.to.farrowing$node
        remove <- sample(pen.id, excess, replace = TRUE)
        culled.row <- residual[residual$node %in% remove,]

        ## define the type for the exit event based on the pentype
        type <- rep(14, length(remove))
        rem.type <- pens.to.farrowing$pentype[match(remove,
                                                    pens.to.farrowing$node)] == "Gilt gestation"
        type[rem.type] <- 15

        event1 <- event(type = type,
                        time = result$time[1],
                        node = remove,
                        dest = 0,
                        n = 1,
                        proportion = 0)
        events <- rbind(events, event1)

        ## remove the culled animals from the gestation.pen.id vector
        gestation.pen.id <- unlist(sapply(union(gestation.pen.id, remove), function(i)
          rep(i, each = length(gestation.pen.id[gestation.pen.id == i]) - length(remove[remove == i]))))
        total.to.farrow <- length(gestation.pen.id)
    }

    ## Crop the vector of empty pens to the length of pens we need if
    ## we need less pens than the section capacity
    empty.farrow.pens <- empty.farrow.pens[seq_len(total.to.farrow)]

    ## Sort out what event type we want for each movement
    type <- rep(2, length(gestation.pen.id))
    gt <- pens.to.farrowing$pentype[match(gestation.pen.id,
                                    pens.to.farrowing$node)] == "Gilt gestation"
    type[gt] <- 6

    event1 <- event(type = type,
                    time = result$time[1],
                    node = gestation.pen.id,
                    dest = empty.farrow.pens,
                    n = 1,
                    proportion = 0)

    ## Piglets are born, currently on a same day as sows arrive
    event2<- event(type = 7,
                   time = result$time[1],
                   node = empty.farrow.pens,
                   dest = 0,
                   n = littersize(total.to.farrow),
                   proportion = 0)

    events <- rbind(events, event1, event2)

    ## setting timers for the destination pens
    result$countdown[result$node %in% empty.farrow.pens] <- timers("Farrowing")
    ## timer for the pens to be empty before new animals
    result$countdown[result$node %in% gestation.pen.id] <- timers("Gestation dt")

    residual <- residual[!index, ]

    return(list(result = result,
                residual = residual,
                events = events))
}

##' weaning
##'
##' In this function we both wean pigs to the grower and to the gilt
##' grower and deal with potentially moving animals to buffer grower
##' pens if that is possible at this stage. We also wean sows to the
##' breeding unit and cull sows that are to be culled. The reason we
##' take all of this together is that all feedback on one one
##' another
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
weaning <- function(result,
                    residual = NULL,
                    events = NULL,
                    fraction_to_gilt = 0.05) {

    if (is.null(residual)) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }

    index <- (residual$pentype == "Farrowing") & (residual$npigs !=0)

    if (sum(index) == 0) {
      return(list(result = result,
                  residual = residual,
                  events = events))
    }

    ##  Weaning piglets and sampling for new gilt growers
    #####################################################

    ## how many piglet pens need to be weaned today
    total.pens.to.weaning <- sum(residual$pentype == "Farrowing" & (residual$Spiglets !=0 | residual$Ipiglets !=0))

    ## how many piglets we have to wean in total
    ## can't use "index" in here because we don't want to count breeding buffers
    pens <- residual[residual$pentype == "Farrowing" & (residual$Spiglets !=0 | residual$Ipiglets !=0), c("Spiglets", "Ipiglets")]
    pigs.per.pen <- rowSums(pens)
    total.piglets.to.wean <-sum(pigs.per.pen)

    ## pen ids of the departure pens
    weaning.pen.id <- residual[residual$pentype == "Farrowing" & (residual$Spiglets !=0 | residual$Ipiglets !=0),"node"]

    ## expected number to gilt growing
    piglet.to.gilt <- round(total.piglets.to.wean * fraction_to_gilt)

    ## Number of gilt growing pens required given the capacity
    npens.gilt <- ceiling(piglet.to.gilt / pen_capacity("Gilt growing"))

    ## finding destination pens, only move animals if we have enough destination pens
    empty.gilt.growing <- free.pens(result, "Gilt growing", residual)[seq_len(npens.gilt)]

    ## if there is no space in the gilt grower then don't move there
    if(any(is.na(empty.gilt.growing)) | is.null(empty.gilt.growing)) {
        piglet.to.gilt <- 0
    }

    source_pens <- sample(unlist(mapply(rep,
                                        times = pigs.per.pen,
                                        x = names(pigs.per.pen))),
                          piglet.to.gilt,
                          replace = FALSE)

    ## finding empty growing pens from one section for the ones that are weaned "normally
    empty.growing.pens <- sectioning_multi(result, "Growing", 0)[seq_len(total.pens.to.weaning)]

    ## Sort out which pigs from the source pens are going which way
    weaned_pigs <- table(source_pens)
    weaned_pigs <- weaned_pigs[match(names(pigs.per.pen), names(weaned_pigs))]
    weaned_pigs[is.na(weaned_pigs)] <- 0
    weaned_pigs <- data.frame(node = as.numeric(rep(names(pigs.per.pen), 2)),
                              n = c(as.integer(weaned_pigs), pigs.per.pen - as.integer(weaned_pigs)),
                              type = rep(c(11, 8), each = length(pigs.per.pen)),
                              stringsAsFactors = FALSE)
    weaned_pigs <- weaned_pigs[weaned_pigs$n != 0, ]
    ## pick pens for gilts
    weaned_pigs$dest <-NA
    ## If one gilt destination pen then don't sample. (See help on sample())
    if(length(empty.gilt.growing) > 1) {
        weaned_pigs$dest[weaned_pigs$type == 11] <- sample(empty.gilt.growing, sum(weaned_pigs$type == 11), replace = TRUE)
    }
    if(length(empty.gilt.growing) == 1) {
        weaned_pigs$dest[weaned_pigs$type == 11] <- empty.gilt.growing
    }
    ## pick pens for growers
    weaned_pigs$dest[weaned_pigs$type == 8] <- empty.growing.pens[seq_len(sum(weaned_pigs$type == 8))]

    ## check that you don't index too many pens in the
    ## sectioning_multi vector
    stopifnot(all(!is.na(empty.growing.pens)))

    ## create the events:
    event1 <- event(type = weaned_pigs$type,
                    time = result$time[1],
                    node = weaned_pigs$node,
                    dest = weaned_pigs$dest,
                    n = weaned_pigs$n,
                    proportion = 0)
    events <- rbind(events, event1)
    ## Set the timers for the destination pens
    result$countdown[result$node %in% weaned_pigs$dest[weaned_pigs$type == 11]] <- timers("Gilt growing")
    result$countdown[result$node %in% weaned_pigs$dest[weaned_pigs$type == 8]] <- timers("Growing")

    ## Moving sows from farrowing to breeding (incl. removal)
    ## merging animals from breeding buffers
    ## filling the missing amount of animals from growing gilts
    ##########################################################
    total.sows.to.breeding <- SI.pentype(residual, "Sows", "Farrowing")
    total.sows.from.buffer <- SI.pentype(result, "Sows", "Sow breeding buffer")
    total.gilts.from.buffer <- SI.pentype(result, "Gilts", "Gilt breeding buffer")
    total.buffer <- total.sows.from.buffer + total.gilts.from.buffer
    sows.in.buffer <- result[result$pentype == "Sow breeding buffer" & (result$Ssows != 0 | result$Isows != 0), ]
    sows.in.buffer <- sows.in.buffer[order(sows.in.buffer$countdown, decreasing = FALSE),]
    gilts.in.buffer <- result[result$pentype == "Gilt breeding buffer" & (result$Sgilts != 0 | result$Igilts != 0), ]
    gilts.in.buffer <- gilts.in.buffer[order(gilts.in.buffer$countdown, decreasing = FALSE),]

    ## priorities 1. sows from farrowing 2. sows/gilts from buffer units 3. refill if space

    if (total.sows.to.breeding > 0) {
      ## Part of the sows are going to be removed from the herd at weaning
      ## sampling the amount of animals and their departure pens
      farrowing.pen.ids <- residual[residual$pentype == "Farrowing" & (residual$Ssows !=0 | residual$Isows !=0),"node"]
      remove.sows <- sow.culling(farrowing.pen.ids)

      ## generating events for sows to be removed
      ## it doesn't matter if it's null, then we don't create anything
      event1 <- event(type = 14,
                      time = result$time[1],
                      node = remove.sows,
                      dest = 0,
                      n = 1,
                      proportion = 0)
      events <- rbind(events, event1)

      result$countdown[result$node %in% remove.sows] <- timers("Farrowing dt")

      ## calculating how many sows we have left to move after the removal
      total.sows.to.breeding <- total.sows.to.breeding - length(remove.sows)

      ## removing the removed sows from the pen list that will be weaned
      sow.pens.to.breeding <- farrowing.pen.ids[!(farrowing.pen.ids %in% remove.sows)]

      ## calculating how many pens we need for the weaned sows
      npens.sow <- total.sows.to.breeding %/% pen_capacity("Sow breeding")
      if(total.sows.to.breeding %% pen_capacity("Sow breeding") != 0) {
        npens.sow <- npens.sow + 1
      }
      if(total.sows.to.breeding < pen_capacity("Sow breeding")){
        npens.sow <- 1
      }

      ## using free pens function instead of sectioning because we group up many animals into one pen
      empty.breeding.pen <- free.pens(result, "Sow breeding", residual)[seq_len(npens.sow)]

      ## creating a destination pen vector for the sows
      breeding.div <- total.sows.to.breeding %/% pen_capacity("Sow breeding")
      breeding.remain <- total.sows.to.breeding %% pen_capacity("Sow breeding")
      div.vec <- NULL
      if (breeding.div == 0){
        dest.pens.vec <- rep(empty.breeding.pen, total.sows.to.breeding)
      }
      else {
        for (i in 1:breeding.div){
          temp <- rep(empty.breeding.pen[i], pen_capacity("Sow breeding"))
          div.vec <- c(div.vec, temp)
        }
        breeding.remain.vec <- rep(tail(empty.breeding.pen, 1), breeding.remain)
        dest.pens.vec <- c(div.vec, breeding.remain.vec)
      }

      ## check for NA in empty pens using check.pens()-function
      checkforspace <- check.pens(npens.sow, empty.breeding.pen)

      if (checkforspace==FALSE){
        sow.breeding.nospace <- rbind(sow.breeding.nospace, result$time[1])
      }

      if (checkforspace) {
        ## create the events
        event1 <- event(type = 3,
                        time = result$time[1],
                        node = sow.pens.to.breeding,
                        dest = dest.pens.vec,
                        n = 1,
                        proportion = 0)

        events <- rbind(events, event1)

        ## setting timers for the destination pens
        result$countdown[result$node %in% empty.breeding.pen] <- timers("Breeding")

        ## timer for the pen to be empty before new animals.
        result$countdown[result$node %in% sow.pens.to.breeding] <- timers("Farrowing dt")

        ## calculating the amount of space and fillings we need from other pentypes
        #########
        ## how much room we have in a pen after sows come from farrowing
        filling <- total.breeding()-total.sows.to.breeding

        ## will we get enough animals from buffer to fill the total breeding rate
        buffer.balance <- filling - total.buffer

        new.in <- integer(0)
        ## do we need to get new animals from gilt growing and how many
        if (buffer.balance >= 0){
          new.in <- buffer.balance
        }
        else {
          new.in <- 0
        }

        ## if we have excess amount of animals in buffer, only take what is needed
        ## prioritize animals with smaller timer
        if (total.buffer >= filling){
          if (total.sows.from.buffer >= filling){
            sow.buffer.to.breed <- filling
            gilt.buffer.to.breed <- 0
          }
          else{
            sow.buffer.to.breed <- total.sows.from.buffer
            gilt.buffer.to.breed <- filling - total.sows.from.buffer
          }
        }
        ## if buffer has less animals than we want, take all we can from buffers
        else {
          sow.buffer.to.breed <- total.sows.from.buffer
          gilt.buffer.to.breed <- total.gilts.from.buffer
        }

        ## if we also have animals in breeding buffers, move them to the same pen

        if(sow.buffer.to.breed > 0){

          ## creating a vector of free space in the pens that is left after sows are weaned
          if (breeding.remain != 0){
            last.pen <- tail(empty.breeding.pen, 1)
            free.space <- pen_capacity("Sow breeding") - breeding.remain
            free.space.vec <- rep(last.pen, free.space)
          }
          else {
            free.space.vec <- integer(0)
          }

          pens <- sows.in.buffer[, c("Ssows", "Isows")]

          if(nrow(pens) > 0){
            ## animals per pen
            pigs.per.pen <- rowSums(pens)
            ## sum of all to see if we have enough
            total.sows <- sum(pigs.per.pen)
            ## pen ids of the departure pens
            sow.pen.id <- sows.in.buffer$node
            ## creating vector of departure pens and crop it to the needed length.
            ## Also save the unused ones for determining later which departure pens have gone empty
              all_dep <- rep(sow.pen.id, pigs.per.pen)
              dep <- all_dep[1:sow.buffer.to.breed]
              stopifnot(length(dep) <= length(all_dep))
              if(length(dep) < length(all_dep)){
                rest <- all_dep[seq((sow.buffer.to.breed +1), length(all_dep))]
              }
              else{
                rest <- NULL
              }


            ## find the destination pens and generate the events.
            ## If we have space left after weaning sows, combine the buffer sows with them
            if (total.sows.to.breeding > 0){
              if ((length(free.space.vec) > 0) & (length(free.space.vec) >= sow.buffer.to.breed)){
                dest.pens.vec <- free.space.vec[seq_len(sow.buffer.to.breed)]
              }

              ## when we don't have enough room in the same pens with weaned sows, we need to find new ones
              else {
                npens.sow <- sow.buffer.to.breed %/% pen_capacity("Sow breeding")
                if(sow.buffer.to.breed %% pen_capacity("Sow breeding") != 0) {
                  npens.sow <- npens.sow + 1
                }
                if(sow.buffer.to.breed < pen_capacity("Sow breeding")){
                  npens.sow <- 1
                }

                ## using empty pens function so that we won't try to put sows into the already filled pen(s)
                empty.breeding.pen <- empty.pens(result, "Sow breeding", events)[seq_len(npens.sow)]

                ## creating a destination pen vector for the sows
                breeding.div <- sow.buffer.to.breed %/% pen_capacity("Sow breeding")
                breeding.remain <- sow.buffer.to.breed %% pen_capacity("Sow breeding")

                # creating a vector of destination pens
                div.vec <- NULL

                if (breeding.div == 0){
                  dest.pens.vec <- rep(empty.breeding.pen, sow.buffer.to.breed)
                }
                else {
                  for (i in 1:breeding.div){
                    temp <- rep(empty.breeding.pen[i], pen_capacity("Sow breeding"))
                    div.vec <- c(div.vec, temp)
                  }
                  breeding.remain.vec <- rep(tail(empty.breeding.pen, 1), breeding.remain)
                  dest.pens.vec <- c(div.vec, breeding.remain.vec)
                }
              }

              ## check for NA in empty pens using check.pens()-function
              checkforspace <- check.pens(npens.sow, empty.breeding.pen)

              if (checkforspace==FALSE){
                sow.breeding.nospace <- rbind(sow.breeding.nospace, result$time[1])
              }

              if (checkforspace){
                ## generate the events
                event1 <- event(type = 1,
                                time = result$time[1],
                                node = dep,
                                dest = dest.pens.vec,
                                n = 1,
                                proportion = 0)

                events <- rbind(events, event1)

                result$countdown[result$node %in% empty.breeding.pen] <- timers("Breeding")

                ## set timer only for those breeding pens that end up empty
                emptied_pens <- unique(all_dep[!(all_dep %in% rest)])
                if(length(emptied_pens)>0){
                  result$countdown[result$node %in% emptied_pens] <- timers("Breeding dt")
                }
              }
            }
          }
        }


        ## if we have gilts to inseminate, move them from buffer to gilt breeding pen
        ## (separate from sows)
        if (gilt.buffer.to.breed > 0){
          pens <- gilts.in.buffer[, c("Sgilts", "Igilts")]
          if(nrow(pens) > 0){
            ## animals per pen
            pigs.per.pen <- rowSums(pens)
            ## sum of all to see if we have enough
            total.gilts <- sum(pigs.per.pen)
            ## pen ids of the departure pens
            gilt.pen.id <- gilts.in.buffer$node

            ## creating vector of departure pens and crop it to the needed length.
            ## Also save the unused ones for determining later which departure pens have gone empty
            all_dep <- rep(gilt.pen.id, pigs.per.pen)
            dep <- all_dep[1:gilt.buffer.to.breed]
            stopifnot(length(dep) <= length(all_dep))
            if(length(dep) < length(all_dep)){
              rest <- all_dep[seq((gilt.buffer.to.breed +1), length(all_dep))]
            }
            else{
              rest <- NULL
            }
            ## creating vector of departure pens, crop it to the amoung of animals we want to move
            dep <- rep(gilt.pen.id, pigs.per.pen)[1:gilt.buffer.to.breed]

            ## how many gilt breeding pens are needed
            npens.gilt <- gilt.buffer.to.breed %/% pen_capacity("Gilt breeding")
            if(gilt.buffer.to.breed %% pen_capacity("Gilt breeding") != 0) {
              npens.gilt <- npens.gilt + 1
            }
            if(gilt.buffer.to.breed < pen_capacity("Gilt breeding")){
              npens.gilt <- 1
            }
            ## find empty destination pen(s)
            empty.gilt.breeding <- empty.pens(result, "Gilt breeding", events)[seq_len(npens.gilt)]

            ## creating a destination pen vector for the gilts
            breeding.div <- gilt.buffer.to.breed %/% pen_capacity("Gilt breeding")
            breeding.remain <- gilt.buffer.to.breed %% pen_capacity("Gilt breeding")

            # creating a vector of destination pens
            div.vec <- NULL

            if (breeding.div == 0){
              dest.pens.vec <- rep(empty.gilt.breeding, gilt.buffer.to.breed)
            }
            else {
              for (i in 1:breeding.div){
                temp <- rep(empty.gilt.breeding[i], pen_capacity("Gilt breeding"))
                div.vec <- c(div.vec, temp)
              }
              breeding.remain.vec <- rep(tail(empty.gilt.breeding, 1), breeding.remain)
              dest.pens.vec <- c(div.vec, breeding.remain.vec)
            }

            ## check for NA in empty pens using check.pens()-function
            checkforspace <- check.pens(1, empty.gilt.breeding)

            if (checkforspace==FALSE){
              gilt.breeding.nospace <- rbind(gilt.breeding.nospace, result$time[1])
            }

            if(checkforspace) {
              ## create the events
              event1 <- event(type = 5,
                              time = result$time[1],
                              node = dep,
                              dest = dest.pens.vec,
                              n = 1,
                              proportion = 0)

              events <- rbind(events, event1)

              ## setting timers for the destination pens
              result$countdown[result$node %in% empty.gilt.breeding] <- timers("Breeding")

              ###  only set timer for departure pen if it has gotten empty
              emptied_pens <- unique(all_dep[!(all_dep %in% rest)])
              if(length(emptied_pens) > 0){
                result$countdown[result$node %in% emptied_pens] <- timers("Breeding dt")
              }

              ## if we also need to add animals from gilt growing to breeding buffer,


              ## if we also need to add animals from gilt growing to breeding buffer,
              ## save the booked pens and how many animals we are moving and calculate room left

              if (new.in > 0){
                if (breeding.remain != 0){
                  last.pen <- tail(empty.gilt.breeding, 1)
                  free.space <- pen_capacity("Gilt breeding") - breeding.remain
                  free.space.vec <- rep(last.pen, free.space)
                }
                else {
                  free.space.vec <- integer(0)
                }
                ## create a vector for the free space in gilt breeding after buffer movements
                free.space.vec <- rep(empty.gilt.breeding, free.space)
              }
            }
          }
        }

        ## Try to put new gilts into the same pen with gilts from buffer
        ## do nothing if we don't have animals left to move
        if (new.in > 0){
          ## make a temporary dataframe for growing gilts
          growing.gilts <- result[result$pentype == "Gilt growing" & (result$Sgilts !=0 | result$Igilts !=0), ]
          pens <- growing.gilts[, c("Sgilts", "Igilts")]

          if(nrow(pens) > 0){
            ## animals per pen
            pigs.per.pen <- rowSums(pens)
            ## sum of all to see if we have enough
            total.gilts <- sum(pigs.per.pen)
            ## pen ids of the departure pens
            gilt.pen.id <- growing.gilts$node

            ## check if we are capable of giving as many animals as needed, if not, give what is possible
            if(new.in > total.gilts){
              new.in <- total.gilts
            }

            ## creating vector of departure pens and crop it to the needed length.
            ## Also save the unused ones for determining later which departure pens have gone empty

            all_dep <- rep(gilt.pen.id, pigs.per.pen)
            dep <- all_dep[1:new.in]
            stopifnot(length(dep) <= length(all_dep))
            if(length(dep) < length(all_dep)){
              rest <- all_dep[seq((new.in+1), length(all_dep))]
            }
            else{
              rest <- NULL
            }

            ## creating vector of departure pens and crop it to the length we need
            dep <- rep(gilt.pen.id, pigs.per.pen)[1:new.in]

            ## find the destination pens and generate the events.
            ## different ways of generating destination pens depending on if we are moving gilts from breeding buffer at the same time
            ## if we have also animals coming from gilt breeding buffer, fill in also the extra space left in their destination pens
            if (gilt.buffer.to.breed > 0){
              ## use previously generated destination pen list cropped to the length of need
              ## if we have animals coming from breeding buffer at the same time
              dest.pens.vec <- free.space.vec[new.in]

              ## generate events for new breeding gilts
              event1 <- event(type = 4,
                              time = events$time[1],
                              node = dep,
                              dest = dest.pens.vec,
                              n = 1,
                              proportion = 0)

              events <- rbind(events, event1)

              result$countdown[result$node %in% empty.gilt.breeding] <- timers("Breeding")

            }

            ## if we have no gilts coming from buffer, just move the growing gilts into new pens
            else{
              ## how many pens we need in total
              npens.gilt <- new.in %/% pen_capacity("Gilt breeding")
              if(new.in %% pen_capacity("Gilt breeding") != 0) {
                npens.gilt <- npens.gilt + 1
              }
              if(new.in < pen_capacity("Gilt breeding")){
                npens.gilt <- 1
              }

              dest.pens <- free.pens(result, "Gilt breeding", residual)[seq_len(npens.gilt)]

              ## checking if we got NAs into destination pens
              checkforspace.gilt <- check.pens(npens.gilt, dest.pens)

              ## generate events if we have empty destination pens
              if(checkforspace.gilt) {

                ## create a vector for the destination pens
                dest.pens.vec <- integer(0)
                for(i in 1:length(dest.pens)){
                  temp <- c(rep(dest.pens[i], pen_capacity("Gilt breeding")))
                  dest.pens.vec <- c(dest.pens.vec, temp)
                }
                dest.pens.vec <- dest.pens.vec[seq_len(new.in)]

                ## create events
                event1 <- event(type = 4,
                                time = result$time[1],
                                node = dep,
                                dest = dest.pens.vec,
                                n = 1,
                                proportion = 0)

                events <- rbind(events, event1)
                result$countdown[result$node %in% dest.pens] <- timers("Breeding")
              }
            }
            ## setting a timer for the gilt growing departure pens if they end up empty
            emptied_pens <- unique(all_dep[!(all_dep %in% rest)])
            if(length(emptied_pens) > 0){
              result$countdown[result$node %in% emptied_pens] <- timers("Gilt growing dt")
            }
          }
        }
      }
    }

    residual <- residual[!index, ]
    return(list(result = result,
                residual = residual,
                events = events))
}

##' Abortion
##'
##' In this function we move bred animals from the gestation section
##' to a breeding buffer section or cull them.
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @param abortion.chance the chance of abortion per day per animal
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
abortion <- function(result,
                     residual = NULL,
                     events = NULL,
                     abortion.chance = 0.0029) {

    gest.timer <- timers("Gestation")

    ## collecting the rows that have sows in the first 28 days of their gestation
    early.gestation <- NULL
    early.gestation <- result[(result$pentype == "Sow gestation"| result$pentype == "Gilt gestation") &
                              (result$countdown <= gest.timer & result$countdown >= (gest.timer-28)), ]

    if (nrow(early.gestation) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }
    for (i in seq_len(nrow(early.gestation))) {
        ## go through each row and have a chance of it being pregnant
        row <- early.gestation[i, ]
        type <- ifelse(row$pentype=="Sow gestation", "Sows", "Gilts")
        animal.sum <- row$npigs
        not.pregnant <- sum(rbinom(animal.sum, 1, abortion.chance))
        ## part of the non-pregnant animals go back to breeding and part are removed from herd
        removed <- sum(rbinom(not.pregnant, 1, 0.5))
        back <- not.pregnant - removed
        ## generating events for the animals that are going to be removed
        if (removed > 0) {
            if (type == "Sows") {
                event1 <- event(type = 14,
                                time = result$time[1],
                                node = row$node,
                                dest = 0,
                                n = removed,
                                proportion = 0)

                events <- rbind(events, event1)
            }
            if (type == "Gilts") {
                event1 <- event(type = 15,
                                time = result$time[1],
                                node = row$node,
                                dest = 0,
                                n = removed,
                                proportion = 0)

                events <- rbind(events, event1)
            }
        }
          ##  rest will returnn back from gestation to breeding (type 1 for sows, 5 for gilts)
          if (back > 0){
              timer.set <- FALSE
              ## attempt to merge the animals with highest timer breeding buffer pen if there is any with animals
              if (type == "Sows"){
                  ## collecting buffer pens with animals in them but are not full
                  buffer.pens <- result[result$pentype == "Sow breeding buffer" &
                                          result$npigs != 0 &
                                          result$npigs < result$capacity,]
                  if (nrow(buffer.pens) > 0) {
                      ## order the pens with highest countdown first
                      buffer.pens <- buffer.pens[order(buffer.pens$countdown, decreasing = TRUE),]
                      ## which pens will have enough space for the sows which aborted
                      space <- (pen_capacity("Sow breeding buffer") - (buffer.pens$npigs + back)) >= 0
                      ## nodes of the pens with space
                      space.pen.id <- buffer.pens[space, "node"]
                      ## take the first pen
                      dest <- space.pen.id[1]
                  }
                  ## if there are no breeding buffer pens with animals in them, take an empty pen
                  else{
                      total.buffer.pens <- ceiling(back / pen_capacity("Sow breeding buffer"))
                      empty.sow.buffer.pens <- free.pens(result, "Sow breeding buffer", residual)
                      ## Check that there are enough pens in the buffer and throw error if not
                      stopifnot(length(empty.sow.buffer.pens) >= total.buffer.pens)
                      dest <- empty.sow.buffer.pens[seq_len(total.buffer.pens)]
                      timer.set <- TRUE
                  }
                  ## schedule events
                  event1 <- event(type = 1,
                                  time = result$time[1],
                                  node = row$node,
                                  dest = dest,
                                  n = back,
                                  proportion = 0)

                  events <- rbind(events, event1)
                  ## only set timer for pens that were previously empty
                  if(timer.set){
                      result$countdown[result$node %in% dest] <- timers("Breeding buffer")
                  }

              }
            ## hangle the gilts in similar way: try to merge ~half of the aborted gilts with gilts in breeding buffer with highest timer
            ## if there are no animals within than pentype, find a new pen
            if (type == "Gilts"){
                buffer.pens <- result[result$pentype == "Gilt breeding buffer" &
                                        result$npigs != 0 &
                                        result$npigs < result$capacity,]
                if (nrow(buffer.pens) > 0) {
                    ## collecting buffer pens with animals in them
                    buffer.pens <- result[result$pentype == "Gilt breeding buffer" & (result$npigs!= 0),]
                    ## order the pens with highest countdown first
                    buffer.pens <- buffer.pens[order(buffer.pens$countdown, decreasing = TRUE),]
                    ## which pens will have enough space for the sows which aborted
                    space <- (pen_capacity("Gilt breeding buffer") - (buffer.pens$npigs + back)) >= 0
                    ## nodes of the pens with space
                    space.pen.id <- buffer.pens[space, "node"]
                    ## take the first pen
                    dest <- space.pen.id[1]
                }
                else{
                    total.buffer.pens <- ceiling(back / pen_capacity("Gilt breeding buffer"))
                    empty.gilt.buffer.pens <- free.pens(result, "Gilt breeding buffer", residual)
                    ## Check that there are enough pens in the buffer and throw error if not
                    stopifnot(length(empty.gilt.buffer.pens) >= total.buffer.pens)
                    dest <- empty.gilt.buffer.pens[seq_len(total.buffer.pens)]
                    timer.set <- TRUE
                }
                ## schedule the events
                event1 <- event(type = 5,
                                time = result$time[1],
                                node = row$node,
                                dest = dest,
                                n = back,
                                proportion = 0)

                events <- rbind(events, event1)

                if(timer.set){
                  result$countdown[result$node %in% dest] <- timers("Breeding buffer")
                }
            }
        }
    }

    return(list(result = result,
                residual = residual,
                events = events))
}

##' growing_to_finishing
##'
##' In this function we move pigs from the grower to the finisher
##' section and possibly to a from buffer sections in those two
##' units.
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @param fraction_to_buffer the proportion of animals to be removed to buffer instead of finishing
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
growing_to_finishing <- function(result,
                                 residual = NULL,
                                 events = NULL,
                                 fraction_to_buffer = 0.1) {

    ## Index of growing pens that have pigs
    index <- residual$pentype == "Growing"

    ## If there is nothing to do then just get out now
    if(sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }

    ## how many animals we have to move in total
    pigs.per.pen <- residual[index, "npigs"]
    total.piglets.to.move <-sum(pigs.per.pen)

    ## expected number to buffer
    weaner.to.buffer <- round(total.piglets.to.move * fraction_to_buffer)

    ## Number of buffer pens required given the capacity
    total.buffer.pens <- ceiling(weaner.to.buffer / pen_capacity("Growing buffer"))

    ## Finding empty buffer pens
    empty.growing.buffer.pens <- free.pens(result,
                                           "Growing buffer",
                                           residual)[seq_len(total.buffer.pens)]

    ## check that you don't index too many pens in the freepens vector
    stopifnot(all(!is.na(empty.growing.buffer.pens)))

    ## generating events for moving pigs from growing to finishing
    ## finding empty finishing pens from one section
    empty.finishing.pens <- sectioning_multi(result, "Finishing", 0)[seq_len(sum(index))]

    ## check that you don't index too many pens in the
    ## sectioning_multi vector
    stopifnot(all(!is.na(empty.finishing.pens)))

    ## pen ids of the departure pens
    growing.pen.id <- residual[index, "node"]

    for(i in seq_len(sum(index))) {

        ## A single source pen
        node <- growing.pen.id[i]

        ## n buffer pens and one finisher pen:
        dest <- c(empty.growing.buffer.pens,
                  empty.finishing.pens[i])

        ## split the fraction to buffer over the buffer pens and the
        ## rest in the finisher pen:
        p <- c(rep(fraction_to_buffer / total.buffer.pens, total.buffer.pens),
               1 - fraction_to_buffer)

        ## Scale the p to the order to the dest pen ID
        scaled_p <- scale_p(dest, p)

        ## set the type for the events:
        type <- c(rep(18, total.buffer.pens), 9)

        ## create the events:
        event1 <- event(type = type,
                        time = result$time[1],
                        node = node,
                        dest = dest,
                        n = 0,
                        proportion = scaled_p)
        events <- rbind(events, event1)

        ## set the emptied pen timer
        result$countdown[result$node %in% node] <- timers("Growing dt")
    }

    ## Set the timers for the destination pens
    result$countdown[result$node %in% empty.growing.buffer.pens] <- timers("Growing buffer")
    result$countdown[result$node %in% empty.finishing.pens] <- timers("Finishing")

    ## Drop processed rows
    residual <- residual[!index, ]

    return(list(result = result,
                residual = residual,
                events = events))
}
##' Moving animals from growing buffer to finishing unit
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
growing_buffer_to_finishing <- function(result,
                                        residual = NULL,
                                        events = NULL) {

    ## Index of growing pens that have pigs
    index <- residual$pentype == "Growing buffer"

    ## If there is nothing to do then just get out now
    if(sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }

    ## generating events for moving pigs from growing to finishing
    ## finding empty finishing pens from one section
    empty.finishing.pens <- sectioning_multi(result, "Finishing", sum(index))

    ## what sections are these
    efp <- result[result$node %in% empty.finishing.pens, "section"]

    ## what pens are in those sections
    efp <- result[result$pentype == "Finishing" &
                  result$section %in% efp, ]

    unempty <- efp[efp$npigs > 0,]

    ## What are the timers in these sections?
    cdown <- tapply(unempty$countdown,
                    unempty$section,
                    max, na.rm = TRUE)

    ## prioritize sections that have some animals in them. If we only have empty sections,
    ## then move buffer animals to an empty section
    if (any(cdown != 99999)){
    cdown <- cdown[!cdown == 99999]
    }

    ## which is the highest?
    maxdown_index <- which.max(cdown)
    cdown_section <- names(cdown)[maxdown_index]
    cdown_timer <- cdown[maxdown_index]

    ## The nodes in this max downtime section:
    efp <- efp[efp$section == cdown_section, "node"]

    ## Just keep those in the empty vector that were in this section:
    empty.finishing.pens <- empty.finishing.pens[empty.finishing.pens %in% efp]

    ##sanity check to make sure we are not trying to move more
    ## pens that fit into finishing unit.
    if (!(sum(index) <= length(empty.finishing.pens))){
      browser()
    }
    stopifnot(sum(index) <= length(empty.finishing.pens))

    ## Move to finisher
    source_pens <- residual[index, "node"]
    dest.pens.vec <- empty.finishing.pens[seq_len(length(source_pens))]

    ## create the events
    event1 <- event(type = 9,
                    time = result$time[1],
                    node = source_pens,
                    dest = dest.pens.vec,
                    n = 0,
                    proportion = 1)

    events <- rbind(events, event1)

    ## timers for destination pens to same as their new section
    result$countdown[result$node %in% dest.pens.vec] <- cdown_timer

    ## timer for the pen to be empty before new animals.
    result$countdown[result$node %in% source_pens] <- timers("Growing buffer dt")
    residual <- residual[!index, ]

    return(list(result = result,
                residual = residual,
                events = events))
}

##' culling animals from sow and gilt breeding buffers if they have
##' been there for too long (timer hits 0)
##'
##' @param result The result from the previous day or modified result
##'     from the current day.
##' @param residual The rows with timer == 0 that still need to be
##'     processed. The default is NULL if there is nothing to be done
##' @param events The events that have already been generated. The
##'     default is NULL in the case that no events have yet been
##'     generated for this day
##' @return A list of the (modified) result, the new residual and the
##'     new events appended to the events you fed to the function
buffer_culling <- function(result,
                           residual = NULL,
                           events = NULL) {

    index <- (residual$pentype == "Sow breeding buffer" | residual$pentype == "Gilt breeding buffer")

    ## If there is nothing to do exit
    if (sum(index) == 0) {
        return(list(result = result,
                    residual = residual,
                    events = events))
    }
    pens.to.remove <- residual[index,]
    pen.id <- pens.to.remove$node
    ## define the type for the exit event based on the pentype
    type <- rep(14, length(pen.id))
    rem.type <- pens.to.remove$pentype[match(pen.id,
                                                pens.to.remove$node)] == "Gilt breeding buffer"
    type[rem.type] <- 15

    ## Remove all animals from these pens
    event1 <- event(type = type,
                    time = result$time[1],
                    node = pen.id,
                    dest = 0,
                    n = 0,
                    proportion = 1)
    events <- rbind(events, event1)

    result$countdown[result$node %in% pen.id] <- timers("Breeding dt")

    ## Drop rows from residual
    residual <- residual[!index, ]

    ## return result
    return(list(result = result,
                residual = residual,
                events = events))
}
