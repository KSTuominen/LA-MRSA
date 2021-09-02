#' sow.count
#'
#' A functiont to count the total number of sows and gilts in the herd on a
#' given day
#'
#'
#' @param x A dataframe which is the output from the trajectory
#'     of the MRSA model from 1 day.
#' @return An integer
sow.count <- function(x){
    ## Check that the dataframe is the correct class

  ## check below commented out for. Doesn't work with growing gilt movement because result doesn't have "MRSA_single_step" at that point
  ## how to fix this?
  stopifnot(identical(class(x), c("data.frame", "MRSA_single_step_trajectory")))

    sum(x$Ssows, x$Isows, x$Sgilts[x$pentype != "Gilt growing"], x$Igilts[x$pentype != "Gilt growing"])
}

#' comp.sum
#'
#' how many animals there are per compartment.
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day.
#' @param comptype A string
#' @return sum of animals in that compartment in all units

comp.sum <- function(x, comptype = c("Ssows",
                                    "Isows",
                                    "Sgilts",
                                    "Igilts",
                                    "Spiglets",
                                    "Ipiglets",
                                    "Sgrowers",
                                    "Igrowers",
                                    "Sfinish",
                                    "Ifinish")){
  comptype <- match.arg(comptype)
  sum(x[ , comptype])


}

## different functions for finding empty pens
##################
## free.pens -function
########
# when we move animals one by one from several different pens to one pen (moving many animals one by one)
# from farrowing to breeding
# from farrowing to gilt growing
# from growing to growing buffer
# from gilt growing to breeding
# does not check if there are already events scheduled for the destination pens (opposite to empty.pens -function which wants to limit this)
# does not check for sections as currently the target pens are continuous flow units

## empty.pens -function
#########
# when we move many animals one by one but not (necessarily) to same pen
# from sow/gilt breeding buffer to breeding
# initializing herd with gilts to gilt breeding
# does not check for sections as currently the target pens are continuous flow units
# checks if there if there are already events scheduled for the destioation pens and removes these from the returned pen list


## sectioning -function
########
# animals moved one by one (n =1) from one pen moved into several pens
# can determine if we want the destination section to be fully empty (all-in-all-out) or just have x number of pens
# sows/gilts from gestation to farrowing
# sows/gilts from breeding to gestation
# checks if there if there are already events scheduled for the destionation pens and removes these from the returned pen list

## sectioning_multi -function
#######
# animals move as group (either n > 1 or proportion) from several pens to several destination pens.
# can determine if we want the destination section to be fully empty (all-in-all-out) or just have x number of pens
# farrowing to growing
# growing to finishing
# growing buffer to finishing
# unlike sectioning function, does not check if there are already events scheduled for the destination pens as this caused issues when moving animals from several pens)

#' empty.pens
#'
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day.
#' @param pentype A string'
#' @param events The queued events for this timestep
#' @return a vector of pen ids
empty.pens <- function(x,
                       pentype = c("Sow breeding",
                                   "Gilt breeding",
                                   "Sow breeding buffer",
                                   "Gilt breeding buffer",
                                   "Sow gestation",
                                   "Gilt gestation",
                                   "Farrowing",
                                   "Growing",
                                   "Growing buffer",
                                   "Finishing",
                                   "Finishing buffer",
                                   "Gilt growing"),
                       events = NULL){
    ## Check arguments
    pentype <- match.arg(pentype)
    stopifnot(identical(class(x), c("data.frame", "MRSA_single_step_trajectory")))

    ## Which pens are the class you are interested in:
    pens_of_class <- x[, "pentype"] == pentype

    ## Which pens are empty:
    emptypens <- x$npigs == 0

    ## which pens have timer 99999 meaning that they are truly available
    counter.check <- x[, "countdown"] == "99999"

    ## Which node ids are both empty and your class of interest. Note
    ## the use of & not &&. The first returns a vector equal to the
    ## length of the vectors, which && only uses the first value in
    ## each vector.
    x <- x[ ,"node"][pens_of_class & emptypens & counter.check]

    ## Remove those ids that are in the booked list
    ## how to return those that are in the same section? Should we assume there are enough pens?

    x[!(x %in% events$dest)]
}

#' free.pens
#'
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day.
#' @param y dataframe with the pens which need events
#' @param pentype A string'
#' @return a vector of pen ids
free.pens <- function(x,
                       pentype = c("Sow breeding",
                                   "Gilt breeding",
                                   "Sow breeding buffer",
                                   "Gilt breeding buffer",
                                   "Sow gestation",
                                   "Gilt gestation",
                                   "Farrowing",
                                   "Growing",
                                   "Growing buffer",
                                   "Finishing",
                                   "Finishing buffer",
                                   "Gilt growing"),
                      y){
  ## Check arguments
  pentype <- match.arg(pentype)
  stopifnot(identical(class(x), c("data.frame", "MRSA_single_step_trajectory")))

  ## Which pens are the class you are interested in:
  pens_of_class <- x[, "pentype"] == pentype

  ## Which node ids have room and your class of interest.

  emptypens <- x$npigs == 0

  if(all(!(emptypens & pens_of_class))) return(NULL)

  if(pentype == "Gilt growing"){
      counter.check <- x[, "countdown"] == "99999"
      if(all(!(emptypens & pens_of_class & counter.check))) return(NULL)
      x <- x[ ,"node"][pens_of_class & emptypens & counter.check]
      return(x)
  }

  x[ ,"node"][pens_of_class & emptypens]
}

#' SI
#'
#' function for returning back both susceptible and infected compartments in spesific node
#' instead of having to type then separately
#' subsetting needs to be added
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day.
#' @param SI.comp A string
#' @param node value that can be coersed to an integer greater than
#'     or equal to 0
#' @return an integer, sum of animals in that compartment in all units

SI <- function(x,
               SI.comp = c("Sows", "Gilts", "Piglets", "Growers", "Finish"),
               node) {

  ## check that the compartment names and pentype are default values
  SI.comp <- match.arg(SI.comp)

  ## changing the general compartment title into the susceptible and infected one
  SI.comp <- switch(SI.comp,
                        'Sows' = list("Ssows", "Isows"),
                        'Gilts' = list("Sgilts", "Igilts"),
                        'Piglets' = list("Spiglets", "Ipiglets"),
                        'Growers' = list("Sgrowers", "Igrowers"),
                        'Finish' = list("Sfinish", "Ifinish"))

  SI.sum(x, SI.comp[[1]], SI.comp[[2]], node) ## send the parameters to SI sum to get the sum of the animals in these compartments

}

#' SI.sum
#'
#'function to sum up animals from the 2 compartments combined in SI
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day, passed from SI-function.
#' @param Scomp A string
#' @param Icomp A string
#' @param node value that can be coersed to an integer greater than
#'     or equal to 0
#' @return an integer, sum of animals in that compartment in all units

SI.sum <- function(x, Scomp, Icomp, node){
  sum(x[ ,Scomp][node], x[ ,Icomp][node])
}

#' SI.pentype
#'
#' function for returning back the sum of of animals both susceptible and infected compartments
#' from all pens with given pentype
#' instead of having to type then separately
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day.
#' @param SI.comp A string
#' @param pentype A string'
#' @return an integer, sum of animals
SI.pentype <- function(x,
                       SI.comp = c("Sows",
                                   "Gilts",
                                   "Piglets",
                                   "Growers",
                                   "Finish"),
                       pentype = c("Sow breeding",
                                   "Gilt breeding",
                                   "Sow breeding buffer",
                                   "Gilt breeding buffer",
                                   "Sow gestation",
                                   "Gilt gestation",
                                   "Farrowing",
                                   "Growing",
                                   "Growing buffer",
                                   "Finishing",
                                   "Finishing buffer",
                                   "Gilt growing")){

    ## check that the compartment names and pentype are default values
    SI.comp <- match.arg(SI.comp)
    pentype <- match.arg(pentype)

    ## changing the general compartment title into the susceptible and infected one
    SI.comp <- switch(SI.comp,
                      'Sows' = list("Ssows", "Isows"),
                      'Gilts' = list("Sgilts", "Igilts"),
                      'Piglets' = list("Spiglets", "Ipiglets"),
                      'Growers' = list("Sgrowers", "Igrowers"),
                      'Finish' = list("Sfinish", "Ifinish"))

    ## indexing only the rows that match the desired pentype
    rows <- x[x$pentype == pentype,]

    SI.sum.pentype(rows, SI.comp[[1]], SI.comp[[2]]) ## send the parameters to SI sum to get the sum of the animals in these compartments

}

#' SI.sum.pentype
#'
#' function to sum up animals from the 2 compartments combined in SI.pentype
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day, passed from SI-function.
#' @param Scomp A string
#' @param Icomp A string
#' @param type A string
#' @return an integer, sum of animals

SI.sum.pentype <- function(rows, Scomp, Icomp){
    sum(rows[,Scomp], rows[, Icomp])
}

#' littersize
#'
#' sampling the littersizes with poisson distribution
#' poisson is used instead of binomial distribution to better reflect the chance of getting littersizes outsife of SD

#' @param litters an integer for the amount of litters needed
#' @return a vector of integers

littersize <- function(litters,
                       size = 14.8){
    rpois(litters, size)
}

#' pen_capacity
#'
#' unction to return pensize of a pentype
#' @param pentype A character vector of pentypes
#' @return an numeric vector of the capacities of these pentypes
pen_capacity <- function(pentype) {
    ## Check arguments
    capacity <- data.frame(pentype = c("Sow breeding",
                                       "Gilt breeding",
                                       "Sow breeding buffer",
                                       "Gilt breeding buffer",
                                       "Sow gestation",
                                       "Gilt gestation",
                                       "Farrowing",
                                       "Growing",
                                       "Growing buffer",
                                       "Finishing",
                                       "Finishing buffer",
                                       "Gilt growing"),
                           capacity = c(22,
                                        22,
                                        22,
                                        22,
                                        11,
                                        11,
                                        100,
                                        12,
                                        12,
                                        10,
                                        10,
                                        10))
    stopifnot(all(pentype %in% capacity$pentype))
    capacity$capacity[match(pentype, capacity$pentype)]
}

#' timers
#'
#' Function to return appropriate countdown timer for the given pentype
#' @param pentype A string, dt = downtime
#' @return an integer for setting countdown timer

timers <- function(pentype = c("Breeding",
                               "Breeding dt",
                               "Breeding buffer",
                               "Gestation",
                               "Gestation dt",
                               "Farrowing",
                               "Farrowing dt",
                               "Growing",
                               "Growing dt",
                               "Growing buffer",
                               "Growing buffer dt",
                               "Finishing",
                               "Finishing dt",
                               "Finishing buffer",
                               "Finishing buffer dt",
                               "Gilt growing",
                               ## "Gilt growing finish",
                               "Gilt growing dt")){
  ## Check arguments
  pentype <- match.arg(pentype)

  ## FIXME buffer timers and gilt growing dt need to be thought through
  ## breeding buffer times are now arbitrary set into numbers after when the animals are slaughtered if they are still in the pen when timer hits 0
  ## Gilt growing finish timer, if it's exceeded the animals go to slaughter calculated by:
  ## first gestation around 9 mo age, so 270d - 56d - 35d and rounded up for extra days
  timerlist <- switch(pentype,
                        'Breeding' = 32L,
                        'Breeding dt' = 2L,
                        'Breeding buffer' = 56L,
                        'Gestation' = 88L,
                        'Gestation dt' = 1L,
                        'Farrowing' = 35L,
                        'Farrowing dt' = 2L,
                        'Growing' = 56L,
                        'Growing dt' = 1L,
                        'Growing buffer' = 23L,
                        'Growing buffer dt' = 4L,
                        'Finishing' = 99L,
                        'Finishing dt' = 6L,
                        'Finishing buffer' = 27L,
                        'Finishing buffer dt' = 6L,
                        'Gilt growing' = 216L, ## previously 56 when separate age stages and 'Gilt growing finish' = 160L,
                        'Gilt growing dt' = 6L)
    timerlist
}

#' mortality
#'
#' Mortality function to handle animal mortality outside siminf transitions if needed
#' using poisson distribution to determine how many animal are dead per pen per day
#' return null if all are alive
#' currently has only pentypes for weaning piglets (to gilt grower and growing units)
#' @param x A dataframe which is the output from the trajectory of
#'     the MRSA model from 1 day.
#' @param pentype A string, what type of pens we want to look at
#' @param comp string, compartments which animals we to calculate, used to call SI-function
#' @param prob double, probability of death for one pig per day
#'
#' @return exit events for dead animals

mortality <- function(x,
                     pentype = c("Sow breeding",
                                 "Gilt breeding",
                                 "Sow breeding buffer",
                                 "Gilt breeding buffer",
                                 "Sow gestation",
                                 "Gilt gestation",
                                 "Farrowing",
                                 "Growing",
                                 "Growing buffer",
                                 "Finishing",
                                 "Finishing buffer",
                                 "Gilt growing"),
                     comp = c("Sows", "Gilts", "Piglets", "Growers", "Finish"),
                     prob){

    comp <- match.arg(comp)
    pentype <- match.arg(pentype)
    event_fun <- NULL

    if(pentype == "Farrowing"){
      index <- (x$pentype == pentype) & (x$npigs !=0) & (x$countdown > 0)
    }
    else{
      index <- (x$pentype == pentype) & (x$npigs !=0) & (x$countdown > 0)
    }

    pen.id <- x$node[index]
    for (i in 1:length(pen.id)){
      animals <- SI(x, comp, pen.id[i])
      remove.sum <- sum(rpois(animals, prob))

      if (remove.sum > 0){
        dep <- x$node[index]
        event1 <- event(type = 12,
                        time = x$time[1],
                        node = pen.id[i],
                        dest = 0,
                        n = remove.sum,
                        proportion = 0)

        event_fun <- rbind(event_fun, event1)
      }

      else{
        event_fun <- rbind(event_fun, NULL)
      }
    }
    return(event_fun)
}

#' check.pens
#'
#' special function to check if there is space in farrowing pens

#' @param animals integer of the amount of animals in departure pen
#' @param dest.pen vector of pen ids
#' @return boolean

check.pens<- function(animals, dest.pen){

  ## check if you have enough empty pens
  empty.pens <- !is.na(dest.pen)
  empty.pens <- empty.pens && length(dest.pen) >= animals
  return(empty.pens)
}

#' sectioning
#'
#' function for returning a vector of pen ids in given pentype that have enough empty pens in one section
#' aimed for cases where we move only one pen to new empty pen(s) (breeding -> gestation, gestation -> farrowing)
#' also an option for making sure that the section is all in - all out by checking it's empty
#' @param x A dataframe which is the output from the trajectory of the MRSA model from 1 day.
#' @param pentype A string
#' @param pen.total integer for required amount of pens, if it's 0, we want to have a whole empty section
#' @param events The queued events for this timestep
#' @return vector of pen ids in the given pentype that have enough empty pens
sectioning <- function(x,
                       pentype = c("Sow breeding",
                                   "Gilt breeding",
                                   "Sow breeding buffer",
                                   "Gilt breeding buffer",
                                   "Sow gestation",
                                   "Gilt gestation",
                                   "Farrowing",
                                   "Growing",
                                   "Growing buffer",
                                   "Finishing",
                                   "Finishing buffer",
                                   "Gilt growing"),
                       pen.total = 0,
                       events = NULL){
  ## Check arguments
  pentype <- match.arg(pentype)
  stopifnot(identical(class(x), c("data.frame", "MRSA_single_step_trajectory")))

  section.pen.sum <- NULL
  free.sections <- NULL
  section.space <- NULL
  emptysection <- NULL

  ## when we require the whole section to be empty
  if(pen.total == 0){
      ## going through each individual section number in given pentype
      ## and checking if all pens in that section are empty by using all & rowsums
      for(i in unique(x$section[x$pentype == pentype])){
          emptycheck <- all(x[x$pentype == pentype & x$section == i, "npigs"] == 0)
          if(emptycheck == TRUE){
              emptysection <- c(emptysection, i)
          }
      }
      ## which pen ids match to the first section that is empty
      section.space <- x$section %in% emptysection[1]
      ## which pens match the given pentype
      pens_of_class <- x[, "pentype"] == pentype

      ## which pens have the desired pentype and belong to a section that is empty
      x <- x[ ,"node"][section.space & pens_of_class]

      ## Remove those ids that are in the booked list
      x[!(x %in% events$dest)]

      return(x)
  }

  ## if we don't require the whole section to be empty, we go through which section have enough free pens
  else{
      ## how many free pens we have in each section for given pentype
      ## going through each unique section number per pentype
      for(i in unique(x$section[x$pentype==pentype])){
          empty.pen.counter <- 0L

          ## loop over each row per section for the given pentype and check if the pen is empty
          ## if the pen is empty, add one to counter
          for(j in seq_len(nrow(x))){
              row <- x[j, ]
              if(row$section == i & row$pentype == pentype){
                  if(row[, "npigs"] == 0){
                      empty.pen.counter <- empty.pen.counter + 1
                  }
              }
          }
          ## save the counter for each section
          section.pen.sum[i] <- empty.pen.counter
      }
      ## section numbers that have enough empty pens
      free.sections <- which(section.pen.sum >= pen.total)
      ## which pen ids match the sections that have enough pens
      section.space <- x$section %in% free.sections
  }

  ## which pens match the given pentype
  pens_of_class <- x[, "pentype"] == pentype

  ## which pens have timer 99999 meaning that they are truly available
  counter.check <- x[, "countdown"] == "99999"

  ## all empty pens in the farm
  emptypens <- x[, "npigs"] == 0

  ##  which pen ids have the appropriate section, have the desired pentype and are empty
  x <- x[ ,"node"][section.space & pens_of_class & emptypens & counter.check]

  ## Remove those ids that are in the booked list
  x[!(x %in% events$dest)]
}

#' sectioning_multi
#'
#'Function for returning a vector of pen ids in given pentype that have enough empty pens in one section,
#' aimed for cases where we move several pens to several destination pens (farrowing -> growing, growing -> finishing)
#' Doesn't remove booked pens inside the function (would cause problems when moving several pens from same section)
#' also an option for making sure that the section is all in - all out by checking it's empty
#' @param x A dataframe which is the output from the trajectory of the MRSA model from 1 day.
#' @param pentype A string
#' @param pen.total integer for required amount of pens
#' @return vector of pen ids in the given pentype that have enough empty pens
sectioning_multi <- function(x,
                       pentype = c("Sow breeding",
                                   "Gilt breeding",
                                   "Sow breeding buffer",
                                   "Gilt breeding buffer",
                                   "Sow gestation",
                                   "Gilt gestation",
                                   "Farrowing",
                                   "Growing",
                                   "Growing buffer",
                                   "Finishing",
                                   "Finishing buffer",
                                   "Gilt growing"),
                       pen.total=0){
  ## Check arguments
  pentype <- match.arg(pentype)
  stopifnot(identical(class(x), c("data.frame", "MRSA_single_step_trajectory")))

  section.pen.sum <- NULL
  free.sections <- NULL
  section.space <- NULL
  emptysection <- NULL


  ## when we require the whole section to be empty
  if(pen.total == 0){
      ## going through each individual section number in given pentype
      ## and checking if all pens in that section are empty by using all & rowsums
      for(i in unique(x$section[x$pentype == pentype])){
          emptycheck <- all(x[x$pentype == pentype & x$section == i, "npigs"] == 0)
          if(emptycheck == TRUE){
            emptysection <- c(emptysection, i)
          }
      }
      ## which pen ids match to the first section that is empty
      section.space <- x$section %in% emptysection[1]
      ## which pens match the given pentype
      pens_of_class <- x[, "pentype"] == pentype

      ## which pens have the desired pentype and belong to a section that is empty
      x <- x[ ,"node"][section.space & pens_of_class]

      return(x)
  }

  ## if we don't require the whole section to be empty, we go through which section have enough free pens
  else{
      ## how many free pens we have in each section for given pentype
      ## going through each unique section number per pentype
      for(i in unique(x$section[x$pentype==pentype])){
        empty.pen.counter <- 0L

        ## loop over each row per section for the given pentype and check if the pen is empty
        ## if the pen is empty, add one to counter
        for(j in seq_len(nrow(x))){
          row <- x[j, ]
          if(row$section == i & row$pentype == pentype){
            if(row[, "npigs"] == 0){
              empty.pen.counter <- empty.pen.counter + 1
            }
          }
        }
        ## save the counter for each section
        section.pen.sum[i] <- empty.pen.counter
      }

      ## section numbers that have enough empty pens
      free.sections <- which(section.pen.sum >= pen.total)

      ## which pen ids match the sections that have enough pens
      section.space <- x$section %in% free.sections
  }


  ## which pens match the given pentype
  pens_of_class <- x[, "pentype"] == pentype

  ## which pens have timer 99999 meaning that they are truly available
  counter.check <- x[, "countdown"] == "99999"

  ## all empty pens in the farm
  emptypens <- x[, "npigs"] == 0

  ##  which pen ids have the appropriate section, have the desired pentype and are empty
  x <- x[ ,"node"][section.space & pens_of_class & emptypens & counter.check]

}

#' pregnancy
#'
#' Function to sample how many animals in given pen are not pregnant. Sums the amount of animals
#' of given category from given node
#'
#' @param x A dataframe which is the output from the trajectory of the MRSA model from 1 day.
#' @param pigs an integer of the amount of animals in pen
#' @param node value that can be coersed to an integer greater than or equal to 0
#' @return an integer representing the amount of animals that are not pregnant in the given pen

pregnancy <- function(pigs){

    ## mean omlöpningsprocent is 5.5 % in Winpig at 2019
    sum(rbinom(pigs, 1, 0.055))
}

#'  sow.culling
#'
#'  Function to calculate how many sows to remove at farrowing
#'  and which pens to take them from
#'
#' @param pen.id integer, vector of pen ids that are being weaned
#'
#' @return vector of integers (pen ids)
sow.culling <- function(pen.id,
                        culling_rate = 0.184) {

    ## Find the mean of the poisson distribution that you will sample
    ## from to determine the number of sows to cull
    removed.mean <- culling_rate * length(pen.id)

    ## Sample from that distribution to determine the number of sows
    removed.tot <- rpois(1, removed.mean)
    if(removed.tot > length(pen.id)) {
        removed.tot <- length(pen.id)
    }
    if(removed.tot == 0) return(NULL)

    ## Return the pens
    sample(pen.id, removed.tot)
}

#' total.breeding
#'
#' Calculate the farrowing rate
#' the average farrowing rate for a sow/piglet according to Winpig
#' 2019 is 86.8 %
#'
#' @param farrowing.rate the proportion of sows that are bred that
#'     successfully farrow
#' @param farrowing.pen.capacity the number of pens in a farrowing
#'     room
#' @return The number to breed in a week
total.breeding <- function(farrowing.rate = 0.868,
                           farrowing.pen.capacity = 22) {
    stopifnot(farrowing.rate <= 1 & farrowing.rate > 0)
    round(farrowing.pen.capacity/farrowing.rate)
}


#' max.size
#'
#' Maximum number of PENS in section
#' @param section String, name of the section
#' @return integer

max.size <- function(section){

  ## Check arguments
  capacity <- data.frame(section = c("Sow breeding",
                                     "Gilt breeding",
                                     "Sow breeding buffer",
                                     "Gilt breeding buffer",
                                     "Sow gestation",
                                     "Gilt gestation",
                                     "Farrowing",
                                     "Growing",
                                     "Growing buffer",
                                     "Finishing",
                                     "Finishing buffer",
                                     "Gilt growing"),
                         capacity = c(8,
                                      8,
                                      5,
                                      20,
                                      30,
                                      30,
                                      25,
                                      25,
                                      25,
                                      25,
                                      25,
                                      25),
                         stringsAsFactors = FALSE)
    stopifnot(all(section %in% capacity$section))
  capacity$capacity[match(section, capacity$section)]
}

#' clean trajectory
#'
#' format the daily trajectory into more readable
#' @param x dataframe
#' @return dataframe

clean_trajectory <- function(x) {

    ## Check that the trajectory is just from one day:
    stopifnot(all(table(x$node) == 1))

    ## Label the pentype
    x$pentype <- factor(x$pentype,
                        labels = c("Sow breeding",
                                   "Sow breeding buffer",
                                   "Gilt breeding",
                                   "Gilt breeding buffer",
                                   "Sow gestation",
                                   "Gilt gestation",
                                   "Farrowing",
                                   "Growing",
                                   "Growing buffer",
                                   "Finishing",
                                   "Finishing buffer",
                                   "Gilt growing"))

    ## Add the number of animals and capacity to the result object
    ## then we don't have to keep calculating it.
    x$capacity <- pen_capacity(x$pentype)
    x$npigs <- nanimals(x)

    class(x) <- c(class(x), "MRSA_single_step_trajectory")
    x
}

## Utility function to coerce the data.frame to a transposed u matrix.
as_u_matrix <- function(x) {
    stopifnot("MRSA_single_step_trajectory" %in% class(x))
    x  <- x[, c("Ssows", "Isows", "Sgilts", "Igilts", "Spiglets",
                "Ipiglets", "Sgrowers", "Igrowers", "Sfinish",
                "Ifinish", "pentype", "section", "countdown")]
    n_col <- ncol(x)
    n_row <- nrow(x)
    lbl <- colnames(x)
    x <- as.integer(t(data.matrix(x)))
    attributes(x) <- NULL
    dim(x) <- c(n_col, n_row)
    rownames(x) <- lbl
    x
}

## Utility function to coerce the data.frame to a transposed v matrix.
as_v_matrix <- function(x) {
    stopifnot("MRSA_single_step_trajectory" %in% class(x))
    x  <- x[, c("phi"), drop = FALSE]
    n_col <- ncol(x)
    n_row <- nrow(x)
    lbl <- colnames(x)
    x <- as.double(t(data.matrix(x)))
    attributes(x) <- NULL
    dim(x) <- c(n_col, n_row)
    rownames(x) <- lbl
    x
}

## initialize the farm structure. This generates also finishing buffer pens which are not currently in use.
initialize_u0 <- function(sow.breeding.pens.persection = 10,
                          sow.breeding.pens.section = 1,

                          gilt.breeding.pens.persection = 8,
                          gilt.breeding.pens.section = 1,

                          sow.breeding.buffer.pens.persection = 5,
                          sow.breeding.buffer.pens.section = 1,

                          gilt.breeding.buffer.pens.persection = 10,
                          gilt.breeding.buffer.pens.section = 1,

                          sow.gestation.pens.persection = 35,
                          sow.gestation.pens.section = 1,

                          gilt.gestation.pens.persection = 30,
                          gilt.gestation.pens.section = 1,

                          farrowing.pens.persection = 26,
                          farrowing.pens.sections = 6,

                          growing.pens.persection = 26,
                          growing.pens.sections = 10,

                          growing.buffer.pens.persection = 26,
                          growing.buffer.pens.sections = 1,

                          finishing.pens.persection = 30,
                          finishing.pens.sections = 18,

                          finishing.buffer.pens.persection = 25,
                          finishing.buffer.pens.sections = 2,

                          growing.gilt.pens.persection = 25,
                          growing.gilt.pens.sections = 1
                          ) {
    ## usually there is only 1 breeding and 1 gestation section
    sow.breed.section.vec <- rep(seq_len(sow.breeding.pens.section),
                                 each = sow.breeding.pens.persection)

    gilt.breed.section.vec <- rep(seq_len(gilt.breeding.pens.section),
                                 each = gilt.breeding.pens.persection)

    sow.breed.buffer.section.vec <- rep(seq_len(sow.breeding.buffer.pens.section),
                                        each = sow.breeding.buffer.pens.persection)

    gilt.breed.buffer.section.vec <- rep(seq_len(gilt.breeding.buffer.pens.section),
                                         each = gilt.breeding.buffer.pens.persection)

    ## gestation pens are normally smaller than breeding pens

    sow.gest.section.vec <- rep(seq_len(sow.gestation.pens.section),
                                each = sow.gestation.pens.persection)

    gilt.gest.section.vec <- rep(seq_len(gilt.gestation.pens.section),
                                 each = gilt.gestation.pens.persection)

    ## Farrowing pens
    farrowing.section.vec <- rep(seq_len(farrowing.pens.sections),
                                 each = farrowing.pens.persection)

    ## Growing pens
    growing.section.vec <- rep(seq_len(growing.pens.sections),
                               each = growing.pens.persection)

    ## Growing buffer pens
    growing.buffer.section.vec <- rep(seq_len(growing.buffer.pens.sections),
                                      each = growing.buffer.pens.persection)

    ## Finishing pens
    finishing.section.vec <- rep(seq_len(finishing.pens.sections),
                                 each = finishing.pens.persection)

    ## Finishing buffer pens
    finishing.buffer.section.vec <- rep(seq_len(finishing.buffer.pens.sections),
                                        each = finishing.buffer.pens.persection)

    ## Growing gilt section
    growing.gilt.vec <- rep(seq_len(growing.gilt.pens.sections),
                            each = growing.gilt.pens.persection)

    ## concatenating all the section number vectors together
    section<- c(sow.breed.section.vec,
                sow.breed.buffer.section.vec,
                gilt.breed.section.vec,
                gilt.breed.buffer.section.vec,
                sow.gest.section.vec,
                gilt.gest.section.vec,
                farrowing.section.vec,
                growing.section.vec,
                growing.buffer.section.vec,
                finishing.section.vec,
                finishing.buffer.section.vec,
                growing.gilt.vec)

    ## Produce the pentypes
    pentype <- factor(c(rep(1,  length(sow.breed.section.vec)),
                        rep(2,  length(sow.breed.buffer.section.vec)),
                        rep(3,  length(gilt.breed.section.vec)),
                        rep(4,  length(gilt.breed.buffer.section.vec)),
                        rep(5,  length(sow.gest.section.vec)),
                        rep(6,  length(gilt.gest.section.vec)),
                        rep(7,  length(farrowing.section.vec)),
                        rep(8,  length(growing.section.vec)),
                        rep(9,  length(growing.buffer.section.vec)),
                        rep(10, length(finishing.section.vec)),
                        rep(11, length(finishing.buffer.section.vec)),
                        rep(12, length(growing.gilt.vec))),
                      labels = c("Sow breeding",
                                 "Sow breeding buffer",
                                 "Gilt breeding",
                                 "Gilt breeding buffer",
                                 "Sow gestation",
                                 "Gilt gestation",
                                 "Farrowing",
                                 "Growing",
                                 "Growing buffer",
                                 "Finishing",
                                 "Finishing buffer",
                                 "Gilt growing"))

    data.frame(Ssows = 0,
               Isows = 0,
               Sgilts = 0,
               Igilts = 0,
               Spiglets = 0,
               Ipiglets = 0,
               Sgrowers = 0,
               Igrowers = 0,
               Sfinish = 0,
               Ifinish = 0,
               pentype = pentype,
               section = section,
               countdown = 99999)
}

#' nanimals
#'
#' Number of animals in a pen
#' @param x A result or residual object
#' @return A numeric vector
nanimals <- function(x) {

    index <- names(x) %in% c("Ssows", "Isows",
                             "Sgilts", "Igilts",
                             "Spiglets", "Ipiglets",
                             "Sgrowers", "Igrowers",
                             "Sfinish", "Ifinish")

        return(rowSums(x[ ,index]))
}

#' scale p
#'
#' scaling proportions when we are sampling for production phase that works with proportions instead of n
#' @param dest The ID(s) of destination node(s)
#' @param p A vector of the proportion to put in each destination node
#' @return A vector of proportions to be passed to the events
scale_p <- function(dest, p) {

    ## We only accept that p has the same length as dest.
    stopifnot(identical(length(dest), length(p)))
    den <- sum(p)

    ## the events in SimInf are sort by node then by dest prior to
    ## execution. We are only operatoring on 1 node and therefore we
    ## just need to sort on dest.
    srt <- order(dest)
    dest <- dest[srt]

    ## We also need to sort the fractions going to each of these dest
    ## pens and scale the ratios to sum to 1:
    p <- p[srt]

    ## If the total probability is close to 1 or greater then we will
    ## scale it to total 1. The reason we do this above 0.999 is to
    ## avoid weird floating point problems.
    if(den > 0.999) {
        p <- p/den
    }

    ## Then the trick is to rescale these fractions so they can be
    ## applied in series. If, for example we need to move one node to
    ## 2 dest pens, then we need to have the first event have p = 0.5
    ## and the second event p = 1.
    scale <- 1 - c(0, cumsum(p))
    scale <- 1 - c(scale, NA) / c(NA, scale)
    scale <- scale[!is.na(scale)]

    ## Now we can reorder the scale back to the order the pens were
    ## submitted to the function:
    scale[order(srt)]
}
