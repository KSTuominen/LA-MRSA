##' event_dataframe
##'
##' A function to generate a single event. This is a dataframe with
##' one line in the correct structure that we expect for events.
##'
##' @param event The event type; one of "extTrans", "enter", "exit".
##' @param time A value that can be coersed to an integer greater than
##'     or equal to 1.
##' @param node A value that can be coersed to an integer greater than
##'     or equal to 1.
##' @param dest A value that can be coersed to an integer greater than
##'     or equal to 0.
##' @param n A value that can be coersed to an integer greater than or
##'     equal to 0.
##' @param proportion A numeric between 0 and 1
##' @param select The select value; an value that can be coersed to an
##'     integer between 1 and 6
##' @param shift The shift value; an value that can be coersed to an
##'     integer between 0 and 4
##' @return A dataframe
event_dataframe <- function(event = c("extTrans", "intTrans", "enter", "exit"),
                            time = 1,
                            node,
                            dest,
                            n,
                            proportion = 0,
                            select = c(1, 2, 3, 4, 5, 6, 7),
                            shift = c(1, 2, 3, 4)) {

    ## Check the arguments:

    ## event should be one of the default values
    event <- match.arg(event, several.ok = TRUE)

    ## We should be able to coerce 'time', 'node', 'dest', and 'n' to
    ## integers. Throw an error if we get NA when converting to
    ## integer:
    stopifnot(identical(length(time), 1L))
    time <- as.integer(time)
    stopifnot(!is.na(time))
    stopifnot(time > 0)

    ## Allow for more than one node
    node <- as.integer(node)
    stopifnot(all(!is.na(node)))
    stopifnot(all(node > 0))

    ## Allow for more than one destination pen
    dest <- as.integer(dest)
    stopifnot(all(!is.na(dest)))
    stopifnot(all(dest >= 0))

    ## Allow for more than one number of animals
    n <- as.integer(n)
    stopifnot(all(!is.na(n)))

    ## Proportion should be between 0 and 1
    proportion <- as.numeric(proportion)
    stopifnot(all(!is.na(proportion)))
    stopifnot(all(proportion >= 0 & proportion <= 1))

    ## Check that we don't both use proportion and number in the same
    ## event
    if(any(n > 0 & proportion > 0)) {
        stop("You should declare proportion 0 if n > 0")
    }

    ## Check the potential select values
    select <- as.integer(select)
    stopifnot(all(select %in% c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)))

    ## Check the potential select values
    shift <- as.integer(shift)
    stopifnot(all(shift %in% c(0L, 1L, 2L, 3L, 4L)))

    df <- data.frame(event = event,
                     time = time,
                     node = node,
                     dest = dest,
                     n = n,
                     proportion = proportion,
                     select = select,
                     shift = shift,
                     stringsAsFactors = TRUE)
    df[order(df$node, df$dest), ]
}
##' event
##'
##' A function that give back an event but doesn't require us to
##' understand what shift and select should be.
##'
##' @param type A value that can be coersed to an integer between 1
##'     and 20. The following values are available:
##'
##' 1 Move a sow from breeding to gestation/the opposite direction: select = 1 shift = 0
##' 2 Move a sow from gestation to farrowing: select = 1 shift = 0
##' 3 Move a sow from farrowing to breeding: select = 1 shift = 0
##' 4 Move a gilt grower to breeding: select = 2 shift = 0
##' 5 Move a gilt from breeding to gestation/the opposite direction: select = 2 shift = 0
##' 6 Move a gilt from gestation to farrowing: select = 2 shift = 1
##' 7 Piglets born: select = 3 shift = 0
##' 8 Piglets weaned to grower: select = 4 shift = 2
##' 9 Growers moved to finisher: select = 5 shift = 3
##' 10 Finishers to slaughter: select= 6 shift = 0
##' 11 piglets move to gilt grower (-> gilts): select = 4, shift = 4
##' 12 piglet mortality: select = 4, shift = 0
##' 13 grower mortality: select = 5, shift = 0
##' 14 sow mortality: select = 1, shift = 0
##' 15 gilt mortality: select = 2, shift = 0
##' 16 Purchase gilt: select = 7, shift = 0
##' 17 Gilt growers to finishing (intTrans): select = 5, shift = 3 ##NOT IN USE 
##' 18 Growers moved to growing buffer: select = 5, shift = 0
##' 19 Finishers moved to finishing buffer: select = 6, shift = 0
##' 20 Purchase infected gilt: select = 8, shift = 0
##' 21 Mingle piglets: select = 4, shift = 0 
##' 22 Mingle finishers: select = 6 shift = 0
##'
##'
##' @param time The time
##' @param node The node
##' @param dest The destination node
##' @param n The number of animals
##' @param proportion The proportion of animals
##' @return A data.frame with one line
event <- function(type,
                  time,
                  node,
                  dest,
                  n,
                  proportion) {

    ## reference dataframe to generate select, shift and event type
    df <- data.frame(type = c(1L,  2L,  3L,  4L,  5L,  6L,  7L,  8L,  9L,  10L, 
                              11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 
                              21L, 22L),
                     select = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 4L, 5L, 6L,
                                4L, 4L, 5L, 1L, 2L, 7L, 5L, 5L, 6L, 8L, 
                                4L, 6L),
                     shift =  c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 2L, 3L, 0L,
                                4L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 
                                0L, 0L),
                     event = c("extTrans", "extTrans", "extTrans",
                               "extTrans", "extTrans", "extTrans",
                               "enter", "extTrans", "extTrans", "exit",
                               "extTrans", "exit", "exit", "exit", "exit",
                               "enter", "intTrans", "extTrans",
                               "extTrans", "enter", "extTrans", "extTrans"),
                     stringsAsFactors = FALSE)

    ## Check the arguments:
    type <- as.integer(type)
    stopifnot(all(!is.na(type)))
    stopifnot(all(type %in% df$type))
    if(is.null(node)) return(NULL)
    if(is.null(dest)) return(NULL)

    index <- match(type, df$type)

    event_dataframe(event = df$event[index],
                    time = time,
                    node = node,
                    dest = dest,
                    n = n,
                    proportion = proportion,
                    select = df$select[index],
                    shift = df$shift[index])
}
