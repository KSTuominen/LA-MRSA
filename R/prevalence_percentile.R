## calculating the median, 10 % and 90 % percentiles prevalence from Broens et al. 2012
## mX = sampling moment
## s = sows, p = piglets

p <- c(0.1, 0.5, 0.9)
pos_m1s <- c(0.063, 0.333, 0.500, 0.000, 0.000, 1.000, 0.700)
d_m1s <- density(pos_m1s, from = 0, to = 1)
d_m1s <- sample(d_m1s$x,
                size = 1000000,
                prob = d_m1s$y,
                replace = TRUE)
q_m1s <- quantile(d_m1s, probs = p)

pos_m2s <- c(0.250, 1.000, 0.625, 0.000, 0.000, 1.000, 0.929)
d_m2s <- density(pos_m2s, from = 0, to = 1)
d_m2s <- sample(d_m2s$x,
                size = 1000000,
                prob = d_m2s$y,
                replace = TRUE)
q_m2s <- quantile(d_m2s, probs = p)

pos_m3s <-c(0.938, 1.000, 1.000, 0.000, 0.000, 1.000, 0.714)
d_m3s <- density(pos_m3s, from = 0, to = 1)
d_m3s <- sample(d_m3s$x,
                size = 1000000,
                prob = d_m3s$y,
                replace = TRUE)
q_m3s <- quantile(d_m3s, probs = p)

pos_m2p <- c(0.427, 1.000, 0.956, 0.044, 0.013, 1.000, 0.683)
d_m2p <- density(pos_m2p, from = 0, to = 1)
d_m2p <- sample(d_m2p$x,
                size = 1000000,
                prob = d_m2p$y,
                replace = TRUE)
q_m2p <- quantile(d_m2p, probs = p)


pos_m3p <- c(0.930, 0.993, 0.989, 0.031, 0.000, 1.000, 0.872)
d_m3p <- density(pos_m3p, from = 0, to = 1)
d_m3p <- sample(d_m3p$x,
                size = 1000000,
                prob = d_m3p$y,
                replace = TRUE)
q_m3p <- quantile(d_m3p, probs = p)

pos_m4p <- c(1.000, 0.559, 1.000, 0.523, 0.209, 1.000, 1.000)
d_m4p <- density(pos_m4p, from = 0, to = 1)
d_m4p <- sample(d_m4p$x,
                size = 1000000,
                prob = d_m4p$y,
                replace = TRUE)
q_m4p <- quantile(d_m4p, probs = p)

pos_m5p <- c(NA, 1.000, 1.000, 0.891, 0.162, 1.000, 1.000)
d_m5p <- density(pos_m5p, from = 0, to = 1, na.rm = TRUE)
d_m5p <- sample(d_m5p$x,
                size = 1000000,
                prob = d_m5p$y,
                replace = TRUE)
q_m5p <- quantile(d_m5p, probs = p)

pos_m6p <- c(NA, 0.862, 0.634, 0.645, 0.132, 1.000, 0.985)
d_m6p <- density(pos_m6p, from = 0, to = 1, na.rm = TRUE)
d_m6p <- sample(d_m6p$x,
                size = 1000000,
                prob = d_m6p$y,
                replace = TRUE)
q_m6p <- quantile(d_m6p, probs = p)


prev_percentiles <- setNames(data.frame(rbind(q_m1s, q_m2s, q_m2p, q_m3s, q_m3p, q_m4p, q_m5p, q_m6p)), c("10%", "50%", "90%"))


## rotate to another form
# prev_percentiles <- data.frame("percentile" = c("10%", "50%", "90%"))
# prev_percentiles <- cbind(prev_percentiles, q_m1s, q_m2s, q_m2p, q_m3s, q_m3p, q_m4p, q_m5p, q_m6p)


