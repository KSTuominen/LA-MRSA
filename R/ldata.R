## Unit and section structure for between pen leakage. Leakage not in use currently.

## mapping the left and night neigbours for ldata and continuous state vector
## this is based on separate farm pen plan

## breeding pens, 1 section with two rows
## row 1 has 15 pens, row 2 has 18 pens
## nodes:
## sow breeding 1-10
## sow breeding buffer 11-15
## gilt breeding 16-23
## gilt breeding buffer 24-33
breeding_r1 <- data.frame(left = c(0, rep(1, 14)),
                right = c(rep(1, 14), 0))
breeding_r2 <- data.frame(left = c(0, rep(1, 17)),
                          right = c(rep(1, 17), 0))
breeding <- NULL
breeding <- rbind(breeding, breeding_r1, breeding_r2)

## gestation pens, 1 section with two rows
## rows 1-3 have 16 pens and row 4 has 17 pens
## nodes:
## sows 34-68
## gilts 69-98
gestation_r13 <- data.frame(left = c(0, rep(1, 15)),
                          right = c(rep(1, 15), 0))
gestation_r4 <- data.frame(left = c(0, rep(1, 16)),
                           right = c(rep(1, 16), 0))
gestation <- NULL
gestation <- rbind(gestation, gestation_r13, gestation_r13, gestation_r13, gestation_r4)

## 6 sections in total, each section has 2 rows, each row has 13 pens (all sections are identical)
## farrowing pens are divided into sections, nodes:
## 1 99-124
## 2 125-150
## 3 151-176
## 4 177-202
## 5 203-228
## 6 229-254

farrowing <- data.frame(left = rep(c(0, rep(1, 12)), 12),
                        right = rep(c(rep(1, 12), 0), 12))

## 10 regular growing sections in total, 1 growing buffer section, each section has 2 rows, each row has 13 pens (all sections are identical)
## growing pens are divided into sections, nodes:
## 1 255-280
## 2 281-306
## 3 307-332
## 4 333-358
## 5 359-384
## 6 385-410
## 7 411-436
## 8 437-462
## 9 463-488
## 10 489-514
## growing buffer pens 515-540

growing <- data.frame(left = rep(c(0, rep(1, 12)), 22),
                      right = rep(c(rep(1, 12), 0), 22))

## 18 regular finishing sections in total, each section has 2 rows, each row has 15 pens (all sections are identical)
## finishing pens are divided into sections, nodes:
## 1 541-570
## 2 571-560
## 3 561-630
## 4 631-660
## 5 661-690
## 6 691-720
## 7 721-750
## 8 751-780
## 9 781-810
## 10 811-840
## 11 841-870
## 12 871-900
## 13 901-930
## 14 931-960
## 15 961-990
## 16 991-1020
## 17 1021-1050
## 18 1051-1080

finishing <- data.frame(left = rep(c(0, rep(1, 14)), 36),
                      right = rep(c(rep(1, 14), 0), 36))

## 2 finishing buffer sections in total, 25 pens each, row 1 has 13 pens and row 2 has 12 pens
## finishing buffer pens are divided into sections, nodes:
## 1 1079-1103
## 2 1104-1128

finishing_b_r1 <- data.frame(left = c(0, rep(1, 12)),
                          right = c(rep(1, 12), 0))
finishing_b_r2 <- data.frame(left = c(0, rep(1, 11)),
                          right = c(rep(1, 11), 0))
finishing_b_sec <- NULL
finishing_b_sec <- rbind(finishing_b_sec, finishing_b_r1, finishing_b_r2)
finishing_b <- NULL
finishing_b <- rbind(finishing_b, finishing_b_sec, finishing_b_sec)

## Growing gilt pens, 1 section with 25 pens (13 in row 1, 12 in row 2)
## nodes 1129-1142

gilt_growing_r1 <- data.frame(left = c(0, rep(1, 12)),
                             right = c(rep(1, 12), 0))
gilt_growing_r2 <- data.frame(left = c(0, rep(1, 11)),
                             right = c(rep(1, 11), 0))
gilt_growing <- NULL
gilt_growing <- rbind(gilt_growing, gilt_growing_r1, gilt_growing_r2)
## combine dataframe
ldata <- NULL
ldata <- rbind(ldata, breeding, gestation, farrowing, growing, finishing, finishing_b, gilt_growing)
