## Setup
install.packages(c("RODBC", "scales", "ggplot2", "tidyr", "dplyr", "gridExtra", "data.table"))
lib <- c("RODBC", "scales", "ggplot2", "tidyr", "dplyr", "gridExtra", "data.table")
lapply(lib, require, character.only = TRUE)

db <- odbcConnectAccess2007("", uid="", pwd="")
LA <- sqlFetch(db, "Final_data")
LA$LA_Date <- as.Date(LA$LA_Date, '%m/%d/%Y')
LA <- LA[LA$LA_Date <= as.Date('2017-10-31'),]
TaqMan.pending <- sqlFetch(db, "52Taqman_Pending")
TaqMan.pending$LA_Date <- as.Date(TaqMan.pending$LA_Date, '%m/%d/%Y')
Sample.Received <- sqlFetch(db, "Samples Received")
Sample.Received$"Date Received" <- as.Date(Sample.Received$"Date Received", '%m/%d/%Y')
Pathology <- sqlFetch(db, "Pathology")
Pathology$"Date received" <- as.Date(Pathology$"Date received", '%m/%d/%Y')
Pathology$"Return Date" <- as.Date(Pathology$"Return Date", '%m/%d/%Y')
LA.LiPA <- sqlFetch(db, "LA_LiPA_Final_Data")
LA.LiPA$LA_Date <- as.Date(LA.LiPA$LA_Date, '%m/%d/%Y')
LA.LiPA$`Detection Date.1` <- as.Date(LA.LiPA$`Detection Date.1`, "%m/%d/%Y")
LiPA.pending <- sqlFetch(db, "LiPA_Pending")
LiPA.pending$LA_Date <- as.Date(LiPA.pending$LA_Date, '%m/%d/%Y')
INNO <- sqlFetch(db, "INNO_LiPA_Final_Data")
INNO$`Detection Date` <- as.Date(INNO$`Detection Date`, '%m/%d/%Y')
LBP.LiPA <- sqlFetch(db, "LBP_LiPA_Final_Data")
LBP.LiPA$`Detection Date` <- as.Date(LBP.LiPA$`Detection Date`, '%m/%d/%Y')

LA.t <- paste0("HPV", c(6, 11, 16, 18, 26, 31, 33, 35, 39, 40, 42, 45, 51, 52, 53, 54, 55, 56, 58, 59, 61, 62, 64, 66, 67, 68, 69, 70, 71, 72, 73, 81, 82, 83, 84, "IS39", 89, "X.1")) #Contains HPVX.1
INNO.t <- paste0("HPV", c(6.1, 11.1, 16.1, 18.1, 26.1, 31.1, 33.1, 35.1, 39.1, 40.1, 43, 44, 45.1, 51.1, 52.1, 53.1, 54.1, 56.1, 58.1, 59.1, 66.1, 68.1, 69.1, 70.1, 71.1, 73.1, 74, 82.1))
LBP.t <- paste0("HPV", c(6.2, 11.2, 16.2, 18.2, 31.2, 33.2, 34, 35.2, 39.2, 40.2, 42.1, 43.1, 44.1, 45.2, 51.2, 52.2, 53.2, 54.2, 56.2, 58.2, 59.2, 66.2, 68.2, 70.2, 73.2, 74.1))

hrhpv1 <-  paste0("HPV", c(16, 18, 26, 31, 33, 34, 35, 39, 45, 51, 52, 53, 56, 58, 59,66, 67, 68, 69, 70, 73, 82, "IS39",
                           "16.1", "16.2", "18.1", "18.2", "26.1", "31.1", "31.2", "33.1", "33.2", "35.1", "35.2",
                           "39.1", "39.2", "45.1", "45.2", "51.1", "51.2", "52.1", "52.2", "53.1", "53.2", "56.1", "56.2", "58.1", "58.2",
                           "59.1", "59.2", "66.1", "66.2", "68.1", "68.2", "69.1", "70.1", "70.2", "73.1", "73.2", "82.1"))

ohpv1 <- paste0("HPV", c(6, 11, 40, 42, 43, 44, 54, 55, 61, 62, 64, 71, 72, 74, 81, 83, 84, 89, "6.1", "6.2", "11.1", "11.2", "40.1",
                         "40.2", "42.1", "43.1","44.1","54.1", "54.2", "71.1", "74.1"))

##1) Pathology data - detailed explanation of path.completed can be found in "HPV IMPACT Final Code"
#Total data
T.path <- Pathology[Pathology$`Date received` <= as.Date('2017-10-31'),]
T.path.fun <- function(){
    T.path.completed <- T.path[!(T.path$`Study ID` %in% Pathology[Pathology$"Return Date" >= as.Date('2017-10-31'),][!is.na(Pathology[Pathology$"Return Date" >= as.Date('2017-10-31'),]$"Return Date"),]$"Study ID"),]

    t.path.table <- rbind(nrow(T.path), nrow(T.path.completed), nrow(T.path) - nrow(T.path.completed), length(grep("^5|,5|/5|/^5|, 5|^5$", T.path.completed$Limitations)),
                          length(grep("^4|,4|/^4|, 4|4.|^4$", T.path.completed$Limitations)), length(grep("^3|,3|/3|/^3|, 3|^3$", T.path.completed$Limitations)),
                          length(grep("^2|,2|/2|/^2|, 2|^2$", T.path.completed$Limitations)), length(grep("^1|,1|/1|/^1|, 1|^1$", T.path.completed$Limitations)))
    t.path.table
  }

T.path.fun()

T.LA <- LA.LiPA[LA.LiPA$LA_Date <= as.Date('2017-10-31'),]
T.LA$LAStatus[rowSums(T.LA[,c(5:42)]) > 0] <- "positive"
T.LA$LAStatus[rowSums(T.LA[,c(5:42)]) == 0] <- "negative"
T.LA$LAStatus[rowSums(T.LA[,c(5:44)]) == 0] <- "inadequate"

T.LBP <- LBP.LiPA[LBP.LiPA$`Detection Date` >= as.Date('2017-10-31'),]
T.LBP <- T.LBP[!grepl("P", T.LBP$`Sample ID`),] #controls removed
T.LBP <- T.LBP[!is.na(T.LBP$`Sample ID`),] # Two NA removed

##2) Genotyping
T.Genotyping <- T.LA[!(T.LA$`Sample ID` %in% T.LBP$`Sample ID`),][!(T.LA[!(T.LA$`Sample ID` %in% T.LBP$`Sample ID`),]$"Sample ID" %in% LiPA.pending[LiPA.pending$LA_Date <= as.Date('2017-10-31'),]$"Sample ID"),]
Repeat <- T.Genotyping[T.Genotyping$LAStatus == "positive",][!is.na(T.Genotyping[T.Genotyping$LAStatus == "positive",]$p1),] #30 LA positive samples were LiPA tested
T.Genotyping[T.Genotyping$`Sample ID` %in% Repeat$`Sample ID`,][,c(INNO.t)] <- NA #NA'ed INNO LiPA data of thse sample
T.Genotyping[T.Genotyping$`Sample ID` %in% Repeat$`Sample ID`,][,c(114,115,117:144)] <- NA

T.positive <- T.Genotyping[rowSums(T.Genotyping[,c(5:18,20:42,114,115,117:144,205,209:236)], na.rm = TRUE) > 0,] #working total positive dataframe
T.inadequate <- T.Genotyping[T.Genotyping$LAStatus == "inadequate",][rowSums(T.Genotyping[T.Genotyping$LAStatus == "inadequate",][,c(111,114,115,117:144,209:236)], na.rm = TRUE) == 0,] #working inadequate with correct number
T.negative <- T.Genotyping[!(T.Genotyping$`Sample ID` %in% T.inadequate$`Sample ID`),][rowSums(T.Genotyping[!(T.Genotyping$`Sample ID` %in% T.inadequate$`Sample ID`),][,c(5:18,20:42,114,115,117:144,205,209:236)], na.rm = TRUE) == 0,]

INNO.HPVX <- as.data.frame(which(rowSums(T.positive[,c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(T.positive[,c(INNO.t)], na.rm = TRUE) == 0))
LBP.HPVX <- as.data.frame(which(rowSums(T.positive[,c(209:236)], na.rm = TRUE) > 0 & rowSums(T.positive[,c(LBP.t)], na.rm = TRUE) == 0))
F.HPVX <- rbindlist(list(INNO.HPVX,LBP.HPVX))
F.HPVX <- T.positive[1:NROW(T.positive) %in% F.HPVX$`which(rowSums(T.positive[, c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(T.positive[, c(INNO.t)], na.rm = TRUE) == 0)`,]
T.positive$HPVX.1[T.positive$`Sample ID` %in% F.HPVX$`Sample ID`] <- 1

T.positive <- T.positive %>%
  mutate(HR = rowSums(T.positive[,hrhpv1], na.rm = TRUE) > 0) %>%
  mutate(OR = rowSums(T.positive[,ohpv1], na.rm = TRUE) > 0) %>%
  mutate(HPVSUM = rowSums(T.positive[,c(LA.t,LBP.t,INNO.t)], na.rm = TRUE))

T.Genotyping.t <-
  t(T.positive %>%
      summarise(complete = nrow(T.Genotyping),
                positive = nrow(T.positive),
                negative = nrow(T.negative),
                inadequate = nrow(T.inadequate),
                single = sum(T.positive$HPVSUM == 1),
                multiple = sum(T.positive$HPVSUM > 1),
                types.per = sum(T.positive$HPVSUM)/nrow(T.positive),
                hpv16 = sum(rowSums(T.positive[,c(7,59,153)], na.rm = TRUE) > 0),
                hpv18 = sum(rowSums(T.positive[,c(8,61,155)], na.rm = TRUE) > 0),
                hpv16.18 = nrow(subset(T.positive, T.positive$HPV16 == 1 & T.positive$HPV18 == 1|T.positive$HPV16.1 & T.positive$HPV18.1|T.positive$HPV16.2 & T.positive$HPV18.2)),
                hr = sum(HR == TRUE & OR == FALSE),
                or = sum(OR == TRUE & HR == FALSE),
                hr.or = sum(OR == TRUE & HR == TRUE),
                HPVX = sum(HPVX.1, na.rm = TRUE)))

##3) Control
T.Control <- rbind(
  T.EW.LA = sum(grepl("E", LA$`Sample ID`)),
  T.PC.LA = sum(grepl("C", LA$`Sample ID`)) + sum(grepl("S", LA$`Sample ID`)),
  T.PW.LA = sum(grepl("PW", LA$`Sample ID`)),
  T.PC.LiPA = sum(grepl("C", INNO$`Sample ID`)) + sum(grepl("S", INNO$`Sample ID`)) + sum(grepl("C", LBP.LiPA$`Sample ID`)),
  T.PW.LiPA = sum(grepl("W", INNO$`Sample ID`)) + sum(grepl("W", LBP.LiPA$`Sample ID`)))

##4) LiPA  log

T.LiPA.t <- function(){

  LA.neg.LiPA.pos <- T.Genotyping[T.Genotyping$LAStatus == "negative",][rowSums(T.Genotyping[T.Genotyping$LAStatus == "negative",][,c(114,115,117:144,205,209:236)], na.rm = TRUE) > 0,]
  LA.inad.LiPA.pos <- T.Genotyping[T.Genotyping$LAStatus == "inadequate",][rowSums(T.Genotyping[T.Genotyping$LAStatus == "inadequate",][,c(114,115,117:144,205,209:236)], na.rm = TRUE) > 0,]
  LA.inad.LiPA.neg <- T.Genotyping[T.Genotyping$LAStatus == "inadequate",][rowSums(T.Genotyping[T.Genotyping$LAStatus == "inadequate",][,c(117:144)], na.rm = TRUE) == 0,]
  overlap <- sum(LA.inad.LiPA.pos$`Sample ID` %in% LA.inad.LiPA.neg$`Sample ID`)

  T.LiPA.pending <- sum(LiPA.pending$LA_Date <= as.Date('2017-10-31')) + sum(T.path$`Study ID` %in% T.LBP$`Sample ID`)
  T.LiPA.completed <- sum(T.LA$LAStatus == "negative") + sum(T.LA$LAStatus == "inadequate") - T.LiPA.pending
  Results.as.LiPA <- (nrow(LA.neg.LiPA.pos) + nrow(LA.inad.LiPA.pos) + nrow(LA.inad.LiPA.neg) - overlap)
  rbind(T.LiPA.completed, T.LiPA.pending, Results.as.LiPA)
}

T.LiPA.t()

##5) TaqMan

T.LA$TAQMAN <- as.Date(T.LA$TAQMAN, "%m/%d/%Y")

TaqMan <-
  rbind(sum(T.LA$TAQMAN <= as.Date('2017-10-31'), na.rm = TRUE),
        sum(T.LA$HPV52[T.LA$TAQMAN <= as.Date('2017-10-31')], na.rm = TRUE),
        nrow(TaqMan.pending[TaqMan.pending$LA_Date <= as.Date('2017-10-31'),]) + sum(T.LA$TAQMAN >= as.Date('2017-10-31'), na.rm = TRUE))

##1] Pathology data
#Quarter data

Q.path <- Pathology[complete.cases(Pathology$`Return Date`),][Pathology[complete.cases(Pathology$`Return Date`),]$`Return Date` >= as.Date('2017-08-01') & Pathology[complete.cases(Pathology$`Return Date`),]$`Return Date` <= as.Date('2017-10-31'),]
Q.path.fun <- function(){
  Q.path.completed <- Q.path[!(Q.path$`Study ID` %in% Pathology[Pathology$`Return Date` >= as.Date('2017-08-01') & Pathology$"Return Date" >= as.Date('2017-10-31'),][!is.na(Pathology[Pathology$`Return Date` >= as.Date('2017-08-01') & Pathology$"Return Date" >= as.Date('2017-10-31'),]$"Return Date"),]$"Study ID"),]

  Q.path.table <- rbind(nrow(Pathology[Pathology$`Date received` >= as.Date('2017-08-01') & Pathology$`Date received` <= as.Date('2017-10-31'),]), nrow(Q.path.completed), nrow(Q.path) - nrow(Q.path.completed), length(grep("^5|,5|/5|/^5|, 5|^5$", Q.path.completed$Limitations)),
                        length(grep("^4|,4|/^4|, 4|4.|^4$", Q.path.completed$Limitations)), length(grep("^3|,3|/3|/^3|, 3|^3$", Q.path.completed$Limitations)),
                        length(grep("^2|,2|/2|/^2|, 2|^2$", Q.path.completed$Limitations)), length(grep("^1|,1|/1|/^1|, 1|^1$", Q.path.completed$Limitations)))
  Q.path.table
}

Q.path.fun()


##2] Genotyping
Q.LA <- T.LA[T.LA$LA_Date >= as.Date('2017-08-01'),]
Q.LBP <- T.LBP[T.LBP$`Detection Date` >= as.Date('2017-08-01'),]

Q.Genotyping <- T.Genotyping[T.Genotyping$LA_Date >= as.Date('2017-08-01'),]
Q.positive <- Q.Genotyping[rowSums(Q.Genotyping[,c(5:18,20:42,114,115,117:144,205,209:236)], na.rm = TRUE) > 0,] #working total positive dataframe
Q.inadequate <- Q.Genotyping[Q.Genotyping$LAStatus == "inadequate",][rowSums(Q.Genotyping[Q.Genotyping$LAStatus == "inadequate",][,c(111,114,115,117:144,209:236)], na.rm = TRUE) == 0,] #working inadequate with correct number
Q.negative <- Q.Genotyping[!(Q.Genotyping$`Sample ID` %in% Q.inadequate$`Sample ID`),][rowSums(Q.Genotyping[!(Q.Genotyping$`Sample ID` %in% Q.inadequate$`Sample ID`),][,c(5:18,20:42,114,115,117:144,205,209:236)], na.rm = TRUE) == 0,]

INNO.HPVX <- as.data.frame(which(rowSums(Q.positive[,c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(Q.positive[,c(INNO.t)], na.rm = TRUE) == 0))
LBP.HPVX <- as.data.frame(which(rowSums(Q.positive[,c(209:236)], na.rm = TRUE) > 0 & rowSums(Q.positive[,c(LBP.t)], na.rm = TRUE) == 0))
F.HPVX <- rbindlist(list(INNO.HPVX,LBP.HPVX))
F.HPVX <- Q.positive[1:NROW(Q.positive) %in% F.HPVX$`which(rowSums(Q.positive[, c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(Q.positive[, c(INNO.t)], na.rm = TRUE) == 0)`,]
Q.positive$HPVX.1[Q.positive$`Sample ID` %in% F.HPVX$`Sample ID`] <- 1

Q.positive <- Q.positive %>%
  mutate(HR = rowSums(Q.positive[,hrhpv1], na.rm = TRUE) > 0) %>%
  mutate(OR = rowSums(Q.positive[,ohpv1], na.rm = TRUE) > 0) %>%
  mutate(HPVSUM = rowSums(Q.positive[,c(LA.t,LBP.t,INNO.t)], na.rm = TRUE))

Q.Genotyping.t <-
  t(Q.positive %>%
      summarise(complete = nrow(Q.Genotyping),
                positive = nrow(Q.positive),
                negative = nrow(Q.negative),
                inadequate = nrow(Q.inadequate),
                single = sum(Q.positive$HPVSUM == 1),
                multiple = sum(Q.positive$HPVSUM > 1),
                types.per = sum(Q.positive$HPVSUM)/nrow(Q.positive),
                hpv16 = sum(rowSums(Q.positive[,c(7,59,153)], na.rm = TRUE) > 0),
                hpv18 = sum(rowSums(Q.positive[,c(8,61,155)], na.rm = TRUE) > 0),
                hpv16.18 = nrow(subset(Q.positive, Q.positive$HPV16 == 1 & Q.positive$HPV18 == 1|Q.positive$HPV16.1 & Q.positive$HPV18.1|Q.positive$HPV16.2 & Q.positive$HPV18.2)),
                hr = sum(HR == TRUE & OR == FALSE),
                or = sum(OR == TRUE & HR == FALSE),
                hr.or = sum(OR == TRUE & HR == TRUE),
                HPVX = sum(HPVX.1, na.rm = TRUE)))

##3) Control
Q.Control <- rbind(
  Q.EW.LA = sum(grepl("E", LA[LA$LA_Date >= as.Date('2017-08-01'),]$`Sample ID`)),
  Q.PC.LA = sum(grepl("C", LA[LA$LA_Date >= as.Date('2017-08-01'),]$`Sample ID`)) + sum(grepl("S", LA[LA$LA_Date >= as.Date('2017-08-01'),]$`Sample ID`)),
  Q.PW.LA = sum(grepl("PW", LA[LA$LA_Date >= as.Date('2017-08-01'),]$`Sample ID`)),
  Q.PC.LiPA = sum(grepl("C", INNO[INNO$`Detection Date` >= as.Date('2017-08-01'),]$`Sample ID`)) + sum(grepl("S", INNO[INNO$`Detection Date` >= as.Date('2017-08-01'),]$`Sample ID`)) + sum(grepl("C", LBP.LiPA[LBP.LiPA$`Detection Date` >= as.Date('2017-08-01'),]$`Sample ID`)),
  Q.PW.LiPA = sum(grepl("W", INNO[INNO$`Detection Date` >= as.Date('2017-08-01'),]$`Sample ID`)) + sum(grepl("W", LBP.LiPA[LBP.LiPA$`Detection Date` >= as.Date('2017-08-01'),]$`Sample ID`)))

##4] LiPA- Log
Q.LiPA.t <- function(){

  LA.neg.LiPA.pos <- Q.Genotyping[Q.Genotyping$LAStatus == "negative",][rowSums(Q.Genotyping[Q.Genotyping$LAStatus == "negative",][,c(114,115,117:144,205,209:236)], na.rm = TRUE) > 0,]
  LA.inad.LiPA.pos <- Q.Genotyping[Q.Genotyping$LAStatus == "inadequate",][rowSums(Q.Genotyping[Q.Genotyping$LAStatus == "inadequate",][,c(114,115,117:144,205,209:236)], na.rm = TRUE) > 0,]
  LA.inad.LiPA.neg <- Q.Genotyping[Q.Genotyping$LAStatus == "inadequate",][rowSums(Q.Genotyping[Q.Genotyping$LAStatus == "inadequate",][,c(117:144)], na.rm = TRUE) == 0,]
  overlap <- sum(LA.inad.LiPA.pos$`Sample ID` %in% LA.inad.LiPA.neg$`Sample ID`)

  Results.as.LiPA <- (nrow(LA.neg.LiPA.pos) + nrow(LA.inad.LiPA.pos) + nrow(LA.inad.LiPA.neg) - overlap)
  Q.LiPA.pending <- nrow(Q.LA[!(Q.LA$`Sample ID` %in% Q.Genotyping$`Sample ID`),])
  Q.LiPA.completed <- sum(LA.LiPA$`Detection Date.1` >= as.Date('2017-08-01') & LA.LiPA$`Detection Date.1` <= as.Date('2017-10-31'), na.rm = TRUE)

  rbind(Q.LiPA.completed, Q.LiPA.pending, Results.as.LiPA)
}


Q.LiPA.t() #with results reported as LiPA combined

##5) TaqMan

Q.TaqMan <-
  rbind(sum(T.LA$TAQMAN >= as.Date('2017-08-01') & T.LA$TAQMAN <= as.Date('2017-10-31'), na.rm = TRUE),
        sum(T.LA$HPV52[T.LA$TAQMAN >= as.Date('2017-08-01') & T.LA$TAQMAN <= as.Date('2017-10-31')], na.rm = TRUE),
        nrow(TaqMan.pending[TaqMan.pending$LA_Date >= as.Date('2017-08-01') & TaqMan.pending$LA_Date <= as.Date('2017-10-31'),]) + sum(Q.LA$TAQMAN >= as.Date('2017-10-31'), na.rm = TRUE))

###3. Table

  tablefun <- function(){

    F.df <- data.frame(cbind(rbind(T.path.fun(), TaqMan, T.LiPA.t(), T.Genotyping.t, T.Control), rbind(Q.path.fun(), Q.TaqMan, Q.LiPA.t(), Q.Genotyping.t,Q.Control)))

    F.df[10,1] <- sprintf("%d  %.0f%%", as.numeric(F.df[10,1]), as.numeric(F.df[10,1])/as.numeric(F.df[9,1]) * 100)
    F.df[10,2] <- sprintf("%d  %.0f%%", as.numeric(F.df[10,2]), as.numeric(F.df[10,2])/as.numeric(F.df[9,2]) * 100)

    F.df[c(4:8),1] <- sprintf("%d  %.0f%%", as.numeric(F.df[c(4:8),1]), as.numeric(F.df[c(4:8),1])/as.numeric(F.df[2,1]) * 100)
    F.df[c(4:8),2] <- sprintf("%d  %.0f%%", as.numeric(F.df[c(4:8),2]), as.numeric(F.df[c(4:8),2])/as.numeric(F.df[2,2]) * 100)

    F.df[c(16:17,19:20,22:28),1] <- sprintf("%d  %.0f%%", as.numeric(F.df[c(16:17,19:20,22:28),1]), as.numeric(F.df[c(16:17,19:20,22:28),1])/(as.numeric(F.df[16,1]) + as.numeric(F.df[17,1])) * 100)
    F.df[c(16:17,19:20,22:28),2] <- sprintf("%d  %.0f%%", as.numeric(F.df[c(16:17,19:20,22:28),2]), as.numeric(F.df[c(16:17,19:20,22:28),2])/(as.numeric(F.df[16,2]) + as.numeric(F.df[17,2])) * 100)

    F.df[18,1] <- sprintf("%d  %.0f%%", as.numeric(F.df[18,1]), as.numeric(F.df[18,1])/as.numeric(F.df[15,1]) * 100)
    F.df[18,2] <- sprintf("%d  %.0f%%", as.numeric(F.df[18,2]), as.numeric(F.df[18,2])/as.numeric(F.df[15,2]) * 100)

    F.df[21,c(1,2)] <- as.numeric(sprintf("%.2f", F.df[21,c(1,2)]))
    F.df[21,1] <- paste(F.df[21,1] ,paste(range(Q.positive[Q.positive$LAStatus == "positive",]$HPVSUM)[1], range(Q.positive[Q.positive$LAStatus == "positive",]$HPVSUM)[2], sep = "-"))
    F.df[21,2] <- paste(F.df[21,2] ,paste(range(Q.positive[Q.positive$LAStatus == "positive",]$HPVSUM)[1], range(Q.positive[Q.positive$LAStatus == "positive",]$HPVSUM)[2], sep = "-"))

    F.df[NROW(F.df)+7,] <- NA
    F.df
  }

table <- tablefun()
table <- table[c(1:8,34,35,36,37,38,39, 9,11,10,12:14,15:18,22:28,19:21,40,29:33), c(2,1)]

rownames(table) <- c("Samples received", "Pathology completed", "Pathology pending", "Excluded", "Insufficient", "Poor preservation", "Minute", "Small lesion", "DNA extraction",
                    "Extraction completed", "Extraction pending", "Genotyping", "LA completed", "LA pending", "52 completed", "52 pending", "52 positive", "LiPA completed",
                    "LiPA pending", "Results as LiPA", "Genotyping completed", "HPV positive", "HPV negative", "Inadequate", "HPV16", "HPV18", "HPV16 and HPV18", "HR", "OR",
                    "HR and OR", "HPVX", "single type", "multiple types", "types per positve","Controls", "H2O Extracted", "LA H2O PCR controls", "LA pHPV16 PCR controls",
                    "LiPA H2O PCR controls", "LiPA pHPV16 PCR controls")

colnames(table) <- c("Quarter data", "Total data")
table

write.csv(table, "impact.csv")

###4. Plot
T.INNO.HPVX <- as.data.frame(which(rowSums(T.Genotyping[,c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(T.Genotyping[,c(INNO.t)], na.rm = TRUE) == 0))
T.LBP.HPVX <- as.data.frame(which(rowSums(T.Genotyping[,c(208:235)], na.rm = TRUE) > 0 & rowSums(T.Genotyping[,c(LBP.t)], na.rm = TRUE) == 0))
T.F.HPVX <- rbindlist(list(T.INNO.HPVX,T.LBP.HPVX))
T.F.HPVX <- T.Genotyping[1:NROW(T.Genotyping) %in% T.F.HPVX$`which(rowSums(T.Genotyping[, c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(T.Genotyping[, c(INNO.t)], na.rm = TRUE) == 0)`,]
T.Genotyping$HPVX.1[T.Genotyping$`Sample ID` %in% T.F.HPVX$`Sample ID`] <- 1 #HPVX to be appended
Q.Genotyping <- T.Genotyping[T.Genotyping$LA_Date >= as.Date('2017-08-01'),]
molten <- T.Genotyping[,c("Sample ID" ,LA.t, INNO.t, LBP.t)]
Q.molten <- Q.Genotyping[,c("Sample ID", LA.t, INNO.t, LBP.t)]


molten <- molten %>%
  summarise(HPV6 = sum(HPV6, HPV6.1, HPV6.2, na.rm = TRUE),
            HPV11 = sum(HPV11, HPV11.1, HPV11.2, na.rm = TRUE),
            HPV16 = sum(HPV16, HPV16.1, HPV16.2, na.rm = TRUE),
            HPV18 = sum(HPV18, HPV18.1, HPV18.2, na.rm = TRUE),
            HPV26 = sum(HPV26, HPV26.1, na.rm = TRUE),
            HPV31 = sum(HPV31, HPV31.1, HPV31.2, na.rm = TRUE),
            HPV33 = sum(HPV33, HPV33.1, HPV33.2, na.rm = TRUE),
            HPV34 = sum(HPV34, na.rm = TRUE),
            HPV35 = sum(HPV35, HPV35.1, HPV35.2, na.rm = TRUE),
            HPV39 = sum(HPV39, HPV39.1, HPV39.2, na.rm = TRUE),
            HPV40 = sum(HPV40, HPV40.1, HPV40.2, na.rm = TRUE),
            HPV42 = sum(HPV42, HPV42.1, na.rm = TRUE),
            HPV43 = sum(HPV43, HPV43.1,  na.rm = TRUE),
            HPV44 = sum(HPV44, HPV44.1,  na.rm = TRUE),
            HPV45 = sum(HPV45, HPV45.1, HPV45.2, na.rm = TRUE),
            HPV51 = sum(HPV51, HPV51.1, HPV51.2, na.rm = TRUE),
            HPV52 = sum(HPV52, HPV52.1, HPV52.2, na.rm = TRUE),
            HPV53 = sum(HPV53, HPV53.1, HPV53.2, na.rm = TRUE),
            HPV54 = sum(HPV54, HPV54.1, HPV54.2, na.rm = TRUE),
            HPV55 = sum(HPV55,  na.rm = TRUE),
            HPV56 = sum(HPV56, HPV56.1, HPV56.2, na.rm = TRUE),
            HPV58 = sum(HPV58, HPV58.1, HPV58.2, na.rm = TRUE),
            HPV59 = sum(HPV59, HPV59.1, HPV59.2, na.rm = TRUE),
            HPV61 = sum(HPV61,  na.rm = TRUE),
            HPV62 = sum(HPV62,  na.rm = TRUE),
            HPV64 = sum(HPV64,  na.rm = TRUE),
            HPV66 = sum(HPV66, HPV66.1, HPV66.2, na.rm = TRUE),
            HPV67 = sum(HPV67, na.rm = TRUE),
            HPV68 = sum(HPV68, HPV68.1, HPV68.2, na.rm = TRUE),
            HPV69 = sum(HPV69, HPV69.1, na.rm = TRUE),
            HPV70 = sum(HPV70, HPV70.1, HPV70.2, na.rm = TRUE),
            HPV71 = sum(HPV71, HPV71.1, na.rm = TRUE),
            HPV72 = sum(HPV72,  na.rm = TRUE),
            HPV73 = sum(HPV73, HPV73.1, HPV73.2, na.rm = TRUE),
            HPV74 = sum(HPV74, HPV74.1,  na.rm = TRUE),
            HPV81 = sum(HPV81,  na.rm = TRUE),
            HPV82 = sum(HPV82, HPV82.1, na.rm = TRUE),
            HPV83 = sum(HPV83,  na.rm = TRUE),
            HPV84 = sum(HPV84,  na.rm = TRUE),
            HPVIS39 = sum(HPVIS39, na.rm = TRUE),
            HPV89 = sum(HPV89, na.rm = TRUE),
            HPVX = sum(HPVX.1, na.rm = TRUE))

Q.molten <- Q.molten %>%
  summarise(HPV6 = sum(HPV6, HPV6.1, HPV6.2, na.rm = TRUE),
            HPV11 = sum(HPV11, HPV11.1, HPV11.2, na.rm = TRUE),
            HPV16 = sum(HPV16, HPV16.1, HPV16.2, na.rm = TRUE),
            HPV18 = sum(HPV18, HPV18.1, HPV18.2, na.rm = TRUE),
            HPV26 = sum(HPV26, HPV26.1, na.rm = TRUE),
            HPV31 = sum(HPV31, HPV31.1, HPV31.2, na.rm = TRUE),
            HPV33 = sum(HPV33, HPV33.1, HPV33.2, na.rm = TRUE),
            HPV34 = sum(HPV34, na.rm = TRUE),
            HPV35 = sum(HPV35, HPV35.1, HPV35.2, na.rm = TRUE),
            HPV39 = sum(HPV39, HPV39.1, HPV39.2, na.rm = TRUE),
            HPV40 = sum(HPV40, HPV40.1, HPV40.2, na.rm = TRUE),
            HPV42 = sum(HPV42, HPV42.1, na.rm = TRUE),
            HPV43 = sum(HPV43, HPV43.1,  na.rm = TRUE),
            HPV44 = sum(HPV44, HPV44.1,  na.rm = TRUE),
            HPV45 = sum(HPV45, HPV45.1, HPV45.2, na.rm = TRUE),
            HPV51 = sum(HPV51, HPV51.1, HPV51.2, na.rm = TRUE),
            HPV52 = sum(HPV52, HPV52.1, HPV52.2, na.rm = TRUE),
            HPV53 = sum(HPV53, HPV53.1, HPV53.2, na.rm = TRUE),
            HPV54 = sum(HPV54, HPV54.1, HPV54.2, na.rm = TRUE),
            HPV55 = sum(HPV55,  na.rm = TRUE),
            HPV56 = sum(HPV56, HPV56.1, HPV56.2, na.rm = TRUE),
            HPV58 = sum(HPV58, HPV58.1, HPV58.2, na.rm = TRUE),
            HPV59 = sum(HPV59, HPV59.1, HPV59.2, na.rm = TRUE),
            HPV61 = sum(HPV61,  na.rm = TRUE),
            HPV62 = sum(HPV62,  na.rm = TRUE),
            HPV64 = sum(HPV64,  na.rm = TRUE),
            HPV66 = sum(HPV66, HPV66.1, HPV66.2, na.rm = TRUE),
            HPV67 = sum(HPV67, na.rm = TRUE),
            HPV68 = sum(HPV68, HPV68.1, HPV68.2, na.rm = TRUE),
            HPV69 = sum(HPV69, HPV69.1, na.rm = TRUE),
            HPV70 = sum(HPV70, HPV70.1, HPV70.2, na.rm = TRUE),
            HPV71 = sum(HPV71, HPV71.1, na.rm = TRUE),
            HPV72 = sum(HPV72,  na.rm = TRUE),
            HPV73 = sum(HPV73, HPV73.1, HPV73.2, na.rm = TRUE),
            HPV74 = sum(HPV74, HPV74.1,  na.rm = TRUE),
            HPV81 = sum(HPV81,  na.rm = TRUE),
            HPV82 = sum(HPV82, HPV82.1, na.rm = TRUE),
            HPV83 = sum(HPV83,  na.rm = TRUE),
            HPV84 = sum(HPV84,  na.rm = TRUE),
            HPVIS39 = sum(HPVIS39, na.rm = TRUE),
            HPV89 = sum(HPV89, na.rm = TRUE),
            HPVX = sum(HPVX.1, na.rm = TRUE))


plotfun <- function(){
  theme_set(theme_classic())

  T.frequency = as.data.frame(t(molten/nrow(T.Genotyping)))
  Q.frequency = as.data.frame(t(Q.molten/nrow(Q.Genotyping)))
  Q.frequency <- cbind(HPV_Types = rownames(Q.frequency), Q.frequency)
  rownames(Q.frequency) <- NULL
  Q.frequency$group = 1
  T.frequency <- cbind(HPV_Types = rownames(T.frequency), T.frequency)
  rownames(T.frequency) <- NULL
  T.frequency$group = 2
  bind <- rbind(Q.frequency, T.frequency)

  labs <- c(6, 11, 16, 18, 26, 31, 33, 34, 35, 39, 40, 42, 43, 44, 45, 51, 52, 53, 54, 55, 56, 58, 59, 61, 62, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 81, 82, 83, 84, "IS39", 89, "X")
  labs <- setNames(labs, paste0("HPV", labs))

  Plot <- ggplot(bind, aes(x=factor(bind$HPV_Types, names(labs)), y= V1)) + geom_bar(position = "dodge", stat = "identity", aes(fill=as.factor(group))) + labs(title = "HPV Type Distribution") + xlab("HPV Types") + ylab("% Positive") + scale_y_continuous(labels = scales::percent, expand = c(0,0)) +theme(axis.text.x = element_text(size = 7, vjust = 0.5, hjust = 1, angle = 90), plot.title = element_text(hjust = 0.5)) + scale_fill_manual("Dataset", labels=c("Quarterly", "Cycle Total"), values = c("blue", "red")) + scale_x_discrete(labels = labs)

  Plot
}

plotfun()
ggsave("plot.tiff", plotfun(), width = 13, height = 4.5, dpi = 300, limitsize = TRUE ,device = "tiff", scale = 0.8)
