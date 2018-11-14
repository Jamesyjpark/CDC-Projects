#########################################
####            Setup                ####
#########################################
install.packages(c("RODBC", "reshape2", "scales", "ggplot2", "dplyr", "data.table", "sas7bdat", "xlsx", "ggpmisc"))
lib <- c("RODBC", "reshape2", "scales", "ggplot2", "dplyr", "data.table", "sas7bdat", "xlsx", "ggpmisc")
lapply(lib, require, character.only = TRUE)

db <- odbcConnectAccess2007("", uid="", pwd="")
LA.LiPA.db <- sqlFetch(db, "")
SASdb <- sas7bdat::read.sas7bdat("")
write.csv(SASdb, file = "")
df.csv <- read.csv("", header = TRUE, stringsAsFactors = FALSE)

#########################################
####         Data wrangling          ####
#########################################
LA.LiPA.Epi <- LA.LiPA.db[LA.LiPA.db$`Sample ID` %in% df.csv$sample_id,] 
Epi.LA.LiPA <- df.csv[df.csv$sample_id %in% LA.LiPA.db$`Sample ID`,]

LA.LiPA.Epi$Collection_Date <- Epi.LA.LiPA$ColDate[match(LA.LiPA.Epi$`Sample ID`, Epi.LA.LiPA$sample_id)]
LA.LiPA.Epi$Specimen_Age <- as.double(difftime(as.POSIXct(LA.LiPA.Epi$PCR_Date), as.POSIXct(LA.LiPA.Epi$Collection_Date, format = '%m/%d/%Y'), units ="weeks"))
LA.LiPA.Epi$Patient_ID <- Epi.LA.LiPA$PtID[match(LA.LiPA.Epi$`Sample ID`, Epi.LA.LiPA$sample_id)] #All necessary columns from epi database appended

LA.LiPA.Epi$Specimen_Age_Quarter <- LA.LiPA.Epi$Specimen_Age %/% 13
LA.LiPA.Epi$Specimen_Age_6mo <- LA.LiPA.Epi$Specimen_Age %/% 26
LA.LiPA.Epi$Specimen_Age_yr <- LA.LiPA.Epi$Specimen_Age %/% 52

INNO <- paste0("HPV", c(6.1, 11.1, 16.1, 18.1, 26.1, 31.1, 33.1, 35.1, 39.1, 40.1, 43, 44, 45.1, 51.1, 52.1, 53.1, 54.1, 56.1, 58.1, 59.1, 66.1, 68.1, 69.1, 70.1, 71.1, 73.1, 74, 82.1))
LBP <- paste0("HPV", c(6.2, 11.2, 16.2, 18.2, 31.2, 33.2, 34, 35.2, 39.2, 40.2, 42.1, 43.1, 44.1, 45.2, 51.2, 52.2, 53.2, 54.2, 56.2, 58.2, 59.2, 66.2, 68.2, 70.2, 73.2, 74.1))
INNO_HPVX_DF <- as.data.frame(which(rowSums(LA.LiPA.Epi[,c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(LA.LiPA.Epi[,c(INNO)], na.rm = TRUE) == 0))
LBP_HPVX_DF <- as.data.frame(which(rowSums(LA.LiPA.Epi[,c(209:236)], na.rm = TRUE) > 0 & rowSums(LA.LiPA.Epi[,c(LBP)], na.rm = TRUE) == 0))
Final_HPVX_DF <- rbindlist(list(INNO_HPVX_DF, LBP_HPVX_DF))
Final_HPVX_DF_From_Total_Positive <- (LA.LiPA.Epi[1:NROW(LA.LiPA.Epi) %in% Final_HPVX_DF$LA.LiPA.Epi[, c(114, 115, 117:144)], na.rm = TRUE) > 0 & rowSums(LA.LiPA.Epi[,c(INNO)], na.rm = TRUE) == 0)
LA.LiPA.Epi$HPVX.1[LA.LiPA.Epi$`Sample ID` %in% Final_HPVX_DF_From_Total_Positive$`Sample ID`] <- 1
LA.LiPA.Epi$HPVX.1[!(LA.LiPA.Epi$`Sample ID` %in% Final_HPVX_DF_From_Total_Positive$`Sample ID`)] <- 0 #72 HPVX

LA.LiPA.Epi$Inadequate <- rowSums(LA.LiPA.Epi[,c(5:44,111,114,115,117:144,209:236,243)], na.rm = TRUE) == 0 #68 inadequate
LA.LiPA.Epi$Positive <- rowSums(LA.LiPA.Epi[,c(5:42,114,115,117:144,209:236,243)], na.rm = TRUE) > 0 #9923 positive
LA.LiPA.Epi$Negative <- rowSums(LA.LiPA.Epi[,c(5:42,114,115,117:144,209:236,243)], na.rm = TRUE) == 0 | rowSums(LA.LiPA.Epi[,c(5:42,111,114,115,117:144,243)], na.rm = TRUE) == 0 | rowSums(LA.LiPA.Epi[,c(5:44,114,115,117:144,243)], na.rm = TRUE) == 0 #643 total with overlap

LA.LiPA.Epi$AssayStatus[LA.LiPA.Epi$Positive == TRUE & LA.LiPA.Epi$Negative == FALSE & LA.LiPA.Epi$Inadequate == FALSE] <- "positive"
LA.LiPA.Epi$AssayStatus[LA.LiPA.Epi$Positive == FALSE & LA.LiPA.Epi$Negative == TRUE & LA.LiPA.Epi$Inadequate == FALSE] <- "negative"
LA.LiPA.Epi$AssayStatus[LA.LiPA.Epi$Inadequate == TRUE] <- "inadequate"

Box <- ggplot(LA.LiPA.Epi, aes(x = AssayStatus, y = Specimen_Age_yr)) + geom_boxplot() +  ggtitle("IMPACT Box Plot with Both LA and LiPA Data")
Histogram <- qplot(LA.LiPA.Epi$Specimen_Age_yr, geom="histogram", main = "Sample Distribution by Specimen Age", xlab = "Specimen Age by Years", ylab = "Number of Samples", fill = I("blue"), col=I("red"), alpha=I(.2))

#########################################
#######        Table/Plot         #######
#########################################
# Proporiton
Quarter.Proportion <- by(LA.LiPA.Epi$AssayStatus, LA.LiPA.Epi$Specimen_Age_Quarter, function(x) table(x)/length(x))
Semiyear.Proportion <- by(LA.LiPA.Epi$AssayStatus, LA.LiPA.Epi$Specimen_Age_6mo, function(x) table(x)/length(x))
Yearly.Proportion <- by(LA.LiPA.Epi$AssayStatus, LA.LiPA.Epi$Specimen_Age_yr, function(x) table(x)/length(x))

Quarter.Proportion.df <- as.data.frame(do.call("rbind", Quarter.Proportion))
Semiyear.Proportion.df <- as.data.frame(do.call("rbind", Semiyear.Proportion))
Yearly.Proportion.df <- as.data.frame(do.call("rbind", Yearly.Proportion))

# Export to xlsx
write.xlsx(Quarter.Proportion.df, file = "Quarter.xlsx")
write.xlsx(Semiyear.Proportion.df, file = "Semi.xlsx")
write.xlsx(Yearly.Proportion.df, file = "Year.xlsx")



# Bar Plot
Yearly.Proportion.df$Year <- 0:9 
Year_Molten <- melt(Yearly.Proportion.df, id = c("Year"), variable.name = "Status")
Year.Proportion.Plot <- ggplot(Year_Molten, aes(x = Year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "Change in Assay Status Proportion over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Semiyear.Proportion.df$Year <- 0:19
SemiYear_Molten <- melt(Semiyear.Proportion.df, id = c("Year"), variable.name = "Status")
Semiyear.Proportion.Plot <- ggplot(SemiYear_Molten, aes(x = Year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "Change in Assay Status Proportion over Time", x = "Specimen Age per 6 months", y = "Assay Status Proportion")

Quarter.Proportion.df$Year <- 0:35
Quarter_Molten <- melt(Quarter.Proportion.df, id = c("Year"), variable.name = "Status")
Quarter.Proportion.Plot <- ggplot(Quarter_Molten, aes(x = Year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "Change in Assay Status Proportion over Time", x = "Specimen Age per each quarter", y = "Assay Status Proportion")

# Dot Plot 
Year.Distribution <- by(LA.LiPA.Epi$AssayStatus, LA.LiPA.Epi$Specimen_Age_yr, table)
Year.Distribution.df <- as.data.frame(do.call("rbind", Year.Distribution))
Year.Distribution.df$Year <- 0:9
YearD_Molten <- melt(Year.Distribution.df, id = c("Year"))
YearD_Molten$Proportion <- Year_Molten$value 
Year.Dotplot <- ggplot(YearD_Molten, aes(x = Year, y = Proportion, fill = variable, size = value)) + geom_point(shape = 21) + ggtitle("Assay proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

Semi.Distribution <- by(LA.LiPA.Epi$AssayStatus, LA.LiPA.Epi$Specimen_Age_6mo, table)
Semi.Distribution.df <- as.data.frame(do.call("rbind", Semi.Distribution))
Semi.Distribution.df$Year <- 0:19
SemiD_Molten <- melt(Semi.Distribution.df, id = c("Year"))
SemiD_Molten$Proportion <- SemiYear_Molten$value 
Semi.Dotplot <- ggplot(SemiD_Molten, aes(x = Year, y = Proportion, fill = variable, size = value)) + geom_point(shape = 21) + ggtitle("Assay proportion and size per 6 months") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per 6 months", y = "Assay Status Proportion")

Quarter.Distribution <- by(LA.LiPA.Epi$AssayStatus, LA.LiPA.Epi$Specimen_Age_Quarter, table)
Quarter.Distribution.df <- as.data.frame(do.call("rbind", Quarter.Distribution))
Quarter.Distribution.df$Year <- 0:35
QuarterD_Molten <- melt(Quarter.Distribution.df, id = c("Year"))
QuarterD_Molten$Proportion <- Quarter_Molten$value #Now I have both proportion and number of samples
Quarter.Dotplot <- ggplot(QuarterD_Molten, aes(x = Year, y = Proportion, fill = variable, size = value)) + geom_point(shape = 21) + ggtitle("Assay proportion and size per quarter") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per quarter", y = "Assay Status Proportion")

#########################################
####         Positive only           ####
#########################################
HPV.Only <- LA.LiPA.Epi[,c(1,5:15,17:42,114,115,117:144,209:236,247)]
LA.LiPA.Epi$HPVSUM <- rowSums(HPV.Only[,2:97], na.rm = TRUE)
LA.LiPA.Epi.NOXR <- LA.LiPA.Epi[,-19]
LA.LiPA.Epi.NOXR$NumType <- LA.LiPA.Epi$HPVSUM

Year.Type.Change <- by(LA.LiPA.Epi.NOXR$NumType, LA.LiPA.Epi.NOXR$Specimen_Age_yr, table)
Year.Type.Change.Prop <- by(LA.LiPA.Epi.NOXR$NumType, LA.LiPA.Epi.NOXR$Specimen_Age_yr, function(x) table(x)/length(x)) 
Semi.Type.Change <- by(LA.LiPA.Epi.NOXR$NumType, LA.LiPA.Epi.NOXR$Specimen_Age_6mo, table)
Semi.Type.Change.Prop <- by(LA.LiPA.Epi.NOXR$NumType, LA.LiPA.Epi.NOXR$Specimen_Age_6mo, function(x) table(x)/length(x))
Quarter.Type.Change <- by(LA.LiPA.Epi.NOXR$NumType, LA.LiPA.Epi.NOXR$Specimen_Age_Quarter, table)
Quarter.Type.Change.Prop <- by(LA.LiPA.Epi.NOXR$NumType, LA.LiPA.Epi.NOXR$Specimen_Age_Quarter, function(x) table(x)/length(x))

# Time slope graphs
library(ggpmisc)
df <- data.frame(x = c(1:100))
df$y <- 2 + 3 * df$x + rnorm(100, sd = 40)
my.formula <- y ~ x

LA.LiPA.Epi$Specimen_Age_yr_formatted <- difftime(as.POSIXct(LA.LiPA.Epi$PCR_Date), as.POSIXct(LA.LiPA.Epi$Collection_Date, format = '%m/%d/%Y'))/365 #This is to format the specimen age by year to better visualize
LA.LiPA.Epi$Collection.Date.formatted <- as.POSIXct(LA.LiPA.Epi$Collection_Date, format = '%m/%d/%Y') #convertto a standard date understood R
LA.LiPA.Epi$PCR.Date.formatted <- as.POSIXct(LA.LiPA.Epi$PCR_Date, format = '%m/%d/%Y')
Collection.Time.Graph <- ggplot(LA.LiPA.Epi, aes(x = Collection.Date.formatted, y = Specimen_Age_yr_formatted)) + geom_point(shape = 1, color = "cyan") + geom_smooth(method = lm, se= FALSE) + ggtitle("IMPACT Collection Date and Specimen Age Graph")+ stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, label.x.npc = "right", label.y.npc = 0.85) + stat_fit_glance(method = "lm", method.args = list(formula = my.formula), geom = "text", aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")), label.x.npc = "right")
Extraction.Time.Graph <- ggplot(LA.LiPA.Epi, aes(x = PCR.Date.formatted, y = Specimen_Age_yr_formatted), parse = TRUE) + geom_point(shape = 1, color = "cyan") + geom_smooth(method = lm, se= FALSE) + ggtitle("IMPACT PCR Date and Specimen Age Graph") + stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, label.x.npc = "right", label.y.npc = 0.85) + stat_fit_glance(method = "lm", method.args = list(formula = my.formula), geom = "text", aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")), label.x.npc = "right")


#########################################
####       Vaccinnation type         ####
#########################################
Only.16 <- paste0("HPV", c(16, 16.1, 16.2))
Only.18 <- paste0("HPV", c(18, 18.1, 18.2))
Both.16.18 <- paste0("HPV", c(16, 16.1, 16.2, 18, 18.1, 18.2))
Four.Plex <- paste0("HPV", c(6, 11, 16, 18, 6.1, 6.2, 11.1, 11.2, 16.1, 16.2, 18.1, 18.2))
Nine.Plex <- paste0("HPV", c(6, 11, 16, 18, 31, 33, 45, 52, 58, 6.1, 6.2, 11.1, 11.2, 16.1, 16.2, 18.1, 18.2, 31.1, 31.2, 33.1, 33.2, 45.1, 45.2, 52.1, 52.2, 58.1, 58.2))

mo <- strftime(LA.LiPA.Epi$Collection.Date.formatted, "%m")
yr <- strftime(LA.LiPA.Epi$Collection.Date.formatted, "%Y")
amt <- runif(10252)
dd <- data.frame(mo, yr, amt, LA.LiPA.Epi)
dd.agg <- aggregate(amt ~ mo + yr, dd, FUN = sum)

dd$HPV16stat <- rowSums(dd[,c("HPV16", "HPV16.1", "HPV16.2")], na.rm = TRUE) > 0
dd$HPV16stat[dd$Inadequate == TRUE] <- "inadequate"
dd$HPV16stat[dd$HPV16stat == TRUE] <- "positive"
dd$HPV16stat[dd$HPV16stat == FALSE] <- "negative"

dd$HPV18stat <- rowSums(dd[,c("HPV18", "HPV18.1", "HPV18.2")], na.rm = TRUE) > 0
dd$HPV18stat[dd$Inadequate == TRUE] <- "inadequate"
dd$HPV18stat[dd$HPV18stat == TRUE] <- "positive"
dd$HPV18stat[dd$HPV18stat == FALSE] <- "negative"

dd$HPV16.18stat <- rowSums(dd[,c("HPV16","HPV16.1","HPV16.2","HPV18", "HPV18.1", "HPV18.2")], na.rm = TRUE) > 0
dd$HPV16.18stat[dd$Inadequate == TRUE] <- "inadequate"
dd$HPV16.18stat[dd$HPV16.18stat == TRUE] <- "positive"
dd$HPV16.18stat[dd$HPV16.18stat == FALSE] <- "negative"

dd$HPV4pstat <- rowSums(dd[,c("HPV6", "HPV6.1", "HPV6.2", "HPV11", "HPV11.1", "HPV11.2","HPV16","HPV16.1","HPV16.2","HPV18", "HPV18.1", "HPV18.2")], na.rm = TRUE) > 0
dd$HPV4pstat[dd$Inadequate == TRUE] <- "inadequate"
dd$HPV4pstat[dd$HPV4pstat == TRUE] <- "positive"
dd$HPV4pstat[dd$HPV4pstat == FALSE] <- "negative"

dd$HPV9pstat <- rowSums(dd[,c("HPV6", "HPV6.1", "HPV6.2", "HPV11", "HPV11.1", "HPV11.2","HPV16","HPV16.1","HPV16.2","HPV18", "HPV18.1", "HPV18.2", "HPV31", "HPV31.1", "HPV31.2", "HPV33", "HPV33.1", "HPV33.2", "HPV45", "HPV45.1", "HPV45.2", "HPV52", "HPV52.1", "HPV52.2", "HPV58", "HPV58.1", "HPV58.2")], na.rm = TRUE) > 0
dd$HPV9pstat[dd$Inadequate == TRUE] <- "inadequate"
dd$HPV9pstat[dd$HPV9pstat == TRUE] <- "positive"
dd$HPV9pstat[dd$HPV9pstat == FALSE] <- "negative"

Only.16.pro <- by(dd$HPV16stat, dd$yr, function(x) table(x)/length(x)) # -pro is proportion and -nu is numerical
Only.16.nu <- by(dd$HPV16stat, dd$yr, function(x) table(x))
Only.18.pro <- by(dd$HPV18stat, dd$yr, function(x) table(x)/length(x))
Only.18.nu <- by(dd$HPV18stat, dd$yr, function(x) table(x))
Only.16.18.pro <- by(dd$HPV16.18stat, dd$yr, function(x) table(x)/length(x))
Only.16.18.nu <- by(dd$HPV16.18stat, dd$yr, function(x) table(x))
Only.4p.pro <- by(dd$HPV4pstat, dd$yr, function(x) table(x)/length(x))
Only.4p.nu <- by(dd$HPV4pstat, dd$yr, function(x) table(x))
Only.9p.pro <- by(dd$HPV9pstat, dd$yr, function(x) table(x)/length(x))
Only.9p.nu <- by(dd$HPV9pstat, dd$yr, function(x) table(x))
#They all show decreasing patterns

##non-vaccinated types analysis
Non.Vaccine <- paste0("HPV", c(26, 35, 39, 40, 42, 51, 53, 54, 55, 56, 59, 61, 62, 64, 66, 67, 68, 69, 70, 71, 72, 73, 81, 82, 83, 84, "IS39", 89, 26.1, 35.1, 39.1, 40.1, 51.1, 53.1, 54.1, 56.1, 59.1, 66.1, 68.1, 69.1, 70.1, 71.1, 73.1, 74, 74.1, 82.1, 34, 43, 43.1, 44, 44.1, 35.2, 39.2, 40.2, 42.1, 51.2, 53.2, 54.2, 56.2, 59.2, 66.2, 68.2, 70.2, 73.2))
dd$non.vaccine.stat <- rowSums(dd[,Non.Vaccine], na.rm = TRUE) > 0
dd$non.vaccine.stat[dd$Inadequate == TRUE] <- "inadequate"
dd$non.vaccine.stat[dd$non.vaccine.stat == TRUE] <- "positive"
dd$non.vaccine.stat[dd$non.vaccine.stat == FALSE] <- "negative"

non.vaccine.pro <- by(dd$non.vaccine.stat, dd$yr, function(x) table(x)/length(x))
non.vaccine.nu <- by(dd$non.vaccine.stat, dd$yr, function(x) table(x))

Only.16.df <- as.data.frame(do.call("rbind", Only.16.pro))
Only.16.df$year <- 0:9
Only.16.Molten <- melt(Only.16.df, id = c("year"), variable.name = "Status")
Only.16.Plot <- ggplot(Only.16.Molten, aes(x = year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "HPV16 Proportion Status over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Only.18.df <- as.data.frame(do.call("rbind", Only.18.pro))
Only.18.df$year <- 0:9
Only.18.Molten <- melt(Only.18.df, id = c("year"), variable.name = "Status")
Only.18.Plot <- ggplot(Only.18.Molten, aes(x = year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "HPV18 Proportion Status over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Only.16.18.df <- as.data.frame(do.call("rbind", Only.16.18.pro))
Only.16.18.df$year <- 0:9
Only.16.18.Molten <- melt(Only.16.18.df, id = c("year"), variable.name = "Status")
Only.16.18.Plot <- ggplot(Only.16.Molten, aes(x = year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "HPV16 and 18 Proportion Status over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Only.4p.df <- as.data.frame(do.call("rbind", Only.4p.pro))
Only.4p.df$year <- 0:9
Only.4p.Molten <- melt(Only.4p.df, id = c("year"), variable.name = "Status")
Only.4p.Plot <- ggplot(Only.4p.Molten, aes(x = year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "Four-Plex Proportion Status over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Only.9p.df <- as.data.frame(do.call("rbind", Only.9p.pro))
Only.9p.df$year <- 0:9
Only.9p.Molten <- melt(Only.9p.df, id = c("year"), variable.name = "Status")
Only.9p.Plot <- ggplot(Only.9p.Molten, aes(x = year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "HPV16 Proportion Status over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Non.df <- as.data.frame(do.call("rbind", non.vaccine.pro))
Non.df$year <- 0:9
Non.Molten <- melt(Non.df, id = c("year"), variable.name = "Status")
Non.Plot <- ggplot(Non.Molten, aes(x = year, y = value, fill = Status)) + geom_bar(stat = "identity") + labs(title = "Unvaccinnated Proportion Status over Time", x = "Specimen Age per 1 year", y = "Assay Status Proportion")

Only.16.nu.df <- as.data.frame(do.call("rbind", Only.16.nu))
Only.16.nu.df$year <- 0:9
Only.16.nu.Molten <- melt(Only.16.nu.df, id = c("year"), variable.name = "Status")
Only.16.nu.Molten$Proportion <- Only.16.nu.Molten$value
Only.16.Dotplot <- ggplot(Only.16.nu.Molten, aes(x = year, y = Proportion, fill = Status, size = value)) + geom_point(shape = 21) + ggtitle("HPV 16 proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

Only.18.nu.df <- as.data.frame(do.call("rbind", Only.18.nu))
Only.18.nu.df$year <- 0:9
Only.18.nu.Molten <- melt(Only.18.nu.df, id = c("year"), variable.name = "Status")
Only.18.nu.Molten$Proportion <- Only.18.nu.Molten$value
Only.18.Dotplot <- ggplot(Only.18.nu.Molten, aes(x = year, y = Proportion, fill = Status, size = value)) + geom_point(shape = 21) + ggtitle("HPV 18 proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

Only.16.18.nu.df <- as.data.frame(do.call("rbind", Only.16.18.nu))
Only.16.18.nu.df$year <- 0:9
Only.16.18.nu.Molten <- melt(Only.16.18.nu.df, id = c("year"), variable.name = "Status")
Only.16.18.nu.Molten$Proportion <- Only.16.18.nu.Molten$value
Only.16.18.Dotplot <- ggplot(Only.16.18.nu.Molten, aes(x = year, y = Proportion, fill = Status, size = value)) + geom_point(shape = 21) + ggtitle("HPV 16 and 18 proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

Only.4p.nu.df <- as.data.frame(do.call("rbind", Only.4p.nu))
Only.4p.nu.df$year <- 0:9
Only.4p.nu.Molten <- melt(Only.4p.nu.df, id = c("year"), variable.name = "Status")
Only.4p.nu.Molten$Proportion <- Only.4p.nu.Molten$value
Only.4p.Dotplot <- ggplot(Only.4p.nu.Molten, aes(x = year, y = Proportion, fill = Status, size = value)) + geom_point(shape = 21) + ggtitle("Four Plex proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

Only.9p.nu.df <- as.data.frame(do.call("rbind", Only.9p.nu))
Only.9p.nu.df$year <- 0:9
Only.9p.nu.Molten <- melt(Only.9p.nu.df, id = c("year"), variable.name = "Status")
Only.9p.nu.Molten$Proportion <- Only.9p.nu.Molten$value
Only.9p.Dotplot <- ggplot(Only.9p.nu.Molten, aes(x = year, y = Proportion, fill = Status, size = value)) + geom_point(shape = 21) + ggtitle("Nine Plex proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

Non.nu.df <- as.data.frame(do.call("rbind", non.vaccine.nu))
Non.nu.df$year <- 0:9
Non.nu.Molten <- melt(Non.nu.df, id = c("year"), variable.name = "Status")
Non.nu.Molten$Proportion <- Non.nu.Molten$value
Non.Dotplot <- ggplot(Non.nu.Molten, aes(x = year, y = Proportion, fill = Status, size = value)) + geom_point(shape = 21) + ggtitle("Unvaccinnated proportion and size per year") + scale_size(range = c(1,7), breaks = c(0, 10, 25, 50, 100, 200, 500, 1000, 2000, 4000)) + labs(x = "Specimen Age per year", y = "Assay Status Proportion")

