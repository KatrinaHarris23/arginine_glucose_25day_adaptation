##this is the analysis of the 25 day evolution experiment whee PA14 was evolved in biofilm and planktonic environments, and in arginine vs glucose media. 

#This experiment was performed because in all studies that I have not been a part of, both in the cooper lab and outside of the cooper lab, mutations in MorA, wspF, pil genes, and lasR are always seen in evolutionary studies of pseudomonas. Both PA14 and PA01. In the PA14 evolution studies that I have been a part of, meaning the 90 day evolution and the invasion assays, I very rarely see wsp and pil genes, and I never see morA and lasR. The major difference between these studies is that I always work with arginine media as opposed to a sugar/mixed amino acids rich media. 

#I hypothesize that through forcing PA to undergo the alternate metabolism that is arginine metabolism, you are changing the selective pressures to the populations. With this rationalle I designed the experiment to directly test the differences seen when PA14 is forced to grow on arginine vs glucose. I expect that we will see the same mutations observed in other studies in the glucose media (MorA, wsp and pil genes, lasR, CDG reguators in general), however we will see very few in the arginne media. I do not know what exactly the arginine grown populations are adapting to, the previous 90 day evolution experiment had 118 cases of gene level parallelism. This expeirment had the two carbon sources, and 20 replicates each in both planktonic and biofilm environments. Of the resulting 80 populations propagated I selected the first 6 of each group for sequencing (24 total populations). 

#24 populations have been sequenced at days 6, 12, and 25. This is the analysis of those seqeuncing results. I first filter the data (first each population individually, and then looking for mutations found in every population), then I quantify the number of mutations in populatins, alpha diversity, and the bray-curtis dissimilarity looking for differences both within and between the different environments. 

#even numbered populations are lac +
#odd numbered populations are lac - 

#breseq files have been edited to eliminate the following special characters as wll as spaces in all columns except the description.
# , % → ← + Δ


library(data.table)
library(reshape2)

##first import all the data and filter out the ancestral mutations:
######

#read in the ancestral mutations file
ancestors <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/Ancestor_Breseq_Output.csv", stringsAsFactors = F)

#need to split up the ancestor file into the lac+ and the lac- files
KBH5 <- ancestors[(ancestors$Sample == "KBH5_ancestor"),] #this one is the lac -
KBH6 <- ancestors[(ancestors$Sample == "KBH6_ancestor"),] #this one is the lac +


#start with filtering the ancestral mutations

##
#import day 6 sequencing
day6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Day6_sequencing/Day6_Breseq_Output.csv", stringsAsFactors = F)

#filter out all of the ancestral mutations
day6_filtered <- day6[!(day6$Position %in% ancestors$Position),]

#automatically filtering out all of the mutations in commonly seen ribosomal genes. These are filtered out because there are multiple copies around the genome that are almost identical to one another, so reads are mismapped. Used either the description or the gene name, which one just depended on which option was easier for me to grab.
day6_filtered <- day6_filtered[!(day6_filtered$Description == "16S ribosomal RNA/hypothetical protein"),]
day6_filtered <- day6_filtered[!(day6_filtered$Description == "23S ribosomal RNA"),]
day6_filtered <- day6_filtered[!(day6_filtered$Gene == "PA14_RS09995>"),]
day6_filtered <- day6_filtered[!(day6_filtered$Description == "30S ribosomal protein S14/30S ribosomal protein S8"),]
day6_filtered <- day6_filtered[!(day6_filtered$Gene == "PA14_RS26645←/←PA14_RS26650"),]
day6_filtered <- day6_filtered[!(day6_filtered$Description == "5S ribosomal RNA/23S ribosomal RNA"),]

#print out the filtered data set
write.csv(day6_filtered, "/Users/katrina/Desktop/month_long_exp/Day6_sequencing/Day6_noancestor.csv")

###
#and do the same for the day 12 seqeuncing results

day12 <- read.csv("/Users/katrina/Desktop/month_long_exp/Day12_sequencing/Day12_Breseq_Output.csv", stringsAsFactors = F)

day12_filtered <- day12[!(day12$Position %in% ancestors$Position),]

#automatically filtering out all of the mutations in commonly seen ribosomal genes. These are filtered out because there are multiple copies around the genome that are almost identical to one another, so reads are mismapped. Used either the description or the gene name, which one just depended on which option was easier for me to grab.
day12_filtered <- day12_filtered[!(day12_filtered$Description == "16S ribosomal RNA/hypothetical protein"),]
day12_filtered <- day12_filtered[!(day12_filtered$Description == "23S ribosomal RNA"),]
day12_filtered <- day12_filtered[!(day12_filtered$Gene == "PA14_RS09995>"),]
day12_filtered <- day12_filtered[!(day12_filtered$Description == "30S ribosomal protein S14/30S ribosomal protein S8"),]
day12_filtered <- day12_filtered[!(day12_filtered$Gene == "PA14_RS26645←/←PA14_RS26650"),]
day12_filtered <- day12_filtered[!(day12_filtered$Description == "5S ribosomal RNA/23S ribosomal RNA"),]

#print out the filtered data set
write.csv(day12_filtered, "/Users/katrina/Desktop/month_long_exp/Day12_sequencing/Day12_noancestor.csv")

##
#and finally do day 25

day25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Day25_sequencing/Day25_Breseq_Output.csv", stringsAsFactors = F)

day25_filtered <- day25[!(day25$Position %in% ancestors$Position),]

#automatically filtering out all of the mutations in commonly seen ribosomal genes. These are filtered out because there are multiple copies around the genome that are almost identical to one another, so reads are mismapped. Used either the description or the gene name, which one just depended on which option was easier for me to grab.
day25_filtered <- day25_filtered[!(day25_filtered$Description == "16S ribosomal RNA/hypothetical protein"),]
day25_filtered <- day25_filtered[!(day25_filtered$Description == "23S ribosomal RNA"),]
day25_filtered <- day25_filtered[!(day25_filtered$Gene == "PA14_RS09995>"),]
day25_filtered <- day25_filtered[!(day25_filtered$Description == "30S ribosomal protein S14/30S ribosomal protein S8"),]
day25_filtered <- day25_filtered[!(day25_filtered$Gene == "PA14_RS26645←/←PA14_RS26650"),]
day25_filtered <- day25_filtered[!(day25_filtered$Description == "5S ribosomal RNA/23S ribosomal RNA"),]

#print out the filtered data set
write.csv(day25_filtered, "/Users/katrina/Desktop/month_long_exp/Day25_sequencing/Day25_noancestor.csv")

#####
#now I am going to split into time series for each individual population, so that I can filter mutations further. 
######

#reimport the correct data sets so that I can start over from here if needed. 

day6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Day6_sequencing/Day6_noancestor.csv", stringsAsFactors = F)
day6 <- day6[,-1]
day12 <- read.csv("/Users/katrina/Desktop/month_long_exp/Day12_sequencing/Day12_noancestor.csv", stringsAsFactors = F)
day12 <- day12[,-1]
day25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Day25_sequencing/Day25_noancestor.csv", stringsAsFactors = F)
day25 <- day25[,-1]


#starting with the planktonic populations

#start by grabbing the day 6 data

#P1-P6 are planktonic arginine evolved populations. 
p1_day6 <- day6[(day6$Sample == "P1"),] #grabs all of the rows with P1 as the sample name in the day 6 data frame
p1_day6$Sample <- "P1_6" #names the sample so that I can combine all of the P1 samples. this tells me that it is the P1 day 6 sample data. 
p2_day6 <- day6[(day6$Sample == "P2"),]
p2_day6$Sample <- "P2_6"
p3_day6 <- day6[(day6$Sample == "P3"),]
p3_day6$Sample <- "P3_6"
p4_day6 <- day6[(day6$Sample == "P4"),]
p4_day6$Sample <- "P4_6"
p5_day6 <- day6[(day6$Sample == "P5"),]
p5_day6$Sample <- "P5_6"
p6_day6 <- day6[(day6$Sample == "P6"),]
p6_day6$Sample <- "P6_6"

#and these are the planktonic glucose populations. P21-P26
p21_day6 <- day6[(day6$Sample == "P21"),]
p21_day6$Sample <- "P21_6"
p22_day6 <- day6[(day6$Sample == "P22"),]
p22_day6$Sample <- "P22_6"
p23_day6 <- day6[(day6$Sample == "P23"),]
p23_day6$Sample <- "P23_6"
p24_day6 <- day6[(day6$Sample == "P24"),]
p24_day6$Sample <- "P24_6"
p25_day6 <- day6[(day6$Sample == "P25"),]
p25_day6$Sample <- "P25_6"
p26_day6 <- day6[(day6$Sample == "P26"),]
p26_day6$Sample <- "P26_6"

#then get the day 12 mutations

p1_day12 <- day12[(day12$Sample == "P1"),]
p1_day12$Sample <- "P1_12"
p2_day12 <- day12[(day12$Sample == "P2"),]
p2_day12$Sample <- "P2_12"
p3_day12 <- day12[(day12$Sample == "P3"),]
p3_day12$Sample <- "P3_12"
p4_day12 <- day12[(day12$Sample == "P4"),]
p4_day12$Sample <- "P4_12"
p5_day12 <- day12[(day12$Sample == "P5"),]
p5_day12$Sample <- "P5_12"
p6_day12 <- day12[(day12$Sample == "P6"),]
p6_day12$Sample <- "P6_12"

p21_day12 <- day12[(day12$Sample == "P21"),]
p21_day12$Sample <- "P21_12"
p22_day12 <- day12[(day12$Sample == "P22"),]
p22_day12$Sample <- "P22_12"
p23_day12 <- day12[(day12$Sample == "P23"),]
p23_day12$Sample <- "P23_12"
p24_day12 <- day12[(day12$Sample == "P24"),]
p24_day12$Sample <- "P24_12"
p25_day12 <- day12[(day12$Sample == "P25"),]
p25_day12$Sample <- "P25_12"
p26_day12 <- day12[(day12$Sample == "P26"),]
p26_day12$Sample <- "P26_12"

#and finally, the day 25 data

p1_day25 <- day25[(day25$Sample == "P1"),]
p1_day25$Sample <- "P1_25"
p2_day25 <- day25[(day25$Sample == "P2"),]
p2_day25$Sample <- "P2_25"
p3_day25 <- day25[(day25$Sample == "P3"),]
p3_day25$Sample <- "P3_25"
p4_day25 <- day25[(day25$Sample == "P4"),]
p4_day25$Sample <- "P4_25"
p5_day25 <- day25[(day25$Sample == "P5"),]
p5_day25$Sample <- "P5_25"
p6_day25 <- day25[(day25$Sample == "P6"),]
p6_day25$Sample <- "P6_25"

p21_day25 <- day25[(day25$Sample == "P21"),]
p21_day25$Sample <- "P21_25"
p22_day25 <- day25[(day25$Sample == "P22"),]
p22_day25$Sample <- "P22_25"
p23_day25 <- day25[(day25$Sample == "P23"),]
p23_day25$Sample <- "P23_25"
p24_day25 <- day25[(day25$Sample == "P24"),]
p24_day25$Sample <- "P24_25"
p25_day25 <- day25[(day25$Sample == "P25"),]
p25_day25$Sample <- "P25_25"
p26_day25 <- day25[(day25$Sample == "P26"),]
p26_day25$Sample <- "P26_25"


#now for the biofilm populations
b1_day6 <- day6[(day6$Sample == "B1"),] 
b1_day6$Sample <- "B1_6"  
b2_day6 <- day6[(day6$Sample == "B2"),]
b2_day6$Sample <- "B2_6"
b3_day6 <- day6[(day6$Sample == "B3"),]
b3_day6$Sample <- "B3_6"
b4_day6 <- day6[(day6$Sample == "B4"),]
b4_day6$Sample <- "B4_6"
b5_day6 <- day6[(day6$Sample == "B5"),]
b5_day6$Sample <- "B5_6"
b6_day6 <- day6[(day6$Sample == "B6"),]
b6_day6$Sample <- "B6_6"

#and these are the biofilm glucose populations. B21-B26
b21_day6 <- day6[(day6$Sample == "B21"),]
b21_day6$Sample <- "B21_6"
b22_day6 <- day6[(day6$Sample == "B22"),]
b22_day6$Sample <- "B22_6"
b23_day6 <- day6[(day6$Sample == "B23"),]
b23_day6$Sample <- "B23_6"
b24_day6 <- day6[(day6$Sample == "B24"),]
b24_day6$Sample <- "B24_6"
b25_day6 <- day6[(day6$Sample == "B25"),]
b25_day6$Sample <- "B25_6"
b26_day6 <- day6[(day6$Sample == "B26"),]
b26_day6$Sample <- "B26_6"

#then get the day 12 mutations

b1_day12 <- day12[(day12$Sample == "B1"),]
b1_day12$Sample <- "B1_12"
b2_day12 <- day12[(day12$Sample == "B2"),]
b2_day12$Sample <- "B2_12"
b3_day12 <- day12[(day12$Sample == "B3"),]
b3_day12$Sample <- "B3_12"
b4_day12 <- day12[(day12$Sample == "B4"),]
b4_day12$Sample <- "B4_12"
b5_day12 <- day12[(day12$Sample == "B5"),]
b5_day12$Sample <- "B5_12"
b6_day12 <- day12[(day12$Sample == "B6"),]
b6_day12$Sample <- "B6_12"

b21_day12 <- day12[(day12$Sample == "B21"),]
b21_day12$Sample <- "B21_12"
b22_day12 <- day12[(day12$Sample == "B22"),]
b22_day12$Sample <- "B22_12"
b23_day12 <- day12[(day12$Sample == "B23"),]
b23_day12$Sample <- "B23_12"
b24_day12 <- day12[(day12$Sample == "B24"),]
b24_day12$Sample <- "B24_12"
b25_day12 <- day12[(day12$Sample == "B25"),]
b25_day12$Sample <- "B25_12"
b26_day12 <- day12[(day12$Sample == "B26"),]
b26_day12$Sample <- "B26_12"

#and finally, the day 25 data

b1_day25 <- day25[(day25$Sample == "B1"),]
b1_day25$Sample <- "B1_25"
b2_day25 <- day25[(day25$Sample == "B2"),]
b2_day25$Sample <- "B2_25"
b3_day25 <- day25[(day25$Sample == "B3"),]
b3_day25$Sample <- "B3_25"
b4_day25 <- day25[(day25$Sample == "B4"),]
b4_day25$Sample <- "B4_25"
b5_day25 <- day25[(day25$Sample == "B5"),]
b5_day25$Sample <- "B5_25"
b6_day25 <- day25[(day25$Sample == "B6"),]
b6_day25$Sample <- "B6_25"

b21_day25 <- day25[(day25$Sample == "B21"),]
b21_day25$Sample <- "B21_25"
b22_day25 <- day25[(day25$Sample == "B22"),]
b22_day25$Sample <- "B22_25"
b23_day25 <- day25[(day25$Sample == "B23"),]
b23_day25$Sample <- "B23_25"
b24_day25 <- day25[(day25$Sample == "B24"),]
b24_day25$Sample <- "B24_25"
b25_day25 <- day25[(day25$Sample == "B25"),]
b25_day25$Sample <- "B25_25"
b26_day25 <- day25[(day25$Sample == "B26"),]
b26_day25$Sample <- "B26_25"



#now, combine the 3 days into the same file, so that you have one data frame for each populaiton

#planktonic arginine populations
P1 <- rbind(p1_day6, p1_day12, p1_day25)
P2 <- rbind(p2_day6, p2_day12, p2_day25)
P3 <- rbind(p3_day6, p3_day12, p3_day25)
P4 <- rbind(p4_day6, p4_day12, p4_day25)
P5 <- rbind(p5_day6, p5_day12, p5_day25)
P6 <- rbind(p6_day6, p6_day12, p6_day25)

#planktonic glucose populations
P21 <- rbind(p21_day6, p21_day12, p21_day25)
P22 <- rbind(p22_day6, p22_day12, p22_day25)
P23 <- rbind(p23_day6, p23_day12, p23_day25)
P24 <- rbind(p24_day6, p24_day12, p24_day25)
P25 <- rbind(p25_day6, p25_day12, p25_day25)
P26 <- rbind(p26_day6, p26_day12, p26_day25)

#biofilm arginine populations
B1 <- rbind(b1_day6, b1_day12, b1_day25)
B2 <- rbind(b2_day6, b2_day12, b2_day25)
B3 <- rbind(b3_day6, b3_day12, b3_day25)
B4 <- rbind(b4_day6, b4_day12, b4_day25)
B5 <- rbind(b5_day6, b5_day12, b5_day25)
B6 <- rbind(b6_day6, b6_day12, b6_day25)

#biofilm glucose populations
B21 <- rbind(b21_day6, b21_day12, b21_day25)
B22 <- rbind(b22_day6, b22_day12, b22_day25)
B23 <- rbind(b23_day6, b23_day12, b23_day25)
B24 <- rbind(b24_day6, b24_day12, b24_day25)
B25 <- rbind(b25_day6, b25_day12, b25_day25)
B26 <- rbind(b26_day6, b26_day12, b26_day25)

#Now organize all into a time series
#first include all information in one info column so that I don't lose data when putting in a table
#plnaktonic arginine populations
P1$info <- paste(P1$Position, P1$Mutation, P1$Annotation, P1$Gene, P1$Description, sep =":::")
P2$info <- paste(P2$Position, P2$Mutation, P2$Annotation, P2$Gene, P2$Description, sep =":::")
P3$info <- paste(P3$Position, P3$Mutation, P3$Annotation, P3$Gene, P3$Description, sep =":::")
P4$info <- paste(P4$Position, P4$Mutation, P4$Annotation, P4$Gene, P4$Description, sep =":::")
P5$info <- paste(P5$Position, P5$Mutation, P5$Annotation, P5$Gene, P5$Description, sep =":::")
P6$info <- paste(P6$Position, P6$Mutation, P6$Annotation, P6$Gene, P6$Description, sep =":::")
#planktonic glucose populations
P21$info <- paste(P21$Position, P21$Mutation, P21$Annotation, P21$Gene, P21$Description, sep =":::")
P22$info <- paste(P22$Position, P22$Mutation, P22$Annotation, P22$Gene, P22$Description, sep =":::")
P23$info <- paste(P23$Position, P23$Mutation, P23$Annotation, P23$Gene, P23$Description, sep =":::")
P24$info <- paste(P24$Position, P24$Mutation, P24$Annotation, P24$Gene, P24$Description, sep =":::")
P25$info <- paste(P25$Position, P25$Mutation, P25$Annotation, P25$Gene, P25$Description, sep =":::")
P26$info <- paste(P26$Position, P26$Mutation, P26$Annotation, P26$Gene, P26$Description, sep =":::")
#biofilm arginine populations
B1$info <- paste(B1$Position, B1$Mutation, B1$Annotation, B1$Gene, B1$Description, sep =":::")
B2$info <- paste(B2$Position, B2$Mutation, B2$Annotation, B2$Gene, B2$Description, sep =":::")
B3$info <- paste(B3$Position, B3$Mutation, B3$Annotation, B3$Gene, B3$Description, sep =":::")
B4$info <- paste(B4$Position, B4$Mutation, B4$Annotation, B4$Gene, B4$Description, sep =":::")
B5$info <- paste(B5$Position, B5$Mutation, B5$Annotation, B5$Gene, B5$Description, sep =":::")
B6$info <- paste(B6$Position, B6$Mutation, B6$Annotation, B6$Gene, B6$Description, sep =":::")
#biofilm glucose populations
B21$info <- paste(B21$Position, B21$Mutation, B21$Annotation, B21$Gene, B21$Description, sep =":::")
B22$info <- paste(B22$Position, B22$Mutation, B22$Annotation, B22$Gene, B22$Description, sep =":::")
B23$info <- paste(B23$Position, B23$Mutation, B23$Annotation, B23$Gene, B23$Description, sep =":::")
B24$info <- paste(B24$Position, B24$Mutation, B24$Annotation, B24$Gene, B24$Description, sep =":::")
B25$info <- paste(B25$Position, B25$Mutation, B25$Annotation, B25$Gene, B25$Description, sep =":::")
B26$info <- paste(B26$Position, B26$Mutation, B26$Annotation, B26$Gene, B26$Description, sep =":::")
#######
#put into time series tables for each population
#####

#Planktonic arginine populations
#starting with the P1 population
#next melt the data
m_P1 <- melt(P1, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
#then put it into a table with frequencies over time
cast_P1 <- t(dcast(m_P1, Sample~info, mean, value.var = "value", fill = 0)) #cast the data in the order that I want it
cast_P1 <- as.data.frame(cast_P1, header=T) #make sure it is in a data frame type
colnames(cast_P1) <- as.character(unlist(cast_P1[1,])) #make sure the column names are what I want them to be
cast_P1$"0" <- 0.0 #add in a time 0 
P1_timeseries <- cast_P1[-1,] # the row names are the same as the first row, so I et rid of the first row
P1_col_order <- c("0", "P1_6", "P1_12", "P1_25")
setcolorder(P1_timeseries, P1_col_order)
# P2 population
m_P2 <- melt(P2, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P2 <- t(dcast(m_P2, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P2 <- as.data.frame(cast_P2, header=T)
colnames(cast_P2) <- as.character(unlist(cast_P2[1,])) 
cast_P2$"0" <- 0.0 
P2_timeseries <- cast_P2[-1,]
P2_col_order <- c("0", "P2_6", "P2_12", "P2_25")
setcolorder(P2_timeseries, P2_col_order)
# P3 population
m_P3 <- melt(P3, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P3 <- t(dcast(m_P3, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P3 <- as.data.frame(cast_P3, header=T)
colnames(cast_P3) <- as.character(unlist(cast_P3[1,])) 
cast_P3$"0" <- 0.0 
P3_timeseries <- cast_P3[-1,]
P3_col_order <- c("0", "P3_6", "P3_12", "P3_25")
setcolorder(P3_timeseries, P3_col_order)
# P4 population
m_P4 <- melt(P4, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P4 <- t(dcast(m_P4, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P4 <- as.data.frame(cast_P4, header=T)
colnames(cast_P4) <- as.character(unlist(cast_P4[1,])) 
cast_P4$"0" <- 0.0 
P4_timeseries <- cast_P4[-1,]
P4_col_order <- c("0", "P4_6", "P4_12", "P4_25")
setcolorder(P4_timeseries, P4_col_order)
# P5 population
m_P5 <- melt(P5, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P5 <- t(dcast(m_P5, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P5 <- as.data.frame(cast_P5, header=T)
colnames(cast_P5) <- as.character(unlist(cast_P5[1,])) 
cast_P5$"0" <- 0.0 
P5_timeseries <- cast_P5[-1,]
P5_col_order <- c("0", "P5_6", "P5_12", "P5_25")
setcolorder(P5_timeseries, P5_col_order)
# P6 population
m_P6 <- melt(P6, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P6 <- t(dcast(m_P6, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P6 <- as.data.frame(cast_P6, header=T)
colnames(cast_P6) <- as.character(unlist(cast_P6[1,])) 
cast_P6$"0" <- 0.0 
P6_timeseries <- cast_P6[-1,]
P6_col_order <- c("0", "P6_6", "P6_12", "P6_25")
setcolorder(P6_timeseries, P6_col_order)


#Planktonic glucose populations
# P21 population
m_P21 <- melt(P21, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P21 <- t(dcast(m_P21, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P21 <- as.data.frame(cast_P21, header=T)
colnames(cast_P21) <- as.character(unlist(cast_P21[1,])) 
cast_P21$"0" <- 0.0 
P21_timeseries <- cast_P21[-1,]
P21_col_order <- c("0", "P21_6", "P21_12", "P21_25")
setcolorder(P21_timeseries, P21_col_order)
# P22 population
m_P22 <- melt(P22, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P22 <- t(dcast(m_P22, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P22 <- as.data.frame(cast_P22, header=T)
colnames(cast_P22) <- as.character(unlist(cast_P22[1,])) 
cast_P22$"0" <- 0.0 
P22_timeseries <- cast_P22[-1,]
P22_col_order <- c("0", "P22_6", "P22_12", "P22_25")
setcolorder(P22_timeseries, P22_col_order)
# P23 population
m_P23 <- melt(P23, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P23 <- t(dcast(m_P23, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P23 <- as.data.frame(cast_P23, header=T)
colnames(cast_P23) <- as.character(unlist(cast_P23[1,])) 
cast_P23$"0" <- 0.0 
P23_timeseries <- cast_P23[-1,]
P23_col_order <- c("0", "P23_6", "P23_12", "P23_25")
setcolorder(P23_timeseries, P23_col_order)
# P24 population
m_P24 <- melt(P24, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P24 <- t(dcast(m_P24, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P24 <- as.data.frame(cast_P24, header=T)
colnames(cast_P24) <- as.character(unlist(cast_P24[1,])) 
cast_P24$"0" <- 0.0 
P24_timeseries <- cast_P24[-1,]
P24_col_order <- c("0", "P24_6", "P24_12", "P24_25")
setcolorder(P24_timeseries, P24_col_order)
# P25 population
m_P25 <- melt(P25, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P25 <- t(dcast(m_P25, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P25 <- as.data.frame(cast_P25, header=T)
colnames(cast_P25) <- as.character(unlist(cast_P25[1,])) 
cast_P25$"0" <- 0.0 
P25_timeseries <- cast_P25[-1,]
P25_col_order <- c("0", "P25_6", "P25_12", "P25_25")
setcolorder(P25_timeseries, P25_col_order)
# P26 population
m_P26 <- melt(P26, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_P26 <- t(dcast(m_P26, Sample~info, mean, value.var = "value", fill = 0)) 
cast_P26 <- as.data.frame(cast_P26, header=T)
colnames(cast_P26) <- as.character(unlist(cast_P26[1,])) 
cast_P26$"0" <- 0.0 
P26_timeseries <- cast_P26[-1,]
P26_col_order <- c("0", "P26_6", "P26_12", "P26_25")
setcolorder(P26_timeseries, P26_col_order)


#Biofilm arginine populations
#B1 population
m_B1 <- melt(B1, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B1 <- t(dcast(m_B1, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B1 <- as.data.frame(cast_B1, header=T) 
colnames(cast_B1) <- as.character(unlist(cast_B1[1,])) 
cast_B1$"0" <- 0.0
B1_timeseries <- cast_B1[-1,]
B1_col_order <- c("0", "B1_6", "B1_12", "B1_25")
setcolorder(B1_timeseries, B1_col_order)
# B2 population
m_B2 <- melt(B2, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B2 <- t(dcast(m_B2, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B2 <- as.data.frame(cast_B2, header=T)
colnames(cast_B2) <- as.character(unlist(cast_B2[1,])) 
cast_B2$"0" <- 0.0 
B2_timeseries <- cast_B2[-1,]
B2_col_order <- c("0", "B2_6", "B2_12", "B2_25")
setcolorder(B2_timeseries, B2_col_order)
# B3 population
m_B3 <- melt(B3, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B3 <- t(dcast(m_B3, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B3 <- as.data.frame(cast_B3, header=T)
colnames(cast_B3) <- as.character(unlist(cast_B3[1,])) 
cast_B3$"0" <- 0.0 
B3_timeseries <- cast_B3[-1,]
B3_col_order <- c("0", "B3_6", "B3_12", "B3_25")
setcolorder(B3_timeseries, B3_col_order)
# B4 population
m_B4 <- melt(B4, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B4 <- t(dcast(m_B4, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B4 <- as.data.frame(cast_B4, header=T)
colnames(cast_B4) <- as.character(unlist(cast_B4[1,])) 
cast_B4$"0" <- 0.0 
B4_timeseries <- cast_B4[-1,]
B4_col_order <- c("0", "B4_6", "B4_12", "B4_25")
setcolorder(B4_timeseries, B4_col_order)
# B5 population
m_B5 <- melt(B5, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B5 <- t(dcast(m_B5, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B5 <- as.data.frame(cast_B5, header=T)
colnames(cast_B5) <- as.character(unlist(cast_B5[1,])) 
cast_B5$"0" <- 0.0 
B5_timeseries <- cast_B5[-1,]
B5_col_order <- c("0", "B5_6", "B5_12", "B5_25")
setcolorder(B5_timeseries, B5_col_order)
# B6 population
m_B6 <- melt(B6, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B6 <- t(dcast(m_B6, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B6 <- as.data.frame(cast_B6, header=T)
colnames(cast_B6) <- as.character(unlist(cast_B6[1,])) 
cast_B6$"0" <- 0.0 
B6_timeseries <- cast_B6[-1,]
B6_col_order <- c("0", "B6_6", "B6_12", "B6_25")
setcolorder(B6_timeseries, B6_col_order)


#Biofilm glucose populations
# B21 population
m_B21 <- melt(B21, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B21 <- t(dcast(m_B21, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B21 <- as.data.frame(cast_B21, header=T)
colnames(cast_B21) <- as.character(unlist(cast_B21[1,])) 
cast_B21$"0" <- 0.0 
B21_timeseries <- cast_B21[-1,]
B21_col_order <- c("0", "B21_6", "B21_12", "B21_25")
setcolorder(B21_timeseries, B21_col_order)
# B22 population
m_B22 <- melt(B22, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B22 <- t(dcast(m_B22, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B22 <- as.data.frame(cast_B22, header=T)
colnames(cast_B22) <- as.character(unlist(cast_B22[1,])) 
cast_B22$"0" <- 0.0 
B22_timeseries <- cast_B22[-1,]
B22_col_order <- c("0", "B22_6", "B22_12", "B22_25")
setcolorder(B22_timeseries, B22_col_order)
# B23 population
m_B23 <- melt(B23, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B23 <- t(dcast(m_B23, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B23 <- as.data.frame(cast_B23, header=T)
colnames(cast_B23) <- as.character(unlist(cast_B23[1,])) 
cast_B23$"0" <- 0.0 
B23_timeseries <- cast_B23[-1,]
B23_col_order <- c("0", "B23_6", "B23_12", "B23_25")
setcolorder(B23_timeseries, B23_col_order)
# B24 population
m_B24 <- melt(B24, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B24 <- t(dcast(m_B24, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B24 <- as.data.frame(cast_B24, header=T)
colnames(cast_B24) <- as.character(unlist(cast_B24[1,])) 
cast_B24$"0" <- 0.0 
B24_timeseries <- cast_B24[-1,]
B24_col_order <- c("0", "B24_6", "B24_12", "B24_25")
setcolorder(B24_timeseries, B24_col_order)
# B25 population
m_B25 <- melt(B25, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B25 <- t(dcast(m_B25, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B25 <- as.data.frame(cast_B25, header=T)
colnames(cast_B25) <- as.character(unlist(cast_B25[1,])) 
cast_B25$"0" <- 0.0 
B25_timeseries <- cast_B25[-1,]
B25_col_order <- c("0", "B25_6", "B25_12", "B25_25")
setcolorder(B25_timeseries, B25_col_order)
# B26 population
m_B26 <- melt(B26, id=c("Sample","Evidence","Position","Mutation","Annotation","Gene","Description","info"), measure.vars = c("Frequency"))
cast_B26 <- t(dcast(m_B26, Sample~info, mean, value.var = "value", fill = 0)) 
cast_B26 <- as.data.frame(cast_B26, header=T)
colnames(cast_B26) <- as.character(unlist(cast_B26[1,])) 
cast_B26$"0" <- 0.0 
B26_timeseries <- cast_B26[-1,]
B26_col_order <- c("0", "B26_6", "B26_12", "B26_25")
setcolorder(B26_timeseries, B26_col_order)

######
#need to extract the info from the row names
######
P1_info <- colsplit(rownames(P1_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P1_timeseries$Position <- P1_info$Position
P1_timeseries$Mutation <- P1_info$Mutation
P1_timeseries$Annotation <- P1_info$Annotation
P1_timeseries$Gene <- P1_info$Gene
P1_timeseries$Description <- P1_info$Description

P2_info <- colsplit(rownames(P2_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P2_timeseries$Position <- P2_info$Position
P2_timeseries$Mutation <- P2_info$Mutation
P2_timeseries$Annotation <- P2_info$Annotation
P2_timeseries$Gene <- P2_info$Gene
P2_timeseries$Description <- P2_info$Description

P3_info <- colsplit(rownames(P3_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P3_timeseries$Position <- P3_info$Position
P3_timeseries$Mutation <- P3_info$Mutation
P3_timeseries$Annotation <- P3_info$Annotation
P3_timeseries$Gene <- P3_info$Gene
P3_timeseries$Description <- P3_info$Description

P4_info <- colsplit(rownames(P4_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P4_timeseries$Position <- P4_info$Position
P4_timeseries$Mutation <- P4_info$Mutation
P4_timeseries$Annotation <- P4_info$Annotation
P4_timeseries$Gene <- P4_info$Gene
P4_timeseries$Description <- P4_info$Description

P5_info <- colsplit(rownames(P5_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P5_timeseries$Position <- P5_info$Position
P5_timeseries$Mutation <- P5_info$Mutation
P5_timeseries$Annotation <- P5_info$Annotation
P5_timeseries$Gene <- P5_info$Gene
P5_timeseries$Description <- P5_info$Description

P6_info <- colsplit(rownames(P6_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P6_timeseries$Position <- P6_info$Position
P6_timeseries$Mutation <- P6_info$Mutation
P6_timeseries$Annotation <- P6_info$Annotation
P6_timeseries$Gene <- P6_info$Gene
P6_timeseries$Description <- P6_info$Description

P21_info <- colsplit(rownames(P21_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P21_timeseries$Position <- P21_info$Position
P21_timeseries$Mutation <- P21_info$Mutation
P21_timeseries$Annotation <- P21_info$Annotation
P21_timeseries$Gene <- P21_info$Gene
P21_timeseries$Description <- P21_info$Description

P22_info <- colsplit(rownames(P22_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P22_timeseries$Position <- P22_info$Position
P22_timeseries$Mutation <- P22_info$Mutation
P22_timeseries$Annotation <- P22_info$Annotation
P22_timeseries$Gene <- P22_info$Gene
P22_timeseries$Description <- P22_info$Description

P23_info <- colsplit(rownames(P23_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P23_timeseries$Position <- P23_info$Position
P23_timeseries$Mutation <- P23_info$Mutation
P23_timeseries$Annotation <- P23_info$Annotation
P23_timeseries$Gene <- P23_info$Gene
P23_timeseries$Description <- P23_info$Description

P24_info <- colsplit(rownames(P24_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P24_timeseries$Position <- P24_info$Position
P24_timeseries$Mutation <- P24_info$Mutation
P24_timeseries$Annotation <- P24_info$Annotation
P24_timeseries$Gene <- P24_info$Gene
P24_timeseries$Description <- P24_info$Description

P25_info <- colsplit(rownames(P25_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P25_timeseries$Position <- P25_info$Position
P25_timeseries$Mutation <- P25_info$Mutation
P25_timeseries$Annotation <- P25_info$Annotation
P25_timeseries$Gene <- P25_info$Gene
P25_timeseries$Description <- P25_info$Description

P26_info <- colsplit(rownames(P26_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
P26_timeseries$Position <- P26_info$Position
P26_timeseries$Mutation <- P26_info$Mutation
P26_timeseries$Annotation <- P26_info$Annotation
P26_timeseries$Gene <- P26_info$Gene
P26_timeseries$Description <- P26_info$Description

B1_info <- colsplit(rownames(B1_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B1_timeseries$Position <- B1_info$Position
B1_timeseries$Mutation <- B1_info$Mutation
B1_timeseries$Annotation <- B1_info$Annotation
B1_timeseries$Gene <- B1_info$Gene
B1_timeseries$Description <- B1_info$Description

B2_info <- colsplit(rownames(B2_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B2_timeseries$Position <- B2_info$Position
B2_timeseries$Mutation <- B2_info$Mutation
B2_timeseries$Annotation <- B2_info$Annotation
B2_timeseries$Gene <- B2_info$Gene
B2_timeseries$Description <- B2_info$Description

B3_info <- colsplit(rownames(B3_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B3_timeseries$Position <- B3_info$Position
B3_timeseries$Mutation <- B3_info$Mutation
B3_timeseries$Annotation <- B3_info$Annotation
B3_timeseries$Gene <- B3_info$Gene
B3_timeseries$Description <- B3_info$Description

B4_info <- colsplit(rownames(B4_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B4_timeseries$Position <- B4_info$Position
B4_timeseries$Mutation <- B4_info$Mutation
B4_timeseries$Annotation <- B4_info$Annotation
B4_timeseries$Gene <- B4_info$Gene
B4_timeseries$Description <- B4_info$Description

B5_info <- colsplit(rownames(B5_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B5_timeseries$Position <- B5_info$Position
B5_timeseries$Mutation <- B5_info$Mutation
B5_timeseries$Annotation <- B5_info$Annotation
B5_timeseries$Gene <- B5_info$Gene
B5_timeseries$Description <- B5_info$Description

B6_info <- colsplit(rownames(B6_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B6_timeseries$Position <- B6_info$Position
B6_timeseries$Mutation <- B6_info$Mutation
B6_timeseries$Annotation <- B6_info$Annotation
B6_timeseries$Gene <- B6_info$Gene
B6_timeseries$Description <- B6_info$Description

B21_info <- colsplit(rownames(B21_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B21_timeseries$Position <- B21_info$Position
B21_timeseries$Mutation <- B21_info$Mutation
B21_timeseries$Annotation <- B21_info$Annotation
B21_timeseries$Gene <- B21_info$Gene
B21_timeseries$Description <- B21_info$Description

B22_info <- colsplit(rownames(B22_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B22_timeseries$Position <- B22_info$Position
B22_timeseries$Mutation <- B22_info$Mutation
B22_timeseries$Annotation <- B22_info$Annotation
B22_timeseries$Gene <- B22_info$Gene
B22_timeseries$Description <- B22_info$Description

B23_info <- colsplit(rownames(B23_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B23_timeseries$Position <- B23_info$Position
B23_timeseries$Mutation <- B23_info$Mutation
B23_timeseries$Annotation <- B23_info$Annotation


B23_timeseries$Gene <- B23_info$Gene
B23_timeseries$Description <- B23_info$Description

B24_info <- colsplit(rownames(B24_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B24_timeseries$Position <- B24_info$Position
B24_timeseries$Mutation <- B24_info$Mutation
B24_timeseries$Annotation <- B24_info$Annotation
B24_timeseries$Gene <- B24_info$Gene
B24_timeseries$Description <- B24_info$Description

B25_info <- colsplit(rownames(B25_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B25_timeseries$Position <- B25_info$Position
B25_timeseries$Mutation <- B25_info$Mutation
B25_timeseries$Annotation <- B25_info$Annotation
B25_timeseries$Gene <- B25_info$Gene
B25_timeseries$Description <- B25_info$Description

B26_info <- colsplit(rownames(B26_timeseries), ":::", names = c("Position", "Mutation", "Annotation", "Gene", "Description"))
B26_timeseries$Position <- B26_info$Position
B26_timeseries$Mutation <- B26_info$Mutation
B26_timeseries$Annotation <- B26_info$Annotation
B26_timeseries$Gene <- B26_info$Gene
B26_timeseries$Description <- B26_info$Description

######
#print out the unfiltered time series 
######
write.csv(P1_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P1.csv")
write.csv(P2_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P2.csv")
write.csv(P3_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P3.csv")
write.csv(P4_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P4.csv")
write.csv(P5_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P5.csv")
write.csv(P6_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P6.csv")

write.csv(P21_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P1.csv")
write.csv(P22_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P22.csv")
write.csv(P23_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P23.csv")
write.csv(P24_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P24.csv")
write.csv(P25_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P25.csv")
write.csv(P26_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P26.csv")


write.csv(B1_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B1.csv")
write.csv(B2_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B2.csv")
write.csv(B3_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B3.csv")
write.csv(B4_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B4.csv")
write.csv(B5_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B5.csv")
write.csv(B6_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B6.csv")

write.csv(B21_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B21.csv")
write.csv(B22_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B22.csv")
write.csv(B23_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B23.csv")
write.csv(B24_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B24.csv")
write.csv(B25_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B25.csv")
write.csv(B26_timeseries, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B26.csv")

######


#read in the unfiltered time series
######
P1_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P1.csv", stringsAsFactors = F)
P2_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P2.csv", stringsAsFactors = F)
P3_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P3.csv", stringsAsFactors = F)
P4_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P4.csv", stringsAsFactors = F)
P5_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P5.csv", stringsAsFactors = F)
P6_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P6.csv", stringsAsFactors = F)

P21_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P1.csv", stringsAsFactors = F)
P22_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P22.csv", stringsAsFactors = F)
P23_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P23.csv", stringsAsFactors = F)
P24_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P24.csv", stringsAsFactors = F)
P25_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P25.csv", stringsAsFactors = F)
P26_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/P26.csv", stringsAsFactors = F)

B1_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B1.csv", stringsAsFactors = F)
B2_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B2.csv", stringsAsFactors = F)
B3_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B3.csv", stringsAsFactors = F)
B4_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B4.csv", stringsAsFactors = F)
B5_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B5.csv", stringsAsFactors = F)
B6_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B6.csv", stringsAsFactors = F)

B21_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B21.csv", stringsAsFactors = F)
B22_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B22.csv", stringsAsFactors = F)
B23_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B23.csv", stringsAsFactors = F)
B24_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B24.csv", stringsAsFactors = F)
B25_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B25.csv", stringsAsFactors = F)
B26_timeseries <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/timeseries_unfiltered/B26.csv", stringsAsFactors = F)
#######
#make a matrix of the number of mutations in the unfiltered populations (well, only filtered for ancestral file mutations)

#######
unfiltered_mutations <- as.data.frame(matrix(ncol = 2, byrow = 2, 
                                             c("Population", "mutations",
                                               "P1", nrow(P1_timeseries),
                                               "P2", nrow(P2_timeseries),
                                               "P3", nrow(P3_timeseries),
                                               "P4", nrow(P4_timeseries),
                                               "P5", nrow(P5_timeseries),
                                               "P6", nrow(P6_timeseries),
                                               
                                               "P21", nrow(P21_timeseries),
                                               "P22", nrow(P22_timeseries),
                                               "P23", nrow(P23_timeseries),
                                               "P24", nrow(P24_timeseries),
                                               "P25", nrow(P25_timeseries),
                                               "P26", nrow(P26_timeseries),
                                               
                                               "B1", nrow(B1_timeseries),
                                               "B2", nrow(B2_timeseries),
                                               "B3", nrow(B3_timeseries),
                                               "B4", nrow(B4_timeseries),
                                               "B5", nrow(B5_timeseries),
                                               "B6", nrow(B6_timeseries),
                                               
                                               "B21", nrow(B21_timeseries),
                                               "B22", nrow(B22_timeseries),
                                               "B23", nrow(B23_timeseries),
                                               "B24", nrow(B24_timeseries),
                                               "B25", nrow(B25_timeseries),
                                               "B26", nrow(B26_timeseries))))
write.csv(unfiltered_mutations, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/unfiltered_numbers_.csv")

######
#now to filter 

#######
#make all frequencies type numberic
#######
P1_timeseries_numeric <- as.data.frame(apply(P1_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P1_timeseries[,2:5] <- P1_timeseries_numeric
P2_timeseries_numeric <- as.data.frame(apply(P2_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P2_timeseries[,2:5] <- P2_timeseries_numeric
P3_timeseries_numeric <- as.data.frame(apply(P3_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P3_timeseries[,2:5] <- P3_timeseries_numeric
P4_timeseries_numeric <- as.data.frame(apply(P4_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P4_timeseries[,2:5] <- P4_timeseries_numeric
P5_timeseries_numeric <- as.data.frame(apply(P5_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P5_timeseries[,2:5] <- P5_timeseries_numeric
P6_timeseries_numeric <- as.data.frame(apply(P6_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P6_timeseries[,2:5] <- P6_timeseries_numeric

P21_timeseries_numeric <- as.data.frame(apply(P21_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P21_timeseries[,2:5] <- P21_timeseries_numeric
P22_timeseries_numeric <- as.data.frame(apply(P22_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P22_timeseries[,2:5] <- P22_timeseries_numeric
P23_timeseries_numeric <- as.data.frame(apply(P23_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P23_timeseries[,2:5] <- P23_timeseries_numeric
P24_timeseries_numeric <- as.data.frame(apply(P24_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P24_timeseries[,2:5] <- P24_timeseries_numeric
P25_timeseries_numeric <- as.data.frame(apply(P25_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P25_timeseries[,2:5] <- P25_timeseries_numeric
P26_timeseries_numeric <- as.data.frame(apply(P26_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
P26_timeseries[,2:5] <- P26_timeseries_numeric


B1_timeseries_numeric <- as.data.frame(apply(B1_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B1_timeseries[,2:5] <- B1_timeseries_numeric
B2_timeseries_numeric <- as.data.frame(apply(B2_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B2_timeseries[,2:5] <- B2_timeseries_numeric
B3_timeseries_numeric <- as.data.frame(apply(B3_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B3_timeseries[,2:5] <- B3_timeseries_numeric
B4_timeseries_numeric <- as.data.frame(apply(B4_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B4_timeseries[,2:5] <- B4_timeseries_numeric
B5_timeseries_numeric <- as.data.frame(apply(B5_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B5_timeseries[,2:5] <- B5_timeseries_numeric
B6_timeseries_numeric <- as.data.frame(apply(B6_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B6_timeseries[,2:5] <- B6_timeseries_numeric

B21_timeseries_numeric <- as.data.frame(apply(B21_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B21_timeseries[,2:5] <- B21_timeseries_numeric
B22_timeseries_numeric <- as.data.frame(apply(B22_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B22_timeseries[,2:5] <- B22_timeseries_numeric
B23_timeseries_numeric <- as.data.frame(apply(B23_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B23_timeseries[,2:5] <- B23_timeseries_numeric
B24_timeseries_numeric <- as.data.frame(apply(B24_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B24_timeseries[,2:5] <- B24_timeseries_numeric
B25_timeseries_numeric <- as.data.frame(apply(B25_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B25_timeseries[,2:5] <- B25_timeseries_numeric
B26_timeseries_numeric <- as.data.frame(apply(B26_timeseries[,2:5],2,function(x) as.numeric(as.character(x))))
B26_timeseries[,2:5] <- B26_timeseries_numeric

#####
#create a count column, which counts all time points that are not 0
#######
P1_timeseries$Count <- rowSums(P1_timeseries[,3:5] != 0.0)
P2_timeseries$Count <- rowSums(P2_timeseries[,3:5] != 0.0)
P3_timeseries$Count <- rowSums(P3_timeseries[,3:5] != 0.0)
P4_timeseries$Count <- rowSums(P4_timeseries[,3:5] != 0.0)
P5_timeseries$Count <- rowSums(P5_timeseries[,3:5] != 0.0)
P6_timeseries$Count <- rowSums(P6_timeseries[,3:5] != 0.0)

P21_timeseries$Count <- rowSums(P21_timeseries[,3:5] != 0.0)
P22_timeseries$Count <- rowSums(P22_timeseries[,3:5] != 0.0)
P23_timeseries$Count <- rowSums(P23_timeseries[,3:5] != 0.0)
P24_timeseries$Count <- rowSums(P24_timeseries[,3:5] != 0.0)
P25_timeseries$Count <- rowSums(P25_timeseries[,3:5] != 0.0)
P26_timeseries$Count <- rowSums(P26_timeseries[,3:5] != 0.0)

B1_timeseries$Count <- rowSums(B1_timeseries[,3:5] != 0.0)
B2_timeseries$Count <- rowSums(B2_timeseries[,3:5] != 0.0)
B3_timeseries$Count <- rowSums(B3_timeseries[,3:5] != 0.0)
B4_timeseries$Count <- rowSums(B4_timeseries[,3:5] != 0.0)
B5_timeseries$Count <- rowSums(B5_timeseries[,3:5] != 0.0)
B6_timeseries$Count <- rowSums(B6_timeseries[,3:5] != 0.0)

B21_timeseries$Count <- rowSums(B21_timeseries[,3:5] != 0.0)
B22_timeseries$Count <- rowSums(B22_timeseries[,3:5] != 0.0)
B23_timeseries$Count <- rowSums(B23_timeseries[,3:5] != 0.0)
B24_timeseries$Count <- rowSums(B24_timeseries[,3:5] != 0.0)
B25_timeseries$Count <- rowSums(B25_timeseries[,3:5] != 0.0)
B26_timeseries$Count <- rowSums(B26_timeseries[,3:5] != 0.0)

#####
#create a column with the maximum frequency each mutation is seen over time
#######
P1_timeseries$Max <- apply(P1_timeseries[,3:5], 1, max)
P2_timeseries$Max <- apply(P2_timeseries[,3:5], 1, max)
P3_timeseries$Max <- apply(P3_timeseries[,3:5], 1, max)
P4_timeseries$Max <- apply(P4_timeseries[,3:5], 1, max)
P5_timeseries$Max <- apply(P5_timeseries[,3:5], 1, max)
P6_timeseries$Max <- apply(P6_timeseries[,3:5], 1, max)

P21_timeseries$Max <- apply(P21_timeseries[,3:5], 1, max)
P22_timeseries$Max <- apply(P22_timeseries[,3:5], 1, max)
P23_timeseries$Max <- apply(P23_timeseries[,3:5], 1, max)
P24_timeseries$Max <- apply(P24_timeseries[,3:5], 1, max)
P25_timeseries$Max <- apply(P25_timeseries[,3:5], 1, max)
P26_timeseries$Max <- apply(P26_timeseries[,3:5], 1, max)

B1_timeseries$Max <- apply(B1_timeseries[,3:5], 1, max)
B2_timeseries$Max <- apply(B2_timeseries[,3:5], 1, max)
B3_timeseries$Max <- apply(B3_timeseries[,3:5], 1, max)
B4_timeseries$Max <- apply(B4_timeseries[,3:5], 1, max)
B5_timeseries$Max <- apply(B5_timeseries[,3:5], 1, max)
B6_timeseries$Max <- apply(B6_timeseries[,3:5], 1, max)

B21_timeseries$Max <- apply(B21_timeseries[,3:5], 1, max)
B22_timeseries$Max <- apply(B22_timeseries[,3:5], 1, max)
B23_timeseries$Max <- apply(B23_timeseries[,3:5], 1, max)
B24_timeseries$Max <- apply(B24_timeseries[,3:5], 1, max)
B25_timeseries$Max <- apply(B25_timeseries[,3:5], 1, max)
B26_timeseries$Max <- apply(B26_timeseries[,3:5], 1, max)
#####
#create a column with the minimum frequency each mutation is seen over time
#######
P1_timeseries$Min <- apply(P1_timeseries[,3:5], 1, min)
P2_timeseries$Min <- apply(P2_timeseries[,3:5], 1, min)
P3_timeseries$Min <- apply(P3_timeseries[,3:5], 1, min)
P4_timeseries$Min <- apply(P4_timeseries[,3:5], 1, min)
P5_timeseries$Min <- apply(P5_timeseries[,3:5], 1, min)
P6_timeseries$Min <- apply(P6_timeseries[,3:5], 1, min)

P21_timeseries$Min <- apply(P21_timeseries[,3:5], 1, min)
P22_timeseries$Min <- apply(P22_timeseries[,3:5], 1, min)
P23_timeseries$Min <- apply(P23_timeseries[,3:5], 1, min)
P24_timeseries$Min <- apply(P24_timeseries[,3:5], 1, min)
P25_timeseries$Min <- apply(P25_timeseries[,3:5], 1, min)
P26_timeseries$Min <- apply(P26_timeseries[,3:5], 1, min)

B1_timeseries$Min <- apply(B1_timeseries[,3:5], 1, min)
B2_timeseries$Min <- apply(B2_timeseries[,3:5], 1, min)
B3_timeseries$Min <- apply(B3_timeseries[,3:5], 1, min)
B4_timeseries$Min <- apply(B4_timeseries[,3:5], 1, min)
B5_timeseries$Min <- apply(B5_timeseries[,3:5], 1, min)
B6_timeseries$Min <- apply(B6_timeseries[,3:5], 1, min)

B21_timeseries$Min <- apply(B21_timeseries[,3:5], 1, min)
B22_timeseries$Min <- apply(B22_timeseries[,3:5], 1, min)
B23_timeseries$Min <- apply(B23_timeseries[,3:5], 1, min)
B24_timeseries$Min <- apply(B24_timeseries[,3:5], 1, min)
B25_timeseries$Min <- apply(B25_timeseries[,3:5], 1, min)
B26_timeseries$Min <- apply(B26_timeseries[,3:5], 1, min)
######
#create the final column which is the total change in frequnecy seen over the course of the experiment.
#####
P1_timeseries$Change <- P1_timeseries$Max - P1_timeseries$Min
P2_timeseries$Change <- P2_timeseries$Max - P2_timeseries$Min
P3_timeseries$Change <- P3_timeseries$Max - P3_timeseries$Min
P4_timeseries$Change <- P4_timeseries$Max - P4_timeseries$Min
P5_timeseries$Change <- P5_timeseries$Max - P5_timeseries$Min
P6_timeseries$Change <- P6_timeseries$Max - P6_timeseries$Min

P21_timeseries$Change <- P21_timeseries$Max - P21_timeseries$Min
P22_timeseries$Change <- P22_timeseries$Max - P22_timeseries$Min
P23_timeseries$Change <- P23_timeseries$Max - P23_timeseries$Min
P24_timeseries$Change <- P24_timeseries$Max - P24_timeseries$Min
P25_timeseries$Change <- P25_timeseries$Max - P25_timeseries$Min
P26_timeseries$Change <- P26_timeseries$Max - P26_timeseries$Min


B1_timeseries$Change <- B1_timeseries$Max - B1_timeseries$Min
B2_timeseries$Change <- B2_timeseries$Max - B2_timeseries$Min
B3_timeseries$Change <- B3_timeseries$Max - B3_timeseries$Min
B4_timeseries$Change <- B4_timeseries$Max - B4_timeseries$Min
B5_timeseries$Change <- B5_timeseries$Max - B5_timeseries$Min
B6_timeseries$Change <- B6_timeseries$Max - B6_timeseries$Min

B21_timeseries$Change <- B21_timeseries$Max - B21_timeseries$Min
B22_timeseries$Change <- B22_timeseries$Max - B22_timeseries$Min
B23_timeseries$Change <- B23_timeseries$Max - B23_timeseries$Min
B24_timeseries$Change <- B24_timeseries$Max - B24_timeseries$Min
B25_timeseries$Change <- B25_timeseries$Max - B25_timeseries$Min
B26_timeseries$Change <- B26_timeseries$Max - B26_timeseries$Min

########
#filter1 mutations have to be identified in at least 2 time points. 
#######
#P1_timeseries <- subset(P1_timeseries, P1_timeseries$Count > 1)
#P2_timeseries <- subset(P2_timeseries, P2_timeseries$Count > 1)
#P3_timeseries <- subset(P3_timeseries, P3_timeseries$Count > 1)
#P4_timeseries <- subset(P4_timeseries, P4_timeseries$Count > 1)
#P5_timeseries <- subset(P5_timeseries, P5_timeseries$Count > 1)
#P6_timeseries <- subset(P6_timeseries, P6_timeseries$Count > 1)

#P21_timeseries <- subset(P21_timeseries, P21_timeseries$Count > 1)
#P22_timeseries <- subset(P22_timeseries, P22_timeseries$Count > 1)
#P23_timeseries <- subset(P23_timeseries, P23_timeseries$Count > 1)
#P24_timeseries <- subset(P24_timeseries, P24_timeseries$Count > 1)
#P25_timeseries <- subset(P25_timeseries, P25_timeseries$Count > 1)
#P26_timeseries <- subset(P26_timeseries, P26_timeseries$Count > 1)

#B1_timeseries <- subset(B1_timeseries, B1_timeseries$Count > 1)
#B2_timeseries <- subset(B2_timeseries, B2_timeseries$Count > 1)
#B3_timeseries <- subset(B3_timeseries, B3_timeseries$Count > 1)
#B4_timeseries <- subset(B4_timeseries, B4_timeseries$Count > 1)
#B5_timeseries <- subset(B5_timeseries, B5_timeseries$Count > 1)
#B6_timeseries <- subset(B6_timeseries, B6_timeseries$Count > 1)

#B21_timeseries <- subset(B21_timeseries, B21_timeseries$Count > 1)
#B22_timeseries <- subset(B22_timeseries, B22_timeseries$Count > 1)
#B23_timeseries <- subset(B23_timeseries, B23_timeseries$Count > 1)
#B24_timeseries <- subset(B24_timeseries, B24_timeseries$Count > 1)
#B25_timeseries <- subset(B25_timeseries, B25_timeseries$Count > 1)
#B26_timeseries <- subset(B26_timeseries, B26_timeseries$Count > 1)

#######
#filter2 mutations have to reach at least 10% frequency at one time point, using the max column we made previously for this. This is a much harsher filter than the normal 10%, becuase I am trying to eliminate as many false positives as I can.
########
P1_timeseries <- subset(P1_timeseries, P1_timeseries$Max >= 10)
P2_timeseries <- subset(P2_timeseries, P2_timeseries$Max >= 10)
P3_timeseries <- subset(P3_timeseries, P3_timeseries$Max >= 10)
P4_timeseries <- subset(P4_timeseries, P4_timeseries$Max >= 10)
P5_timeseries <- subset(P5_timeseries, P5_timeseries$Max >= 10)
P6_timeseries <- subset(P6_timeseries, P6_timeseries$Max >= 10)

P21_timeseries <- subset(P21_timeseries, P21_timeseries$Max >= 10)
P22_timeseries <- subset(P22_timeseries, P22_timeseries$Max >= 10)
P23_timeseries <- subset(P23_timeseries, P23_timeseries$Max >= 10)
P24_timeseries <- subset(P24_timeseries, P24_timeseries$Max >= 10)
P25_timeseries <- subset(P25_timeseries, P25_timeseries$Max >= 10)
P26_timeseries <- subset(P26_timeseries, P26_timeseries$Max >= 10)

B1_timeseries <- subset(B1_timeseries, B1_timeseries$Max >= 10)
B2_timeseries <- subset(B2_timeseries, B2_timeseries$Max >= 10)
B3_timeseries <- subset(B3_timeseries, B3_timeseries$Max >= 10)
B4_timeseries <- subset(B4_timeseries, B4_timeseries$Max >= 10)
B5_timeseries <- subset(B5_timeseries, B5_timeseries$Max >= 10)
B6_timeseries <- subset(B6_timeseries, B6_timeseries$Max >= 10)

B21_timeseries <- subset(B21_timeseries, B21_timeseries$Max >= 10)
B22_timeseries <- subset(B22_timeseries, B22_timeseries$Max >= 10)
B23_timeseries <- subset(B23_timeseries, B23_timeseries$Max >= 10)
B24_timeseries <- subset(B24_timeseries, B24_timeseries$Max >= 10)
B25_timeseries <- subset(B25_timeseries, B25_timeseries$Max >= 10)
B26_timeseries <- subset(B26_timeseries, B26_timeseries$Max >= 10)

######
#filter 3, mutations have to change in frequency by at least 10% 
#######
P1_timeseries <- subset(P1_timeseries, P1_timeseries$Change >= 10)
P2_timeseries <- subset(P2_timeseries, P2_timeseries$Change >= 10)
P3_timeseries <- subset(P3_timeseries, P3_timeseries$Change >= 10)
P4_timeseries <- subset(P4_timeseries, P4_timeseries$Change >= 10)
P5_timeseries <- subset(P5_timeseries, P5_timeseries$Change >= 10)
P6_timeseries <- subset(P6_timeseries, P6_timeseries$Change >= 10)

P21_timeseries <- subset(P21_timeseries, P21_timeseries$Change >= 10)
P22_timeseries <- subset(P22_timeseries, P22_timeseries$Change >= 10)
P23_timeseries <- subset(P23_timeseries, P23_timeseries$Change >= 10)
P24_timeseries <- subset(P24_timeseries, P24_timeseries$Change >= 10)
P25_timeseries <- subset(P25_timeseries, P25_timeseries$Change >= 10)
P26_timeseries <- subset(P26_timeseries, P26_timeseries$Change >= 10)

B1_timeseries <- subset(B1_timeseries, B1_timeseries$Change >= 10)
B2_timeseries <- subset(B2_timeseries, B2_timeseries$Change >= 10)
B3_timeseries <- subset(B3_timeseries, B3_timeseries$Change >= 10)
B4_timeseries <- subset(B4_timeseries, B4_timeseries$Change >= 10)
B5_timeseries <- subset(B5_timeseries, B5_timeseries$Change >= 10)
B6_timeseries <- subset(B6_timeseries, B6_timeseries$Change >= 10)

B21_timeseries <- subset(B21_timeseries, B21_timeseries$Change >= 10)
B22_timeseries <- subset(B22_timeseries, B22_timeseries$Change >= 10)
B23_timeseries <- subset(B23_timeseries, B23_timeseries$Change >= 10)
B24_timeseries <- subset(B24_timeseries, B24_timeseries$Change >= 10)
B25_timeseries <- subset(B25_timeseries, B25_timeseries$Change >= 10)
B26_timeseries <- subset(B26_timeseries, B26_timeseries$Change >= 10)



######
#I am filtering out the intergenic mutations and the synonymous mutations from these data tables. I am only concentrating on nonsynonymous. 

######
#remove the intergenic mutations
P1_filter <- P1_timeseries[!(grepl("intergenic", P1_timeseries$Annotation)),]
P2_filter <- P2_timeseries[!(grepl("intergenic", P2_timeseries$Annotation)),]
P3_filter <- P3_timeseries[!(grepl("intergenic", P3_timeseries$Annotation)),]
P4_filter <- P4_timeseries[!(grepl("intergenic", P4_timeseries$Annotation)),]
P5_filter <- P5_timeseries[!(grepl("intergenic", P5_timeseries$Annotation)),]
P6_filter <- P6_timeseries[!(grepl("intergenic", P6_timeseries$Annotation)),]

P21_filter <- P21_timeseries[!(grepl("intergenic", P21_timeseries$Annotation)),]
P22_filter <- P22_timeseries[!(grepl("intergenic", P22_timeseries$Annotation)),]
P23_filter <- P23_timeseries[!(grepl("intergenic", P23_timeseries$Annotation)),]
P24_filter <- P24_timeseries[!(grepl("intergenic", P24_timeseries$Annotation)),]
P25_filter <- P25_timeseries[!(grepl("intergenic", P25_timeseries$Annotation)),]
P26_filter <- P26_timeseries[!(grepl("intergenic", P26_timeseries$Annotation)),]


B1_filter <- B1_timeseries[!(grepl("intergenic", B1_timeseries$Annotation)),]
B2_filter <- B2_timeseries[!(grepl("intergenic", B2_timeseries$Annotation)),]
B3_filter <- B3_timeseries[!(grepl("intergenic", B3_timeseries$Annotation)),]
B4_filter <- B4_timeseries[!(grepl("intergenic", B4_timeseries$Annotation)),]
B5_filter <- B5_timeseries[!(grepl("intergenic", B5_timeseries$Annotation)),]
B6_filter <- B6_timeseries[!(grepl("intergenic", B6_timeseries$Annotation)),]

B21_filter <- B21_timeseries[!(grepl("intergenic", B21_timeseries$Annotation)),]
B22_filter <- B22_timeseries[!(grepl("intergenic", B22_timeseries$Annotation)),]
B23_filter <- B23_timeseries[!(grepl("intergenic", B23_timeseries$Annotation)),]
B24_filter <- B24_timeseries[!(grepl("intergenic", B24_timeseries$Annotation)),]
B25_filter <- B25_timeseries[!(grepl("intergenic", B25_timeseries$Annotation)),]
B26_filter <- B26_timeseries[!(grepl("intergenic", B26_timeseries$Annotation)),]

#And to take out the synonymous mutations so that I am only looking at the nonsynonymous ones

#split the annotatio colum to extract the amino acid
P1_aas <- colsplit(P1_filter$Annotation,"\\(", names = c("AA","DNA"))
P2_aas <- colsplit(P2_filter$Annotation,"\\(", names = c("AA","DNA"))
P3_aas <- colsplit(P3_filter$Annotation,"\\(", names = c("AA","DNA"))
P4_aas <- colsplit(P4_filter$Annotation,"\\(", names = c("AA","DNA"))
P5_aas <- colsplit(P5_filter$Annotation,"\\(", names = c("AA","DNA"))
P6_aas <- colsplit(P6_filter$Annotation,"\\(", names = c("AA","DNA"))

P21_aas <- colsplit(P21_filter$Annotation,"\\(", names = c("AA","DNA"))
P22_aas <- colsplit(P22_filter$Annotation,"\\(", names = c("AA","DNA"))
P23_aas <- colsplit(P23_filter$Annotation,"\\(", names = c("AA","DNA"))
P24_aas <- colsplit(P24_filter$Annotation,"\\(", names = c("AA","DNA"))
P25_aas <- colsplit(P25_filter$Annotation,"\\(", names = c("AA","DNA"))
P26_aas <- colsplit(P26_filter$Annotation,"\\(", names = c("AA","DNA"))

B1_aas <- colsplit(B1_filter$Annotation,"\\(", names = c("AA","DNA"))
B2_aas <- colsplit(B2_filter$Annotation,"\\(", names = c("AA","DNA"))
B3_aas <- colsplit(B3_filter$Annotation,"\\(", names = c("AA","DNA"))
B4_aas <- colsplit(B4_filter$Annotation,"\\(", names = c("AA","DNA"))
B5_aas <- colsplit(B5_filter$Annotation,"\\(", names = c("AA","DNA"))
B6_aas <- colsplit(B6_filter$Annotation,"\\(", names = c("AA","DNA"))

B21_aas <- colsplit(B21_filter$Annotation,"\\(", names = c("AA","DNA"))
B22_aas <- colsplit(B22_filter$Annotation,"\\(", names = c("AA","DNA"))
B23_aas <- colsplit(B23_filter$Annotation,"\\(", names = c("AA","DNA"))
B24_aas <- colsplit(B24_filter$Annotation,"\\(", names = c("AA","DNA"))
B25_aas <- colsplit(B25_filter$Annotation,"\\(", names = c("AA","DNA"))
B26_aas <- colsplit(B26_filter$Annotation,"\\(", names = c("AA","DNA"))

#split the amino acid annotation to see the starting and final aa
P1_splitaas <- colsplit(P1_aas$AA, "[0-9]+", names = c("first", "last"))
P2_splitaas <- colsplit(P2_aas$AA, "[0-9]+", names = c("first", "last"))
P3_splitaas <- colsplit(P3_aas$AA, "[0-9]+", names = c("first", "last"))
P4_splitaas <- colsplit(P4_aas$AA, "[0-9]+", names = c("first", "last"))
P5_splitaas <- colsplit(P5_aas$AA, "[0-9]+", names = c("first", "last"))
P6_splitaas <- colsplit(P6_aas$AA, "[0-9]+", names = c("first", "last"))

P21_splitaas <- colsplit(P21_aas$AA, "[0-9]+", names = c("first", "last"))
P22_splitaas <- colsplit(P22_aas$AA, "[0-9]+", names = c("first", "last"))
P23_splitaas <- colsplit(P23_aas$AA, "[0-9]+", names = c("first", "last"))
P24_splitaas <- colsplit(P24_aas$AA, "[0-9]+", names = c("first", "last"))
P25_splitaas <- colsplit(P25_aas$AA, "[0-9]+", names = c("first", "last"))
P26_splitaas <- colsplit(P26_aas$AA, "[0-9]+", names = c("first", "last"))

B1_splitaas <- colsplit(B1_aas$AA, "[0-9]+", names = c("first", "last"))
B2_splitaas <- colsplit(B2_aas$AA, "[0-9]+", names = c("first", "last"))
B3_splitaas <- colsplit(B3_aas$AA, "[0-9]+", names = c("first", "last"))
B4_splitaas <- colsplit(B4_aas$AA, "[0-9]+", names = c("first", "last"))
B5_splitaas <- colsplit(B5_aas$AA, "[0-9]+", names = c("first", "last"))
B6_splitaas <- colsplit(B6_aas$AA, "[0-9]+", names = c("first", "last"))

B21_splitaas <- colsplit(B21_aas$AA, "[0-9]+", names = c("first", "last"))
B22_splitaas <- colsplit(B22_aas$AA, "[0-9]+", names = c("first", "last"))
B23_splitaas <- colsplit(B23_aas$AA, "[0-9]+", names = c("first", "last"))
B24_splitaas <- colsplit(B24_aas$AA, "[0-9]+", names = c("first", "last"))
B25_splitaas <- colsplit(B25_aas$AA, "[0-9]+", names = c("first", "last"))
B26_splitaas <- colsplit(B26_aas$AA, "[0-9]+", names = c("first", "last"))

P1_filter$Amino_first <- P1_splitaas$first #save just the starting aa
P2_filter$Amino_first <- P2_splitaas$first
P3_filter$Amino_first <- P3_splitaas$first
P4_filter$Amino_first <- P4_splitaas$first
P5_filter$Amino_first <- P5_splitaas$first
P6_filter$Amino_first <- P6_splitaas$first

P21_filter$Amino_first <- P21_splitaas$first
P22_filter$Amino_first <- P22_splitaas$first
P23_filter$Amino_first <- P23_splitaas$first
P24_filter$Amino_first <- P24_splitaas$first
P25_filter$Amino_first <- P25_splitaas$first
P26_filter$Amino_first <- P26_splitaas$first

B1_filter$Amino_first <- B1_splitaas$first
B2_filter$Amino_first <- B2_splitaas$first
B3_filter$Amino_first <- B3_splitaas$first
B4_filter$Amino_first <- B4_splitaas$first
B5_filter$Amino_first <- B5_splitaas$first
B6_filter$Amino_first <- B6_splitaas$first

B21_filter$Amino_first <- B21_splitaas$first
B22_filter$Amino_first <- B22_splitaas$first
B23_filter$Amino_first <- B23_splitaas$first
B24_filter$Amino_first <- B24_splitaas$first
B25_filter$Amino_first <- B25_splitaas$first
B26_filter$Amino_first <- B26_splitaas$first



P1_filter$Amino_last <- P1_splitaas$last #save just the final aa
#filter the dataset so that only mutations that change the aa are saved
P2_filter$Amino_last <- P2_splitaas$last 
P3_filter$Amino_last <- P3_splitaas$last 
P4_filter$Amino_last <- P4_splitaas$last 
P5_filter$Amino_last <- P5_splitaas$last 
P6_filter$Amino_last <- P6_splitaas$last 

P21_filter$Amino_last <- P21_splitaas$last
P22_filter$Amino_last <- P22_splitaas$last 
P23_filter$Amino_last <- P23_splitaas$last 
P24_filter$Amino_last <- P24_splitaas$last 
P25_filter$Amino_last <- P25_splitaas$last 
P26_filter$Amino_last <- P26_splitaas$last 

B1_filter$Amino_last <- B1_splitaas$last
B2_filter$Amino_last <- B2_splitaas$last 
B3_filter$Amino_last <- B3_splitaas$last 
B4_filter$Amino_last <- B4_splitaas$last 
B5_filter$Amino_last <- B5_splitaas$last 
B6_filter$Amino_last <- B6_splitaas$last 

B21_filter$Amino_last <- B21_splitaas$last
B22_filter$Amino_last <- B22_splitaas$last 
B23_filter$Amino_last <- B23_splitaas$last 
B24_filter$Amino_last <- B24_splitaas$last 
B25_filter$Amino_last <- B25_splitaas$last 
B26_filter$Amino_last <- B26_splitaas$last

#AND FINALLY FILTER so that I only keep NS mutations, those that change the amino acid
P1_NS <- P1_filter[!(P1_filter$Amino_first == P1_filter$Amino_last),]
P2_NS <- P2_filter[!(P2_filter$Amino_first == P2_filter$Amino_last),]
P3_NS <- P3_filter[!(P3_filter$Amino_first == P3_filter$Amino_last),]
P4_NS <- P4_filter[!(P4_filter$Amino_first == P4_filter$Amino_last),]
P5_NS <- P5_filter[!(P5_filter$Amino_first == P5_filter$Amino_last),]
P6_NS <- P6_filter[!(P6_filter$Amino_first == P6_filter$Amino_last),]
 
P21_NS <- P21_filter[!(P21_filter$Amino_first == P21_filter$Amino_last),]
P22_NS <- P22_filter[!(P22_filter$Amino_first == P22_filter$Amino_last),]
P23_NS <- P23_filter[!(P23_filter$Amino_first == P23_filter$Amino_last),]
P24_NS <- P24_filter[!(P24_filter$Amino_first == P24_filter$Amino_last),]
P25_NS <- P25_filter[!(P25_filter$Amino_first == P25_filter$Amino_last),]
P26_NS <- P26_filter[!(P26_filter$Amino_first == P26_filter$Amino_last),]

B1_NS <- B1_filter[!(B1_filter$Amino_first == B1_filter$Amino_last),]
B2_NS <- B2_filter[!(B2_filter$Amino_first == B2_filter$Amino_last),]
B3_NS <- B3_filter[!(B3_filter$Amino_first == B3_filter$Amino_last),]
B4_NS <- B4_filter[!(B4_filter$Amino_first == B4_filter$Amino_last),]
B5_NS <- B5_filter[!(B5_filter$Amino_first == B5_filter$Amino_last),]
B6_NS <- B6_filter[!(B6_filter$Amino_first == B6_filter$Amino_last),]

B21_NS <- B21_filter[!(B21_filter$Amino_first == B21_filter$Amino_last),]
B22_NS <- B22_filter[!(B22_filter$Amino_first == B22_filter$Amino_last),]
B23_NS <- B23_filter[!(B23_filter$Amino_first == B23_filter$Amino_last),]
B24_NS <- B24_filter[!(B24_filter$Amino_first == B24_filter$Amino_last),]
B25_NS <- B25_filter[!(B25_filter$Amino_first == B25_filter$Amino_last),]
B26_NS <- B26_filter[!(B26_filter$Amino_first == B26_filter$Amino_last),]

#NOW SAVE THESE FILTERED DATASETS

write.csv(P1_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P1_NS.csv")
write.csv(P2_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P2_NS.csv")
write.csv(P3_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P3_NS.csv")
write.csv(P4_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P4_NS.csv")
write.csv(P5_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P5_NS.csv")
write.csv(P6_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P6_NS.csv")

write.csv(P21_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P21_NS.csv")
write.csv(P22_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P22_NS.csv")
write.csv(P23_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis//P23_NS.csv")
write.csv(P24_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P24_NS.csv")
write.csv(P25_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P25_NS.csv")
write.csv(P26_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P26_NS.csv")


write.csv(B1_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B1_NS.csv")
write.csv(B2_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B2_NS.csv")
write.csv(B3_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B3_NS.csv")
write.csv(B4_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B4_NS.csv")
write.csv(B5_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B5_NS.csv")
write.csv(B6_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B6_NS.csv")

write.csv(B21_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B21_NS.csv")
write.csv(B22_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B22_NS.csv")
write.csv(B23_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B23_NS.csv")
write.csv(B24_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B24_NS.csv")
write.csv(B25_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B25_NS.csv")
write.csv(B26_NS, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B26_NS.csv")

#######
#READ IN THE FILTERED DATA SETS

P1_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P1_NS.csv", stringsAsFactors = F)
P2_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P2_NS.csv", stringsAsFactors = F)
P3_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P3_NS.csv", stringsAsFactors = F)
P4_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P4_NS.csv", stringsAsFactors = F)
P5_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P5_NS.csv", stringsAsFactors = F)
P6_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P6_NS.csv", stringsAsFactors = F)

P21_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P21_NS.csv", stringsAsFactors = F)
P22_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P22_NS.csv", stringsAsFactors = F)
P23_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis//P23_NS.csv", stringsAsFactors = F)
P24_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P24_NS.csv", stringsAsFactors = F)
P25_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P25_NS.csv", stringsAsFactors = F)
P26_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P26_NS.csv", stringsAsFactors = F)

B1_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B1_NS.csv", stringsAsFactors = F)
B2_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B2_NS.csv", stringsAsFactors = F)
B3_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B3_NS.csv", stringsAsFactors = F)
B4_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B4_NS.csv", stringsAsFactors = F)
B5_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B5_NS.csv", stringsAsFactors = F)
B6_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B6_NS.csv", stringsAsFactors = F)

B21_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B21_NS.csv", stringsAsFactors = F)
B22_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B22_NS.csv", stringsAsFactors = F)
B23_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B23_NS.csv", stringsAsFactors = F)
B24_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B24_NS.csv", stringsAsFactors = F)
B25_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B25_NS.csv", stringsAsFactors = F)
B26_filtered <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B26_NS.csv", stringsAsFactors = F)


######
#read in the initially filtered timeseries tables
#####

##print a table of the initially filtered numbers of mutations
initialfiltered_mutations <- as.data.frame(matrix(ncol = 2, byrow = 2, 
                                             c("Population", "mutations",
                                               "P1", nrow(P1_filtered),
                                               "P2", nrow(P2_filtered),
                                               "P3", nrow(P3_filtered),
                                               "P4", nrow(P4_filtered),
                                               "P5", nrow(P5_filtered),
                                               "P6", nrow(P6_filtered),
                                               
                                               "P21", nrow(P21_filtered),
                                               "P22", nrow(P22_filtered),
                                               "P23", nrow(P23_filtered),
                                               "P24", nrow(P24_filtered),
                                               "P25", nrow(P25_filtered),
                                               "P26", nrow(P26_filtered),
                                               
                                               "B1", nrow(B1_filtered),
                                               "B2", nrow(B2_filtered),
                                               "B3", nrow(B3_filtered),
                                               "B4", nrow(B4_filtered),
                                               "B5", nrow(B5_filtered),
                                               "B6", nrow(B6_filtered),
                                               
                                               "B21", nrow(B21_filtered),
                                               "B22", nrow(B22_filtered),
                                               "B23", nrow(B23_filtered),
                                               "B24", nrow(B24_filtered),
                                               "B25", nrow(B25_filtered),
                                               "B26", nrow(B26_filtered))))
write.csv(initialfiltered_mutations, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/initial_filtered_numbers_.csv")




#then just save the information column and the maximum frequency column. This will allow me to make a big table of all populations and all mutations by position. Saving the maximum frequency will allow me to color code by max frequency at a later time point.

P1_max <- P1_filtered[,c(2,13)] #this saves the maximum frequency column and the column that contains all of the descriptive information
P1_max$sample <- "P1" #this adds a column that will have the population sample name sot hat I can make a table by sample name later.
P2_max <- P2_filtered[,c(2,13)]
P2_max$sample <- "P2"
P3_max <- P3_filtered[,c(2,13)]
P3_max$sample <- "P3"
P4_max <- P4_filtered[,c(2,13)]
P4_max$sample <- "P4"
P5_max <- P5_filtered[,c(2,13)]
P5_max$sample <- "P5"
P6_max <- P6_filtered[,c(2,13)]
P6_max$sample <- "P6"

P21_max <- P21_filtered[,c(2,13)]
P21_max$sample <- "P21"
P22_max <- P22_filtered[,c(2,13)]
P22_max$sample <- "P22"
P23_max <- P23_filtered[,c(2,13)]
P23_max$sample <- "P23"
P24_max <- P24_filtered[,c(2,13)]
P24_max$sample <- "P24"
P25_max <- P25_filtered[,c(2,13)]
P25_max$sample <- "P25"
P26_max <- P26_filtered[,c(2,13)]
P26_max$sample <- "P26"


B1_max <- B1_filtered[,c(2,13)]
B1_max$sample <- "B1"
B2_max <- B2_filtered[,c(2,13)]
B2_max$sample <- "B2"
B3_max <- B3_filtered[,c(2,13)]
B3_max$sample <- "B3"
B4_max <- B4_filtered[,c(2,13)]
B4_max$sample <- "B4"
B5_max <- B5_filtered[,c(2,13)]
B5_max$sample <- "B5"
B6_max <- B6_filtered[,c(2,13)]
B6_max$sample <- "B6"

B21_max <- B21_filtered[,c(2,13)]
B21_max$sample <- "B21"
B22_max <- B22_filtered[,c(2,13)]
B22_max$sample <- "B22"
B23_max <- B23_filtered[,c(2,13)]
B23_max$sample <- "B23"
B24_max <- B24_filtered[,c(2,13)]
B24_max$sample <- "B24"
B25_max <- B25_filtered[,c(2,13)]
B25_max$sample <- "B25"
B26_max <- B26_filtered[,c(2,13)]
B26_max$sample <- "B26"

#####
#group by environment to later be able to filter by common mutations
######

#first I need to bind all of the samples for each of the 4 environmental groupings.

arginine_planktonic <- rbind(P1_max,P2_max,P3_max,P4_max,P5_max,P6_max)
glucose_planktonic <- rbind(P21_max,P22_max,P23_max,P24_max,P25_max,P26_max)

arginine_biofilm <- rbind(B1_max,B2_max,B3_max,B4_max,B5_max,B6_max)
glucose_biofilm <- rbind(B21_max,B22_max,B23_max,B24_max,B25_max,B26_max)


all_populations <- rbind(arginine_planktonic, arginine_biofilm, glucose_biofilm, glucose_planktonic)
#then I will take each of the 4 groups and put them into tables where a column is the max frequency in a given population, and the rows are unique mutations at the position level. This will show me parallelism across all 6 samples in the same environment. Doing this in the 6 environments and then I will be doing this for all 24 of the populations to be able to see all of them. 



#first, the arginine planktonic populations
m_arginine_planktonic <- melt(arginine_planktonic, id=c("X","sample"), measure.vars = "Max") #melt the data so that it can be rearranged later
cast_arginine_planktonic <- as.data.frame(t(dcast(m_arginine_planktonic, sample~X, mean, value.var="value", fill = 0))) #cast the data at a data frame with rows being mutations and columns being the frequency of that mutation in a given population.
colnames(cast_arginine_planktonic) <- as.character(unlist(cast_arginine_planktonic[1,])) #make the column names the 
arginine_planktonic_final <- cast_arginine_planktonic[-1,]
#at this point I am NOT going to be extracting the information from the info in the rownames. It is irrelevant to the filtering process what the identity of the mutation is.


#The glucose planktonic populations
m_glucose_planktonic <- melt(glucose_planktonic, id=c("X","sample"), measure.vars = "Max")
cast_glucose_planktonic <- as.data.frame(t(dcast(m_glucose_planktonic, sample~X, mean, value.var="value", fill = 0)))
colnames(cast_glucose_planktonic) <- as.character(unlist(cast_glucose_planktonic[1,]))
glucose_planktonic_final <- cast_glucose_planktonic[-1,]


#The arginine biofilm populations
m_arginine_biofilm <- melt(arginine_biofilm, id=c("X","sample"), measure.vars = "Max")
cast_arginine_biofilm <- as.data.frame(t(dcast(m_arginine_biofilm, sample~X, mean, value.var="value", fill = 0)))
colnames(cast_arginine_biofilm) <- as.character(unlist(cast_arginine_biofilm[1,]))
arginine_biofilm_final <- cast_arginine_biofilm[-1,]

#The glucose biofilm populations
m_glucose_biofilm <- melt(glucose_biofilm, id=c("X","sample"), measure.vars = "Max")
cast_glucose_biofilm <- as.data.frame(t(dcast(m_glucose_biofilm, sample~X, mean, value.var="value", fill = 0)))
colnames(cast_glucose_biofilm) <- as.character(unlist(cast_glucose_biofilm[1,]))
glucose_biofilm_final <- cast_glucose_biofilm[-1,]

#pint these
write.csv(arginine_planktonic_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/arginine_plankotnic.csv")
write.csv(glucose_planktonic_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/glucose_plankotnic.csv")
write.csv(arginine_biofilm_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/arginine_biofilm.csv")
write.csv(glucose_biofilm_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/glucose_biofilm.csv")


#finally make a table with the number of mutations in each of the populations.
environment_mutations <- as.data.frame(matrix(ncol = 2, byrow = T, 
                                              c("Condition", "Mutation",
                                              "Arginine Planktonic", nrow(arginine_planktonic_final),
                                              "Glucose Planktonic", nrow(glucose_planktonic_final),
                                              "Arginine Biofilm", nrow(arginine_biofilm_final),
                                              "Glucose Biofilm", nrow(glucose_biofilm_final))))
write.csv(environment_mutations, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/environment_mutations.csv")


#now do all of the populations
m_all_populations <- melt(all_populations, id=c("X","sample"), measure.vars = "Max")
cast_all_populations <- as.data.frame(t(dcast(m_all_populations, sample~X, mean, value.var="value", fill = 0)))
colnames(cast_all_populations) <- as.character(unlist(cast_all_populations[1,]))
all_populations_final <- cast_all_populations[-1,]

colorder <- c("P1", "P2","P3","P4","P5","P6", "P21", "P22", "P23", "P24", "P25", "P26","B1", "B2",
              "B3","B4", "B5", "B6", "B21", "B22", "B23", "B24", "B25", "B26")
all_populations_final <- setcolorder(all_populations_final, colorder)

#set rows as type numeric
all_pops <- as.data.frame(apply(all_populations_final[,1:24],2,function(x) as.numeric(as.character(x))))
all_populations_final[,1:24] <- all_pops[,1:24]

#mak an additional column on the all which has a count of how many populations see that mutation.
all_populations_final$count <- rowSums(all_populations_final[1:24] != 0.0)

#######
#final filter, if the mutetion is seen in >90% of the populations it is likely background
#####
final_mutations <- subset(all_populations_final, all_populations_final$count < 12)
final_mutations$info <- rownames(final_mutations)
#I have previously filtered out a bunch of mutations in known regions of variability. going to import that list and filter now. 
filtered_out <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200422_filteredout.csv", stringsAsFactors = F)

final_mutations <- final_mutations[!(final_mutations$info %in% filtered_out$X),]

#print
write.csv(final_mutations,  "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/Paralellism_final_filtering.csv")

#I have filtered mutations by hand one last time. SO import the saved mutations here 
save <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_filter.csv", stringsAsFactors = F)
save_info <- colsplit(save$X, ":::", c("one","two","three","four","five"))
save$position <- save_info$one

#and filter the final mutations so that they only include the positions in the save files
p1_final <- P1_filtered[(P1_filtered$Position %in% save$position),]
p2_final <- P2_filtered[(P2_filtered$Position %in% save$position),]
p3_final <- P3_filtered[(P3_filtered$Position %in% save$position),]
p4_final <- P4_filtered[(P4_filtered$Position %in% save$position),]
p5_final <- P5_filtered[(P5_filtered$Position %in% save$position),]
p6_final <- P6_filtered[(P6_filtered$Position %in% save$position),]

p21_final <- P21_filtered[(P21_filtered$Position %in% save$position),]
p22_final <- P22_filtered[(P22_filtered$Position %in% save$position),]
p23_final <- P23_filtered[(P23_filtered$Position %in% save$position),]
p24_final <- P24_filtered[(P24_filtered$Position %in% save$position),]
p25_final <- P25_filtered[(P25_filtered$Position %in% save$position),]
p26_final <- P26_filtered[(P26_filtered$Position %in% save$position),]

b1_final <- B1_filtered[(B1_filtered$Position %in% save$position),]
b2_final <- B2_filtered[(B2_filtered$Position %in% save$position),]
b3_final <- B3_filtered[(B3_filtered$Position %in% save$position),]
b4_final <- B4_filtered[(B4_filtered$Position %in% save$position),]
b5_final <- B5_filtered[(B5_filtered$Position %in% save$position),]
b6_final <- B6_filtered[(B6_filtered$Position %in% save$position),]

b21_final <- B21_filtered[(B21_filtered$Position %in% save$position),]
b22_final <- B22_filtered[(B22_filtered$Position %in% save$position),]
b23_final <- B23_filtered[(B23_filtered$Position %in% save$position),]
b24_final <- B24_filtered[(B24_filtered$Position %in% save$position),]
b25_final <- B25_filtered[(B25_filtered$Position %in% save$position),]
b26_final <- B26_filtered[(B26_filtered$Position %in% save$position),]

######
#unload all of the information in the X column to make these data tables pretty.
######
p1_info2 <- colsplit(p1_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p1_final$Position <- p1_info2$Position
p1_final$Mutation <- p1_info2$Mutation
p1_final$Annotation <- p1_info2$Annotation
p1_final$Gene <- p1_info2$Gene
p1_final$Description <- p1_info2$Description

p2_info2 <- colsplit(p2_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p2_final$Position <- p2_info2$Position
p2_final$Mutation <- p2_info2$Mutation
p2_final$Annotation <- p2_info2$Annotation
p2_final$Gene <- p2_info2$Gene
p2_final$Description <- p2_info2$Description

p3_info2 <- colsplit(p3_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p3_final$Position <- p3_info2$Position
p3_final$Mutation <- p3_info2$Mutation
p3_final$Annotation <- p3_info2$Annotation
p3_final$Gene <- p3_info2$Gene
p3_final$Description <- p3_info2$Description

p4_info2 <- colsplit(p4_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p4_final$Position <- p4_info2$Position
p4_final$Mutation <- p4_info2$Mutation
p4_final$Annotation <- p4_info2$Annotation
p4_final$Gene <- p4_info2$Gene
p4_final$Description <- p4_info2$Description

p5_info2 <- colsplit(p5_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p5_final$Position <- p5_info2$Position
p5_final$Mutation <- p5_info2$Mutation
p5_final$Annotation <- p5_info2$Annotation
p5_final$Gene <- p5_info2$Gene
p5_final$Description <- p5_info2$Description


p6_info2 <- colsplit(p6_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p6_final$Position <- p6_info2$Position
p6_final$Mutation <- p6_info2$Mutation
p6_final$Annotation <- p6_info2$Annotation
p6_final$Gene <- p6_info2$Gene
p6_final$Description <- p6_info2$Description


p21_info2 <- colsplit(p21_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p21_final$Position <- p21_info2$Position
p21_final$Mutation <- p21_info2$Mutation
p21_final$Annotation <- p21_info2$Annotation
p21_final$Gene <- p21_info2$Gene
p21_final$Description <- p21_info2$Description

p22_info2 <- colsplit(p22_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p22_final$Position <- p22_info2$Position
p22_final$Mutation <- p22_info2$Mutation
p22_final$Annotation <- p22_info2$Annotation
p22_final$Gene <- p22_info2$Gene
p22_final$Description <- p22_info2$Description

p23_info2 <- colsplit(p23_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p23_final$Position <- p23_info2$Position
p23_final$Mutation <- p23_info2$Mutation
p23_final$Annotation <- p23_info2$Annotation
p23_final$Gene <- p23_info2$Gene
p23_final$Description <- p23_info2$Description

p24_info2 <- colsplit(p24_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p24_final$Position <- p24_info2$Position
p24_final$Mutation <- p24_info2$Mutation
p24_final$Annotation <- p24_info2$Annotation
p24_final$Gene <- p24_info2$Gene
p24_final$Description <- p24_info2$Description

p25_info2 <- colsplit(p25_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p25_final$Position <- p25_info2$Position
p25_final$Mutation <- p25_info2$Mutation
p25_final$Annotation <- p25_info2$Annotation
p25_final$Gene <- p25_info2$Gene
p25_final$Description <- p25_info2$Description

p26_info2 <- colsplit(p26_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Description"))
p26_final$Position <- p26_info2$Position
p26_final$Mutation <- p26_info2$Mutation
p26_final$Annotation <- p26_info2$Annotation
p26_final$Gene <- p26_info2$Gene
p26_final$Description <- p26_info2$Description


b1_info2 <- colsplit(b1_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b1_final$Position <- b1_info2$Position
b1_final$Mutation <- b1_info2$Mutation
b1_final$Annotation <- b1_info2$Annotation
b1_final$Gene <- b1_info2$Gene
b1_final$Describtion <- b1_info2$Describtion

b2_info2 <- colsplit(b2_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b2_final$Position <- b2_info2$Position
b2_final$Mutation <- b2_info2$Mutation
b2_final$Annotation <- b2_info2$Annotation
b2_final$Gene <- b2_info2$Gene
b2_final$Describtion <- b2_info2$Describtion

b3_info2 <- colsplit(b3_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b3_final$Position <- b3_info2$Position
b3_final$Mutation <- b3_info2$Mutation
b3_final$Annotation <- b3_info2$Annotation
b3_final$Gene <- b3_info2$Gene
b3_final$Describtion <- b3_info2$Describtion

b4_info2 <- colsplit(b4_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b4_final$Position <- b4_info2$Position
b4_final$Mutation <- b4_info2$Mutation
b4_final$Annotation <- b4_info2$Annotation
b4_final$Gene <- b4_info2$Gene
b4_final$Describtion <- b4_info2$Describtion

b5_info2 <- colsplit(b5_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b5_final$Position <- b5_info2$Position
b5_final$Mutation <- b5_info2$Mutation
b5_final$Annotation <- b5_info2$Annotation
b5_final$Gene <- b5_info2$Gene
b5_final$Describtion <- b5_info2$Describtion


b6_info2 <- colsplit(b6_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b6_final$Position <- b6_info2$Position
b6_final$Mutation <- b6_info2$Mutation
b6_final$Annotation <- b6_info2$Annotation
b6_final$Gene <- b6_info2$Gene
b6_final$Describtion <- b6_info2$Describtion


b21_info2 <- colsplit(b21_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b21_final$Position <- b21_info2$Position
b21_final$Mutation <- b21_info2$Mutation
b21_final$Annotation <- b21_info2$Annotation
b21_final$Gene <- b21_info2$Gene
b21_final$Describtion <- b21_info2$Describtion

b22_info2 <- colsplit(b22_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b22_final$Position <- b22_info2$Position
b22_final$Mutation <- b22_info2$Mutation
b22_final$Annotation <- b22_info2$Annotation
b22_final$Gene <- b22_info2$Gene
b22_final$Describtion <- b22_info2$Describtion

b23_info2 <- colsplit(b23_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b23_final$Position <- b23_info2$Position
b23_final$Mutation <- b23_info2$Mutation
b23_final$Annotation <- b23_info2$Annotation
b23_final$Gene <- b23_info2$Gene
b23_final$Describtion <- b23_info2$Describtion

b24_info2 <- colsplit(b24_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b24_final$Position <- b24_info2$Position
b24_final$Mutation <- b24_info2$Mutation
b24_final$Annotation <- b24_info2$Annotation
b24_final$Gene <- b24_info2$Gene
b24_final$Describtion <- b24_info2$Describtion

b25_info2 <- colsplit(b25_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b25_final$Position <- b25_info2$Position
b25_final$Mutation <- b25_info2$Mutation
b25_final$Annotation <- b25_info2$Annotation
b25_final$Gene <- b25_info2$Gene
b25_final$Describtion <- b25_info2$Describtion

b26_info2 <- colsplit(b26_final$X, ":::", names=c("Position", "Mutation", "Annotation", "Gene", "Describtion"))
b26_final$Position <- b26_info2$Position
b26_final$Mutation <- b26_info2$Mutation
b26_final$Annotation <- b26_info2$Annotation
b26_final$Gene <- b26_info2$Gene
b26_final$Describtion <- b26_info2$Describtion



#####
#PRINT OUT THE FINAL MUTATIONS
write.csv(p1_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P1_final.csv")
write.csv(p2_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P2_final.csv")
write.csv(p3_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P3_final.csv")
write.csv(p4_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P4_final.csv")
write.csv(p5_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P5_final.csv")
write.csv(p6_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P6_final.csv")

write.csv(p21_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P21_final.csv")
write.csv(p22_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P22_final.csv")
write.csv(p23_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P23_final.csv")
write.csv(p24_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P24_final.csv")
write.csv(p25_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P25_final.csv")
write.csv(p26_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/P26_final.csv")

write.csv(b1_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B1_final.csv")
write.csv(b2_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B2_final.csv")
write.csv(b3_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B3_final.csv")
write.csv(b4_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B4_final.csv")
write.csv(b5_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B5_final.csv")
write.csv(b6_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B6_final.csv")

write.csv(b21_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B21_final.csv")
write.csv(b22_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B22_final.csv")
write.csv(b23_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B23_final.csv")
write.csv(b24_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B24_final.csv")
write.csv(b25_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B25_final.csv")
write.csv(b26_final, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/B26_final.csv")


#####
#combine into one large table
#######
#give them all a population designator 
p1_final$Population <- "P1"
p1_final <- p1_final[,-16]
p2_final$Population <- "P2"
p3_final$Population <- "P3"
p4_final$Population <- "P4"
p5_final$Population <- "P5"
p6_final$Population <- "P6"

p21_final$Population <- "P21"
p22_final$Population <- "P22"
p23_final$Population <- "P23"
p24_final$Population <- "P24"
p25_final$Population <- "P25"
p26_final$Population <- "P26"

b1_final$Population <- "B1"
b2_final$Population <- "B2"
b3_final$Population <- "B3"
b4_final$Population <- "B4"
b5_final$Population <- "B5"
b6_final$Population <- "B6"

b21_final$Population <- "B21"
b22_final$Population <- "B22"
b23_final$Population <- "B23"
b24_final$Population <- "B24"
b25_final$Population <- "B25"
b26_final$Population <- "B26"

#change all of the column names to be identical
columnnames <- c("X.1", "X","X0","X6","X12","X25","Position","Mutation","Annotation","Gene","Description","Count","Max","Min","CHange","Amino_first","Amino_last","Population")
colnames(p1_final) <- columnnames
colnames(p2_final) <- columnnames
colnames(p3_final) <- columnnames
colnames(p4_final) <- columnnames
colnames(p5_final) <- columnnames
colnames(p6_final) <- columnnames

colnames(p21_final) <- columnnames
colnames(p22_final) <- columnnames
colnames(p23_final) <- columnnames
colnames(p24_final) <- columnnames
colnames(p25_final) <- columnnames
colnames(p26_final) <- columnnames

colnames(b1_final) <- columnnames
colnames(b2_final) <- columnnames
colnames(b3_final) <- columnnames
colnames(b4_final) <- columnnames
colnames(b5_final) <- columnnames
colnames(b6_final) <- columnnames

colnames(b21_final) <- columnnames
colnames(b22_final) <- columnnames
colnames(b23_final) <- columnnames
colnames(b24_final) <- columnnames
colnames(b25_final) <- columnnames
colnames(b26_final) <- columnnames

AP <- rbind(p1_final, p2_final, p3_final, p4_final, p5_final, p6_final)
GP <- rbind(p21_final, p22_final, p23_final, p24_final, p25_final, p26_final)

AB <- rbind(b1_final, b2_final, b3_final, b4_final, b5_final, b6_final)
AB <- AB[,-18]
colnames(AB) <- columnnames
GB <- rbind(b21_final, b22_final, b23_final, b24_final, b25_final, b26_final)
GB <- GB[,-18]
colnames(GB) <- columnnames
All <- rbind(AP, GP, AB, GB)

#print the table of all mutations to be able to look at them
write.csv(All, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/all_final_mutations.csv")


######
#import the final mutation master table so that I can restart from here at any point.
Final_mutations <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/all_final_mutations.csv", stringsAsFactors = F)

######
#make a table of the number of mutations seen in each population

mutations <- table(Final_mutations$Population)

write.csv(mutations, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_mutations_numbers.csv")

######
#going to make a table of the gene level parallelism for the populations. Going to just save the population, maximum frequency, and gene name for all. 
GeneParallelism <- Final_mutations[,c(8,11,16)]

Cast_geneparallel <- as.data.frame(t(dcast(GeneParallelism, Population~Gene, mean, value.var = "Max", fill = 0)))
colnames(Cast_geneparallel) <- as.character(unlist(Cast_geneparallel[1,]))
GeneParallelism <- Cast_geneparallel[-1,]

#print this table
write.csv(GeneParallelism, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_geneparallelism.csv")



#now do for day 25 only
day25_distribution <- Final_mutations[,c(4,8,16)]
cast_day25 <- as.data.frame(t(dcast(day25_distribution, Population~Gene, mean, value.var = "X25", fill = 0)))
colnames(cast_day25) <- as.character(unlist(cast_day25[1,]))
day25_distribution <- cast_day25[-1,]
write.csv(day25_distribution, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_day25_geneparallelism.csv")
######
#separate into individual populations

P1 <- Final_mutations[(Final_mutations$Population == "P1"),]
P2 <- Final_mutations[(Final_mutations$Population == "P2"),]
P3 <- Final_mutations[(Final_mutations$Population == "P3"),]
P4 <- Final_mutations[(Final_mutations$Population == "P4"),]
P5 <- Final_mutations[(Final_mutations$Population == "P5"),]
P6 <- Final_mutations[(Final_mutations$Population == "P6"),]

P21 <- Final_mutations[(Final_mutations$Population == "P21"),]
P22 <- Final_mutations[(Final_mutations$Population == "P22"),]
P23 <- Final_mutations[(Final_mutations$Population == "P23"),]
P24 <- Final_mutations[(Final_mutations$Population == "P24"),]
P25 <- Final_mutations[(Final_mutations$Population == "P25"),]
P26 <- Final_mutations[(Final_mutations$Population == "P26"),]

B1 <- Final_mutations[(Final_mutations$Population == "B1"),]
B2 <- Final_mutations[(Final_mutations$Population == "B2"),]
B3 <- Final_mutations[(Final_mutations$Population == "B3"),]
B4 <- Final_mutations[(Final_mutations$Population == "B4"),]
B5 <- Final_mutations[(Final_mutations$Population == "B5"),]
B6 <- Final_mutations[(Final_mutations$Population == "B6"),]

B21 <- Final_mutations[(Final_mutations$Population == "B21"),]
B22 <- Final_mutations[(Final_mutations$Population == "B22"),]
B23 <- Final_mutations[(Final_mutations$Population == "B23"),]
B24 <- Final_mutations[(Final_mutations$Population == "B24"),]
B25 <- Final_mutations[(Final_mutations$Population == "B25"),]
B26 <- Final_mutations[(Final_mutations$Population == "B26"),]
######
#Calculate the alpha diversity in all populations at all time points sampled. Doing Shannon
library(vegan)
######
#going to create a table with the amount of alpha diversity seen in all of the populations at all of the time points sampled. 

#starting with shannon's index
shannon_diversity <- as.data.frame(matrix(ncol = 5, byrow=T, c("Population", "Day 0", "Day 6", "Day 12", "Day 25", "P1", diversity(P1$X0, index="shannon"), 
                  diversity(P1$X6, index="shannon"), 
                  diversity(P1$X12, index="shannon"), 
                  diversity(P1$X25, index="shannon"),
                  "P2", diversity(P2$X0, index="shannon"), 
                  diversity(P2$X6, index="shannon"), 
                  diversity(P2$X12, index="shannon"), 
                  diversity(P2$X25, index="shannon"),
                  "P3", diversity(P3$X0, index="shannon"),
                  diversity(P3$X6, index="shannon"), 
                  diversity(P3$X12, index="shannon"), 
                  diversity(P3$X25, index="shannon"),
                  "P4", diversity(P4$X0, index="shannon"), 
                  diversity(P4$X6, index="shannon"), 
                  diversity(P4$X12, index="shannon"), 
                  diversity(P4$X25, index="shannon"), 
                  "P5", diversity(P5$X0, index="shannon"), 
                  diversity(P5$X6, index="shannon"), 
                  diversity(P5$X12, index="shannon"), 
                  diversity(P5$X25, index="shannon"),
                  "P6", diversity(P6$X0, index="shannon"), 
                  diversity(P6$X6, index="shannon"), 
                  diversity(P6$X12, index="shannon"), 
                  diversity(P6$X25, index="shannon"),
                  "P21", diversity(P21$X0, index="shannon"), 
                  diversity(P21$X6, index="shannon"), 
                  diversity(P21$X12, index="shannon"), 
                  diversity(P21$X25, index="shannon"),
                  "P22", diversity(P22$X0, index="shannon"), 
                  diversity(P22$X6, index="shannon"), 
                  diversity(P22$X12, index="shannon"),
                  diversity(P22$X25, index="shannon"),
                  "P23", diversity(P23$X0, index="shannon"), 
                  diversity(P23$X6, index="shannon"), 
                  diversity(P23$X12, index="shannon"), 
                  diversity(P23$X25, index="shannon"),
                  "P24", diversity(P24$X0, index="shannon"), 
                  diversity(P24$X6, index="shannon"), 
                  diversity(P24$X12, index="shannon"), 
                  diversity(P24$X25, index="shannon"),
                  "P25", diversity(P25$X0, index="shannon"), 
                  diversity(P25$X6, index="shannon"), 
                  diversity(P25$X12, index="shannon"), 
                  diversity(P25$X25, index="shannon"),
                  "P26", diversity(P26$X0, index="shannon"), 
                  diversity(P26$X6, index="shannon"), 
                  diversity(P26$X12, index="shannon"), 
                  diversity(P26$X25, index="shannon"),
                  "B1", diversity(B1$X0, index="shannon"), 
                  diversity(B1$X6, index="shannon"), 
                  diversity(B1$X12, index="shannon"), 
                  diversity(B1$X25, index="shannon"),
                  "B2", diversity(B2$X0, index="shannon"), 
                  diversity(B2$X6, index="shannon"), 
                  diversity(B2$X12, index="shannon"), 
                  diversity(B2$X25, index="shannon"),
                  "B3", diversity(B3$X0, index="shannon"), 
                  diversity(B3$X6, index="shannon"), 
                  diversity(B3$X12, index="shannon"), 
                  diversity(B3$X25, index="shannon"),
                  "B4", diversity(B4$X0, index="shannon"), 
                  diversity(B4$X6, index="shannon"), 
                  diversity(B4$X12, index="shannon"), 
                  diversity(B4$X25, index="shannon"),
                  "B5", diversity(B5$X0, index="shannon"), 
                  diversity(B5$X6, index="shannon"), 
                  diversity(B5$X12, index="shannon"), 
                  diversity(B5$X25, index="shannon"),
                  "B6", diversity(B6$X0, index="shannon"), 
                  diversity(B6$X6, index="shannon"),
                  diversity(B6$X12, index="shannon"), 
                  diversity(B6$X25, index="shannon"),
                  "B21", diversity(B21$X0, index="shannon"), 
                  diversity(B21$X6, index="shannon"), 
                  diversity(B21$X12, index="shannon"), 
                  diversity(B21$X25, index="shannon"),
                  "B22", diversity(B22$X0, index="shannon"), 
                  diversity(B22$X6, index="shannon"), 
                  diversity(B22$X12, index="shannon"),
                  diversity(B22$X25, index="shannon"),
                  "B23", diversity(B23$X0, index="shannon"), 
                  diversity(B23$X6, index="shannon"), 
                  diversity(B23$X12, index="shannon"), 
                  diversity(B23$X25, index="shannon"),
                  "B24", diversity(B24$X0, index="shannon"), 
                  diversity(B24$X6, index="shannon"), 
                  diversity(B24$X12, index="shannon"), 
                  diversity(B24$X25, index="shannon"),
                  "B25", diversity(B25$X0, index="shannon"), 
                  diversity(B25$X6, index="shannon"), 
                  diversity(B25$X12, index="shannon"), 
                  diversity(B25$X25, index="shannon"),
                  "B26", diversity(B26$X0, index="shannon"), 
                  diversity(B26$X6, index="shannon"), 
                  diversity(B26$X12, index="shannon"), 
                  diversity(B26$X25, index="shannon"))))

##and print it off 
write.csv(shannon_diversity, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/Shannon_diversity.csv")


#####
#and write the individual population files so that I can use lolipop to make muller plots.
write.csv(P1, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P1.csv")
write.csv(P2, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P2.csv")
write.csv(P3, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P3.csv")
write.csv(P4, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P4.csv")
write.csv(P5, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P5.csv")
write.csv(P6, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P6.csv")

write.csv(P21, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P21.csv")
write.csv(P22, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P22.csv")
write.csv(P23, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P23.csv")
write.csv(P24, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P24.csv")
write.csv(P25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P25.csv")
write.csv(P26, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P26.csv")

write.csv(B1, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B1.csv")
write.csv(B2, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B2.csv")
write.csv(B3, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B3.csv")
write.csv(B4, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B4.csv")
write.csv(B5, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B5.csv")
write.csv(B6, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B6.csv")

write.csv(B21, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B21.csv")
write.csv(B22, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B22.csv")
write.csv(B23, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B23.csv")
write.csv(B24, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B24.csv")
write.csv(B25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B25.csv")
write.csv(B26, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B26.csv")
#####

#now look at the types of mutations that evolve in the different populations


#####
#first make the function that will determine the different types of mutations
Mutations_analysis <- function(Mutation_Data, AminoAcid="description", Bases = "Mutation") {
  #Mutation_Data is the input CSV file that is a breseq output with rows containing different mutations and columns containing the various characteristics for those mutations.
  #AminoAcid is the column name that holds the breseq information on amino acid changes. this column will look like "*342C(TGA>TGC)" or say coding, pseudogene, etc. The default that will be looked for is "Description". This is case sensitive!
  #Bases is the column name that holds the breseq information for the nucleotide mutations. These are things like A>C, or T>G in the breseq output. This is case sensitive!!!
  
  ##############
  #first i am going to deal with the nucleotide level - the mutation column. This uses the Mutation_Data that you put in and grabs the information in the column that you specified under Bases. It looks for the > symbol, because that is how breseq separates the two bases, and splits the data. It then creates a new data set called Nucleotides that has 2 columns containing the original base (from the ancestor) and the mutated base.
  Nucleotides <- colsplit(Mutation_Data[,Bases], ">", names = c("original", "mutant"))
  #View(Nucleotides)
  #I want to calculate the total number of mutations present in the sample.
  allmutations <- nrow(Nucleotides)
  #I want to determine the number of mutations that are not just substituting one base from another. These are indels, because this is how breseq represents them.
  indel <- sum(grepl("[^ACGT]", Nucleotides$original))
  
  
  #I need to find all of the different combinations for base replacements. To do this I am going to find the index for each base, and then I will look at that index in the second column and see what the new base is.
  C <- grep("C", Nucleotides$original) #find placeswhere the original was C
  CT <- sum(ifelse(Nucleotides$mutant[C]=="T",1,0)) # find when there was a C to T transition
  CA <- sum(ifelse(Nucleotides$mutant[C]=="A",1,0)) # find when there was a C to A transition
  CG <- sum(ifelse(Nucleotides$mutant[C]=="G",1,0)) # find when there was a C to G transition
  
  Ts <- grep("T", Nucleotides$original) #find when the original was a T
  TC <- sum(ifelse(Nucleotides$mutant[Ts]=="C",1,0)) #find when there were T to C transitions
  TG <- sum(ifelse(Nucleotides$mutant[Ts]=="G",1,0)) #find when there were T to G transitions
  TA <- sum(ifelse(Nucleotides$mutant[Ts]=="A",1,0)) #find when there were T to A transitions
  
  G <- grep("G", Nucleotides$original) #find placeswhere the original was G
  GA <- sum(ifelse(Nucleotides$mutant[G]=="A",1,0)) # find when there was a G to A transition
  GT <- sum(ifelse(Nucleotides$mutant[G]=="T",1,0)) # find when there was a G to T transition
  GC <- sum(ifelse(Nucleotides$mutant[G]=="C",1,0)) # find when there was a G to C transition
  
  A <- grep("A", Nucleotides$original) #find placeswhere the original was A
  AG <- sum(ifelse(Nucleotides$mutant[A]=="G",1,0)) # find when there was a A to G transition
  AC <- sum(ifelse(Nucleotides$mutant[A]=="C",1,0)) # find when there was a A to C transition
  AT <- sum(ifelse(Nucleotides$mutant[A]=="T",1,0)) # find when there was a A to T transition
  
  # Now that I have the numbers of all of the possible base changes, I can look for the
  transitions <- sum(c(CT,TC,GA,AG)) #there are 4 options for transitions. C>T, T>C, G>A, A>G. this adds up all of those changes
  transversions <- AT+AC+GC+GT+CA+CG+TG+TA # need to do to check that the sums of the transition categories actually do add up to the number of transitions that there should be (assuming transitions and indel numbers are correct) when I turn this into a function I need to stop if transversions != trans -- should be fine but just an extra error checking step.
  
  ###############
  
  ### now at the Amino acid level
  #have to get the amino acid column that the user specifies out of the input data
  Protein <- colsplit(Mutation_Data[,AminoAcid], "\\(", names = c("AA", "DNA"))
  
  #there are a few options that you can get for this one and I can't just split the column. I need to look for all of them in the strings FOr this I need to use regular expressions. I will have to use gregexpr which returns a list of positions and then I will have to find the length of that to determine the number of them. The regular expressios that I will use for each option are as follows.
  
  # reg expressions to use "coding", "intergenic", "pseudogene",
  #"[A-Z][0-9]+[A-Z]" #this looks for a base, followed by a number that can be any size 1 or more, followed by a base.
  #"\\*[0-9]+[A-Z]" #this looks for an asterisk, followed by at least 1 number, followed by a base
  #"[A-Z][0-9]*\\*" #this looks for a base, followed by at least 1 number, followed by an asterisk
  
  coding = sum(grepl("coding",Protein$AA)) #Breseq's coding region
  intergenic = sum(grepl("intergenic",Protein$AA)) #intergenic breseq designation
  pseudogene = sum(grepl("pseudogene",Protein$AA)) #pseudogene breseq designation
  prematurestop = sum(lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]*\\*", Protein$AA)))) #these are when you have a coding amino acid that gets mutated into a stop codon
  elongating = sum(lengths(regmatches(Protein$AA,gregexpr("\\*[0-9]*[A-Z]", Protein$AA)))) #these are stop codons that get mutated into coding amino acids that elongate your protein.
  aamutation = sum(lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]+[A-Z]", Protein$AA)))) # these are all of the mutations that dont fit other categories. so these mutations change one amino acid to another with no other breseq designation.
  
  #I now need to determine if the amino acid mutations category are synonymous or nonsynonymous. The above just determines the number of leter to leter strings exist. Now I need to actually look at that subset of mutations and determine if
  aas <- lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]+[A-Z]", Protein$AA))) #this returns a list of logicals that tell you if there is an aamutation at that spot (1) or not (0)
  #aas
  aminos <- as.matrix(Protein$AA[aas==1]) #this find the amino acid changes at the previously found indexes, so these are the actual identities of the aamutation mutations.
  aminos2 <- colsplit(aminos, "[0-9]+", names = c("first", "last")) # splitting the breseq annotation. I am taking out the number, because the position is irrelevent, and separating the two letters into two columns containing the first, original aa, and the last, or mutated, aa.
  
  synonymous <- sum(ifelse(as.character(aminos2$first) == as.character(aminos2$last),1,0)) #if the letters are the same before and after the mutation it is synonymous
  nonsynonymous <- sum(ifelse(as.character(aminos2$first) == as.character(aminos2$last),0,1)) #if the letters are different then it is nonsynonymous
  dnds <- (nonsynonymous/synonymous)/2.60 
  #2.60 is the pseudomonas genomic dN/dS ratio. So this reported value will be gorrected for the basal genomic distribution.
  
  # I am now making a table of all of the mutational types that I can print out later. The other thing that I would like to do is give it a name specific to the data set that you put in, but I don't know how to do that. For the moment it will just always return the same named file each time, so you have to change the name before you use the code again. 
  table<- matrix(c("Mutations: ", allmutations,
                   "Nucleotide level mutations", "",
                   "Indels: ", indel,
                   "Transitions: ",transitions,
                   "C>T: ", CT,
                   "T>C: ", TC,
                   "A>G: ", AG,
                   "G>A: ", GA,
                   "Transversions: ", transversions,
                   "A>T: ", AT,
                   "A>C: ", AC,
                   "G>C: ", GC,
                   "G>T: ", GT,
                   "C>A: ", CA,
                   "C>G: ", CG,
                   "T>G: ", TG,
                   "T>A: ", TA,
                   "Amino acid level mutations", "",
                   "Coding: ", coding,
                   "Intergenic: ", intergenic,
                   "Pseudogene: ", pseudogene,
                   "Premature stop: ", prematurestop,
                   "Elongating: ", elongating,
                   "Synonymous: ", synonymous,
                   "Non Synonymous: ", nonsynonymous,
                   "dN/dS: ", dnds), ncol = 2, byrow=T)
  
  #write out the file. Adding col.names won't actually give the columns names but it will keep things in two different columns instead of compiling it all.
  write.csv(table, file = "Mutations_table.csv", col.names = T)
  
}

library("reshape2")
#set the working directory so that the output goes to the right place
setwd("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown")

#did this for all populations, but because the function writes out the same file each name it will just overwrite the same file if I dont change the names. This is just to make sure that I don't mix any samples up. 
Mutations_analysis(P1_final, "Annotation", "Mutation")

#now import all of the tables to be able to combine them. 
b1_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B1.csv", stringsAsFactors = F)
b2_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B2.csv", stringsAsFactors = F)
b3_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B3.csv", stringsAsFactors = F)
b4_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B4.csv", stringsAsFactors = F)
b5_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B5.csv", stringsAsFactors = F)
b6_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B6.csv", stringsAsFactors = F)

b21_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B21.csv", stringsAsFactors = F)
b22_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B22.csv", stringsAsFactors = F)
b23_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B23.csv", stringsAsFactors = F)
b24_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B24.csv", stringsAsFactors = F)
b25_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B25.csv", stringsAsFactors = F)
b26_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/B26.csv", stringsAsFactors = F)



p1_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P1.csv", stringsAsFactors = F)
p2_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P2.csv", stringsAsFactors = F)
p3_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P3.csv", stringsAsFactors = F)
p4_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P4.csv", stringsAsFactors = F)
p5_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P5.csv", stringsAsFactors = F)
p6_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P6.csv", stringsAsFactors = F)

p21_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P21.csv", stringsAsFactors = F)
p22_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P22.csv", stringsAsFactors = F)
p23_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P23.csv", stringsAsFactors = F)
p24_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P24.csv", stringsAsFactors = F)
p25_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P25.csv", stringsAsFactors = F)
p26_breakdown <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/P26.csv", stringsAsFactors = F)


#And combine all rows to make one large table 
mutation_breakdown <- cbind(p1_breakdown[,2:3], p2_breakdown[,3], p3_breakdown[,3], p4_breakdown[,3], p5_breakdown[,3], p6_breakdown[,3],p21_breakdown[,3], p22_breakdown[,3], p23_breakdown[,3], p24_breakdown[,3], p25_breakdown[,3], p26_breakdown[,3],b1_breakdown[,3], b2_breakdown[,3], b3_breakdown[,3], b4_breakdown[,3], b5_breakdown[,3], b6_breakdown[,3],b21_breakdown[,3], b22_breakdown[,3], b23_breakdown[,3], b24_breakdown[,3], b25_breakdown[,3], b26_breakdown[,3])

colnames(mutation_breakdown) <- c("Mutation type","P1", "P2", "P3","P4","P5","P6","P21","P22","P23","P24","P25","P26","B1", "B2", "B3","B4","B5","B6","B21","B22","B23","B24","B25","B26")


#and print the file

write.csv(mutation_breakdown,"/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_breakdown/all_breakdown.csv")

######


#going to be doing the bray curtis dissimilarity for the day 25 population samples.


######
#actually, I need to do a few things at the three different time points. I need to calculate bray curtis, plot mutations, and do mutational number analyses for all three of the time points. the code will be very similar, so it will all be here.
######


#####
#read in the files
P1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p1_filtered_final.csv", stringsAsFactors = F)
P2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p2_filtered_final.csv", stringsAsFactors = F)
P3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p3_filtered_final.csv", stringsAsFactors = F)
P4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p4_filtered_final.csv", stringsAsFactors = F)
P5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p5_filtered_final.csv", stringsAsFactors = F)
P6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p6_filtered_final.csv", stringsAsFactors = F)

P21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p21_filtered_final.csv", stringsAsFactors = F)
P22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p22_filtered_final.csv", stringsAsFactors = F)
P23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p23_filtered_final.csv", stringsAsFactors = F)
P24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p24_filtered_final.csv", stringsAsFactors = F)
P25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p25_filtered_final.csv", stringsAsFactors = F)
P26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p26_filtered_final.csv", stringsAsFactors = F)

B1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b1_filtered_final.csv", stringsAsFactors = F)
B2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b2_filtered_final.csv", stringsAsFactors = F)
B3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b3_filtered_final.csv", stringsAsFactors = F)
B4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b4_filtered_final.csv", stringsAsFactors = F)
B5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b5_filtered_final.csv", stringsAsFactors = F)
B6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b6_filtered_final.csv", stringsAsFactors = F)

B21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b21_filtered_final.csv", stringsAsFactors = F)
B22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b22_filtered_final.csv", stringsAsFactors = F)
B23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b23_filtered_final.csv", stringsAsFactors = F)
B24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b24_filtered_final.csv", stringsAsFactors = F)
B25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b25_filtered_final.csv", stringsAsFactors = F)
B26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b26_filtered_final.csv", stringsAsFactors = F)


######
#I need to make day 6 tables for each population, where the columns are Gene, frequency, and sample ID. rows will be genes.

#####
columnnames_6 <- c("Day_6", "Gene", "Sample")

P1_6 <- P1[,c(5,11)]
P1_6$Sample <- "P1"
P2_6 <- P2[,c(5,11)]
P2_6$Sample <- "P2"
P3_6 <- P3[,c(5,11)]
P3_6$Sample <- "P3"
P4_6 <- P4[,c(5,11)]
P4_6$Sample <- "P4"
P5_6 <- P5[,c(5,11)]
P5_6$Sample <- "P5"
P6_6 <- P6[,c(5,11)]
P6_6$Sample <- "P6"

P21_6 <- P21[,c(5,11)]
P21_6$Sample <- "P21"
P22_6 <- P22[,c(5,11)]
P22_6$Sample <- "P22"
P23_6 <- P23[,c(5,11)]
P23_6$Sample <- "P23"
P24_6 <- P24[,c(5,11)]
P24_6$Sample <- "P24"
P25_6 <- P25[,c(5,11)]
P25_6$Sample <- "P25"
P26_6 <- P26[,c(5,11)]
P26_6$Sample <- "P26"

B1_6 <- B1[,c(5,11)]
B1_6$Sample <- "B1"
B2_6 <- B2[,c(5,11)]
B2_6$Sample <- "B2"
B3_6 <- B3[,c(5,11)]
B3_6$Sample <- "B3"
B4_6 <- B4[,c(5,11)]
B4_6$Sample <- "B4"
B5_6 <- B5[,c(5,11)]
B5_6$Sample <- "B5"
B6_6 <- B6[,c(5,11)]
B6_6$Sample <- "B6"

B21_6 <- B21[,c(5,11)]
B21_6$Sample <- "B21"
B22_6 <- B22[,c(5,11)]
B22_6$Sample <- "B22"
B23_6 <- B23[,c(5,11)]
B23_6$Sample <- "B23"
B24_6 <- B24[,c(5,11)]
B24_6$Sample <- "B24"
B25_6 <- B25[,c(5,11)]
B25_6$Sample <- "B25"
B26_6 <- B26[,c(5,11)]
B26_6$Sample <- "B26"

colnames(P1_6) <- columnnames_6
colnames(P2_6) <- columnnames_6
colnames(P3_6) <- columnnames_6
colnames(P4_6) <- columnnames_6
colnames(P5_6) <- columnnames_6
colnames(P6_6) <- columnnames_6

colnames(P21_6) <- columnnames_6
colnames(P22_6) <- columnnames_6
colnames(P23_6) <- columnnames_6
colnames(P24_6) <- columnnames_6
colnames(P25_6) <- columnnames_6
colnames(P26_6) <- columnnames_6

colnames(B1_6) <- columnnames_6
colnames(B2_6) <- columnnames_6
colnames(B3_6) <- columnnames_6
colnames(B4_6) <- columnnames_6
colnames(B5_6) <- columnnames_6
colnames(B6_6) <- columnnames_6

colnames(B21_6) <- columnnames_6
colnames(B22_6) <- columnnames_6
colnames(B23_6) <- columnnames_6
colnames(B24_6) <- columnnames_6
colnames(B25_6) <- columnnames_6
colnames(B26_6) <- columnnames_6

#now group into the different environments.
arg_plank_6 <- rbind(P1_6, P2_6, P3_6, P4_6, P5_6, P6_6)
glu_plank_6 <- rbind(P21_6, P22_6, P23_6, P24_6, P25_6, P26_6)
arg_biof_6 <- rbind(B1_6, B2_6, B3_6, B4_6, B5_6, B6_6)
glu_biof_6 <- rbind(B21_6, B22_6, B23_6, B24_6, B25_6, B26_6)

glucose_6 <- rbind(glu_plank_6, glu_biof_6)
arginine_6 <- rbind(arg_biof_6, arg_plank_6)
biofilm_6 <- rbind(arg_biof_6, glu_biof_6)
planktonic_6 <- rbind(glu_plank_6, arg_plank_6)

all_6 <- rbind(arg_biof_6, arg_plank_6, glu_biof_6, glu_plank_6)

#####remove rows that are zero

arg_plank_6 <- arg_plank_6[!(arg_plank_6$Day_6 == 0),]
arg_biof_6 <- arg_biof_6[!(arg_biof_6$Day_6 == 0),]

glu_plank_6 <- glu_plank_6[!(glu_plank_6$Day_6 == 0),]
glu_biof_6 <- glu_biof_6[!(glu_biof_6$Day_6 == 0),]

glucose_6 <- glucose_6[!(glucose_6$Day_6 == 0),]
arginine_6 <- arginine_6[!(arginine_6$Day_6 == 0),]

biofilm_6 <- biofilm_6[!(biofilm_6$Day_6 == 0),]
planktonic_6 <- planktonic_6[!(planktonic_6$Day_6 == 0),]


all_6 <- all_6[!(all_6$Day_6 == 0),]

#write the table of all
write.csv(all_6, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/Day6_bray_curtis_formatted.csv")

#######
#I need to make day 12 tables for each population, where the columns are Gene, frequency, and sample ID. rows will be genes.

#####
columnnames_12 <- c("Day_12", "Gene", "Sample")

P1_12 <- P1[,c(6,11)]
P1_12$Sample <- "P1"
P2_12 <- P2[,c(6,11)]
P2_12$Sample <- "P2"
P3_12 <- P3[,c(6,11)]
P3_12$Sample <- "P3"
P4_12 <- P4[,c(6,11)]
P4_12$Sample <- "P4"
P5_12 <- P5[,c(6,11)]
P5_12$Sample <- "P5"
P6_12 <- P6[,c(6,11)]
P6_12$Sample <- "P6"

P21_12 <- P21[,c(6,11)]
P21_12$Sample <- "P21"
P22_12 <- P22[,c(6,11)]
P22_12$Sample <- "P22"
P23_12 <- P23[,c(6,11)]
P23_12$Sample <- "P23"
P24_12 <- P24[,c(6,11)]
P24_12$Sample <- "P24"
P25_12 <- P25[,c(6,11)]
P25_12$Sample <- "P25"
P26_12 <- P26[,c(6,11)]
P26_12$Sample <- "P26"

B1_12 <- B1[,c(6,11)]
B1_12$Sample <- "B1"
B2_12 <- B2[,c(6,11)]
B2_12$Sample <- "B2"
B3_12 <- B3[,c(6,11)]
B3_12$Sample <- "B3"
B4_12 <- B4[,c(6,11)]
B4_12$Sample <- "B4"
B5_12 <- B5[,c(6,11)]
B5_12$Sample <- "B5"
B6_12 <- B6[,c(6,11)]
B6_12$Sample <- "B6"

B21_12 <- B21[,c(6,11)]
B21_12$Sample <- "B21"
B22_12 <- B22[,c(6,11)]
B22_12$Sample <- "B22"
B23_12 <- B23[,c(6,11)]
B23_12$Sample <- "B23"
B24_12 <- B24[,c(6,11)]
B24_12$Sample <- "B24"
B25_12 <- B25[,c(6,11)]
B25_12$Sample <- "B25"
B26_12 <- B26[,c(6,11)]
B26_12$Sample <- "B26"

colnames(P1_12) <- columnnames_12
colnames(P2_12) <- columnnames_12
colnames(P3_12) <- columnnames_12
colnames(P4_12) <- columnnames_12
colnames(P5_12) <- columnnames_12
colnames(P6_12) <- columnnames_12

colnames(P21_12) <- columnnames_12
colnames(P22_12) <- columnnames_12
colnames(P23_12) <- columnnames_12
colnames(P24_12) <- columnnames_12
colnames(P25_12) <- columnnames_12
colnames(P26_12) <- columnnames_12

colnames(B1_12) <- columnnames_12
colnames(B2_12) <- columnnames_12
colnames(B3_12) <- columnnames_12
colnames(B4_12) <- columnnames_12
colnames(B5_12) <- columnnames_12
colnames(B6_12) <- columnnames_12

colnames(B21_12) <- columnnames_12
colnames(B22_12) <- columnnames_12
colnames(B23_12) <- columnnames_12
colnames(B24_12) <- columnnames_12
colnames(B25_12) <- columnnames_12
colnames(B26_12) <- columnnames_12

#now group into the different environments.
arg_plank_12 <- rbind(P1_12, P2_12, P3_12, P4_12, P5_12, P6_12)
glu_plank_12 <- rbind(P21_12, P22_12, P23_12, P24_12, P25_12, P26_12)
arg_biof_12 <- rbind(B1_12, B2_12, B3_12, B4_12, B5_12, B6_12)
glu_biof_12 <- rbind(B21_12, B22_12, B23_12, B24_12, B25_12, B26_12)

glucose_12 <- rbind(glu_plank_12, glu_biof_12)
arginine_12 <- rbind(arg_biof_12, arg_plank_12)
biofilm_12 <- rbind(arg_biof_12, glu_biof_12)
planktonic_12 <- rbind(glu_plank_12, arg_plank_12)

all_12 <- rbind(arg_biof_12, arg_plank_12, glu_biof_12, glu_plank_12)


#####remove rows that are zero

arg_plank_12 <- arg_plank_12[!(arg_plank_12$Day_12 == 0),]
arg_biof_12 <- arg_biof_12[!(arg_biof_12$Day_12 == 0),]

glu_plank_12 <- glu_plank_12[!(glu_plank_12$Day_12 == 0),]
glu_biof_12 <- glu_biof_12[!(glu_biof_12$Day_12 == 0),]

glucose_12 <- glucose_12[!(glucose_12$Day_12 == 0),]
arginine_12 <- arginine_12[!(arginine_12$Day_12 == 0),]

biofilm_12 <- biofilm_12[!(biofilm_12$Day_12 == 0),]
planktonic_12 <- planktonic_12[!(planktonic_12$Day_12 == 0),]


all_12 <- all_12[!(all_12$Day_12 == 0),]

write.csv(all_12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/Day12_bray_curtis_formatted.csv")


######
#I need to make day 25 tables for each population, where the columns are Gene, frequency, and sample ID. rows will be genes.

#####
columnnames_25 <- c("Day_25", "Gene", "Sample")

P1_25 <- P1[,c(7,11)]
P1_25$Sample <- "P1"
P2_25 <- P2[,c(7,11)]
P2_25$Sample <- "P2"
P3_25 <- P3[,c(7,11)]
P3_25$Sample <- "P3" #no mutations
P4_25 <- P4[,c(7,11)]
P4_25$Sample <- "P4"
P5_25 <- P5[,c(7,11)]
P5_25$Sample <- "P5"
P6_25 <- P6[,c(7,11)]
P6_25$Sample <- "P6"

P21_25 <- P21[,c(7,11)]
P21_25$Sample <- "P21"
P22_25 <- P22[,c(7,11)]
P22_25$Sample <- "P22"
P23_25 <- P23[,c(7,11)]
P23_25$Sample <- "P23"
P24_25 <- P24[,c(7,11)]
P24_25$Sample <- "P24"
P25_25 <- P25[,c(7,11)]
P25_25$Sample <- "P25"
P26_25 <- P26[,c(7,11)]
P26_25$Sample <- "P26" #no mutations

B1_25 <- B1[,c(7,11)]
B1_25$Sample <- "B1"
B2_25 <- B2[,c(7,11)]
B2_25$Sample <- "B2"
B3_25 <- B3[,c(7,11)]
B3_25$Sample <- "B3"
B4_25 <- B4[,c(7,11)]
B4_25$Sample <- "B4"
B5_25 <- B5[,c(7,11)]
B5_25$Sample <- "B5"
B6_25 <- B6[,c(7,11)]
B6_25$Sample <- "B6"

B21_25 <- B21[,c(7,11)]
B21_25$Sample <- "B21"
B22_25 <- B22[,c(7,11)]
B22_25$Sample <- "B22"
B23_25 <- B23[,c(7,11)]
B23_25$Sample <- "B23"
B24_25 <- B24[,c(7,11)]
B24_25$Sample <- "B24"
B25_25 <- B25[,c(7,11)]
B25_25$Sample <- "B25"
B26_25 <- B26[,c(7,11)]
B26_25$Sample <- "B26" #no mutations

colnames(P1_25) <- columnnames_25
colnames(P2_25) <- columnnames_25
colnames(P3_25) <- columnnames_25 #no mutations
colnames(P4_25) <- columnnames_25
colnames(P5_25) <- columnnames_25
colnames(P6_25) <- columnnames_25

colnames(P21_25) <- columnnames_25
colnames(P22_25) <- columnnames_25
colnames(P23_25) <- columnnames_25
colnames(P24_25) <- columnnames_25
colnames(P25_25) <- columnnames_25
colnames(P26_25) <- columnnames_25 #no mutations

colnames(B1_25) <- columnnames_25
colnames(B2_25) <- columnnames_25
colnames(B3_25) <- columnnames_25
colnames(B4_25) <- columnnames_25
colnames(B5_25) <- columnnames_25
colnames(B6_25) <- columnnames_25

colnames(B21_25) <- columnnames_25
colnames(B22_25) <- columnnames_25
colnames(B23_25) <- columnnames_25
colnames(B24_25) <- columnnames_25
colnames(B25_25) <- columnnames_25
colnames(B26_25) <- columnnames_25 #no mutations


#now group into the different environments.
arg_plank_25 <- rbind(P1_25, P2_25, P4_25, P5_25, P6_25)
glu_plank_25 <- rbind(P21_25, P22_25, P23_25, P24_25, P25_25)
arg_biof_25 <- rbind(B1_25, B2_25, B3_25, B4_25, B5_25, B6_25)
glu_biof_25 <- rbind(B21_25, B22_25, B23_25, B24_25, B25_25)

glucose_25 <- rbind(glu_plank_25, glu_biof_25)
arginine_25 <- rbind(arg_biof_25, arg_plank_25)
biofilm_25 <- rbind(arg_biof_25, glu_biof_25)
planktonic_25 <- rbind(glu_plank_25, arg_plank_25)

all_25 <- rbind(arg_biof_25, arg_plank_25, glu_biof_25, glu_plank_25)

#####remove rows that are zero

arg_plank_25 <- arg_plank_25[!(arg_plank_25$Day_25 == 0),]
arg_biof_25 <- arg_biof_25[!(arg_biof_25$Day_25 == 0),]

glu_plank_25 <- glu_plank_25[!(glu_plank_25$Day_25 == 0),]
glu_biof_25 <- glu_biof_25[!(glu_biof_25$Day_25 == 0),]

glucose_25 <- glucose_25[!(glucose_25$Day_25 == 0),]
arginine_25 <- arginine_25[!(arginine_25$Day_25 == 0),]

biofilm_25 <- biofilm_25[!(biofilm_25$Day_25 == 0),]
planktonic_25 <- planktonic_25[!(planktonic_25$Day_25 == 0),]


all_25 <- all_25[!(all_25$Day_25 == 0),]

write.csv(all_25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/Day25_bray_curtis_formatted.csv")


######
#I am now going to make tables of numbers of mutations identified for each population in the different groupings. I will be doing this for all three of the days. The tables include only those mutations that are non zero frequencies at the time point. So I will be going in and putting 0,s into the data tables for populations that do not have any observed mutations. 

#####
#day 6 first
arg_plank_6_table <- as.data.frame(table(arg_plank_6$Sample)) #P6 does not have any
arg_biof_6_table <- as.data.frame(table(arg_biof_6$Sample)) #
glu_plank_6_table <- as.data.frame(table(glu_plank_6$Sample))#P22 and P24 do not have any
glu_biof_6_table <- as.data.frame(table(glu_biof_6$Sample))
glucose_6_table <- as.data.frame(table(glucose_6$Sample))
arginine_6_table <- as.data.frame(table(arginine_6$Sample))
biofilm_6_table <- as.data.frame(table(biofilm_6$Sample))
planktonic_6_table <- as.data.frame(table(planktonic_6$Sample))
all_6_table <- as.data.frame(table(all_6$Sample))

#day 12
arg_plank_12_table <- as.data.frame(table(arg_plank_12$Sample))
arg_biof_12_table <- as.data.frame(table(arg_biof_12$Sample))
glu_plank_12_table <- as.data.frame(table(glu_plank_12$Sample))
glu_biof_12_table <- as.data.frame(table(glu_biof_12$Sample))
glucose_12_table <- as.data.frame(table(glucose_12$Sample))
arginine_12_table <- as.data.frame(table(arginine_12$Sample))
biofilm_12_table <- as.data.frame(table(biofilm_12$Sample))
planktonic_12_table <- as.data.frame(table(planktonic_12$Sample))
all_12_table <- as.data.frame(table(all_12$Sample))

#day 25
arg_plank_25_table <- as.data.frame(table(arg_plank_25$Sample))
arg_biof_25_table <- as.data.frame(table(arg_biof_25$Sample)) 
glu_plank_25_table <- as.data.frame(table(glu_plank_25$Sample))
glu_biof_25_table <- as.data.frame(table(glu_biof_25$Sample))
glucose_25_table <- as.data.frame(table(glucose_25$Sample))
arginine_25_table <- as.data.frame(table(arginine_25$Sample))
biofilm_25_table <- as.data.frame(table(biofilm_25$Sample))
planktonic_25_table <- as.data.frame(table(planktonic_25$Sample))
all_25_table <- as.data.frame(table(all_25$Sample))

#print off the tables that include all of the populations for each time point.
write.csv(all_6_table, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_numbers/all_6.csv")
write.csv(all_12_table, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_numbers/all_12.csv")
write.csv(all_25_table, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/mutation_numbers/all_25.csv")










########
#Now I want to be able to do bray curtis. I am going to do day 25. I think that I am going to have to do all three time points, but I am starting with Day 25 to see how similar the populations are at the end of the experiment.

#######
#going to first calculate bray curtis at the position level

#######
library(reshape2)
library(vegan)
library(gtools)
library(combinat)
library(data.table)

#read in the individual population time series data. 
#####
#read in the files
#####
#read in the files
P1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p1_final.csv", stringsAsFactors = F)
P2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p2_final.csv", stringsAsFactors = F)
P2 <- P2[1:8,]
P3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p3_final.csv", stringsAsFactors = F)
P4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p4_final.csv", stringsAsFactors = F)
P4 <- P4[1:7,]
P5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p5_final.csv", stringsAsFactors = F)
P5 <- P5[1:4,]
P6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p6_final.csv", stringsAsFactors = F)

P21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p21_final.csv", stringsAsFactors = F)
P22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p22_final.csv", stringsAsFactors = F)
P23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p23_final.csv", stringsAsFactors = F)
P24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p24_final.csv", stringsAsFactors = F)
P24 <- P24[1:6,]
P25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p25_final.csv", stringsAsFactors = F)
P26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p26_final.csv", stringsAsFactors = F)

B1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b1_final.csv", stringsAsFactors = F)
B2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b2_final.csv", stringsAsFactors = F)
B3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b3_final.csv", stringsAsFactors = F)
B4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b4_final.csv", stringsAsFactors = F)
B5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b5_final.csv", stringsAsFactors = F)
B6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b6_final.csv", stringsAsFactors = F)

B21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b21_final.csv", stringsAsFactors = F)
B22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b22_final.csv", stringsAsFactors = F)
B23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b23_final.csv", stringsAsFactors = F)
B24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b24_final.csv", stringsAsFactors = F)
B25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b25_final.csv", stringsAsFactors = F)
B26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b26_final.csv", stringsAsFactors = F)




#and the day 12 information
#####
columnnames_bc <- c("Frequency","Gene","Sample")

P1_bc_12 <- P1[,c(6,3)]
P1_bc_12$Sample <- "P1"
colnames(P1_bc_12) <- columnnames_bc
P2_bc_12 <- P2[,c(6,3)]
P2_bc_12$Sample <- "P2"
colnames(P2_bc_12) <- columnnames_bc
P3_bc_12 <- P3[,c(6,3)]
P3_bc_12$Sample <- "P3"
colnames(P3_bc_12) <- columnnames_bc
P4_bc_12 <- P4[,c(6,3)]
P4_bc_12$Sample <- "P4"
colnames(P4_bc_12) <- columnnames_bc
P5_bc_12 <- P5[,c(6,3)]
P5_bc_12$Sample <- "P5"
colnames(P5_bc_12) <- columnnames_bc
P6_bc_12 <- P6[,c(6,3)]
P6_bc_12$Sample <- "P6"
colnames(P6_bc_12) <- columnnames_bc

P21_bc_12 <- P21[,c(7,4)]
P21_bc_12$Sample <- "P21"
colnames(P21_bc_12) <- columnnames_bc
P22_bc_12 <- P22[,c(6,3)]
P22_bc_12$Sample <- "P22"
colnames(P22_bc_12) <- columnnames_bc
P23_bc_12 <- P23[,c(6,3)]
P23_bc_12$Sample <- "P23"
colnames(P23_bc_12) <- columnnames_bc
P24_bc_12 <- P24[,c(6,3)]
P24_bc_12$Sample <- "P24"
colnames(P24_bc_12) <- columnnames_bc
P25_bc_12 <- P25[,c(6,3)]
P25_bc_12$Sample <- "P25"
colnames(P25_bc_12) <- columnnames_bc
P26_bc_12 <- P26[,c(6,3)]
P26_bc_12$Sample <- "P26"
colnames(P26_bc_12) <- columnnames_bc

B1_bc_12 <- B1[,c(6,3)]
B1_bc_12$Sample <- "B1"
colnames(B1_bc_12) <- columnnames_bc
B2_bc_12 <- B2[,c(6,3)]
B2_bc_12$Sample <- "B2"
colnames(B2_bc_12) <- columnnames_bc
B3_bc_12 <- B3[,c(6,3)]
B3_bc_12$Sample <- "B3"
colnames(B3_bc_12) <- columnnames_bc
B4_bc_12 <- B4[,c(6,3)]
B4_bc_12$Sample <- "B4"
colnames(B4_bc_12) <- columnnames_bc
B5_bc_12 <- B5[,c(6,3)]
B5_bc_12$Sample <- "B5"
colnames(B5_bc_12) <- columnnames_bc
B6_bc_12 <- B6[,c(6,3)]
B6_bc_12$Sample <- "B6"
colnames(B6_bc_12) <- columnnames_bc

B21_bc_12 <- B21[,c(6,3)]
B21_bc_12$Sample <- "B21"
colnames(B21_bc_12) <- columnnames_bc
B22_bc_12 <- B22[,c(6,3)]
B22_bc_12$Sample <- "B22"
colnames(B22_bc_12) <- columnnames_bc
B23_bc_12 <- B23[,c(6,3)]
B23_bc_12$Sample <- "B23"
colnames(B23_bc_12) <- columnnames_bc
B24_bc_12 <- B24[,c(6,3)]
B24_bc_12$Sample <- "B24"
colnames(B24_bc_12) <- columnnames_bc
B25_bc_12 <- B25[,c(6,3)]
B25_bc_12$Sample <- "B25"
colnames(B25_bc_12) <- columnnames_bc
B26_bc_12 <- B26[,c(6,3)]
B26_bc_12$Sample <- "B26"
colnames(B26_bc_12) <- columnnames_bc


######
#and day 25
######

P1_bc_25 <- P1[,c(8,3)]
P1_bc_25$Sample <- "P1"
colnames(P1_bc_25) <- columnnames_bc
P2_bc_25 <- P2[,c(8,3)]
P2_bc_25$Sample <- "P2"
colnames(P2_bc_25) <- columnnames_bc
P3_bc_25 <- P3[,c(8,3)]
P3_bc_25$Sample <- "P3"
colnames(P3_bc_25) <- columnnames_bc
P4_bc_25 <- P4[,c(8,3)]
P4_bc_25$Sample <- "P4"
colnames(P4_bc_25) <- columnnames_bc
P5_bc_25 <- P5[,c(8,3)]
P5_bc_25$Sample <- "P5"
colnames(P5_bc_25) <- columnnames_bc
P6_bc_25 <- P6[,c(8,3)]
P6_bc_25$Sample <- "P6"
colnames(P6_bc_25) <- columnnames_bc

P21_bc_25 <- P21[,c(9,4)]
P21_bc_25$Sample <- "P21"
colnames(P21_bc_25) <- columnnames_bc
P22_bc_25 <- P22[,c(8,3)]
P22_bc_25$Sample <- "P22"
colnames(P22_bc_25) <- columnnames_bc
P23_bc_25 <- P23[,c(8,3)]
P23_bc_25$Sample <- "P23"
colnames(P23_bc_25) <- columnnames_bc
P24_bc_25 <- P24[,c(8,3)]
P24_bc_25$Sample <- "P24"
colnames(P24_bc_25) <- columnnames_bc
P25_bc_25 <- P25[,c(8,3)]
P25_bc_25$Sample <- "P25"
colnames(P25_bc_25) <- columnnames_bc
P26_bc_25 <- P26[,c(8,3)]
P26_bc_25$Sample <- "P26"
colnames(P26_bc_25) <- columnnames_bc

B1_bc_25 <- B1[,c(8,3)]
B1_bc_25$Sample <- "B1"
colnames(B1_bc_25) <- columnnames_bc
B2_bc_25 <- B2[,c(8,3)]
B2_bc_25$Sample <- "B2"
colnames(B2_bc_25) <- columnnames_bc
B3_bc_25 <- B3[,c(8,3)]
B3_bc_25$Sample <- "B3"
colnames(B3_bc_25) <- columnnames_bc
B4_bc_25 <- B4[,c(8,3)]
B4_bc_25$Sample <- "B4"
colnames(B4_bc_25) <- columnnames_bc
B5_bc_25 <- B5[,c(8,3)]
B5_bc_25$Sample <- "B5"
colnames(B5_bc_25) <- columnnames_bc
B6_bc_25 <- B6[,c(8,3)]
B6_bc_25$Sample <- "B6"
colnames(B6_bc_25) <- columnnames_bc

B21_bc_25 <- B21[,c(8,3)]
B21_bc_25$Sample <- "B21"
colnames(B21_bc_25) <- columnnames_bc
B22_bc_25 <- B22[,c(8,3)]
B22_bc_25$Sample <- "B22"
colnames(B22_bc_25) <- columnnames_bc
B23_bc_25 <- B23[,c(8,3)]
B23_bc_25$Sample <- "B23"
colnames(B23_bc_25) <- columnnames_bc
B24_bc_25 <- B24[,c(8,3)]
B24_bc_25$Sample <- "B24"
colnames(B24_bc_25) <- columnnames_bc
B25_bc_25 <- B25[,c(8,3)]
B25_bc_25$Sample <- "B25"
colnames(B25_bc_25) <- columnnames_bc
B26_bc_25 <- B26[,c(8,3)]
B26_bc_25$Sample <- "B26"
colnames(B26_bc_25) <- columnnames_bc
#####
#now create one data frame for each day and print to a file
#####
bc_6 <- rbind(P1_bc_6,
              P2_bc_6,
              #P3_bc_6,
              P4_bc_6,
              P5_bc_6,
              P6_bc_6,
              
              P21_bc_6,
              P22_bc_6,
              P23_bc_6,
              P24_bc_6,
              P25_bc_6,
              #P26_bc_6,
              
              B1_bc_6,
              B2_bc_6,
              B3_bc_6,
              B4_bc_6,
              B5_bc_6,
              B6_bc_6,
              
              B21_bc_6,
              B22_bc_6,
              B23_bc_6,
              B24_bc_6,
              B25_bc_6,
              B26_bc_6)

bc_12 <- rbind(P1_bc_12,
               P2_bc_12,
               P3_bc_12,
               P4_bc_12,
               P5_bc_12,
               P6_bc_12,
               
               P21_bc_12,
               P22_bc_12,
               P23_bc_12,
               P24_bc_12,
               P25_bc_12,
               P26_bc_12,
               
               B1_bc_12,
               B2_bc_12,
               B3_bc_12,
               B4_bc_12,
               B5_bc_12,
               B6_bc_12,
               
               B21_bc_12,
               B22_bc_12,
               B23_bc_12,
               B24_bc_12,
               B25_bc_12,
               B26_bc_12)


bc_25 <- rbind(P1_bc_25,
               P2_bc_25,
               P4_bc_25,
               P5_bc_25,
               P6_bc_25,
               
               P21_bc_25,
               P22_bc_25,
               P23_bc_25,
               P24_bc_25,
               P25_bc_25,
               
               B1_bc_25,
               B2_bc_25,
               B3_bc_25,
               B4_bc_25,
               B5_bc_25,
               B6_bc_25,
               
               B21_bc_25,
               B22_bc_25,
               B23_bc_25,
               B24_bc_25,
               B25_bc_25)
#write this to a file
#write.csv(bc_6, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day6_bray_curtis_formatted.csv")
write.csv(bc_12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day12_bray_curtis_formatted.csv")
write.csv(bc_25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day25_bray_curtis_formatted.csv")




######
#read in the data formatted for bray crtis. Doing this so that It is always in the same format and I don't have to run all of the lines of code each time I do this.

#####
#bc_6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day6_bray_curtis_formatted.csv", stringsAsFactors = F)
bc_12 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day12_bray_curtis_formatted.csv", stringsAsFactors = F)
bc_25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day25_bray_curtis_formatted.csv", stringsAsFactors = F)
######
#reformat data frame to be able to perform bray curtis analysis
#####


#Start with Day 6
m_bc_6 <- melt(bc_6, id=c("Gene", "Sample"), measure.vars = "Frequency")
cast_bc_6 <- as.data.frame(t(dcast(m_bc_6, Sample~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_bc_6) <- as.character(unlist(cast_bc_6[1,]))
braycurtis_6 <- cast_bc_6[-1,]

#now reorder the columns
colorder <- c("P1","P2","P3","P4","P5","P6",
              "P21","P22","P23","P24","P25","P26",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25","B26")
setcolorder(braycurtis_6, colorder)
#need to make all of the entries numeric
braycurtis_6 <- apply(braycurtis_6, 2, function(x) as.numeric(as.character(x)))
#transform the data so that each row is a population and the columns are the diffeent mutations
t_braycurtis_6 <- t(braycurtis_6)

#####Then day 12
m_bc_12 <- melt(bc_12, id=c("Gene", "Sample"), measure.vars = "Frequency")
cast_bc_12 <- as.data.frame(t(dcast(m_bc_12, Sample~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_bc_12) <- as.character(unlist(cast_bc_12[1,]))
braycurtis_12 <- cast_bc_12[-1,]

#now reorder the columns
colorder <- c("P1","P2","P3","P4","P5","P6",
              "P21","P22","P23","P24","P25","P26",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25","B26")
setcolorder(braycurtis_12, colorder)
#need to make all of the entries numeric
braycurtis_12 <- apply(braycurtis_12, 2, function(x) as.numeric(as.character(x)))
#transform the data so that each row is a population and the columns are the diffeent mutations
t_braycurtis_12 <- t(braycurtis_12)

#finaly day 25

m_bc_25 <- melt(bc_25, id=c("Gene", "Sample"), measure.vars = "Frequency")
cast_bc_25 <- as.data.frame(t(dcast(m_bc_25, Sample~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_bc_25) <- as.character(unlist(cast_bc_25[1,]))
braycurtis_25 <- cast_bc_25[-1,]

#now reorder the columns
colorder <- c("P1","P2","P4","P5","P6",
              "P21","P22","P23","P24","P25",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25")
setcolorder(braycurtis_25, colorder)
braycurtis_25 <- braycurtis_25[c(12,13,19,20,21,14,15,16,17,18,1,2,8,9,10,11,3,4,5,6,7),]
#need to make all of the entries numeric
braycurtis_25 <- apply(braycurtis_25, 2, function(x) as.numeric(as.character(x)))
#transform the data so that each row is a population and the columns are the diffeent mutations
t_braycurtis_25 <- t(braycurtis_25)

########
#actually calculate the bray curtis dissimilarity

#using Caroline's code for Bray curtis dissimilarity, but I have modified for use with my populations 

######


#day 6
#now going to calculate bray curtis and turn it into a data frame so that I can visualize the data.
braycurtis6 <- vegdist(t_braycurtis_6, method = "bray") #actually calculating the bray curtis dissimilarity for each population vs. each other population. 
bc_matrix6 <- as.matrix(braycurtis6) #make this as a matrix to be able to look at the data.
bcsim_matrix6 <- (1-bc_matrix6) #reverse so that I'm looking at the similarity
treat_list <- c(rep("AP",6), rep("GP",6), rep ("AB",6), rep("GB",6)) #list that will tell me what treatments the various populations are evolved in. AP = Arginine Plankotnic, GP = Glucose Planktonic, AB = Arginine Biofilm, GB = Glucose biofilm
rownames(bcsim_matrix6) <- treat_list #make the row names the treatments
colnames(bcsim_matrix6) <-treat_list #also make the column names the treatments
braycurtis_frame6 <- as.data.frame(as.matrix(braycurtis6)) #make it a data frame


#write out the matrix
write.csv(braycurtis_frame6, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day6.csv")

#day 12
#now going to calculate bray curtis and turn it into a data frame so that I can visualize the data.
braycurtis12 <- vegdist(t_braycurtis_12, method = "bray") #actually calculating the bray curtis dissimilarity for each population vs. each other population. 
bc_matrix12 <- as.matrix(braycurtis12) #make this as a matrix to be able to look at the data.
bcsim_matrix12 <- (1-bc_matrix12) #reverse so that I'm looking at the similarity
treat_list <- c(rep("AP",6), rep("GP",6), rep ("AB",6), rep("GB",6)) #list that will tell me what treatments the various populations are evolved in. AP = Arginine Plankotnic, GP = Glucose Planktonic, AB = Arginine Biofilm, GB = Glucose biofilm
rownames(bcsim_matrix12) <- treat_list #make the row names the treatments
colnames(bcsim_matrix12) <-treat_list #also make the column names the treatments
braycurtis_frame12 <- as.data.frame(as.matrix(braycurtis12)) #make it a data frame


#write out the matrix
write.csv(braycurtis_frame12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day12.csv")
write.csv(bcsim_matrix12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day12_similarity.csv")



#day 25
#now going to calculate bray curtis and turn it into a data frame so that I can visualize the data.
braycurtis25 <- vegdist(t_braycurtis_25, method = "bray") #actually calculating the bray curtis dissimilarity for each population vs. each other population. 
bc_matrix25 <- as.matrix(braycurtis25) #make this as a matrix to be able to look at the data.
bcsim_matrix25 <- (1-bc_matrix25) #reverse so that I'm looking at the similarity
treat_list <- c(rep("AP",5), rep("GP",5), rep ("AB",6), rep("GB",5)) #list that will tell me what treatments the various populations are evolved in. AP = Arginine Plankotnic, GP = Glucose Planktonic, AB = Arginine Biofilm, GB = Glucose biofilm
rownames(bcsim_matrix25) <- treat_list #make the row names the treatments
colnames(bcsim_matrix25) <-treat_list #also make the column names the treatments
braycurtis_frame25 <- as.data.frame(as.matrix(braycurtis25)) #make it a data frame


#write out the matrix
write.csv(braycurtis_frame25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day25.csv")
write.csv(bcsim_matrix25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day25_similarity.csv")


######

#now I am creating functions that I will then be able to apply to the data frame. I was smart and kept Caroline's same names, so this should work.....
######

#A function that tests for the bray curtis dissimilarity within a treatment group.
bc_within <- function(within){
  #returns mean within treatment values in a matrix
  sum <- 0
  count <- 0
  for(i in 1:(nrow(within)-1)){
    for(j in (i+1):ncol(within)){
      if(rownames(within)[i] == colnames(within)[j]){
        sum <- sum+within[i,j]
        count <- count+1
      }
    }
  }
  return(sum/count)
}
#within <- as.matrix(bc_within(bcsim_matrix))

#and now a function that tests the bray curtis dissimilarity between two different treatment groups.
bc_between <- function(between){
  #returns mean between treatment values in a matrix
  sum <- 0
  count <- 0
  for(i in 1:(nrow(between)-1)){
    for(j in (i+1):ncol(between)){
      if(rownames(between)[i] != colnames(between)[j]){
        sum <- sum+between[i,j]
        count <- count+1
      }
    }
  }
  return(sum/count)
}

#and this is a function that will find the difference between the same and the different treatments.
bc_similarity <- function(similar){
  #returns difference between mean bc of similar and non-similar (different) treatments
  sum_sim <- 0
  count_sim <- 0
  sum_diff <- 0
  count_diff <- 0
  for(i in 1:(nrow(similar)-1)){
    for(j in (i+1):ncol(similar)){
      if(rownames(similar)[i] != colnames(similar)[j]){
        if(rownames(similar)[i]==("GP")||rownames(similar)[i]=="AB"){ 
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("AB") && colnames(similar)[j]=="GP"){ 
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("GB") && colnames(similar)[j]=="AP"){ #testing carbon source in biofilm
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("AP") && colnames(similar)[j]=="GB"){ #testing carbon source differences in planktonic
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
          
        } else {
          sum_sim <- sum_sim+similar[i,j]
          count_sim <- count_sim+1
        }
      }
    }
  }
  return(sum_sim/count_sim-sum_diff/count_diff)
}



bc_within(bcsim_matrix12) #0.3087198
bc_between(bcsim_matrix12) #0.2142211
bc_similarity(bcsim_matrix12) #-0.003767487

bc_within(bcsim_matrix25) #0.2986108
bc_between(bcsim_matrix25) #0.1883547
bc_similarity(bcsim_matrix25) #0.04395233

#difference6 <- bc_within(bcsim_matrix6)-bc_between(bcsim_matrix6) #NA
difference12 <- bc_within(bcsim_matrix12)-bc_between(bcsim_matrix12) #0.09449874
difference25 <- bc_within(bcsim_matrix25)-bc_between(bcsim_matrix25) #0.11025615243



#Testing whether bray-curtis values are more similar within treatments than between treatments

#day 12

difference_rand <- character(length(10^5))
envir_similarity <- character(length(10^5))
randomizing_matrix <-bcsim_matrix12

for(i in 1:10^5)
{
  treat_list_random <- sample(treat_list)
  rownames(randomizing_matrix) <- treat_list_random
  colnames(randomizing_matrix) <- treat_list_random
  between <- bc_between(randomizing_matrix)
  within <- bc_within(randomizing_matrix)
  similarity <- bc_similarity(randomizing_matrix)
  difference_rand[i] <- (between-within)
  envir_similarity[i] <- similarity
}

mean(as.numeric(difference_rand)>=difference12) #0
mean(as.numeric(envir_similarity)>=bc_similarity(bcsim_matrix12)) #0.054981

upermn <- function(x) {
  #calculates number of unique permutations of x
  #copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  n <- length(x)
  duplicates <- as.numeric(table(x))
  factorial(n) / prod(factorial(duplicates))
}

uperm <- function(x) {
  #returns a list of all the unique permutations of the elements of vector x
  # copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  u <- sort(unique(x))
  l <- length(u)
  if (l == length(x)) {
    return(do.call(rbind,permn(x)))
  }
  if (l == 1) return(x)
  result <- matrix(NA, upermn(x), length(x))
  index <- 1
  for (i in 1:l) {
    v <- x[-which(x==u[i])[1]]
    newindex <- upermn(v)
    if (table(x)[i] == 1) {
      result[index:(index+newindex-1),] <- cbind(u[i], do.call(rbind, unique(permn(v))))
    } else {
      result[index:(index+newindex-1),] <- cbind(u[i], uperm(v))
    }
    index <- index+newindex
  }
  return(result)
}

permutation_bc <- function(x){
  #
  p <- uperm(rownames(x))
  temp <- x
  diff_rand <- numeric(length=nrow(p))
  for(row in 1:nrow(p)){
    rownames(temp) <- p[row,]
    colnames(temp) <- p[row,]
    diff_rand[row] <- (bc_within(temp) - bc_between(temp))
  }
  return(diff_rand)
}

#Testing whether pairs of treatments differ significantly more between treatments than within treatments


# Loop through each treatment pair
treatments <- c("AP","GP","AB", "GB")
for(m in 1:(length(treatments)-1)){
  for(n in (m+1):length(treatments)){
    a <- treatments[m]
    b <- treatments[n]
    
    #make a matrix which only includes 2 selested treatments
    include_list <- c(a,b)
    subcomm <- bcsim_matrix12[include_list,include_list]
    subcommunity<- bcsim_matrix12[(rownames(bcsim_matrix12)== a)|(rownames(bcsim_matrix12)== b),(colnames(bcsim_matrix12)==a | colnames(bcsim_matrix12)==b)]
    
    #Calculate difference between between and within treatment bray curtis
    between_ab <- bc_between(subcommunity)
    within_ab <- bc_within(subcommunity)
    #For all permutations of treatment assignemnts, calculate within-between
    permutation_diff <- permutation_bc(subcommunity)
    #Calculate significance
    significance <- mean(as.numeric(permutation_diff)>=(within_ab-between_ab))
    #Record outcome
    print(c(a,b,significance))
  }
}
#[1] "AP"                "GP"                "0.720779220779221"
#[1] "AP"                  "AB"                  "0.00432900432900433"
#[1] "AP"                "GB"                "0.175324675324675"
#[1] "GP"                  "AB"                  "0.00432900432900433"
#[1] "GP"                "GB"                "0.151515151515152"
#[1] "AB"                  "GB"                  "0.00432900432900433"




#for day 25
difference_rand <- character(length(10^5))
envir_similarity <- character(length(10^5))
randomizing_matrix <-bcsim_matrix25

for(i in 1:10^5)
{
  treat_list_random <- sample(treat_list)
  rownames(randomizing_matrix) <- treat_list_random
  colnames(randomizing_matrix) <- treat_list_random
  between <- bc_between(randomizing_matrix)
  within <- bc_within(randomizing_matrix)
  similarity <- bc_similarity(randomizing_matrix)
  difference_rand[i] <- (between-within)
  envir_similarity[i] <- similarity
}

mean(as.numeric(difference_rand)>=difference25) #0
mean(as.numeric(envir_similarity)>=bc_similarity(bcsim_matrix25)) #0.03089

upermn <- function(x) {
  #calculates number of unique permutations of x
  #copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  n <- length(x)
  duplicates <- as.numeric(table(x))
  factorial(n) / prod(factorial(duplicates))
}

uperm <- function(x) {
  #returns a list of all the unique permutations of the elements of vector x
  # copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  u <- sort(unique(x))
  l <- length(u)
  if (l == length(x)) {
    return(do.call(rbind,permn(x)))
  }
  if (l == 1) return(x)
  result <- matrix(NA, upermn(x), length(x))
  index <- 1
  for (i in 1:l) {
    v <- x[-which(x==u[i])[1]]
    newindex <- upermn(v)
    if (table(x)[i] == 1) {
      result[index:(index+newindex-1),] <- cbind(u[i], do.call(rbind, unique(permn(v))))
    } else {
      result[index:(index+newindex-1),] <- cbind(u[i], uperm(v))
    }
    index <- index+newindex
  }
  return(result)
}

permutation_bc <- function(x){
  #
  p <- uperm(rownames(x))
  temp <- x
  diff_rand <- numeric(length=nrow(p))
  for(row in 1:nrow(p)){
    rownames(temp) <- p[row,]
    colnames(temp) <- p[row,]
    diff_rand[row] <- (bc_within(temp) - bc_between(temp))
  }
  return(diff_rand)
}

#Testing whether pairs of treatments differ significantly more between treatments than within treatments


# Loop through each treatment pair
treatments <- c("AP","GP","AB", "GB")
for(m in 1:(length(treatments)-1)){
  for(n in (m+1):length(treatments)){
    a <- treatments[m]
    b <- treatments[n]
    
    #make a matrix which only includes 2 selested treatments
    include_list <- c(a,b)
    subcomm <- bcsim_matrix25[include_list,include_list]
    subcommunity<- bcsim_matrix25[(rownames(bcsim_matrix25)== a)|(rownames(bcsim_matrix25)== b),(colnames(bcsim_matrix25)==a | colnames(bcsim_matrix25)==b)]
    
    #Calculate difference between between and within treatment bray curtis
    between_ab <- bc_between(subcommunity)
    within_ab <- bc_within(subcommunity)
    #For all permutations of treatment assignemnts, calculate within-between
    permutation_diff <- permutation_bc(subcommunity)
    #Calculate significance
    significance <- mean(as.numeric(permutation_diff)>=(within_ab-between_ab))
    #Record outcome
    print(c(a,b,significance))
  }
}
#[1] "AP"              "GP"              "0.4004329004329"
#[1] "AP"                  "AB"                  "0.00216450216450216"
#[1] "AP"                  "GB"                  "0.00865800865800866"
#[1] "GP"                  "AB"                  "0.00216450216450216"
#[1] "GP"                 "GB"                 "0.0281385281385281"
#[1] "AB"                  "GB"                  "0.00432900432900433"



######
#now I need to calculate bray curtis at the gene level. 

####


#####
#read in the files
#####
P1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p1_final.csv", stringsAsFactors = F)
P2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p2_final.csv", stringsAsFactors = F)
P2 <- P2[1:8,]
P3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p3_final.csv", stringsAsFactors = F)
P4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p4_final.csv", stringsAsFactors = F)
P4 <- P4[1:7,]
P5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p5_final.csv", stringsAsFactors = F)
P5 <- P5[1:4,]
P6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p6_final.csv", stringsAsFactors = F)

P21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p21_final.csv", stringsAsFactors = F)
P22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p22_final.csv", stringsAsFactors = F)
P23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p23_final.csv", stringsAsFactors = F)
P24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p24_final.csv", stringsAsFactors = F)
P24 <- P24[1:6,]
P25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p25_final.csv", stringsAsFactors = F)
P26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/p26_final.csv", stringsAsFactors = F)

B1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b1_final.csv", stringsAsFactors = F)
B2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b2_final.csv", stringsAsFactors = F)
B3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b3_final.csv", stringsAsFactors = F)
B4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b4_final.csv", stringsAsFactors = F)
B5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b5_final.csv", stringsAsFactors = F)
B6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b6_final.csv", stringsAsFactors = F)

B21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b21_final.csv", stringsAsFactors = F)
B22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b22_final.csv", stringsAsFactors = F)
B23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b23_final.csv", stringsAsFactors = F)
B24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b24_final.csv", stringsAsFactors = F)
B25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b25_final.csv", stringsAsFactors = F)
B26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b26_final.csv", stringsAsFactors = F)

#Grab the day 12 frequency information and gene names
columnnames <- c("Frequency", "Gene", "Population")

P1_gene_12 <- P1[,c(6,11)]
P1_gene_12$pop <- "P1"
colnames(P1_gene_12) <- columnnames

P2_gene_12 <- P2[,c(6,11)]
P2_gene_12$pop <- "P2"
colnames(P2_gene_12) <- columnnames
#there are 2 mutations in the pitA gene, so I am only going to keep the highest frequency one.
P2_gene_12 <- P2_gene_12[-3,] #removes a pitA mutation at 0 frequency

P3_gene_12 <- P3[,c(6,11)]
P3_gene_12$pop <- "P3"
colnames(P3_gene_12) <- columnnames

P4_gene_12 <- P4[,c(6,11)]
P4_gene_12$pop <- "P4"
colnames(P4_gene_12) <- columnnames

P5_gene_12 <- P5[,c(6,11)]
P5_gene_12$pop <- "P5"
colnames(P5_gene_12) <- columnnames

P6_gene_12 <- P6[,c(6,11)]
P6_gene_12$pop <- "P6"
colnames(P6_gene_12) <- columnnames


P21_gene_12 <- P21[,c(7,12)]
P21_gene_12$pop <- "P21"
colnames(P21_gene_12) <- columnnames

P22_gene_12 <- P22[,c(6,11)]
P22_gene_12$pop <- "P22"
colnames(P22_gene_12) <- columnnames

P23_gene_12 <- P23[,c(6,11)]
P23_gene_12$pop <- "P23"
colnames(P23_gene_12) <- columnnames
#there are 2 fleQ mutations in this population so I am, again, choosing the highest frequency at this time point and getting rid of the lower frequency mutation. 
P23_gene_12 <- P23_gene_12[-4,]

P24_gene_12 <- P24[,c(6,11)]
P24_gene_12$pop <- "P24"
colnames(P24_gene_12) <- columnnames
#again has 2 fleQ mutations
P24_gene_12 <- P24_gene_12[-3,]

P25_gene_12 <- P25[,c(6,11)]
P25_gene_12$pop <- "P25"
colnames(P25_gene_12) <- columnnames

P26_gene_12 <- P26[,c(6,11)]
P26_gene_12$pop <- "P26"
colnames(P26_gene_12) <- columnnames


B1_gene_12 <- B1[,c(6,11)]
B1_gene_12$pop <- "B1"
colnames(B1_gene_12) <- columnnames

B2_gene_12 <- B2[,c(6,11)]
B2_gene_12$pop <- "B2"
colnames(B2_gene_12) <- columnnames
B2_gene_12 <- B2_gene_12[-c(2,7),] #had 2 01160 and intergenic upstream lasA mutations that I deleted the lower frequencies of.

B3_gene_12 <- B3[,c(6,11)]
B3_gene_12$pop <- "B3"
colnames(B3_gene_12) <- columnnames
B3_gene_12 <- B3_gene_12[-c(8,9),]

B4_gene_12 <- B4[,c(6,11)]
B4_gene_12$pop <- "B4"
colnames(B4_gene_12) <- columnnames
B4_gene_12 <- B4_gene_12[-c(3,9),]

B5_gene_12 <- B5[,c(6,11)]
B5_gene_12$pop <- "B5"
colnames(B5_gene_12) <- columnnames
B5_gene_12 <- B5_gene_12[-13,]

B6_gene_12 <- B6[,c(6,11)]
B6_gene_12$pop <- "B6"
colnames(B6_gene_12) <- columnnames


B21_gene_12 <- B21[,c(6,11)]
B21_gene_12$pop <- "B21"
colnames(B21_gene_12) <- columnnames
B21_gene_12 <- B21_gene_12[-c(3,5),]

B22_gene_12 <- B22[,c(6,11)]
B22_gene_12$pop <- "B22"
colnames(B22_gene_12) <- columnnames

B23_gene_12 <- B23[,c(6,11)]
B23_gene_12$pop <- "B23"
colnames(B23_gene_12) <- columnnames

B24_gene_12 <- B24[,c(6,11)]
B24_gene_12$pop <- "B24"
colnames(B24_gene_12) <- columnnames

B25_gene_12 <- B25[,c(6,11)]
B25_gene_12$pop <- "B25"
colnames(B25_gene_12) <- columnnames

B26_gene_12 <- B26[,c(6,11)]
B26_gene_12$pop <- "B26"
colnames(B26_gene_12) <- columnnames
B26_gene_12 <- B26_gene_12[-4,] #has 2 argJ mutations 

#####
#now to create a data frame with all populations in it for day 12

bc_gene_12 <- rbind(P1_gene_12,
                    P2_gene_12,
                    P3_gene_12,
                    P4_gene_12,
                    P5_gene_12,
                    P6_gene_12,
                    
                    P21_gene_12,
                    P22_gene_12,
                    P23_gene_12,
                    P24_gene_12,
                    P25_gene_12,
                    P26_gene_12,
                    
                    B1_gene_12,
                    B2_gene_12,
                    B3_gene_12,
                    B4_gene_12,
                    B5_gene_12,
                    B6_gene_12,
                    
                    B21_gene_12,
                    B22_gene_12,
                    B23_gene_12,
                    B24_gene_12,
                    B25_gene_12,
                    B26_gene_12)
######

write.csv(bc_gene_12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day12_bray_curtis_formatted_by_gene.csv")
####
bc_12 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day12_bray_curtis_formatted_by_gene.csv", stringsAsFactors = F)


m_bc_12 <- melt(bc_12, id=c("Gene", "Population"), measure.vars = "Frequency")
cast_bc_12 <- as.data.frame(t(dcast(m_bc_12, Population~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_bc_12) <- as.character(unlist(cast_bc_12[1,]))
braycurtis_12 <- cast_bc_12[-1,]

#now reorder the columns
colorder <- c("P1","P2","P3","P4","P5","P6",
              "P21","P22","P23","P24","P25","P26",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25","B26")
setcolorder(braycurtis_12, colorder)
#need to make all of the entries numeric
braycurtis_12 <- apply(braycurtis_12, 2, function(x) as.numeric(as.character(x)))
#transform the data so that each row is a population and the columns are the diffeent mutations
t_braycurtis_12 <- t(braycurtis_12)

####
#and, finally, going to calculate the bray curtis dissimilarity for each combination

####
#now going to calculate bray curtis and turn it into a data frame so that I can visualize the data.
braycurtis12 <- vegdist(t_braycurtis_12, method = "bray") #actually calculating the bray curtis dissimilarity for each population vs. each other population. 
bc_matrix12 <- as.matrix(braycurtis12) #make this as a matrix to be able to look at the data.
bcsim_matrix12 <- (1-bc_matrix12) #reverse so that I'm looking at the similarity
treat_list <- c(rep("AP",6), rep("GP",6), rep ("AB",6), rep("GB",6)) #list that will tell me what treatments the various populations are evolved in. AP = Arginine Plankotnic, GP = Glucose Planktonic, AB = Arginine Biofilm, GB = Glucose biofilm
rownames(bcsim_matrix12) <- c("AP","AP","AP","AP","AP","AP","GP","GP","GP","GP","GP","GP","AB","AB","AB","AB","AB","AB","GB","GB","GB","GB","GB","GB")  #treat_list #make the row names the treatments
colnames(bcsim_matrix12) <-treat_list #also make the column names the treatments
braycurtis_frame12 <- as.data.frame(as.matrix(braycurtis12)) #make it a data frame


#write out the matrix
write.csv(braycurtis_frame12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day12_byGene.csv")
write.csv(bcsim_matrix12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/reverse_braycurtis_day12_byGene.csv")


#now do the day 25
columnnames <- c("Frequency", "Gene", "Population")

P1_gene_25 <- P1[,c(7,11)]
P1_gene_25$pop <- "P1"
colnames(P1_gene_25) <- columnnames

P2_gene_25 <- P2[,c(7,11)]
P2_gene_25$pop <- "P2"
colnames(P2_gene_25) <- columnnames
#there are 2 mutations in the pitA gene, so I am only going to keep the highest frequency one.
P2_gene_25 <- P2_gene_25[-2,] #removes a pitA mutation at 0 frequency

P3_gene_25 <- P3[,c(7,11)]
P3_gene_25$pop <- "P3"
colnames(P3_gene_25) <- columnnames

P4_gene_25 <- P4[,c(7,11)]
P4_gene_25$pop <- "P4"
colnames(P4_gene_25) <- columnnames

P5_gene_25 <- P5[,c(7,11)]
P5_gene_25$pop <- "P5"
colnames(P5_gene_25) <- columnnames

P6_gene_25 <- P6[,c(7,11)]
P6_gene_25$pop <- "P6"
colnames(P6_gene_25) <- columnnames


P21_gene_25 <- P21[,c(8,12)]
P21_gene_25$pop <- "P21"
colnames(P21_gene_25) <- columnnames

P22_gene_25 <- P22[,c(7,11)]
P22_gene_25$pop <- "P22"
colnames(P22_gene_25) <- columnnames

P23_gene_25 <- P23[,c(7,11)]
P23_gene_25$pop <- "P23"
colnames(P23_gene_25) <- columnnames
#there are 2 fleQ mutations in this population so I am, again, choosing the highest frequency at this time point and getting rid of the lower frequency mutation. 
P23_gene_25 <- P23_gene_25[-3,]

P24_gene_25 <- P24[,c(7,11)]
P24_gene_25$pop <- "P24"
colnames(P24_gene_25) <- columnnames
#again has 2 fleQ mutations
P24_gene_25 <- P24_gene_25[-2,]

P25_gene_25 <- P25[,c(7,11)]
P25_gene_25$pop <- "P25"
colnames(P25_gene_25) <- columnnames

P26_gene_25 <- P26[,c(7,11)]
P26_gene_25$pop <- "P26"
colnames(P26_gene_25) <- columnnames


B1_gene_25 <- B1[,c(7,11)]
B1_gene_25$pop <- "B1"
colnames(B1_gene_25) <- columnnames

B2_gene_25 <- B2[,c(7,11)]
B2_gene_25$pop <- "B2"
colnames(B2_gene_25) <- columnnames
B2_gene_25 <- B2_gene_25[-c(2,7),] #had 2 01160 and intergenic upstream lasA mutations that I deleted the lower frequencies of.

B3_gene_25 <- B3[,c(7,11)]
B3_gene_25$pop <- "B3"
colnames(B3_gene_25) <- columnnames
B3_gene_25 <- B3_gene_25[-c(9,10),]

B4_gene_25 <- B4[,c(7,11)]
B4_gene_25$pop <- "B4"
colnames(B4_gene_25) <- columnnames
B4_gene_25 <- B4_gene_25[-c(3,8),]

B5_gene_25 <- B5[,c(7,11)]
B5_gene_25$pop <- "B5"
colnames(B5_gene_25) <- columnnames
B5_gene_25 <- B5_gene_25[-12,]

B6_gene_25 <- B6[,c(7,11)]
B6_gene_25$pop <- "B6"
colnames(B6_gene_25) <- columnnames


B21_gene_25 <- B21[,c(7,11)]
B21_gene_25$pop <- "B21"
colnames(B21_gene_25) <- columnnames
B21_gene_25 <- B21_gene_25[-c(4,5),]

B22_gene_25 <- B22[,c(7,11)]
B22_gene_25$pop <- "B22"
colnames(B22_gene_25) <- columnnames

B23_gene_25 <- B23[,c(7,11)]
B23_gene_25$pop <- "B23"
colnames(B23_gene_25) <- columnnames

B24_gene_25 <- B24[,c(7,11)]
B24_gene_25$pop <- "B24"
colnames(B24_gene_25) <- columnnames

B25_gene_25 <- B25[,c(7,11)]
B25_gene_25$pop <- "B25"
colnames(B25_gene_25) <- columnnames

B26_gene_25 <- B26[,c(7,11)]
B26_gene_25$pop <- "B26"
colnames(B26_gene_25) <- columnnames
B26_gene_12 <- B26_gene_12[-4,] #has 2 argJ mutations 

#####
#now to create a data frame with all populations in it for day 12

bc_gene_25 <- rbind(P1_gene_25,
                    P2_gene_25,
                    P3_gene_25,
                    P4_gene_25,
                    P5_gene_25,
                    P6_gene_25,
                    
                    P21_gene_25,
                    P22_gene_25,
                    P23_gene_25,
                    P24_gene_25,
                    P25_gene_25,
                    P26_gene_25,
                    
                    B1_gene_25,
                    B2_gene_25,
                    B3_gene_25,
                    B4_gene_25,
                    B5_gene_25,
                    B6_gene_25,
                    
                    B21_gene_25,
                    B22_gene_25,
                    B23_gene_25,
                    B24_gene_25,
                    B25_gene_25,
                    B26_gene_25)
######

write.csv(bc_gene_25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day25_bray_curtis_formatted_by_gene.csv")
####
bc_25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day25_bray_curtis_formatted_by_gene.csv", stringsAsFactors = F)


m_bc_25 <- melt(bc_25, id=c("Gene", "Population"), measure.vars = "Frequency")
cast_bc_25 <- as.data.frame(t(dcast(m_bc_25, Population~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_bc_25) <- as.character(unlist(cast_bc_25[1,]))
braycurtis_25 <- cast_bc_25[-1,]

#now reorder the columns
colorder <- c("P1","P2","P3","P4","P5","P6",
              "P21","P22","P23","P24","P25","P26",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25","B26")
setcolorder(braycurtis_25, colorder)
#need to make all of the entries numeric
braycurtis_25 <- apply(braycurtis_25, 2, function(x) as.numeric(as.character(x)))
#transform the data so that each row is a population and the columns are the diffeent mutations
t_braycurtis_25 <- t(braycurtis_25)

####
#and, finally, going to calculate the bray curtis dissimilarity for each combination

####
#now going to calculate bray curtis and turn it into a data frame so that I can visualize the data.
braycurtis25 <- vegdist(t_braycurtis_25, method = "bray") #actually calculating the bray curtis dissimilarity for each population vs. each other population. 
bc_matrix25 <- as.matrix(braycurtis25) #make this as a matrix to be able to look at the data.
bcsim_matrix25 <- (1-bc_matrix25) #reverse so that I'm looking at the similarity
treat_list <- c(rep("AP",6), rep("GP",6), rep ("AB",6), rep("GB",6)) #list that will tell me what treatments the various populations are evolved in. AP = Arginine Plankotnic, GP = Glucose Planktonic, AB = Arginine Biofilm, GB = Glucose biofilm
rownames(bcsim_matrix25) <- c("AP","AP","AP","AP","AP","AP","GP","GP","GP","GP","GP","GP","AB","AB","AB","AB","AB","AB","GB","GB","GB","GB","GB","GB")  #treat_list #make the row names the treatments
colnames(bcsim_matrix25) <-treat_list #also make the column names the treatments
braycurtis_frame25 <- as.data.frame(as.matrix(braycurtis25)) #make it a data frame


#write out the matrix
write.csv(braycurtis_frame25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/braycurtis_day25_byGene.csv")
write.csv(bcsim_matrix25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/reverse_braycurtis_day25_byGene.csv")



































######

#now I am creating functions that I will then be able to apply to the data frame. I was smart and kept Caroline's same names, so this should work.....
######

#A function that tests for the bray curtis dissimilarity within a treatment group.
bc_within <- function(within){
  #returns mean within treatment values in a matrix
  sum <- 0
  count <- 0
  for(i in 1:(nrow(within)-1)){
    for(j in (i+1):ncol(within)){
      if(rownames(within)[i] == colnames(within)[j]){
        sum <- sum+within[i,j]
        count <- count+1
      }
    }
  }
  return(sum/count)
}

bc_within(bcsim_matrix12) #0.3180288

########

#and now a function that tests the bray curtis dissimilarity between two different treatment groups.
bc_between <- function(between){
  #returns mean between treatment values in a matrix
  sum <- 0
  count <- 0
  for(i in 1:(nrow(between)-1)){
    for(j in (i+1):ncol(between)){
      if(rownames(between)[i] != colnames(between)[j]){
        sum <- sum+between[i,j]
        count <- count+1
      }
    }
  }
  return(sum/count)
}

bc_between(bcsim_matrix12) #0.2277652




######


bc_similarity <- function(similar){
  #returns difference between mean bc of similar and non-similar (different) treatments
  sum_sim <- 0
  count_sim <- 0
  sum_diff <- 0
  count_diff <- 0
  for(i in 1:(nrow(similar)-1)){
    for(j in (i+1):ncol(similar)){
      if(rownames(similar)[i] != colnames(similar)[j]){
        if(rownames(similar)[i]==("GP")||rownames(similar)[i]=="AB"){ 
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("AB") && colnames(similar)[j]=="GP"){ 
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("GB") && colnames(similar)[j]=="AP"){ #testing carbon source in biofilm
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("AP") && colnames(similar)[j]=="GB"){ #testing carbon source differences in planktonic
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
          
        } else {
          sum_sim <- sum_sim+similar[i,j]
          count_sim <- count_sim+1
        }
      }
    }
  }
  return(sum_sim/count_sim-sum_diff/count_diff)
}

bc_similarity(bcsim_matrix12) #-0.0182
difference12 <- bc_within(bcsim_matrix12)-bc_between(bcsim_matrix12) #0.09026355


###testing weather there is more difference within treatment groups or between groups
difference_rand <- character(length(10^5))
envir_similarity <- character(length(10^5))
randomizing_matrix <-bcsim_matrix12

for(i in 1:10^5)
{
  treat_list_random <- sample(treat_list) #randomize the treatment list
  rownames(randomizing_matrix) <- treat_list_random #this will assign a random treatment to a row
  colnames(randomizing_matrix) <- treat_list_random #this will assign a random treatment to a column.
  between <- bc_between(randomizing_matrix) #now you calculate the average between treatment value in this random sample
  within <- bc_within(randomizing_matrix) #within treatment average for this randomly assigned treatment matrix
  similarity <- bc_similarity(randomizing_matrix) #how similar are the between and withing groups
  difference_rand[i] <- (between-within) #get the average
  envir_similarity[i] <- similarity
}

mean(as.numeric(difference_rand)>=difference12) #0
mean(as.numeric(envir_similarity)>=bc_similarity(bcsim_matrix12)) #0.749839


upermn <- function(x) {
  #calculates number of unique permutations of x
  #copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  n <- length(x)
  duplicates <- as.numeric(table(x))
  factorial(n) / prod(factorial(duplicates))
}

uperm <- function(x) {
  #returns a list of all the unique permutations of the elements of vector x
  # copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  u <- sort(unique(x))
  l <- length(u)
  if (l == length(x)) {
    return(do.call(rbind,permn(x)))
  }
  if (l == 1) return(x)
  result <- matrix(NA, upermn(x), length(x))
  index <- 1
  for (i in 1:l) {
    v <- x[-which(x==u[i])[1]]
    newindex <- upermn(v)
    if (table(x)[i] == 1) {
      result[index:(index+newindex-1),] <- cbind(u[i], do.call(rbind, unique(permn(v))))
    } else {
      result[index:(index+newindex-1),] <- cbind(u[i], uperm(v))
    }
    index <- index+newindex
  }
  return(result)
}

permutation_bc <- function(x){
  #
  p <- uperm(rownames(x))
  temp <- x
  diff_rand <- numeric(length=nrow(p))
  for(row in 1:nrow(p)){
    rownames(temp) <- p[row,]
    colnames(temp) <- p[row,]
    diff_rand[row] <- (bc_within(temp) - bc_between(temp))
  }
  return(diff_rand)
}





#Testing whether pairs of treatments differ significantly more between treatments than within treatments


# Loop through each treatment pair
treatments <- c("AP","GP","AB", "GB")
for(m in 1:(length(treatments)-1)){
  for(n in (m+1):length(treatments)){
    a <- treatments[m]
    b <- treatments[n]
    
    #make a matrix which only includes 2 selected treatments
    include_list <- c(a,b)
    subcomm <- bcsim_matrix12[include_list,include_list]
    subcommunity<- bcsim_matrix12[(rownames(bcsim_matrix12)== a)|(rownames(bcsim_matrix12)== b),(colnames(bcsim_matrix12)==a | colnames(bcsim_matrix12)==b)]
    
    #Calculate difference between and within treatment bray curtis
    between_ab <- bc_between(subcommunity)
    within_ab <- bc_within(subcommunity)
    #For all permutations of treatment assignemnts, calculate within-between
    permutation_diff <- permutation_bc(subcommunity)
    #Calculate significance
    significance <- mean(as.numeric(permutation_diff)>=(within_ab-between_ab))
    #Record outcome
    print(c(a,b,significance))
  }
}

#[1] "AP"                "GP"                "0.805194805194805"
#[1] "AP"                  "AB"                  "0.00432900432900433"
#[1] "AP"                "GB"                "0.350649350649351"
#[1] "GP"                  "AB"                  "0.00216450216450216"
#[1] "GP"                "GB"                "0.155844155844156"
#[1] "AB"                  "GB"                  "0.00432900432900433"








###### 
#Day 25

#now I am creating functions that I will then be able to apply to the data frame. I was smart and kept Caroline's same names, so this should work.....
######

#A function that tests for the bray curtis dissimilarity within a treatment group.
bc_within <- function(within){
  #returns mean within treatment values in a matrix
  sum <- 0
  count <- 0
  for(i in 1:(nrow(within)-1)){
    for(j in (i+1):ncol(within)){
      if(rownames(within)[i] == colnames(within)[j]){
        sum <- sum+within[i,j]
        count <- count+1
      }
    }
  }
  return(sum/count)
}

bc_within(bcsim_matrix25) #0.3451426


#and now a function that tests the bray curtis dissimilarity between two different treatment groups.
bc_between <- function(between){
  #returns mean between treatment values in a matrix
  sum <- 0
  count <- 0
  for(i in 1:(nrow(between)-1)){
    for(j in (i+1):ncol(between)){
      if(rownames(between)[i] != colnames(between)[j]){
        sum <- sum+between[i,j]
        count <- count+1
      }
    }
  }
  return(sum/count)
}

bc_between(bcsim_matrix25) #0.2039803




bc_similarity <- function(similar){
  #returns difference between mean bc of similar and non-similar (different) treatments
  sum_sim <- 0
  count_sim <- 0
  sum_diff <- 0
  count_diff <- 0
  for(i in 1:(nrow(similar)-1)){
    for(j in (i+1):ncol(similar)){
      if(rownames(similar)[i] != colnames(similar)[j]){
        if(rownames(similar)[i]==("GP")||rownames(similar)[i]=="AB"){ 
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("AB") && colnames(similar)[j]=="GP"){ 
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("GB") && colnames(similar)[j]=="AP"){ #testing carbon source in biofilm
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
        } else if(rownames(similar)[i]==("AP") && colnames(similar)[j]=="GB"){ #testing carbon source differences in planktonic
          sum_diff <- sum_diff+similar[i,j]
          count_diff <- count_diff+1
          
        } else {
          sum_sim <- sum_sim+similar[i,j]
          count_sim <- count_sim+1
        }
      }
    }
  }
  return(sum_sim/count_sim-sum_diff/count_diff)
}

bc_similarity(bcsim_matrix25) #0.05001143
difference25 <- bc_within(bcsim_matrix25)-bc_between(bcsim_matrix25) #0.1411623


###testing weather there is more difference within treatment groups or between groups
difference_rand <- character(length(10^5))
envir_similarity <- character(length(10^5))
randomizing_matrix <-bcsim_matrix25

for(i in 1:10^5)
{
  treat_list_random <- sample(treat_list) #randomize the treatment list
  rownames(randomizing_matrix) <- treat_list_random #this will assign a random treatment to a row
  colnames(randomizing_matrix) <- treat_list_random #this will assign a random treatment to a column.
  between <- bc_between(randomizing_matrix) #now you calculate the average between treatment value in this random sample
  within <- bc_within(randomizing_matrix) #within treatment average for this randomly assigned treatment matrix
  similarity <- bc_similarity(randomizing_matrix) #how similar are the between and withing groups
  difference_rand[i] <- (between-within) #get the average
  envir_similarity[i] <- similarity
}

mean(as.numeric(difference_rand)>=difference25) #0
mean(as.numeric(envir_similarity)>=bc_similarity(bcsim_matrix25)) #0.01892


upermn <- function(x) {
  #calculates number of unique permutations of x
  #copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  n <- length(x)
  duplicates <- as.numeric(table(x))
  factorial(n) / prod(factorial(duplicates))
}

uperm <- function(x) {
  #returns a list of all the unique permutations of the elements of vector x
  # copied from http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r
  u <- sort(unique(x))
  l <- length(u)
  if (l == length(x)) {
    return(do.call(rbind,permn(x)))
  }
  if (l == 1) return(x)
  result <- matrix(NA, upermn(x), length(x))
  index <- 1
  for (i in 1:l) {
    v <- x[-which(x==u[i])[1]]
    newindex <- upermn(v)
    if (table(x)[i] == 1) {
      result[index:(index+newindex-1),] <- cbind(u[i], do.call(rbind, unique(permn(v))))
    } else {
      result[index:(index+newindex-1),] <- cbind(u[i], uperm(v))
    }
    index <- index+newindex
  }
  return(result)
}

permutation_bc <- function(x){
  #
  p <- uperm(rownames(x))
  temp <- x
  diff_rand <- numeric(length=nrow(p))
  for(row in 1:nrow(p)){
    rownames(temp) <- p[row,]
    colnames(temp) <- p[row,]
    diff_rand[row] <- (bc_within(temp) - bc_between(temp))
  }
  return(diff_rand)
}





#Testing whether pairs of treatments differ significantly more between treatments than within treatments


# Loop through each treatment pair
treatments <- c("AP","GP","AB", "GB")
for(m in 1:(length(treatments)-1)){
  for(n in (m+1):length(treatments)){
    a <- treatments[m]
    b <- treatments[n]
    
    #make a matrix which only includes 2 selected treatments
    include_list <- c(a,b)
    subcomm <- bcsim_matrix25[include_list,include_list]
    subcommunity<- bcsim_matrix25[(rownames(bcsim_matrix25)== a)|(rownames(bcsim_matrix25)== b),(colnames(bcsim_matrix25)==a | colnames(bcsim_matrix25)==b)]
    
    #Calculate difference between and within treatment bray curtis
    between_ab <- bc_between(subcommunity)
    within_ab <- bc_within(subcommunity)
    #For all permutations of treatment assignemnts, calculate within-between
    permutation_diff <- permutation_bc(subcommunity)
    #Calculate significance
    significance <- mean(as.numeric(permutation_diff)>=(within_ab-between_ab))
    #Record outcome
    print(c(a,b,significance))
  }
}

#[1] "AP"                 "GP"                 "0.0541125541125541"
#[1] "AP"                  "AB"                  "0.00649350649350649"
#[1] "AP"                  "GB"                  "0.00649350649350649"
#[1] "GP"                  "AB"                  "0.00216450216450216"
#[1] "GP"                 "GB"                 "0.0151515151515152"
#[1] "AB"                  "GB"                  "0.00432900432900433"

#######




########
#make tables of parallelism at the gene level for DFE examination. 

######
#using tables made in bray curtis. 
day12 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day12_bray_curtis_formatted_by_gene.csv", stringsAsFactors = F)

day25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day25_bray_curtis_formatted_by_gene.csv", stringsAsFactors = F)





#now I need to melt the data frame and cast it again
m_day12 <- melt(day12, id=c("Gene", "Population"), measure.vars = "Frequency")
cast_day12 <- as.data.frame(t(dcast(m_day12, Population~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_day12) <- as.character(unlist(cast_day12[1,]))
day_12 <- cast_day12[-1,]

#now reorder the columns
colorder <- c("P1","P2","P3","P4","P5","P6",
              "P21","P22","P23","P24","P25","P26",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25","B26")
setcolorder(day_12, colorder)
#print the tables 
write.csv(day_12, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day12_gene_parallelism.csv")




m_day25 <- melt(day25, id=c("Gene", "Population"), measure.vars = "Frequency")
cast_day25 <- as.data.frame(t(dcast(m_day25, Population~Gene, mean, value.var = "value", fill = 0)))
colnames(cast_day25) <- as.character(unlist(cast_day25[1,]))
day_25 <- cast_day25[-1,]

#now reorder the columns
colorder <- c("P1","P2","P3","P4","P5","P6",
              "P21","P22","P23","P24","P25","P26",
              "B1","B2","B3","B4","B5","B6",
              "B21","B22","B23","B24","B25","B26")
setcolorder(day_25, colorder)
#print the tables 
write.csv(day_25, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Bray_curtis/Day25_gene_parallelism.csv")







#######
#example command for running lolipop on my computer 

#######
# lolipop lineage -i /Users/katrina/Desktop/month_long_exp/Sequencing_analysis/final_timeseries/final/Lolipop/b24_final.xlsx -o /Users/katrina/Desktop/month_long_exp/Sequencing_analysis/Muller_output/B24

#####

#need to do the fishers exact test for all of the cases of parallelism
######
a <-  matrix (c(145, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(100, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(98, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(86, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(756, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(170, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1239, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1509, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(2451, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1464, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(2172, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(543, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(95, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1713, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(2277, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(4254, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(382, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(2493, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1140, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(585, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(345, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(169, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(804, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(981, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(2925, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(639, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(73, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(798, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1500, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(144, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(127, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(138, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(927, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(172, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1278, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1443, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(165, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(102, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(2058, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1911, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(2778, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(101, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1044, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(945, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(336, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(275, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(2229, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1239, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(855, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(1617, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(71, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(73, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(936, 2, 6537648, 624), nrow = 2)

a <-  matrix (c(1143, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(1371, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(129, 2, 6537648, 624), nrow = 2)
a <-  matrix (c(121, 2, 6537648, 624), nrow = 2)

b <-  matrix (c(97, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(1218, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(1206, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(1461, 3, 6537648, 624), nrow = 2)

b <-  matrix (c(729, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(660, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(2568, 3, 6537648, 624), nrow = 2)

b <-  matrix (c(609, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(2718, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(2226, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(99, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(90, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(2721, 3, 6537648, 624), nrow = 2)

b <-  matrix (c(244, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(1254, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(615, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(4254, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(549, 3, 6537648, 624), nrow = 2)

b <-  matrix (c(102, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(5430, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(4074, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(346, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(3018, 3, 6537648, 624), nrow = 2)
b <-  matrix (c(154, 3, 6537648, 624), nrow = 2)

c <-  matrix (c(217, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(960, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(167, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(2040, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(975, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(405, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(1470, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(189, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(474, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(164, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(1437, 4, 6537648, 624), nrow = 2)
c <-  matrix (c(516, 4, 6537648, 624), nrow = 2)

d <-  matrix (c(1287, 5, 6537648, 624), nrow = 2)
d <-  matrix (c(915, 5, 6537648, 624), nrow = 2)
d <-  matrix (c(555, 5, 6537648, 624), nrow = 2)
d <-  matrix (c(1104, 5, 6537648, 624), nrow = 2)
d <-  matrix (c(616, 5, 6537648, 624), nrow = 2)

e <-  matrix (c(921, 6, 6537648, 624), nrow = 2)
e <-  matrix (c(2670, 6, 6537648, 624), nrow = 2)
e <-  matrix (c(1392, 6, 6537648, 624), nrow = 2)
e <-  matrix (c(1443, 6, 6537648, 624), nrow = 2)

f <-  matrix (c(645, 7, 6537648, 624), nrow = 2)
f <-  matrix (c(48, 7, 6537648, 624), nrow = 2)
f <-  matrix (c(960, 7, 6537648, 624), nrow = 2)

g <-  matrix (c(1161, 10, 6537648, 624), nrow = 2)
g <-  matrix (c(2361, 10, 6537648, 624), nrow = 2)



fisher.test(one)

#####


#read in all of the final mutation tables

B1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B1.csv", stringsAsFactors = F)
B2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B2.csv", stringsAsFactors = F)
B3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B3.csv", stringsAsFactors = F)
B4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B4.csv", stringsAsFactors = F)
B5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B5.csv", stringsAsFactors = F)
B6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B6.csv", stringsAsFactors = F)

B21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B21.csv", stringsAsFactors = F)
B22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B22.csv", stringsAsFactors = F)
B23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B23.csv", stringsAsFactors = F)
B24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B24.csv", stringsAsFactors = F)
B25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B25.csv", stringsAsFactors = F)
B26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/B26.csv", stringsAsFactors = F)


P1 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P1.csv", stringsAsFactors = F)
P2 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P2.csv", stringsAsFactors = F)
P3 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P3.csv", stringsAsFactors = F)
P4 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P4.csv", stringsAsFactors = F)
P5 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P5.csv", stringsAsFactors = F)
P6 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P6.csv", stringsAsFactors = F)

P21 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P21.csv", stringsAsFactors = F)
P22 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P22.csv", stringsAsFactors = F)
P23 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P23.csv", stringsAsFactors = F)
P24 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P24.csv", stringsAsFactors = F)
P25 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P25.csv", stringsAsFactors = F)
P26 <- read.csv("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/P26.csv", stringsAsFactors = F)



#figuring out the mutational breakdown for all of the populations

library(reshape2)

setwd("/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/mutation_breakdown")
Mutations_analysis <- function(Mutation_Data, AminoAcid="description", Bases = "Mutation") {
  #Mutation_Data is the input CSV file that is a breseq output with rows containing different mutations and columns containing the various characteristics for those mutations.
  #AminoAcid is the column name that holds the breseq information on amino acid changes. this column will look like "*342C(TGA>TGC)" or say coding, pseudogene, etc. The default that will be looked for is "Description". This is case sensitive!
  #Bases is the column name that holds the breseq information for the nucleotide mutations. These are things like A>C, or T>G in the breseq output. This is case sensitive!!!
  
  ##############
  #first i am going to deal with the nucleotide level - the mutation column. This uses the Mutation_Data that you put in and grabs the information in the column that you specified under Bases. It looks for the > symbol, because that is how breseq separates the two bases, and splits the data. It then creates a new data set called Nucleotides that has 2 columns containing the original base (from the ancestor) and the mutated base.
  Nucleotides <- colsplit(Mutation_Data[,Bases], ">", names = c("original", "mutant"))
  #View(Nucleotides)
  #I want to calculate the total number of mutations present in the sample.
  allmutations <- nrow(Nucleotides)
  #I want to determine the number of mutations that are not just substituting one base from another. These are indels, because this is how breseq represents them.
  indel <- sum(grepl("[^ACGT]", Nucleotides$original))
  
  
  #I need to find all of the different combinations for base replacements. To do this I am going to find the index for each base, and then I will look at that index in the second column and see what the new base is.
  C <- grep("C", Nucleotides$original) #find placeswhere the original was C
  CT <- sum(ifelse(Nucleotides$mutant[C]=="T",1,0)) # find when there was a C to T transition
  CA <- sum(ifelse(Nucleotides$mutant[C]=="A",1,0)) # find when there was a C to A transition
  CG <- sum(ifelse(Nucleotides$mutant[C]=="G",1,0)) # find when there was a C to G transition
  
  Ts <- grep("T", Nucleotides$original) #find when the original was a T
  TC <- sum(ifelse(Nucleotides$mutant[Ts]=="C",1,0)) #find when there were T to C transitions
  TG <- sum(ifelse(Nucleotides$mutant[Ts]=="G",1,0)) #find when there were T to G transitions
  TA <- sum(ifelse(Nucleotides$mutant[Ts]=="A",1,0)) #find when there were T to A transitions
  
  G <- grep("G", Nucleotides$original) #find placeswhere the original was G
  GA <- sum(ifelse(Nucleotides$mutant[G]=="A",1,0)) # find when there was a G to A transition
  GT <- sum(ifelse(Nucleotides$mutant[G]=="T",1,0)) # find when there was a G to T transition
  GC <- sum(ifelse(Nucleotides$mutant[G]=="C",1,0)) # find when there was a G to C transition
  
  A <- grep("A", Nucleotides$original) #find placeswhere the original was A
  AG <- sum(ifelse(Nucleotides$mutant[A]=="G",1,0)) # find when there was a A to G transition
  AC <- sum(ifelse(Nucleotides$mutant[A]=="C",1,0)) # find when there was a A to C transition
  AT <- sum(ifelse(Nucleotides$mutant[A]=="T",1,0)) # find when there was a A to T transition
  
  # Now that I have the numbers of all of the possible base changes, I can look for the
  transitions <- sum(c(CT,TC,GA,AG)) #there are 4 options for transitions. C>T, T>C, G>A, A>G. this adds up all of those changes
  transversions <- AT+AC+GC+GT+CA+CG+TG+TA # need to do to check that the sums of the transition categories actually do add up to the number of transitions that there should be (assuming transitions and indel numbers are correct) when I turn this into a function I need to stop if transversions != trans -- should be fine but just an extra error checking step.
  
  ###############
  
  ### now at the Amino acid level
  #have to get the amino acid column that the user specifies out of the input data
  Protein <- colsplit(Mutation_Data[,AminoAcid], "\\(", names = c("AA", "DNA"))
  
  #there are a few options that you can get for this one and I can't just split the column. I need to look for all of them in the strings FOr this I need to use regular expressions. I will have to use gregexpr which returns a list of positions and then I will have to find the length of that to determine the number of them. The regular expressios that I will use for each option are as follows.
  
  # reg expressions to use "coding", "intergenic", "pseudogene",
  #"[A-Z][0-9]+[A-Z]" #this looks for a base, followed by a number that can be any size 1 or more, followed by a base.
  #"\\*[0-9]+[A-Z]" #this looks for an asterisk, followed by at least 1 number, followed by a base
  #"[A-Z][0-9]*\\*" #this looks for a base, followed by at least 1 number, followed by an asterisk
  
  coding = sum(grepl("coding",Protein$AA)) #Breseq's coding region
  intergenic = sum(grepl("intergenic",Protein$AA)) #intergenic breseq designation
  pseudogene = sum(grepl("pseudogene",Protein$AA)) #pseudogene breseq designation
  prematurestop = sum(lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]*\\*", Protein$AA)))) #these are when you have a coding amino acid that gets mutated into a stop codon
  elongating = sum(lengths(regmatches(Protein$AA,gregexpr("\\*[0-9]*[A-Z]", Protein$AA)))) #these are stop codons that get mutated into coding amino acids that elongate your protein.
  aamutation = sum(lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]+[A-Z]", Protein$AA)))) # these are all of the mutations that dont fit other categories. so these mutations change one amino acid to another with no other breseq designation.
  
  #I now need to determine if the amino acid mutations category are synonymous or nonsynonymous. The above just determines the number of leter to leter strings exist. Now I need to actually look at that subset of mutations and determine if
  aas <- lengths(regmatches(Protein$AA,gregexpr("[A-Z][0-9]+[A-Z]", Protein$AA))) #this returns a list of logicals that tell you if there is an aamutation at that spot (1) or not (0)
  #aas
  aminos <- as.matrix(Protein$AA[aas==1]) #this find the amino acid changes at the previously found indexes, so these are the actual identities of the aamutation mutations.
  aminos2 <- colsplit(aminos, "[0-9]+", names = c("first", "last")) # splitting the breseq annotation. I am taking out the number, because the position is irrelevent, and separating the two letters into two columns containing the first, original aa, and the last, or mutated, aa.
  
  synonymous <- sum(ifelse(as.character(aminos2$first) == as.character(aminos2$last),1,0)) #if the letters are the same before and after the mutation it is synonymous
  nonsynonymous <- sum(ifelse(as.character(aminos2$first) == as.character(aminos2$last),0,1)) #if the letters are different then it is nonsynonymous
  dnds <- (nonsynonymous/synonymous)/2.60 
  #2.60 is the pseudomonas genomic dN/dS ratio. So this reported value will be gorrected for the basal genomic distribution.
  
  # I am now making a table of all of the mutational types that I can print out later. The other thing that I would like to do is give it a name specific to the data set that you put in, but I don't know how to do that. For the moment it will just always return the same named file each time, so you have to change the name before you use the code again. 
  table<- matrix(c("Mutations: ", allmutations,
                   "Nucleotide level mutations", "",
                   "Indels: ", indel,
                   "Transitions: ",transitions,
                   "C>T: ", CT,
                   "T>C: ", TC,
                   "A>G: ", AG,
                   "G>A: ", GA,
                   "Transversions: ", transversions,
                   "A>T: ", AT,
                   "A>C: ", AC,
                   "G>C: ", GC,
                   "G>T: ", GT,
                   "C>A: ", CA,
                   "C>G: ", CG,
                   "T>G: ", TG,
                   "T>A: ", TA,
                   "Amino acid level mutations", "",
                   "Coding: ", coding,
                   "Intergenic: ", intergenic,
                   "Pseudogene: ", pseudogene,
                   "Premature stop: ", prematurestop,
                   "Elongating: ", elongating,
                   "Synonymous: ", synonymous,
                   "Non Synonymous: ", nonsynonymous,
                   "dN/dS: ", dnds), ncol = 2, byrow=T)
  
  #write out the file. Adding col.names won't actually give the columns names but it will keep things in two different columns instead of compiling it all.
  write.csv(table, file = "Mutations_table.csv", col.names = T)
  
}

planktonic_arginine <- rbind(P1, P2, P3, P4, P5, P6)
planktonic_glucose <- rbind(P21, P22, P23, P24, P25, P26)
biofilm_arginine <- rbind(B1, B2, B3, B4, B5, B6)
biofilm_glucose <- rbind(B21, B22, B23, B24, B25, B26)

planktonic <- rbind(planktonic_arginine, planktonic_glucose)
biofilm <- rbind(biofilm_arginine, biofilm_glucose)
glucose <- rbind(planktonic_glucose, biofilm_glucose)
arginine <- rbind(planktonic_arginine, biofilm_arginine)

all <- rbind(glucose, arginine)

Mutations_analysis(all, AminoAcid = "Annotation", Bases = "Mutation")


#want to write files for the final sets of mutations
write.csv(all, "/Users/katrina/Desktop/month_long_exp/Sequencing_analysis/200520_finalmutations/all_final_mutations_210213.csv")

#putting together a table of the different breakdowns by environment. 
#need to read in the data

