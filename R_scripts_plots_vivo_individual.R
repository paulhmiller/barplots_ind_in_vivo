# Figures for plerixafor experiments
# In vivo data, individual samples
# Barplots and 3-way comparison statistics
# Paul Miller, paulhmiller@gmail.com

library(reshape2)
library(ggplot2)
library(grid)
library(plyr)

source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")
#setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures')
#codes<-read.csv("plerixafor_individual_sample_codes.csv",header=T,check.names=F) #,row.names=1)
setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/in vivo individual')

dat<-read.csv("plerixafor_vivo_individual5.csv",header=T,check.names=F,colClasses=c("Experiment"="factor")) #,row.names=1)
#dat<-merge(codes,vivo,"ID",all=T,na.rm=T)

# Take assay names
assays<-names(dat[9:19])
# Change low values in PB leukocyte columns to detection threshold:
dat[,9:12][dat[,9:12] < 0.25] <- 0.25
# Change low values in platelet column to detection threshold:
dat[,13][dat[,13] < 25] <- 25
# Change low values in BM columns to detection threshold:
dat[,14:19][dat[,14:19] < 0.01] <- 0.01

# Anonymize donor names
levels(dat$Donor) <- 1:10

# Add fold increase columns
	
# Create single value per donor. There is little inter-mouse variation for each sample, so just use the triplicate geometric means. 
w3 <- NULL
BLm <- 10^(colMeans(log10(dat[dat$Timepoint=="BL" & dat$Week==3, c(9:19)]), na.rm=T))
for(i in 1:10){
	#take geometric means for donor i BL:
	BL <- dat[dat$Donor==i & dat$Timepoint=="BL" & dat$Week==3, c(1:6,8,9:19)]
	BL[,8:18] <- log10(BL[8:18])
	BLl <- BL[1,] #just to take the row names
	foldBL <- BL[1,] #just to take the row names
	BLl[,8:18] <- colMeans(BL[,8:18],na.rm=T) 
	BL <- BLl
	BL[,8:18] <- 10^(BL[8:18]) #convert back to arithmetic
	#take geometric means for donor i D-1:
	Dm1 <- dat[dat$Donor==i & dat$Timepoint=="-24hrs" & dat$Week==3, c(1:6,8,9:19)]
	Dm1[,8:18] <- log10(Dm1[8:18])
	Dm1l <- Dm1[1,] #just to take the row names
	Dm1l[,8:18] <- colMeans(Dm1[,8:18],na.rm=T) 
	Dm1 <- Dm1l
	Dm1[,8:18] <- 10^(Dm1[8:18])
	foldDm1 <- Dm1
	#take geometric means for donor i D0:
	D0 <- dat[dat$Donor==i & dat$Timepoint=="4hrs" & dat$Week==3, c(1:6,8,9:19)]
	D0[,8:18] <- log10(D0[8:18])
	D0l <- D0[1,] #just to take the row names
	D0l[,8:18] <- colMeans(D0[,8:18],na.rm=T) 
	D0 <- D0l
	D0[,8:18] <- 10^(D0[8:18])
	foldD0 <- D0
	# calculate fold increase of D-1 over BL:
	foldDm1[,8:18] <- Dm1[,8:18] / BL[,8:18]
	# calculate fold increase of D0 over BL:
	foldD0[,8:18] <- D0[,8:18] / BL[,8:18]
	# calculate fold change of donor i BL over BL mean:
	foldBL[,8:18] <- BL[,8:18] / BLm
	w3 <- rbind(w3,foldDm1,foldD0,foldBL)
	}


w6 <- NULL
BLm <- 10^(colMeans(log10(dat[dat$Timepoint=="BL" & dat$Week==6, c(9:19)]), na.rm=T))
for(i in 1:10){
  #take geometric means for donor i BL:
  BL <- dat[dat$Donor==i & dat$Timepoint=="BL" & dat$Week==6, c(1:6,8,9:19)]
  BL[,8:18] <- log10(BL[8:18])
  BLl <- BL[1,] #just to take the row names
  foldBL <- BL[1,] #just to take the row names
  BLl[,8:18] <- colMeans(BL[,8:18],na.rm=T) 
  BL <- BLl
  BL[,8:18] <- 10^(BL[8:18]) #convert back to arithmetic
  #take geometric means for donor i D-1:
  Dm1 <- dat[dat$Donor==i & dat$Timepoint=="-24hrs" & dat$Week==6, c(1:6,8,9:19)]
  Dm1[,8:18] <- log10(Dm1[8:18])
  Dm1l <- Dm1[1,] #just to take the row names
  Dm1l[,8:18] <- colMeans(Dm1[,8:18],na.rm=T) 
  Dm1 <- Dm1l
  Dm1[,8:18] <- 10^(Dm1[8:18])
  foldDm1 <- Dm1
  #take geometric means for donor i D0:
  D0 <- dat[dat$Donor==i & dat$Timepoint=="4hrs" & dat$Week==6, c(1:6,8,9:19)]
  D0[,8:18] <- log10(D0[8:18])
  D0l <- D0[1,] #just to take the row names
  D0l[,8:18] <- colMeans(D0[,8:18],na.rm=T) 
  D0 <- D0l
  D0[,8:18] <- 10^(D0[8:18])
  foldD0 <- D0
  # calculate fold increase of D-1 over BL:
  foldDm1[,8:18] <- Dm1[,8:18] / BL[,8:18]
  # calculate fold increase of D0 over BL:
  foldD0[,8:18] <- D0[,8:18] / BL[,8:18]
  # calculate fold change of donor i BL over BL mean:
  foldBL[,8:18] <- BL[,8:18] / BLm
  w6 <- rbind(w6,foldDm1,foldD0,foldBL)
}

	
w12 <- NULL
BLm <- 10^(colMeans(log10(dat[dat$Timepoint=="BL" & dat$Week==12, c(9:19)]), na.rm=T))
for(i in 1:10){
  #take geometric means for donor i BL:
  BL <- dat[dat$Donor==i & dat$Timepoint=="BL" & dat$Week==12, c(1:6,8,9:19)]
  BL[,8:18] <- log10(BL[8:18])
  BLl <- BL[1,] #just to take the row names
  foldBL <- BL[1,] #just to take the row names
  BLl[,8:18] <- colMeans(BL[,8:18],na.rm=T) 
  BL <- BLl
  BL[,8:18] <- 10^(BL[8:18]) #convert back to arithmetic
  #take geometric means for donor i D-1:
  Dm1 <- dat[dat$Donor==i & dat$Timepoint=="-24hrs" & dat$Week==12, c(1:6,8,9:19)]
  Dm1[,8:18] <- log10(Dm1[8:18])
  Dm1l <- Dm1[1,] #just to take the row names
  Dm1l[,8:18] <- colMeans(Dm1[,8:18],na.rm=T) 
  Dm1 <- Dm1l
  Dm1[,8:18] <- 10^(Dm1[8:18])
  foldDm1 <- Dm1
  #take geometric means for donor i D0:
  D0 <- dat[dat$Donor==i & dat$Timepoint=="4hrs" & dat$Week==12, c(1:6,8,9:19)]
  D0[,8:18] <- log10(D0[8:18])
  D0l <- D0[1,] #just to take the row names
  D0l[,8:18] <- colMeans(D0[,8:18],na.rm=T) 
  D0 <- D0l
  D0[,8:18] <- 10^(D0[8:18])
  foldD0 <- D0
  # calculate fold increase of D-1 over BL:
  foldDm1[,8:18] <- Dm1[,8:18] / BL[,8:18]
  # calculate fold increase of D0 over BL:
  foldD0[,8:18] <- D0[,8:18] / BL[,8:18]
  # calculate fold change of donor i BL over BL mean:
  foldBL[,8:18] <- BL[,8:18] / BLm
  w12 <- rbind(w12,foldDm1,foldD0,foldBL)
}
dat <- rbind(w3,w6,w12)
	
# Remove NA columns (created by above for loop running d-1 fold changes on group A, which doesn't have d-1):
dat<-dat[!is.na(dat$Experiment),]

#Finished fold increase additions


#Log10 transform
dat[,8:18]<-log10(dat[8:18])

# Convert to tall format 
dat<-melt(dat,id.vars=c("Experiment","Donor","ID","Sample","Timepoint","Group","Week"),variable.name="assay",value.name="value",na.rm=TRUE) 

# Rename group factors (had to separate because couldn't rename factor based on 
# two columns. Consequence is that BL values now removed, so don't need to do that 
# separately)
P <- dat[dat$Group=="A" & dat$Timepoint=="4hrs",]
GP <- dat[dat$Group=="B" & dat$Timepoint=="4hrs",]
G <- dat[dat$Group=="B" & dat$Timepoint=="-24hrs",]
levels(G$Group)[levels(G$Group)=="B"] <- "G"
dat <- rbind(P,GP,G)
levels(dat$Group)[levels(dat$Group)=="A"] <- "P"   
levels(dat$Group)[levels(dat$Group)=="B"] <- "G+P"

# Re-order factor levels:
dat$Group<-factor(dat$Group,levels=c("G","P","G+P")) 


########## Fold increase - Means with SEM#############

PB <- dat

#**** N.B.: add 2 log to all values so plots look OK. MUST change y-axis breaks. 

PB[, 9] <- PB[, 9] +2

# Note: donors merged in following:
PB <- summarySE(PB, measurevar="value", groupvars=c("Timepoint","Group","Week","assay"),na.rm=TRUE) # summarySE provides std, SEM, and (default 95%) CI. #measurevar is the x-axis. "Donor",

# Exlude baseline
PB <- PB[PB$Timepoint!="BL",]



# Plot by week

PB <- PB[PB$assay=="CD45" | PB$assay=="CD33" | PB$assay=="CD19" | PB$assay=="Platelets",]
w3 <- PB[PB$Week==3,]
w6 <- PB[PB$Week==6,]
w12 <- PB[PB$Week==12,]


plot<-ggplot(w3,aes(x=as.factor(assay), y=value, fill=Group))+ 
  annotation_logticks(sides = "l", size=0.3) +
  geom_bar(stat="identity", position=position_dodge())+ # colour="black",
  scale_fill_manual(values=p1) + 
  scale_color_manual(values=p1) +	
  geom_errorbar(aes(ymin=value-se, ymax=value+se, color=Group), lty=1, width=0.5, size=0.75, position=position_dodge(width=0.9)) +
  ggtitle("Week 3") + 
  guides(fill=F,color=F)+
  #	scale_y_continuous(breaks=c(-1, 0, 1, 2, 3, 4),labels = c(0.1, 1, 10, 100, 1000, 10000)) + 
  scale_y_continuous(breaks=c(1, 2, 3, 4, 5, 6),labels = c(expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4))) + 
  coord_cartesian(ylim=c(2.105,4.3)) +
  ylab("fold-increase") + 
  scale_x_discrete("")+
  themePM1() + 
  theme(plot.margin = unit(c(0+0.2, 0.2, 0, 0.1), "cm"))
plot
ggsave(filename="fold w3 vivo ind.pdf",width=5.75,height=5.75, units="cm")

plot %+% w6 +	ggtitle("Week 6")
ggsave(filename="fold w6 vivo ind.pdf",width=5.75,height=5.75, units="cm")

plot %+% w12 + ggtitle("Week 12")
ggsave(filename="fold w12 vivo ind.pdf",width=5.75,height=5.75, units="cm")




## Plot by lineage

CD45 <- PB[PB$assay=="CD45",]
CD33 <- PB[PB$assay=="CD33",]
CD19 <- PB[PB$assay=="CD19",]
CD3 <- PB[PB$assay=="CD3",]
PLT <- PB[PB$assay=="Platelets",]


plot<-ggplot(CD45,aes(x=as.factor(Week), y=value, fill=Group))+ 
  annotation_logticks(sides = "l", size=0.3) +
  geom_bar(stat="identity", position=position_dodge())+ # colour="black",
	scale_fill_manual(values=p1) + #use F if also plotting BL
	scale_color_manual(values=p1) +	
	geom_errorbar(aes(ymin=value-se, ymax=value+se, color=Group), lty=1, width=0.5, size=0.75, position=position_dodge(width=0.9)) +
	ggtitle("CD45") + 
	guides(fill=F,color=F)+
#	scale_y_continuous(breaks=c(-1, 0, 1, 2, 3, 4),labels = c(0.1, 1, 10, 100, 1000, 10000)) + 
	scale_y_continuous(breaks=c(1, 2, 3, 4, 5, 6),labels = c(expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4))) + 
#ylab(expression(paste("CD45 x 10"^"3","/mL mouse blood"))) +
	coord_cartesian(ylim=c(2.105,4.3)) +
	ylab("fold-increase") +
	scale_x_discrete("Week")+
	themePM1() +
  theme(plot.margin = unit(c(0+0.2, 0.2, 0, 0.1), "cm"))
plot
ggsave(filename="fold CD45 vivo ind.pdf",width=5.75,height=5.75, units="cm")

plot %+% CD33 +	ggtitle("CD33/15") #+ coord_cartesian(ylim=c(0,1.69897)) 
ggsave(filename="fold GM vivo ind.pdf",width=5.75,height=5.75, units="cm")

plot %+% CD19 + ggtitle("CD19") + coord_cartesian(ylim=c(2.09,3.91)) +  theme(plot.margin = unit(c(0+0.2, 0.2, 0, 0.1-0.12), "cm"))
ggsave(filename="fold CD19 vivo ind.pdf",width=5.75,height=5.75, units="cm")

plot %+% PLT + ggtitle("Platelets") +  coord_cartesian(ylim=c(2.09,3.91)) 
ggsave(filename="fold PLT vivo ind.pdf",width=5.75,height=5.75, units="cm")







#### STATS ####

# t.test (assumes equal variance) - FDR corrected for multiple testing
# Comparison between arms:
tmp1<-dat[dat$Timepoint=="4hrs",]
tmp2<-dat[dat$Timepoint=="-24hrs",]
PB <- rbind(tmp1,tmp2)
   
sink('ttest FDR vivo.txt')
for(i in 1:5){ 
  tmp <- PB[PB$assay==assays[[i]], c(6:9)]
  cat("\n")
  cat("\n")
  cat("\n")
  print(assays[[i]])
  for(f in c(3,6,12)){
    cat("\n")
    print(paste("Week",f))
    tmp2 <- tmp[tmp$Week==f, c(1,4)]
    out <- pairwise.t.test(tmp2$value, tmp2$Group, p.adjust.method="fdr")
    print(out)
  }
}
sink()


