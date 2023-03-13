# This script reads the AddMyPet database
# via Matlab files .mat
# Author: Jaap van der Meer
# Last update: March 6, 2023

# Set working directory and load required packages
# setwd("R/AmPdata")
rm(list=ls())
require("R.matlab")

# Read workspace 
load("M:\\My Documents\\R\\AmPdata\\allAmP.RData")

# If workspace was indeed available, then skip the next reading data sections 
# and go directly to function definitions

# Read allLabel.mat and allUnits.mat
# These files contain all labels and units of DEB variables and constants
# Presently (17-2-2023) they contain 371 character string and 1 list with ecocodes
# Three labels have an empty unit
allLabel <- readMat("allLabel.mat")$allLabel
k <- length(allLabel)
rownames(allLabel)
tmp <- rep(NA,k)
for (i in 1:k) tmp[i] <- length(allLabel[[i]])
table(tmp)
rm(k,tmp)

allUnits <- readMat("allUnits.mat")$allUnits
k <- length(allUnits)
rownames(allUnits)
tmp <- rep(NA,k)
for (i in 1:k) tmp[i] <- length(allUnits[[i]])
table(tmp)
rm(k,tmp)

# Read allStat.mat
# Presently (17-2-2023) the number of species equals 3983
# The number of statistics vary between 239 and 280 per species
allStat <- readMat("allStat.mat")$allStat
k <- length(allStat)
rownames(allStat)
tmp <- rep(NA,k)
for (i in 1:k) tmp[i] <- length(allStat[[i]])
table(tmp)
rm(k,tmp)

# Read popStat.mat
# Presently (17-2-2023) the number of species equals 3983
# The number of statistics vary between 6 and 7 per species
# T, c.T, f0, ff, f1, model, par
# Four of them have sub-statistics at various levels
# The statistic par is also given in allStat.mat
popStat <- readMat("popStat.mat")$popStat
n.species <- k <- length(popStat)
rownames(popStat)
tmp <- rep(NA,k)
for (i in 1:n.species) tmp[i] <- length(popStat[[i]])
table(tmp)
rm(tmp)
# End reading


# Begin function definitions
# Two functions to show content of allStat and popStat
# Arguments are File, which can either be allStat or popStat
# and Species which can either be a number (ranging from 1 to 
# the number of species in the collection) or a species name, e.g.
# Homo.sapiens (note the dot as separator)
# Values are lists of all statistics and of those statistics presented as a list 
 
Show.Statistics <- function(File=allStat,Species="Homo.sapiens"){
	n.species <- dim(File)[1]
	if (mode(Species)=="character") {
		Index <- rownames(File)==Species
		SpeciesNumber <- (1:n.species)[Index]
		if (sum(Index)==0) print("Wrong species names; use a . as seperator")
	} else SpeciesNumber <- Species
	tmp1 <- File[[SpeciesNumber]]
	return(rownames(tmp1))
}

Show.Lists <- function(File=allStat,Species="Homo.sapiens"){
	n.species <- dim(File)[1]
	if (mode(Species)=="character") {
		Index <- rownames(File)==Species
		SpeciesNumber <- (1:n.species)[Index]
		if (sum(Index)==0) print("Wrong species names; use a . as seperator")
	} else SpeciesNumber <- Species
	tmp1 <- File[[SpeciesNumber]]
	n.statistics <- length(tmp1)
	for (i in 1:n.statistics){
			if (mode(tmp1[[i]])=="list") {
			print(paste("Element",i,rownames(tmp1)[[i]],"is a list of length",length(tmp1[[i]]),sep=" "))
		}
	}
}

# Two functions to draw specific(sub)statistic(s) for a specific species

Draw.Statistic <- function(File=allStat,Species="Homo.sapiens",StatisticName="L.p"){
	n.species <- dim(File)[1]
	if (mode(Species)=="character") {
		Index <- rownames(File)==Species
		SpeciesNumber <- (1:n.species)[Index]
	} else SpeciesNumber <- Species
	tmp <- File[[SpeciesNumber]]
	Index <- rownames(tmp)==StatisticName
	n.statistics <- length(tmp)
	out <- NA
	if (sum(Index)==1) {
		StatisticNumber <- (1:n.statistics)[Index]
		out <- tmp[[StatisticNumber]]
		names(out) <- StatisticName
	} else print("StatisticName does not exist")
	if (StatisticName=="date.acc"|StatisticName=="date.subm") {
		out <- as.Date(paste(out[1],out[2],out[3],sep="-"),"%Y-%m-%d")
	}
	if (mode(out)=="list") {
		print("Statistic is a list")
		names(out) <- rownames(out)
	}
	return(out)
}

Draw.Statistic.from.List <- function(File=popStat,Species="Homo.sapiens",ListName="f0",
	SubListName="thin0",Sex="f",StatisticName="Y.CX"){
	out <- NA
	tmp1 <- Draw.Statistic(File,Species,ListName)
	if (ListName=="author"){
		out <- paste(unlist(tmp1),collapse=" and ")
		names(out) <- ListName
	}
	if (ListName=="ecoCode"){
		Index <- rownames(tmp1)==SubListName
		n1 <- length(tmp1)
		tmp2 <- NA
		if (sum(Index)==1) {
			Number1 <- (1:n1)[Index]
			tmp2 <- tmp1[[Number1]]
		} else print("Data for the specific ecocode do not exist")
		out <- paste(unlist(tmp2),collapse="-")
		names(out) <- paste(ListName,SubListName,sep="-")
	}
	if (ListName=="f0"|ListName=="ff"|ListName=="f1"){
		Index <- rownames(tmp1)==SubListName # either thin0 or thin1
		n1 <- length(tmp1)
		tmp2 <- NA
		if (sum(Index)==1) {
			Number1 <- (1:n1)[Index]
			tmp2 <- tmp1[[Number1]]
		} else print("Data for the specific thinning rule do not exist")
		Index <- rownames(tmp2)==Sex # either f or m
		n2 <- length(tmp2)
		tmp3 <- NA
		if (sum(Index)==1) {
			Number2 <- (1:n2)[Index]
			tmp3 <- tmp2[[Number2]]
		} else print("Data for the specific sex do not exist")
		Index <- rownames(tmp3)==StatisticName
		n3 <- length(tmp3)
		out <- NA
		if (sum(Index)==1) {
			Number3 <- (1:n3)[Index]
			out <- tmp3[[Number3]]
		} else print("Datum for the specific statistic does not exist")
		names(out) <- paste(ListName,SubListName,SexName,StatisticName,sep="-")
	}
	return(out)
}
# End function definitions

# Begin examples
# List of eight ecoCodes (without units)
allLabel[[(1:length(allLabel))[rownames(allLabel)=="ecoCode"]]]

# Show statistics
Show.Statistics()

# Draw a particular (sub)statistic for a particular species
# Arguments are file, species (name or number), statistic (and more) 
n.species <- length(allStat)
Draw.Statistic(allStat,n.species,"L.p")
Species <- "Gorilla.gorilla"
Draw.Statistic(popStat,Species,"f0")
Draw.Statistic(popStat,Species,"f0")[[2]]
Draw.Statistic(popStat,Species,"f0")[[2]][[1]]
Draw.Statistic.from.List(popStat,"Gorilla.gorilla","f0","thin0","f","mu.TX")
Draw.Statistic.from.List(allStat,"Homo.sapiens","author")
Draw.Statistic.from.List(allStat,"Homo.sapiens","ecoCode","climate")

# Draw a particular (sub)statistic for all species
#out4 <- vector("list", n.species)
out4 <- vector("character", n.species)
for (i in 1:n.species) out4[i] <- Draw.Statistic.from.List(allStat,i,"ecoCode","climate")
#out4 <- unlist(out4)

# Draw yield factors for f0, thinning and females. For all species.
# Relevant yield factors are Y.PX,Y.VX.d,Y.EX.d,Y.CX, and Y.NX*n.CN
out1 <- out2 <- out3 <- out4 <- out5 <- out6 <- SpeciesList <- vector("numeric", n.species)
for (i in 1:n.species) {
	out1[i] <- Draw.Statistic.from.List(popStat,i,"f0","thin1","f","Y.PX")
	out2[i] <- Draw.Statistic.from.List(popStat,i,"f0","thin1","f","Y.VX.d")
	out3[i] <- Draw.Statistic.from.List(popStat,i,"f0","thin1","f","Y.EX.d")
	out4[i] <- Draw.Statistic.from.List(popStat,i,"f0","thin1","f","Y.CX")
	out5[i] <- Draw.Statistic.from.List(popStat,i,"f0","thin1","f","Y.NX")
	out6[i] <- Draw.Statistic(allStat,i,"n.CN")
	SpeciesList[i] <- Draw.Statistic(allStat,i,"species")
}
out <- cbind(out1,out2,out3,out4,out5*out6)
colnames(out) <- c("Y.PX","Y.VX.d","Y.EX.d","Y.CX","Y.NX*n.CN")
rownames(out) <- SpeciesList
apply(out,c(1),sum)
