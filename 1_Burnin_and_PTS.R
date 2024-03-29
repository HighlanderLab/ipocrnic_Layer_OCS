# Long-term selection in layers (Founders/Burn-in PTS/PTS)
# Ivan Pocrnic

# Clean
rm(list = ls())
# Load packages
library("AlphaSimR")
library("tidyverse")
# Load fje
source("../functions.R")
# Set WD
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Must provide replicate number [1] and number of sires [2] !!!")
}
rep = args[1]
noSires = as.numeric(args[2])
scenarioMain = paste("Additive", noSires, sep="_")
scenarioName = paste("Rep_", rep, "_Additive", sep="")
dire = paste("../../data",scenarioMain,scenarioName, sep = "/")
unlink(dire, recursive = TRUE)
dir.create(path = dire, recursive = TRUE, showWarnings = FALSE)
setwd(dir = dire)

noDams = 1080

# Create blupf90 parameter files and bash scripts
prepare_par()

#### Founder Population ####
BaseNe = 100
ChrSize = (1.2 * 10^9) / 39
MutRate = 5E-8
RecRate = 2.5E-8
MaCSeNFlags = "-eN 1.25 2.5 -eN 2.5 10 -eN 12.5 100 -eN 25 500 -eN 250 2500 -eN 2500 5000"
founderPop = runMacs(nInd = 2500,
                     nChr = 39,
                     segSites = 2500,
                     manualCommand = paste(as.integer(ChrSize),
                                           "-t", MutRate * 4 * BaseNe,
                                           "-r", RecRate * 4 * BaseNe,
                                           MaCSeNFlags),
                     manualGenLen = RecRate * ChrSize)


save.image("founders.RData")

# Trait parameters:
Addh2 = c(0.18, 0.22, 0.25)
GenCor = matrix(c(1.00, 0.75, 0.60, 0.75, 1.00, 0.70, 0.60, 0.70, 1.00), 3, 3)
iweight = c(0.20, 0.35, 0.45)
meangv = c(0, 0, 0)
addvar = c(1, 1, 1)

# Simulation parameters:
SP = SimParam$new(founderPop)

SP$restrSegSites(
  minQtlPerChr = 250,
  minSnpPerChr = 2000,
  overlap = FALSE,
  minSnpFreq = NULL)

SP$addTraitA(nQtlPerChr = 250, corA = GenCor, mean = meangv, var = addvar)

SP$setVarE(h2 = Addh2)

SP$setSexes("yes_sys")

# Creating new SNP-panels from available segregating sites

exclude = SP$invalidSnp

# Pick loci one chromosome at a time
chip1 = chip2 = vector("list", 39) # One element per chr
all_loc = 1:2500 # Indicator for all loci on chr
for(i in 1:39){
  take = sample( all_loc[-exclude[[i]]], 2000)
  # First 1000 for chip 1
  chip1[[i]] = take[1:1000]
  chip1[[i]] = sort(chip1[[i]]) # Must be sorted
  # Last 1000 for chip 1
  chip2[[i]] = c(take[1001:2000])
  chip2[[i]] = sort(chip2[[i]]) # Must be sorted
}
chip1 = do.call("c", chip1) # Collapse list to vector
chip2 =  do.call("c", chip2)

# Create LociMap objects
snpChip1 = new("LociMap",
               nLoci=39000L,
               lociPerChr=rep(1000L,39),
               lociLoc=chip1)
snpChip2 = new("LociMap",
               nLoci=39000L,
               lociPerChr=rep(1000L,39),
               lociLoc=chip2)

# We will trace 3 panels; (1) QTL-panel, (2) SNP-panel used for GS, and (3) SNP-Neutral-panel for the control 
SP$snpChips[[1]] = snpChip1
SP$snpChips[[2]] = snpChip2

#### Burn-In ####

# Generate initial founders
Parents = newPop(founderPop, simParam = SP)

# Select randomly 1080 Females and 40/120 Males from the founder population
StartMales = selectInd(Parents, noSires, use = "rand", sex = "M")
StartFemales = selectInd(Parents, noDams, use = "rand", sex = "F")

# Gather initial frequencies and initial heterozygosity 
M = pullSnpGeno(c(StartMales, StartFemales))
N = pullSnpGeno(c(StartMales, StartFemales), snpChip = 2)
QT1 = pullQtlGeno(c(StartMales, StartFemales), trait = 1) 
QT2 = pullQtlGeno(c(StartMales, StartFemales), trait = 2) 
QT3 = pullQtlGeno(c(StartMales, StartFemales), trait = 3) 
# Calculate Observed and Expected Heterozygosity
## For selected SNP
het0_snp = apply(X = M, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) ) 
het0_snp = ifelse(het0_snp > 0, yes = het0_snp, no = het0_snp + 1e-12)
p0_snp = colMeans(M)/2
## For neutral SNP
het0_neutral = apply(X = N, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
het0_neutral = ifelse(het0_neutral > 0, yes = het0_neutral, no = het0_neutral + 1e-12)
p0_neutral = colMeans(N)/2
## For QTL
het0_qtl1 = apply(X = QT1, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) ) 
het0_qtl1 = ifelse(het0_qtl1 > 0, yes = het0_qtl1, no = het0_qtl1 + 1e-12)
het0_qtl2 = apply(X = QT2, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) ) 
het0_qtl2 = ifelse(het0_qtl2 > 0, yes = het0_qtl2, no = het0_qtl2 + 1e-12)
het0_qtl3 = apply(X = QT3, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
het0_qtl3 = ifelse(het0_qtl3 > 0, yes = het0_qtl3, no = het0_qtl3 + 1e-12)
p0_qtl1 = colMeans(QT1)/2
p0_qtl2 = colMeans(QT2)/2
p0_qtl3 = colMeans(QT3)/2

rm(Parents, M, N, QT1, QT2, QT3)

# Initialize recording database ("RecSys") and generate records for the initial parents:
RecSys = NULL
gen = 0
candidatesgroup = 0
year = 0
Program = "Founders"

RecSys = RecSysMale(RecSys, StartMales)
RecSys = RecSysFemale(RecSys, StartFemales)

# Initialize Candidates Summary (Sex separated and Together)
CSumm = NULL
CSummTogether = NULL
for(gen in 1:5){
    Program = "Burnin"
    cat("Working on the round:",Program, ":", gen,"\n")

    year = year + 1

    # Mate initial sires and dams:
    p_start = selectCross(pop = c(StartMales, StartFemales), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)

    # Select first 4 M & 9 F hatched per mating (per full-sib family):
    p2m = selectWithinFam(p_start,nInd = 4, famType = "B", sex = "M", use = "rand")
    p2f = selectWithinFam(p_start,nInd = 9, famType = "B", sex = "F", use = "rand")

    candidatesgroup = 1
    RecSys = RecSysMale(RecSys, p2m)
    RecSys = RecSysFemale(RecSys, p2f)

    OldDir = getwd()
    Dir = paste("Blup", Program, gen, candidatesgroup, sep = "_")
    unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
    dir.create(path = Dir, showWarnings = FALSE)
    setwd(dir = Dir)

    # Prepare pedigree and datafile for blupf90
    # Removes phenotypes that are not really available at given timepoint
    run_prepare(RecSys)

    # Run blupf90, and update EBV in recording system
    RecSys = run_blup(RecSys)

    setwd(dir = OldDir)

    # Set EBV for selection candidates:
    p2m@ebv = as.matrix(RecSys[RecSys$IId %in% p2m@id, c("EbvT1", "EbvT2","EbvT3")])
    p2f@ebv = as.matrix(RecSys[RecSys$IId %in% p2f@id, c("EbvT1", "EbvT2","EbvT3")])

    CSumm = PullSumm(CSumm,p2m,"M")
    CSumm = PullSumm(CSumm,p2f,"F")
    CSummTogether = PullSummTogether(CSummTogether, c(p2m,p2f))

    # Select next generations F & M based on ebv:
    # Select males:
    p2ms = selectInd(p2m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
    # Select females:
    p2fs = selectInd(p2f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)

    p33=selectCross(pop=c(p2ms,p2fs), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
    
    year = year + 1
    
    # Selected first 4 M & 9 F hatched per mating (per full-sib family):
    p3m = selectWithinFam(p33,nInd = 4, famType = "B", sex = "M", use = "rand")
    p3f = selectWithinFam(p33,nInd = 9, famType = "B", sex = "F", use = "rand")

    candidatesgroup = 2
    RecSys = RecSysMale(RecSys, p3m)
    RecSys = RecSysFemale(RecSys, p3f)

    OldDir = getwd()
    Dir = paste("Blup", Program, gen, candidatesgroup, sep = "_")
    unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
    dir.create(path = Dir, showWarnings = FALSE)
    setwd(dir = Dir)

    # Prepare pedigree and datafile for blupf90
    # Removes phenotypes that are not really available at given timepoint
    run_prepare(RecSys)

    # Run blupf90, and update EBV in recording system
    RecSys = run_blup(RecSys)

    setwd(dir = OldDir)

    # Set EBV for selection candidates:
    p3m@ebv = as.matrix(RecSys[RecSys$IId %in% p3m@id, c("EbvT1", "EbvT2","EbvT3")])
    p3f@ebv = as.matrix(RecSys[RecSys$IId %in% p3f@id, c("EbvT1", "EbvT2","EbvT3")])

    CSumm = PullSumm(CSumm,p3m,"M")
    CSumm = PullSumm(CSumm,p3f,"F")
    CSummTogether = PullSummTogether(CSummTogether, c(p3m,p3f))

    # Select next generations F & M based on ebv:
    # Select males:
    StartMales = selectInd(p3m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
    # Select females:
    StartFemales = selectInd(p3f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
}

# Save last generation start males and females so they won't be overwritten
StartMalesBurnin = StartMales
StartFemalesBurnin = StartFemales

# Save last year of burn-in
year_burnin = year

save.image("burnin.RData")


#### Standard Truncation Selection ####

for(gen in 6:15){
  Program = "Standard"
  cat("Working on the round:",Program, ":", gen,"\n")
  
  year = year + 1
  
  # Mate initial sires and dams:
  p_start = selectCross(pop=c(StartMales,StartFemales), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Select first 4 M & 9 F hatched per mating (per full-sib family):
  p2m = selectWithinFam(p_start,nInd = 4, famType = "B", sex = "M", use = "rand")
  p2f = selectWithinFam(p_start,nInd = 9, famType = "B", sex = "F", use = "rand")
  
  candidatesgroup = 1
  RecSys = RecSysMale(RecSys, p2m)
  RecSys = RecSysFemale(RecSys, p2f)
  
  OldDir = getwd()
  Dir = paste("Blup", Program, gen, candidatesgroup, sep = "_")
  unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
  dir.create(path = Dir, showWarnings = FALSE)
  setwd(dir = Dir)

  # Prepare pedigree and datafile for blupf90
  # Removes phenotypes that are not really available at given timepoint
  run_prepare(RecSys)
  
  # Run blupf90, and update EBV in recording system
  RecSys = run_blup(RecSys)
  
  setwd(dir = OldDir)
  
  # Set EBV for selection candidates:
  p2m@ebv = as.matrix(RecSys[RecSys$IId %in% p2m@id, c("EbvT1", "EbvT2","EbvT3")])
  p2f@ebv = as.matrix(RecSys[RecSys$IId %in% p2f@id, c("EbvT1", "EbvT2","EbvT3")])
  
  CSumm = PullSumm(CSumm,p2m,"M")
  CSumm = PullSumm(CSumm,p2f,"F")
  CSummTogether = PullSummTogether(CSummTogether, c(p2m,p2f))

  # Select next generations F & M based on ebv:
  # Select males:
  p2ms = selectInd(p2m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
  # Select females:
  p2fs = selectInd(p2f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
  
  p33=selectCross(pop=c(p2ms,p2fs), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  year = year + 1
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p3m = selectWithinFam(p33,nInd = 4, famType = "B", sex = "M", use = "rand")
  p3f = selectWithinFam(p33,nInd = 9, famType = "B", sex = "F", use = "rand")
  
  candidatesgroup = 2
  RecSys = RecSysMale(RecSys, p3m)
  RecSys = RecSysFemale(RecSys, p3f)
  
  OldDir = getwd()
  Dir = paste("Blup", Program, gen, candidatesgroup, sep = "_")
  unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
  dir.create(path = Dir, showWarnings = FALSE)
  setwd(dir = Dir)
  
  # Prepare pedigree and datafile for blupf90
  # Removes phenotypes that are not really available at given timepoint
  run_prepare(RecSys)
  
  # Run blupf90, and update EBV in recording system
  RecSys = run_blup(RecSys)
  
  setwd(dir = OldDir)
  
  # Set EBV for selection candidates:
  p3m@ebv = as.matrix(RecSys[RecSys$IId %in% p3m@id, c("EbvT1", "EbvT2","EbvT3")])
  p3f@ebv = as.matrix(RecSys[RecSys$IId %in% p3f@id, c("EbvT1", "EbvT2","EbvT3")])
  
  CSumm = PullSumm(CSumm,p3m,"M")
  CSumm = PullSumm(CSumm,p3f,"F")
  CSummTogether = PullSummTogether(CSummTogether, c(p3m,p3f))
  
  # Select next generations F & M based on ebv:
  # Select males:
  StartMales = selectInd(p3m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
  # Select females:
  StartFemales = selectInd(p3f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
}

save.image("results.RData")

