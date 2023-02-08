# Chick Breeding Scheme / Additive Scenario with Random Selection 2x a Year
# Ivan Pocrnic
# Latest Update: November 2020.

# Clean
rm(list = ls())
# Load libs
library("AlphaSimR")
library("tidyverse")
# Load fje
source("../funkcije_breeding.R")
# Set WD
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Must provide replicate number [1] and number of sires [2] !!!")
}
rep = args[1]
noSires = as.numeric(args[2])
scenarioMain = paste("Additive", noSires, sep="_")
scenarioName = paste("Rep_", rep, "_Additive_RND", sep="")
dire = paste("../../data", scenarioMain, scenarioName, sep = "/")
unlink(dire, recursive = TRUE)
dir.create(path = dire, recursive = TRUE, showWarnings = FALSE)
setwd(dir = dire)
# getwd()

noDams = 1080

# Create blupf90 parameter files and bash scripts
prepare_par()

# Sys.getenv("PATH")
# system("echo $PATH")
# sessionInfo()
# Sys.setenv(PATH=" ")

# For the sake of comparison, use the same Burn-in as for the Standard scenario

bdirname = paste("Rep_", rep, "_Additive", sep="")
bdir = paste("..", bdirname, "burnin.RData", sep = "/")
load(bdir)


# Selection: 10 generations of random selection (10 x 2 years) ----

# Parents in the 1st gen of selection are from the last generation of burnin:
GSStartMales   = StartMalesBurnin
GSStartFemales = StartFemalesBurnin

# Year is the last year of the burn-in:
year = year_burnin

# Separately keep selected/genotyped dams to use for the 1st gen of GS
# Those dams are the same as dams in the PG from the last gen (for the 1st gen of selection)
damg = vector("list", 30)
damg[[6]] = StartFemalesBurnin
# In previous version these were labeled as p3fs

for(gen in 6:15){
  Program = "Random"
  cat("Working on the round:", Program, ":", gen,"\n")
  
  # Q1:
  year = year + 1
  
  # Mate initial sires and dams:
  p_start = selectCross(pop = c(GSStartMales, GSStartFemales), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p2m = selectWithinFam(p_start, nInd = 4, famType = "B", sex = "M", use = "rand")
  p2f = selectWithinFam(p_start, nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # First batch of selection candidates:
  candidatesgroup = 1
  # Generate records:
  RecSys = RecSysMale(RecSys, p2m)
  RecSys = RecSysFemale(RecSys, p2f)

  # Just a dummy "EBV" so function doesn't crash
  p2m@ebv = p2m@gv
  p2f@ebv = p2f@gv
  
  # Summarize Selection Candidates from Batch 1
  CSumm = PullSumm(CSumm, p2m, "M")
  CSumm = PullSumm(CSumm, p2f, "F")
  CSummTogether = PullSummTogether(CSummTogether, c(p2m, p2f))
  
  # Randomly select next generation F & M 
  # Select males:
  p2ms = selectInd(p2m, noSires, use = "rand", sex = "M")
  
  # Mate young selected males with selected females (From Previous Generation)
  # Those females were already selected in previous generation and saved as damg[[gen]]
  p33 = selectCross(pop=c(p2ms, damg[[gen]]), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p3m = selectWithinFam(p33, nInd = 4, famType = "B", sex = "M", use = "rand")
  p3f = selectWithinFam(p33, nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Second batch of selection candidates:
  candidatesgroup = 2
  # Generate records:
  RecSys = RecSysMale(RecSys, p3m)
  RecSys = RecSysFemale(RecSys, p3f)
  
  # Just a dummy "EBV" so function doesn't crash
  p3m@ebv = p3m@gv
  p3f@ebv = p3f@gv
  
  # Summarize Selection Candidates from Batch 2
  CSumm = PullSumm(CSumm, p3m, "M")
  CSumm = PullSumm(CSumm, p3f, "F")
  CSummTogether = PullSummTogether(CSummTogether, c(p3m, p3f))
  
  # Select males:
  p3ms = selectInd(p3m, noSires, use = "rand", sex = "M")
  # Select females:
  p2fs = selectInd(p2f, noDams, use = "rand", sex = "F")

  # Q3:
  year = year + 1
  
  # Mate young selected males with selected females (From Previous Generation)
  p44 = selectCross(pop = c(p3ms, p2fs), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p4m = selectWithinFam(p44, nInd = 4, famType = "B", sex = "M", use = "rand")
  p4f = selectWithinFam(p44, nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Third batch of selection candidates:
  candidatesgroup = 3
  # Generate records:
  RecSys = RecSysMale(RecSys, p4m)
  RecSys = RecSysFemale(RecSys, p4f)
  
  # Just a dummy "EBV" so function doesn't crash
  p4m@ebv = p4m@gv
  p4f@ebv = p4f@gv
  
  # Summarize Selection Candidates from Batch 3
  CSumm = PullSumm(CSumm, p4m, "M")
  CSumm = PullSumm(CSumm, p4f, "F")
  CSummTogether = PullSummTogether(CSummTogether, c(p4m, p4f))
  
  # Select males:
  p4ms = selectInd(p4m, noSires, use = "rand", sex = "M")
  # Select females:
  p3fs = selectInd(p3f, noDams, use = "rand", sex = "F")
  # Actually EBV of these females was already estimated previously.
  
  # Q4:
  
  # Mate young selected males with selected females (From Previous Generation)
  p55 = selectCross(pop = c(p4ms, p3fs), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p5m = selectWithinFam(p55, nInd = 4, famType = "B", sex = "M", use = "rand")
  p5f = selectWithinFam(p55, nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Forth batch of selection candidates:
  candidatesgroup = 4
  
  # Generate records:
  RecSys = RecSysMale(RecSys, p5m)
  RecSys = RecSysFemale(RecSys, p5f)
  
  # Just a dummy "EBV" so function doesn't crash
  p5m@ebv = p5m@gv
  p5f@ebv = p5f@gv
  
  # Summarize Selection Candidates from Batch 4
  CSumm = PullSumm(CSumm, p5m, "M")
  CSumm = PullSumm(CSumm, p5f, "F")
  CSummTogether = PullSummTogether(CSummTogether, c(p5m, p5f))
  
  
  # This are initial parents for the next round:
  # Select males:
  GSStartMales = selectInd(p5m, noSires, use = "rand", sex = "M")
  # In previous schemes called p5ms
  # Select females:
  GSStartFemales = selectInd(p4f, noDams, use = "rand", sex = "F")
  # In previous schemes called p4fs
  # Actually EBV of these females was already estimated previously.
  
  # These are actually selected in Q1 next year after they got T2 Phenotype
  # But same as before, we pretend as we already have their EBV
  # Select females:
  p5fs = selectInd(p5f, noDams, use = "rand", sex = "F")
  damg[[gen+1]] = p5fs
  
}

# Save full image at the end
# save.image("genomic.RData")
save.image("results.RData")
# load("genomic.RData")

