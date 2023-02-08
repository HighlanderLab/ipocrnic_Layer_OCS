# Chick Breeding Scheme / Additive Scenario with truncation GS
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
scenarioName = paste("Rep_", rep, "_Additive_GS", sep="")
dire = paste("../../data",scenarioMain,scenarioName, sep = "/")
unlink(dire, recursive = TRUE)
dir.create(path = dire, recursive = TRUE, showWarnings = FALSE)
setwd(dir = dire)
# getwd()

noDams=1080

# Create blupf90 parameter files and bash scripts
prepare_par()

# Sys.getenv("PATH")
# system("echo $PATH")
# sessionInfo()
# Sys.setenv(PATH=" ")

# For the sake of comparison, use the same Burn-in as for the Standard scenario

bdirname = paste("Rep_", rep, "_Additive", sep="")
bdir = paste("..",bdirname, "burnin.RData", sep = "/")
load(bdir)


# Selection: 10 generations of genomic selection (10 x 2 years) ----

# Parents in the 1st gen of selection are from the last generation of burnin:
GSStartMales   = StartMalesBurnin
GSStartFemales = StartFemalesBurnin

# Year is the last year of the burnin:
year = year_burnin

# Set this dams to something else -> so they are not their mothers?

# Separately keep selected/genotyped dams to use for the 1st gen of GS
# Those dams are the same as dams in the PG from the last gen (for the 1st gen of selection)
damg= vector("list",30)
damg[[6]] = StartFemalesBurnin
# In previous version these were labeled as p3fs

for(gen in 6:15){
  Program = "Genomic"
  cat("Working on the round:",Program, ":", gen,"\n")
  
  # Q1:
  year = year + 1
  
  # Mate initial sires and dams:
  p_start = selectCross(pop=c(GSStartMales,GSStartFemales), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p2m = selectWithinFam(p_start,nInd = 4, famType = "B", sex = "M", use = "rand")
  p2f = selectWithinFam(p_start,nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Prepare genotypes I:
  # Each time select all relevant genotypes (ugly code now!)
  # Also need if/else due to transition phase from burnin to genomics
  if(gen > 6){
    # Selected dams from the G-4 keep genotypes:
    run_prepareGENO(damg[[gen-1]])
    # Selected sires from the G-4 keep genotypes:
    run_prepareGENO(p2ms)
    # Selected dams from the G-3 keep genotypes:
    run_prepareGENO(p2fs)
    # Selected sires from the G-3 keep genotypes:
    run_prepareGENO(p3ms)
    # Selected dams from the G-2 keep genotypes:
    run_prepareGENO(p3fs)
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(p4ms)
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(GSStartMales)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p2m)
  } else {
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(GSStartMales)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p2m)
  }
  
  
  # First batch of selection candidates:
  candidatesgroup = 1
  # Generate records:
  RecSys = RecSysMale(RecSys, p2m)
  RecSys = RecSysFemale(RecSys, p2f)
  
  # Event @ 25 weeks (6m) - p2f get T1 phenotypes
  # Event @ 52 weeks (12m) - damg[[gen]] / p5f get T2 phenotypes (From Previous Generation)
  # Event @ 100 weeks (23m) - p3f get T3 phenotypes (From Previous Generation)
  
  # Q2:
  
  OldDir = getwd()
  Dir = paste("GBlup", Program, gen, candidatesgroup, sep = "_")
  unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
  dir.create(path = Dir, showWarnings = FALSE)
  setwd(dir = Dir)
  #getwd()
  
  # Prepare pedigree and datafile for blupf90
  # Removes phenotypes that are not really available at given timepoint
  run_prepare(RecSys)
  
  # Run blupf90, and update GEBV in recording system
  RecSys = run_gblup(RecSys)
  
  setwd(dir = OldDir)
  
  # Clean marker files from the last round of selection:
  unlink("mrk*")
  unlink("markeri*")
  
  # Set EBV for selection candidates:
  p2m@ebv = as.matrix(RecSys[RecSys$IId %in% p2m@id, c("EbvT1", "EbvT2","EbvT3")])
  p2f@ebv = as.matrix(RecSys[RecSys$IId %in% p2f@id, c("EbvT1", "EbvT2","EbvT3")])
  
  # Summarize Selection Candidates from Batch 1
  CSumm = PullSumm(CSumm,p2m,"M")
  CSumm = PullSumm(CSumm,p2f,"F")
  CSummTogether = PullSummTogether(CSummTogether, c(p2m,p2f))
  
  # WARNING: For the females (p2f) from batch 1 we get EBV on available T1 and T2
  # In reality they dont have T2 at this time point
  # Why doesn't matter: We dont use them till later when they actualy get T2.
  
  # Select next generations F & M based on ebv:
  # Select males:
  p2ms = selectInd(p2m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
  
  # Mate young selected males with selected females (From Previous Generation)
  # Those females were already selected in previous generation and saved as damg[[gen]]
  p33=selectCross(pop=c(p2ms,damg[[gen]]), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p3m = selectWithinFam(p33,nInd = 4, famType = "B", sex = "M", use = "rand")
  p3f = selectWithinFam(p33,nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Prepare genotypes II:
  # Each time select all relevant genotypes (ugly code now!)
  # Also need if/else due to transition phase from burnin to genomics
  if(gen > 6){
    # Selected dams from the G-4 keep genotypes:
    run_prepareGENO(p2fs)
    # Selected sires from the G-4 keep genotypes:
    run_prepareGENO(p3ms)
    # Selected dams from the G-3 keep genotypes:
    run_prepareGENO(p3fs)
    # Selected sires from the G-3 keep genotypes:
    run_prepareGENO(p4ms)
    # Selected dams from the G-2 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(GSStartMales)
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(damg[[gen]])
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(p2ms)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p3m)
  } else {
    # Selected dams from the G-2 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(GSStartMales)
    # Selected dams from the G-1 keep genotypes:
    # run_prepareGENO(damg[[gen]])
    # In first gen of GS these are the same as GSStartFemales
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(p2ms)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p3m)
  }
  
  # Second batch of selection candidates:
  candidatesgroup = 2
  # Generate records:
  RecSys = RecSysMale(RecSys, p3m)
  RecSys = RecSysFemale(RecSys, p3f)
  
  # Event @ 25 weeks (6m) - p3f get T1 phenotypes
  # Event @ 52 weeks (12m) - p2f get T2 phenotypes
  # Event @ 100 weeks (23m) - p4f get T3 phenotypes (From Previous Generation)
  
  OldDir = getwd()
  Dir = paste("GBlup", Program, gen, candidatesgroup, sep = "_")
  unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
  dir.create(path = Dir, showWarnings = FALSE)
  setwd(dir = Dir)
  #getwd()
  
  # Prepare pedigree and datafile for blupf90
  # Removes phenotypes that are not really available at given timepoint
  run_prepare(RecSys)
  
  # Run blupf90, and update GEBV in recording system
  RecSys = run_gblup(RecSys)
  
  setwd(dir = OldDir)
  
  # Clean marker files from the last round of selection:
  unlink("mrk*")
  unlink("markeri*")
  
  # Set EBV for selection candidates:
  p3m@ebv = as.matrix(RecSys[RecSys$IId %in% p3m@id, c("EbvT1", "EbvT2","EbvT3")])
  p3f@ebv = as.matrix(RecSys[RecSys$IId %in% p3f@id, c("EbvT1", "EbvT2","EbvT3")])
  
  # Summarize Selection Candidates from Batch 1
  CSumm = PullSumm(CSumm,p3m,"M")
  CSumm = PullSumm(CSumm,p3f,"F")
  CSummTogether = PullSummTogether(CSummTogether, c(p3m,p3f))
  
  # WARNING: For the females (p3f) from batch 2 we get EBV on available T1 and T2
  # In reality they dont have T2 at this time point
  # Why doesn't matter: We dont use them till later when they actualy get T2.
  
  # Select males:
  p3ms = selectInd(p3m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
  
  # Select females:
  p2fs = selectInd(p2f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
  # Actually EBV of these females was already estimated previously.
  
  # Q3:
  year = year + 1
  
  # Mate young selected males with selected females (From Previous Generation)
  p44=selectCross(pop=c(p3ms,p2fs), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p4m = selectWithinFam(p44,nInd = 4, famType = "B", sex = "M", use = "rand")
  p4f = selectWithinFam(p44,nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Prepare genotypes III:
  # Each time select all relevant genotypes (ugly code now!)
  # Also need if/else due to transition phase from burnin to genomics
  if(gen > 6){
    # Selected dams from the G-4 keep genotypes:
    run_prepareGENO(p3fs)
    # Selected sires from the G-4 keep genotypes:
    run_prepareGENO(p4ms)
    # Selected dams from the G-3 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-3 keep genotypes:
    run_prepareGENO(GSStartMales)
    # Selected dams from the G-2 keep genotypes:
    run_prepareGENO(damg[[gen]])
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(p2ms)
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(p2fs)
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(p3ms)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p4m)
  } else {
    # Selected dams from the G-3 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-3 keep genotypes:
    run_prepareGENO(GSStartMales)
    # Selected dams from the G-2 keep genotypes:
    # run_prepareGENO(damg[[gen]])
    # In first gen of GS these are the same as GSStartFemales
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(p2ms)
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(p2fs)
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(p3ms)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p4m)
  }
  
  # Third batch of selection candidates:
  candidatesgroup = 3
  # Generate records:
  RecSys = RecSysMale(RecSys, p4m)
  RecSys = RecSysFemale(RecSys, p4f)
  
  # Event @ 25 weeks (6m) - p4f get T1 phenotypes
  # Event @ 52 weeks (12m) - p3f get T2 phenotypes
  # Event @ 100 weeks (23m) - damg[[gen]] / p5f get T3 phenotypes (From Previous Generation)
  
  OldDir = getwd()
  Dir = paste("GBlup", Program, gen, candidatesgroup, sep = "_")
  unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
  dir.create(path = Dir, showWarnings = FALSE)
  setwd(dir = Dir)
  #getwd()
  
  # Prepare pedigree and datafile for blupf90
  # Removes phenotypes that are not rally available at given timepoint
  run_prepare(RecSys)
  
  # Run blupf90, and update GEBV in recording system
  RecSys = run_gblup(RecSys)
  
  setwd(dir = OldDir)
  
  # Clean marker files from the last round of selection:
  unlink("mrk*")
  unlink("markeri*")
  
  # Set EBV for selection candidates:
  p4m@ebv = as.matrix(RecSys[RecSys$IId %in% p4m@id, c("EbvT1", "EbvT2","EbvT3")])
  p4f@ebv = as.matrix(RecSys[RecSys$IId %in% p4f@id, c("EbvT1", "EbvT2","EbvT3")])
  
  # Summarize Selection Candidates from Batch 1
  CSumm = PullSumm(CSumm,p4m,"M")
  CSumm = PullSumm(CSumm,p4f,"F")
  CSummTogether = PullSummTogether(CSummTogether, c(p4m,p4f))
  
  # WARNING: For the females (p4f) from batch 3 we get EBV on available T1 and T2
  # In reality they dont have T2 at this time point
  # Why doesn't matter: We dont use them till later when they actualy get T2.
  
  # Select males:
  p4ms = selectInd(p4m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
  
  # Select females:
  p3fs = selectInd(p3f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
  # Actually EBV of these females was already estimated previously.
  
  # Q4:
  
  # Mate young selected males with selected females (From Previous Generation)
  p55=selectCross(pop=c(p4ms,p3fs), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
  
  # Selected first 4 M & 9 F hatched per mating (per full-sib family):
  p5m = selectWithinFam(p55,nInd = 4, famType = "B", sex = "M", use = "rand")
  p5f = selectWithinFam(p55,nInd = 9, famType = "B", sex = "F", use = "rand")
  
  # Prepare genotypes IV:
  # Each time select all relevant genotypes (ugly code now!)
  # Also need if/else due to transition phase from burnin to genomics
  if(gen > 6){
    # Selected dams from the G-4 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-4 keep genotypes:
    run_prepareGENO(GSStartMales)
    # Selected dams from the G-3 keep genotypes:
    run_prepareGENO(damg[[gen]])
    # Selected sires from the G-3 keep genotypes:
    run_prepareGENO(p2ms)
    # Selected dams from the G-2 keep genotypes:
    run_prepareGENO(p2fs)
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(p3ms)
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(p3fs)
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(p4ms)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p5m)
  } else {
    # Selected dams from the G-4 keep genotypes:
    run_prepareGENO(GSStartFemales)
    # Selected sires from the G-4 keep genotypes:
    run_prepareGENO(GSStartMales)
    # Selected dams from the G-3 keep genotypes:
    # run_prepareGENO(damg[[gen]])
    # In first gen of GS these are the same as GSStartFemales
    # Selected sires from the G-3 keep genotypes:
    run_prepareGENO(p2ms)
    # Selected dams from the G-2 keep genotypes:
    run_prepareGENO(p2fs)
    # Selected sires from the G-2 keep genotypes:
    run_prepareGENO(p3ms)
    # Selected dams from the G-1 keep genotypes:
    run_prepareGENO(p3fs)
    # Selected sires from the G-1 keep genotypes:
    run_prepareGENO(p4ms)
    # All male candidates from this generation get genotypes:
    run_prepareGENO(p5m)
  }
  
  # Forth batch of selection candidates:
  candidatesgroup = 4
  # Generate records:
  RecSys = RecSysMale(RecSys, p5m)
  RecSys = RecSysFemale(RecSys, p5f)
  
  # Event @ 25 weeks (6m) - p5f get T1 phenotypes
  # Event @ 52 weeks (12m) - p4f get T2 phenotypes
  # Event @ 100 weeks (23m) - p2f get T3 phenotypes
  
  OldDir = getwd()
  Dir = paste("GBlup", Program, gen, candidatesgroup, sep = "_")
  unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
  dir.create(path = Dir, showWarnings = FALSE)
  setwd(dir = Dir)
  #getwd()
  
  # Prepare pedigree and datafile for blupf90
  # Removes phenotypes that are not rally available at given timepoint
  run_prepare(RecSys)
  
  # Run blupf90, and update GEBV in recording system
  RecSys = run_gblup(RecSys)
  
  setwd(dir = OldDir)
  
  # Clean marker files from the last round of selection:
  unlink("mrk*")
  unlink("markeri*")
  
  # Set EBV for selection candidates:
  p5m@ebv = as.matrix(RecSys[RecSys$IId %in% p5m@id, c("EbvT1", "EbvT2","EbvT3")])
  p5f@ebv = as.matrix(RecSys[RecSys$IId %in% p5f@id, c("EbvT1", "EbvT2","EbvT3")])
  
  # Summarize Selection Candidates from Batch 1
  CSumm = PullSumm(CSumm,p5m,"M")
  CSumm = PullSumm(CSumm,p5f,"F")
  CSummTogether = PullSummTogether(CSummTogether, c(p5m,p5f))
  
  # WARNING: For the females (p5f) from batch 3 we get EBV on available T1 and T2
  # In reality they don't have T2 at this time point
  # Why it doesn't matter: We don't use them till later when they actualy get T2
  
  # This are initial parents for the next ronud:
  # Select males:
  GSStartMales = selectInd(p5m, noSires, use = "ebv", sex = "M", trait=selIndex, b=iweight)
  # In previous schemes called p5ms
  # Select females:
  GSStartFemales = selectInd(p4f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
  # In previous schemes called p4fs
  # Actually EBV of these females was already estimated previously.
  
  # These are actualy selected in Q1 next year after they got T2 Phenotype
  # But same as before, we pretend as we already have their EBV
  # Select females:
  p5fs = selectInd(p5f, noDams, use = "ebv", sex = "F", trait=selIndex, b=iweight)
  damg[[gen+1]] = p5fs
  
  # For testing purposes, save image after each cycle (Remove from the final code)
  # put=paste("genomic", gen, ".RData", sep = "_")
  # save.image(put)
}

# Save full image at the end
# save.image("genomic.RData")
save.image("results.RData")
# load("genomic.RData")

# If sucesseful then clean: 
# rm genomic_*
# rm mr*
# rm mar*
# rm renum*
# rm extractsol.sh
# rm prepgeno.sh
