# Long-term selection in layers (UGOCS)
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
  stop("Must provide replicate number [1], number of sires [2] and OCS degrees [3] !!!")
}
rep = args[1]
noSires = as.numeric(args[2])
degrees = as.numeric(args[3])
scenarioMain = paste("Additive", noSires, sep="_")
scenarioName = paste("Rep_", rep, "_Additive_OCSFREE_", degrees, sep = "")
dire = paste("../../data",scenarioMain,scenarioName, sep = "/")
unlink(dire, recursive = TRUE)
dir.create(path = dire, recursive = TRUE, showWarnings = FALSE)
setwd(dir = dire)

noDams = 1080

# Create blupf90 parameter files and bash scripts
prepare_par()

# For the sake of comparison, use the same Burn-in as for the Standard scenario

bdirname = paste("Rep_", rep, "_Additive", sep = "")
bdir = paste("..",bdirname, "burnin.RData", sep = "/")
load(bdir)


#### Selection: 10 generations of genomic selection (10 x 2 years) ####

# Parents in the 1st gen of GS  are from the last generation of burnin:
GSStartMales   = StartMalesBurnin
GSStartFemales = StartFemalesBurnin

# Year is the last year of the burnin:
year = year_burnin

# Separately keep selection candidates dams to use for the 1st gen of OCS
# For this scenario (OCS for both males and females) these dams are not selected, but are from the last gen. burnin selection candidates (p3f)
# Kind of bizarre situation: Those selected in burn-in, are also starting females for this scenario  
# This has to be done when switching from BLUP burn-in to GS-OCS scenarios 
damg = vector("list",30)
damg[[6]] = p3f

# Initiate the object for saving number of OCS-selected sires, and their contributions
SireContribution = NULL

for(gen in 6:15){
    Program = paste("OCSFREE_", degrees, sep = "")
    cat("Working on the round:", Program, ":", gen,"\n")
    
    # Q1:
    year = year + 1
     
    if(gen > 6){
      p_start = makeCross2(females = GSStartFemales, males = GSStartMales, crossPlan = crossPlanM)
    } else {
      p_start = selectCross(pop = c(GSStartMales, GSStartFemales), nFemale = noDams, nMale = noSires, nCrosses = noDams, nProgeny = 20)
    }
    
    # Selected first 4 M & 9 F hatched per mating (per full-sib family):
    p2m = selectWithinFam(p_start, nInd = 4, famType = "B", sex = "M", use = "rand")
    p2f = selectWithinFam(p_start, nInd = 9, famType = "B", sex = "F", use = "rand")

    # Prepare genotypes I:
    # Each time select all relevant genotypes (ugly code now!)
    # Also need if/else due to transition phase from burnin to genomics
    if(gen > 6){
      # Selected dams from the G-4 keep genotypes:
      run_prepareGENO(p1fs)
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


    #### First batch of selection candidates ####
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
    # getwd()

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
    CSumm = PullSumm(CSumm, p2m, "M")
    CSumm = PullSumm(CSumm, p2f, "F")
    CSummTogether = PullSummTogether(CSummTogether, c(p2m,p2f))

    # WARNING: For the females (p2f) from batch 1 we get EBV on available T1 and T2
    # In reality they don't have T2 at this time point
    # Why doesn't matter: We don't use them till later when they actually get T2.

    # Select next generations F & M based on OCS
    
    # Select populations to be included in OCS
    candidates_plan = c(damg[[gen]], p2m)
    
    # Run OCS:
    run_OCS_EASED(candidates_plan, RecSys, degrees, noDams)
    system(command = paste("cp AlphaMate.log AlphaMate.log",gen,year,candidatesgroup,sep = "_"))
    
    # Read crossPlan
    crossPlan = read.table(file="ContributorsModeOptTarget1.txt", header = TRUE)
    
    # Save summary of sires contribution
    SireContribution = PullSireContribution(SireContribution, crossPlan)

    # Population of selected males based on OCS:
    p2ms = p2m[p2m@id %in% filter(crossPlan, Gender == 1)[[1]],]
    # a1=(filter(crossPlan, Gender == 1) %>% dplyr::select(Id))
    
    # Population of selected females based on OCS:
    p1fs_tmp = damg[[gen]]
    p1fs = p1fs_tmp[p1fs_tmp@id %in% filter(crossPlan, Gender == 2)[[1]],]
    
    # Create mate allocation, and expand
    crossPlanM = run_allocate_matings(20, crossPlan)
    
    # Mate young OCS selected males with OCS selected females (from previous generation)
    p33 = makeCross2(females = p1fs, males = p2ms, crossPlan = crossPlanM)
    
    # Selected first 4 M & 9 F hatched per mating (per full-sib family):
    p3m = selectWithinFam(p33, nInd = 4, famType = "B", sex = "M", use = "rand")
    p3f = selectWithinFam(p33, nInd = 9, famType = "B", sex = "F", use = "rand")

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
      run_prepareGENO(p1fs)
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
      # Some of those can be equal to initial dams so remove duplicates:
      tmp_gd = p1fs[!(p1fs@id) %in% GSStartFemales@id,]
      run_prepareGENO(tmp_gd)
      # Selected sires from the G-1 keep genotypes:
      run_prepareGENO(p2ms)
      # All male candidates from this generation get genotypes:
      run_prepareGENO(p3m)
    }


    ##### Second batch of selection candidates ####
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

    # Summarize Selection Candidates from Batch 2
    CSumm = PullSumm(CSumm, p3m, "M")
    CSumm = PullSumm(CSumm, p3f, "F")
    CSummTogether = PullSummTogether(CSummTogether, c(p3m,p3f))

    # WARNING: For the females (p3f) from batch 2 we get EBV on available T1 and T2
    # In reality they don't have T2 at this time point
    # Why doesn't matter: We don't use them till later when they actually get T2.

    # Select next generations F & M based on OCS
    
    # Select populations to be included in OCS
    candidates_plan = c(p2f, p3m)
    
    # Run OCS
    run_OCS_EASED(candidates_plan, RecSys, degrees, noDams)
    system(command = paste("cp AlphaMate.log AlphaMate.log",gen,year,candidatesgroup,sep = "_"))
    
    # Read crossPlan
    crossPlan = read.table(file = "ContributorsModeOptTarget1.txt", header = TRUE)
    
    # Save summary of sires contribution
    SireContribution = PullSireContribution(SireContribution, crossPlan)
    
    # Population of selected males based on OCS:
    p3ms = p3m[p3m@id %in% filter(crossPlan, Gender == 1)[[1]],]

    # Population of selected females based on OCS:
    p2fs = p2f[p2f@id %in% filter(crossPlan, Gender == 2)[[1]],]
    
    # Create mate allocation, and expand
    crossPlanM = run_allocate_matings(20, crossPlan)
    
    # Q3:
    year = year + 1
    # Mate young OCS selected males with OCS selected females (from previous generation)
    p44 = makeCross2(females = p2fs, males = p3ms, crossPlan = crossPlanM)

    # Selected first 4 M & 9 F hatched per mating (per full-sib family):
    p4m = selectWithinFam(p44, nInd = 4, famType = "B", sex = "M", use = "rand")
    p4f = selectWithinFam(p44, nInd = 9, famType = "B", sex = "F", use = "rand")

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
      run_prepareGENO(p1fs)
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
      # Some of those can be equal to initial dams so remove duplicates:
      # tmp_gd = p1fs[!(p1fs@id) %in% GSStartFemales,]
      run_prepareGENO(tmp_gd)
      # Selected sires from the G-2 keep genotypes:
      run_prepareGENO(p2ms)
      # Selected dams from the G-1 keep genotypes:
      run_prepareGENO(p2fs)
      # Selected sires from the G-1 keep genotypes:
      run_prepareGENO(p3ms)
      # All male candidates from this generation get genotypes:
      run_prepareGENO(p4m)
    }


    #### Third batch of selection candidates ####
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
    # Removes phenotypes that are not really available at given timepoint
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

    # Summarize Selection Candidates from Batch 3
    CSumm = PullSumm(CSumm, p4m, "M")
    CSumm = PullSumm(CSumm, p4f, "F")
    CSummTogether = PullSummTogether(CSummTogether, c(p4m,p4f))

    # WARNING: For the females (p4f) from batch 3 we get EBV on available T1 and T2
    # In reality they dont have T2 at this time point
    # Why doesn't matter: We don't use them till later when they actually get T2.
    
    # Select next generations F & M based on OCS
    
    # Select populations to be included in OCS
    candidates_plan = c(p3f, p4m)
    
    # Run OCS
    run_OCS_EASED(candidates_plan, RecSys, degrees, noDams)
    system(command = paste("cp AlphaMate.log AlphaMate.log",gen,year,candidatesgroup,sep = "_"))
    
    # Read crossPlan
    crossPlan = read.table(file="ContributorsModeOptTarget1.txt", header = TRUE)
    
    # Save summary of sires contribution
    SireContribution = PullSireContribution(SireContribution, crossPlan)
    
    # Population of selected males based on OCS:
    p4ms = p4m[p4m@id %in% filter(crossPlan, Gender == 1)[[1]],]
    
    # Population of selected females based on OCS:
    p3fs = p3f[p3f@id %in% filter(crossPlan, Gender == 2)[[1]],]
    
    # Create mate allocation, and expand
    crossPlanM = run_allocate_matings(20, crossPlan)
    
    # Q4:
    
    # Mate young OCS selected males with OCS selected females (from previous generation)
    p55 = makeCross2(females = p3fs, males = p4ms, crossPlan = crossPlanM)

    # Selected first 4 M & 9 F hatched per mating (per full-sib family):
    p5m = selectWithinFam(p55, nInd = 4, famType = "B", sex = "M", use = "rand")
    p5f = selectWithinFam(p55, nInd = 9, famType = "B", sex = "F", use = "rand")

    # Prepare genotypes IV:
    # Each time select all relevant genotypes (ugly code now!)
    # Also need if/else due to transition phase from burnin to genomics
    if(gen > 6){
      # Selected dams from the G-4 keep genotypes:
      run_prepareGENO(GSStartFemales)
      # Selected sires from the G-4 keep genotypes:
      run_prepareGENO(GSStartMales)
      # Selected dams from the G-3 keep genotypes:
      run_prepareGENO(p1fs)
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
      # Some of those can be equal to initial dams so remove duplicates:
      # tmp_gd = p1fs[!(p1fs@id) %in% GSStartFemales,]
      run_prepareGENO(tmp_gd)
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

    #### Forth batch of selection candidates ####
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
    # Removes phenotypes that are not really available at given timepoint
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

    # Summarize Selection Candidates from Batch 4
    CSumm = PullSumm(CSumm, p5m, "M")
    CSumm = PullSumm(CSumm, p5f, "F")
    CSummTogether = PullSummTogether(CSummTogether, c(p5m,p5f))

    # WARNING: For the females (p5f) from batch 4 we get EBV on available T1 and T2
    # In reality they don't have T2 at this time point
    # Why it doesn't matter: We don't use them till later when they actually get T2

    # This are initial parents for the next round
    # Select next generations F & M based on OCS
    
    # Select populations to be included in OCS
    candidates_plan = c(p4f, p5m)
    
    # Run OCS:
    run_OCS_EASED(candidates_plan, RecSys, degrees, noDams)
    system(command = paste("cp AlphaMate.log AlphaMate.log",gen,year,candidatesgroup,sep = "_"))
    
    # Read crossPlan
    crossPlan = read.table(file="ContributorsModeOptTarget1.txt", header = TRUE)
    
    # Save summary of sires contribution
    SireContribution = PullSireContribution(SireContribution, crossPlan)

    # Population of selected males based on OCS:
    GSStartMales = p5m[p5m@id %in% filter(crossPlan, Gender == 1)[[1]],]

    # Population of selected females based on OCS:
    GSStartFemales = p4f[p4f@id %in% filter(crossPlan, Gender == 2)[[1]],]
    
    # Create mate allocation, and expand
    crossPlanM = run_allocate_matings(20, crossPlan)
    
    # These are actually selected in Q1 next year after they got T2 Phenotype
    # But same as before, we pretend as we already have their EBV
    
    damg[[gen+1]] = p5f

    put = paste("Sires_Free", gen, ".RData", sep = "_")
    save(SireContribution, file = put)
}

save.image("results.RData")

