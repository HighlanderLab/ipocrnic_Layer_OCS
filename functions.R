# Long-term selection in layers (support functions)
# Ivan Pocrnic


RecSysMale <- function(datafile, popname) {
datafile = rbind(datafile,
                 tibble(Generation = gen,
                        IId        = popname@id,
                        FId        = popname@father,
                        MId        = popname@mother,
                        Sex        = popname@sex,
                        Year       = year,
                        Program    = Program,
                        SelCand    = candidatesgroup,
                        PhenoT1    = NA,
                        PhenoT2    = NA,
                        PhenoT3    = NA,
                        EbvT1      = NA,
                        EbvT2      = NA,
                        EbvT3      = NA,
                        TbvT1      = popname@gv[, 1],
                        TbvT2      = popname@gv[, 2],
                        TbvT3      = popname@gv[, 3]))
#return(datafile)
}


RecSysFemale <- function(datafile, popname) {
datafile = rbind(datafile,
                 tibble(Generation = gen,
                        IId        = popname@id,
                        FId        = popname@father,
                        MId        = popname@mother,
                        Sex        = popname@sex,
                        Year       = year,
                        Program    = Program,
                        SelCand    = candidatesgroup,
                        PhenoT1    = popname@pheno[,1],
                        PhenoT2    = popname@pheno[,2],
                        PhenoT3    = popname@pheno[,3],
                        EbvT1      = NA,
                        EbvT2      = NA,
                        EbvT3      = NA,
                        TbvT1      = popname@gv[, 1],
                        TbvT2      = popname@gv[, 2],
                        TbvT3      = popname@gv[, 3]))
#return(datafile)
}


run_blup <- function(datafile) {
  system(command = "cp ../renum.par .")
  system(command = "cp ../extractsol.sh .")
  system(command = "echo renum.par | $HOME/bin/renumf90 | tee renum.log")
  system(command = "echo renf90.par | $HOME/bin/blupf90 | tee blup.log")
  system(command = "sh extractsol.sh")
  # Update breeding values in database:
  sol1 <- read.table("sl3.tmp")
  sol1 <- sol1[sol1$V2 %in% datafile$IId,]
  sol1 <- sol1[order(match(sol1$V2, datafile$IId)),]
  datafile$EbvT1 <- sol1$V3
  datafile$EbvT2 <- sol1$V4
  datafile$EbvT3 <- sol1$V5
  return(datafile)
}


run_gblup <- function(datafile) {
  system(command = "cp ../renumGS.par .")
  system(command = "cp ../extractsol.sh .")
  system(command = "echo renumGS.par | $HOME/bin/renumf90 | tee renum.log")
  system(command = "echo renf90.par | $HOME/bin/blupf90 | tee blup.log")
  system(command = "sh extractsol.sh")
  # Update breeding values in database:
  sol1 <- read_table2("sl3.tmp",col_names=FALSE,
                    col_types = cols(.default = col_character(),
                    X1 = col_character(),
                    X2 = col_character(),
                    X3 = col_double(),
                    X4 = col_double(),
                    X5 = col_double()))
  sol1_1 <- sol1[sol1$X2 %in% datafile$IId,]
  sol1_1 <- sol1_1[order(match(sol1_1$X2, datafile$IId)),]
  datafile$EbvT1 <- sol1_1$X3
  datafile$EbvT2 <- sol1_1$X4
  datafile$EbvT3 <- sol1_1$X5
  return(datafile)
}


# Input: candidates_plan & RecSys 
run_OCS_NOMP <- function(ocspop, datafile, deg, sires, dams) {
  # Create index of EBV for criterion:
  tmp = data.frame(id=ocspop@id, crit=(0.2*ocspop@ebv[,1]+0.35*ocspop@ebv[,2])+0.45*ocspop@ebv[,3])
  write.table(x=tmp, file="ebv.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # Assign sex:
  tmp = data.frame(id=ocspop@id, genderRole=ocspop@sex)
  tmp$genderRole = as.numeric(factor(ocspop@sex, levels=c("M", "F")))
  write.table(x=tmp, file="GenderRole.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # Full (All Animals) Pedigree: 
  tmp <- data.frame(datafile[,c("IId", "FId", "MId")])
  tmp$FId[is.na(tmp$FId)] <- 0
  tmp$MId[is.na(tmp$MId)] <- 0
  write.table(tmp, file="Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # ID of the OCS candidates: 
  idcand=ocspop@id
  write.table(idcand, file="PedigreeSet.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # AlphaRelate parameter file: 
  sink(file = "AlphaRelateSpec.txt")
  cat("PedigreeFile                , Blupf90.ped\n")
  cat("PedInbreeding               , Yes\n")
  cat("PedNrm                      , Yes\n")
  cat("PedNrmSubsetFile            , PedigreeSet.ped\n")
  cat("OutputFormat                , f7.4\n")
  cat("Stop\n")
  sink()
  
  # Run AlphaRelate:
  system(command="$HOME/bin/AlphaRelate | tee AlphaRelate.log")
  
  # AlphaMate parameter file: 
  sink(file = "AlphaMateSpec.txt")
  cat("GenderFile                  , GenderRole.txt\n")
  cat("NrmMatrixFile               , PedigreeNrm.txt\n")
  cat("SelCriterionFile            , ebv.txt\n")
  cat("MateAllocation              , No\n")
  cat("EvolAlgStopTolerance        , 0.01\n")
  cat("NumberOfMatings             , ",dams,"\n", sep="")
  cat("NumberOfMaleParents         ,  ",sires,"\n", sep="")
  cat("EqualizeMaleContributions   , Yes\n")
  cat("NumberOfFemaleParents       , ",dams,"\n", sep="")
  cat("EqualizeFemaleContributions , Yes\n")
  cat("TargetDegree                ,  ",deg,"\n", sep="")
  cat("Stop\n")
  sink()
  
  # Run AlphaMate:
  system(command="$HOME/bin/AlphaMate | tee AlphaMate.log")
  
  #return()
}

# Input: candidates_plan & RecSys 
run_OCS_EASED <- function(ocspop, datafile, deg, dams) {
  # Create index of EBV for criterion:
  tmp = data.frame(id = ocspop@id, crit = (0.2*ocspop@ebv[,1]+0.35*ocspop@ebv[,2])+0.45*ocspop@ebv[,3])
  write.table(x = tmp, file = "ebv.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Assign sex:
  tmp = data.frame(id = ocspop@id, genderRole = ocspop@sex)
  tmp$genderRole = as.numeric(factor(ocspop@sex, levels=c("M", "F")))
  write.table(x = tmp, file = "GenderRole.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Full (All Animals) Pedigree: 
  tmp <- data.frame(datafile[,c("IId", "FId", "MId")])
  tmp$FId[is.na(tmp$FId)] <- 0
  tmp$MId[is.na(tmp$MId)] <- 0
  write.table(tmp, file = "Blupf90.ped", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
  # ID of the OCS candidates: 
  idcand = ocspop@id
  write.table(idcand, file = "PedigreeSet.ped", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0")
  
  # AlphaRelate parameter file: 
  sink(file = "AlphaRelateSpec.txt")
  cat("PedigreeFile                , Blupf90.ped\n")
  cat("PedInbreeding               , Yes\n")
  cat("PedNrm                      , Yes\n")
  cat("PedNrmSubsetFile            , PedigreeSet.ped\n")
  cat("OutputFormat                , f7.4\n")
  cat("Stop\n")
  sink()
  
  # Run AlphaRelate:
  system(command="$HOME/bin/AlphaRelate | tee AlphaRelate.log")
  
  # AlphaMate parameter file: 
  sink(file = "AlphaMateSpec.txt")
  cat("GenderFile                  , GenderRole.txt\n")
  cat("NrmMatrixFile               , PedigreeNrm.txt\n")
  cat("SelCriterionFile            , ebv.txt\n")
  cat("MateAllocation              , No\n")
  cat("EvolAlgStopTolerance        , 0.01\n")
  cat("NumberOfMatings             , ",dams,"\n", sep="")
  cat("NumberOfFemaleParents       , ",dams,"\n", sep="")
  cat("EqualizeFemaleContributions , Yes\n")
  cat("TargetDegree                ,  ",deg,"\n", sep="")
  cat("Stop\n")
  sink()
  
  # Run AlphaMate:
  system(command="$HOME/bin/AlphaMate | tee AlphaMate.log")
  
  #return()
}


# Input: candidates_plan & RecSys 
run_MinF <- function(ocspop, datafile, sires, dams) {
  # Create index of EBV for criterion:
  tmp = data.frame(id=ocspop@id, crit=(0.2*ocspop@ebv[,1]+0.35*ocspop@ebv[,2])+0.45*ocspop@ebv[,3])
  write.table(x=tmp, file="ebv.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # Assign sex:
  tmp = data.frame(id=ocspop@id, genderRole=ocspop@sex)
  tmp$genderRole = as.numeric(factor(ocspop@sex, levels=c("M", "F")))
  write.table(x=tmp, file="GenderRole.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # Full (All Animals) Pedigree: 
  tmp <- data.frame(datafile[,c("IId", "FId", "MId")])
  tmp$FId[is.na(tmp$FId)] <- 0
  tmp$MId[is.na(tmp$MId)] <- 0
  write.table(tmp, file="Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # ID of the OCS candidates: 
  idcand=ocspop@id
  write.table(idcand, file="PedigreeSet.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  
  # AlphaRelate parameter file: 
  sink(file = "AlphaRelateSpec.txt")
  cat("PedigreeFile                , Blupf90.ped\n")
  cat("PedInbreeding               , Yes\n")
  cat("PedNrm                      , Yes\n")
  cat("PedNrmSubsetFile            , PedigreeSet.ped\n")
  cat("OutputFormat                , f7.4\n")
  cat("Stop\n")
  sink()
  
  # Run AlphaRelate:
  system(command="$HOME/bin/AlphaRelate | tee AlphaRelate.log")
  
  # AlphaMate parameter file: 
  sink(file = "AlphaMateSpec.txt")
  cat("GenderFile                  , GenderRole.txt\n")
  cat("NrmMatrixFile               , PedigreeNrm.txt\n")
  cat("SelCriterionFile            , ebv.txt\n")
  cat("NumberOfMatings             , ",dams,"\n", sep="")
  cat("NumberOfMaleParents         ,  ",sires,"\n", sep="")
  cat("EqualizeMaleContributions   , Yes\n")
  cat("NumberOfFemaleParents       , ",dams,"\n", sep="")
  cat("EqualizeFemaleContributions , Yes\n")
  cat("ModeMinInbreeding           , Yes\n")
  cat("Stop\n")
  sink()
  
  # Run AlphaMate:
  system(command="$HOME/bin/AlphaMate | tee AlphaMate.log")
  
  #return()
}

# N.B. This function is redundant with new versions of AlphaSimR, but main code should be updated as well
run_expandCrossPlan <- function(nProgenyPerMat, CrossPlanFile) {
  for (set in 1:length(nProgenyPerMat)) {
    crp = matrix(ncol=2, nrow=nrow(CrossPlanFile) * nProgenyPerMat[set])
    k = 0
    for (i in 1:nrow(CrossPlanFile)) {
      for (j in 1:nProgenyPerMat[set]) {
        k = k + 1
        crp[k, 1] = CrossPlanFile[i, 2]
        crp[k, 2] = CrossPlanFile[i, 3]
      }
    }
  }
  crp=data.frame(crp)
  tmp = crp %>%
    mutate(sm=as.character(crp$X1)) %>%
    mutate(sf=as.character(crp$X2))
  crossPlanM=as.matrix(cbind(tmp$sf,tmp$sm))
  return(crossPlanM)
}


run_allocate_matings <- function(nProgenyPerMat, CrossPlanFile) {
  dfF = CrossPlanFile %>%  
    filter(Gender == 2) %>% 
    select("Id") %>% 
    rename(Id_F = Id) %>%
    dplyr::mutate(Id_F = as.character(Id_F))
  dfM = CrossPlanFile %>%  
    filter(Gender == 1) %>% 
    select("Id", "nContribution") %>% 
    rename(Id_M = Id) %>%  
    slice(rep(1:n(), times = nContribution) ) %>% 
    dplyr::mutate(Id_M = as.character(Id_M))
  crossPlanMF = cbind( slice(dfF, sample(1:n())),
                       slice(dfM, sample(1:n())) )
  crossPlanMF_expanded = crossPlanMF %>% 
    slice(rep(1:n(), each = nProgenyPerMat))
  crossPlanMF_final = as.matrix(cbind(crossPlanMF_expanded$Id_F, crossPlanMF_expanded$Id_M))
  return(crossPlanMF_final)
}


run_prepare <- function(datafile) {
  # Create pedigree file:
  blupPed <- datafile[,c("IId", "FId", "MId")]
  blupPed$FId[is.na(blupPed$FId)] <- 0
  blupPed$MId[is.na(blupPed$MId)] <- 0
  # blupPed <- blupPed[order(blupPed$IId),]
  # This ordering might not be correct, but doesn't matter since we run renumf90 later anyways
  write.table(blupPed, "Blupf90.ped", quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  # Delete pedigree - each time create a new pedigree:
  rm(blupPed)
  # Create phenotype file:
    # Add the overall mean vector for blupf90: rep(1,nrow(Pheno))
    # Only females have phenotypes, and animals in Generation 0 have no phenotypes
    # Females from the current selection group don't have Phenotype T3
      # gen and candidatesgroup
  Pheno = datafile %>%
    filter(Sex == "F" & Generation >= 1)
  Pheno = Pheno %>%
    dplyr::mutate(PhenoT3 = replace(PhenoT3, Generation == gen & SelCand == candidatesgroup, NA))
  write.table(cbind(Pheno[,c("IId", "PhenoT1", "PhenoT2", "PhenoT3", "Generation")], rep(1,nrow(Pheno))), "Blupf901.dat",
              quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ", na = "0")
  # Delete datafile - each time create a new datafile:
  rm(Pheno)
}

# Prepare genotype file in blupf90 format
# Extremely inefficient (copying) !!!
run_prepareGENO <- function(popname) {
  idg = popname@id 
  pr=pullSnpGeno(popname)
  prt = apply(pr[,1:ncol(pr)],1,paste,sep ="",collapse="")
  idgpr=cbind(idg, prt)
  write.table(idgpr,"mrk.tmp",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE,sep=" ")
  system(command = "sh prepgeno.sh")
}


# Create summary report table
PullSumm <- function(datafile, popname, sex) {
  gePa = genParam(popname)
  datafile = rbind(datafile,
                   tibble(Replicate      = rep,
                          Generation     = gen,
                          Program        = Program,
                          SelCand        = candidatesgroup,
                          Sex            = sex,
                          Year           = year,
                          NumberID       = popname@nInd,
                          MeanPhenoT1    = (meanP(popname))[1],
                          MeanPhenoT2    = (meanP(popname))[2],
                          MeanPhenoT3    = (meanP(popname))[3],
                          MeanGenoT1     = gePa$mu[1],
                          MeanGenoT2     = gePa$mu[2],
                          MeanGenoT3     = gePa$mu[3],
                          MeanA_T1       = mean(gePa$gv_a[, 1]),
                          MeanA_T2       = mean(gePa$gv_a[, 2]),
                          MeanA_T3       = mean(gePa$gv_a[, 3]),
                          MeanD_T1       = mean(gePa$gv_d[, 1]),
                          MeanD_T2       = mean(gePa$gv_d[, 2]),
                          MeanD_T3       = mean(gePa$gv_d[, 3]),
                          VarA_T1        = gePa$varA[1,1],
                          VarA_T2        = gePa$varA[2,2],
                          VarA_T3        = gePa$varA[3,3],
                          VarD_T1        = gePa$varD[1,1],
                          VarD_T2        = gePa$varD[2,2],
                          VarD_T3        = gePa$varD[3,3],
                          VarG_T1        = gePa$varG[1,1],
                          VarG_T2        = gePa$varG[2,2],
                          VarG_T3        = gePa$varG[3,3],
                          GenicVA_T1     = gePa$genicVarA[1],
                          GenicVA_T2     = gePa$genicVarA[2],
                          GenicVA_T3     = gePa$genicVarA[3],
                          GenicVD_T1     = gePa$genicVarD[1],
                          GenicVD_T2     = gePa$genicVarD[2],
                          GenicVD_T3     = gePa$genicVarD[3],
                          GenicVG_T1     = gePa$genicVarG[1],
                          GenicVG_T2     = gePa$genicVarG[2],
                          GenicVG_T3     = gePa$genicVarG[3]
                          ))
  #return(datafile)
}

# N.B. Huge function creating several metrics in summary report table
PullSummTogether <- function(datafile, popname) {
  # Calculate population statistics
  gePa = genParam(popname)
  # Get the True Index
  true_index = selIndex(popname@gv, b = iweight)
  ebv_index = selIndex(popname@ebv, b = iweight)
  # Extract Genotypes 
  M = pullSnpGeno(popname)
  N = pullSnpGeno(popname, snpChip = 2)
  QT1 = pullQtlGeno(popname, trait = 1) 
  QT2 = pullQtlGeno(popname, trait = 2) 
  QT3 = pullQtlGeno(popname, trait = 3) 
  # Calculate Observed and Expected Heterozygosity
  # For selected SNP
  # het_snp = sum(apply(X = M, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) ) )
  het_snp = apply(X = M, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
  p_snp = colMeans(M)/2
  # For neutral SNP
  het_neutral = apply(X = N, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
  p_neutral = colMeans(N)/2
  # For QTL
  het_qtl1 = apply(X = QT1, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
  het_qtl2 = apply(X = QT2, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
  het_qtl3 = apply(X = QT3, MARGIN = 2, FUN = function(X) (sum(X==1)/length(X)) )
  p_qtl1 = colMeans(QT1)/2
  p_qtl2 = colMeans(QT2)/2
  p_qtl3 = colMeans(QT3)/2
  # Param.
  nSNP = ncol(M)
  nQTL = ncol(QT1)
  LDB = NULL
  LDW = NULL
  for (traitN in 1:3) {
    ChrN = SP$nChr
    QtlN = unique(SP$traits[[traitN]]@lociPerChr)
    QTL_genos = sweep(pullQtlGeno(popname), 2,
                      SP$traits[[traitN]]@addEff, 
                      "*")
    QTL.matrix = popVar(QTL_genos)
    # chr.genetic = NULL
    # chr.genic = NULL 
    chr.LD.without = NULL
    chr.LD = NULL 
    for(c in 1:ChrN){
      i = ((c * QtlN) - QtlN) + 1
      j = c * QtlN
      # genetic.value <- sum(QTL.matrix[i:j, i:j])
      # genic.value <- sum(diag(QTL.matrix[i:j, i:j]))
      LD.without <- (sum(QTL.matrix[ ,i:j]) - sum(QTL.matrix[i:j, i:j]))
      matrix.chr <- popVar(QTL_genos[ ,i:j])
      LD.chr <- (sum(matrix.chr) - sum(diag(matrix.chr)))
      # chr.genetic[c] <- list(genetic.value)
      # chr.genic[c] <- list(genic.value)
      chr.LD.without[c] <- list(LD.without)
      chr.LD[c] <- list(LD.chr)
    }
    LDB[traitN] = list(sum(unlist(chr.LD)) )
    LDW[traitN] = list(sum(unlist(chr.LD.without)) )
  }
  
  # Accumulate all info into datafile
  datafile = rbind(datafile,
                   tibble(Replicate      = rep,
                          Generation     = gen,
                          Program        = Program,
                          SelCand        = candidatesgroup,
                          Year           = year,
                          NumberID       = popname@nInd,
                          MeanPhenoT1    = (meanP(popname))[1],
                          MeanPhenoT2    = (meanP(popname))[2],
                          MeanPhenoT3    = (meanP(popname))[3],
                          MeanGenoT1     = gePa$mu[1],
                          MeanGenoT2     = gePa$mu[2],
                          MeanGenoT3     = gePa$mu[3],
                          Accuracy_T1    = cor(gv(popname)[,1], ebv(popname)[,1], use = "pairwise.complete.obs"),
                          Accuracy_T2    = cor(gv(popname)[,2], ebv(popname)[,2], use = "pairwise.complete.obs"),
                          Accuracy_T3    = cor(gv(popname)[,3], ebv(popname)[,3], use = "pairwise.complete.obs"),
                          Accuracy_index = cor(true_index, ebv_index, use = "pairwise.complete.obs"),
                          MeanA_T1       = mean(gePa$gv_a[, 1]),
                          MeanA_T2       = mean(gePa$gv_a[, 2]),
                          MeanA_T3       = mean(gePa$gv_a[, 3]),
                          MeanD_T1       = mean(gePa$gv_d[, 1]),
                          MeanD_T2       = mean(gePa$gv_d[, 2]),
                          MeanD_T3       = mean(gePa$gv_d[, 3]),
                          VarA_T1        = gePa$varA[1,1],
                          VarA_T2        = gePa$varA[2,2],
                          VarA_T3        = gePa$varA[3,3],
                          VarD_T1        = gePa$varD[1,1],
                          VarD_T2        = gePa$varD[2,2],
                          VarD_T3        = gePa$varD[3,3],
                          VarG_T1        = gePa$varG[1,1],
                          VarG_T2        = gePa$varG[2,2],
                          VarG_T3        = gePa$varG[3,3],
                          GenicVA_T1     = gePa$genicVarA[1],
                          GenicVA_T2     = gePa$genicVarA[2],
                          GenicVA_T3     = gePa$genicVarA[3],
                          GenicVD_T1     = gePa$genicVarD[1],
                          GenicVD_T2     = gePa$genicVarD[2],
                          GenicVD_T3     = gePa$genicVarD[3],
                          GenicVG_T1     = gePa$genicVarG[1],
                          GenicVG_T2     = gePa$genicVarG[2],
                          GenicVG_T3     = gePa$genicVarG[3],
                          LDB_T1         = LDB[1],
                          LDB_T2         = LDB[2],
                          LDB_T3         = LDB[3],
                          LDW_T1         = LDW[1],
                          LDW_T2         = LDW[2],
                          LDW_T3         = LDW[3],
                          LDT_T1         = LDB[[1]] + LDW[[1]],
                          LDT_T2         = LDB[[2]] + LDW[[2]],
                          LDT_T3         = LDB[[3]] + LDW[[3]],
                          covG_L_T1      = gePa$covG_L[1],
                          covG_L_T2      = gePa$covG_L[2],
                          covG_L_T3      = gePa$covG_L[3],
                          LD_T1          = gePa$varG[1,1] - gePa$genicVarG[1],
                          LD_T2          = gePa$varG[2,2] - gePa$genicVarG[2],
                          LD_T3          = gePa$varG[3,3] - gePa$genicVarG[3],
                          covG_HW_T1     = gePa$covG_HW[1],
                          covG_HW_T2     = gePa$covG_HW[2],
                          covG_HW_T3     = gePa$covG_HW[3],
                          Mean_index     = mean(true_index),
                          Var_index      = var(true_index) * ( (length(true_index) - 1) / length(true_index) ),
                          MarkerHet      = sum(het_snp),
                          NeutralHet     = sum(het_neutral),
                          QTLHet_T1      = sum(het_qtl1),
                          QTLHet_T2      = sum(het_qtl2),
                          QTLHet_T3      = sum(het_qtl3),
                          F_snp_obs      = 1 - ( sum(het_snp) / nSNP) ,
                          F_snp_exp      = 1 - (sum(2*p_snp*(1-p_snp)) / nSNP),
                          F_neutral_obs  = 1 - ( sum(het_neutral) / nSNP) ,
                          F_neutral_exp  = 1 - (sum(2*p_neutral*(1-p_neutral)) / nSNP),
                          F_QTL1_obs     = 1 - ( sum(het_qtl1) / nQTL),
                          F_QTL2_obs     = 1 - ( sum(het_qtl2) / nQTL),
                          F_QTL3_obs     = 1 - ( sum(het_qtl3) / nQTL),
                          F_QTL1_exp     = 1 - (sum(2*p_qtl1*(1-p_qtl1)) / nQTL),
                          F_QTL2_exp     = 1 - (sum(2*p_qtl2*(1-p_qtl2)) / nQTL),
                          F_QTL3_exp     = 1 - (sum(2*p_qtl3*(1-p_qtl3)) / nQTL)
                   ))
  # return(datafile)
}



PullSireContribution <- function(datafile, CrossPlanFile) {
  dfM = CrossPlanFile %>%  
    filter(Gender == 1) %>% 
    select("Id", "nContribution")
  datafile = rbind(datafile,
                   tibble(Replicate      = rep,
                          Generation     = gen,
                          Program        = Program,
                          SelCand        = candidatesgroup,
                          Year           = year,
                          N_Sires        = nrow(dfM),
                          Mean_Contri    = mean(dfM$nContribution),
                          Min_Contri     = min(dfM$nContribution),
                          Max_Contri     = max(dfM$nContribution)
                   ))
  #return(datafile)
}


prepare_par <- function() {

# extractsol.sh can be rather replaced by R-function
sink("extractsol.sh", type="output")
writeLines("#!/bin/bash
awk '{if ($1==1 && $2==2) print $3,$4}' solutions | sort +0 -1 > sol1.tmp
awk '{if ($1==2 && $2==2) print $3,$4}' solutions | sort +0 -1 > sol2.tmp
awk '{if ($1==3 && $2==2) print $3,$4}' solutions | sort +0 -1 > sol3.tmp
awk '{print $1, $10}' renadd02.ped | sort +0 -1 > id.tmp
join -1 +1 -2 +1 id.tmp sol1.tmp > sl1.tmp
join -1 +1 -2 +1 sl1.tmp sol2.tmp > sl2.tmp
join -1 +1 -2 +1 sl2.tmp sol3.tmp > sl3.tmp
# rm *.tmp
")
sink()


sink("prepgeno.sh", type="output")
writeLines("#!/bin/bash
# Prepare Marker Genotypes file for blupf90 format
# And automatically cut old markers for blupf90 limit (25k)

# Don't forget escaping quotes and backslashes when in R :)

# Blupf90 fixed format
awk '{printf(\"%-8s %s\\n\", $1,i$2)}' mrk.tmp > mrk2.tmp

gen=$(< mrk2.tmp wc -l)
# hard limit = 25000
limit=10200

if [ \"$gen\" -gt \"$limit\" ]
 then
   cold=(`expr $gen - $limit + 1`)
   tail -n +$cold mrk2.tmp > markeri.dat
 else
   cp mrk2.tmp markeri.dat
fi
")
sink()


sink("renum.par", type="output")
writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2 3 4
#Phen
FIELDS_PASSED TO OUTPUT
1 5
#id generation
WEIGHT(S)

RESIDUAL_VARIANCE
4.5556 0 0
0 3.5455 0
0 0 3.0
EFFECT
6 6 6 cross alpha
#mean
EFFECT
1 1 1 cross alpha
#animal
RANDOM
animal
OPTIONAL

FILE
Blupf90.ped
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
1.00 0.75 0.60
0.75 1.00 0.70
0.60 0.70 1.00
OPTION blksize 3
")
sink()


sink("renumGS.par", type="output")
writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2 3 4
#Phen
FIELDS_PASSED TO OUTPUT
1 5
#id generation
WEIGHT(S)

RESIDUAL_VARIANCE
4.5556 0 0
0 3.5455 0
0 0 3.0
EFFECT
6 6 6 cross alpha
#mean
EFFECT
1 1 1 cross alpha
#animal
RANDOM
animal
OPTIONAL

FILE
Blupf90.ped
SNP_FILE
../markeri.dat
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
1.00 0.75 0.60
0.75 1.00 0.70
0.60 0.70 1.00
OPTION blksize 3
OPTION msg 1
OPTION no_quality_control
OPTION use_yams
OPTION thrStopCorAG 0
")
sink()

}
