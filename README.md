# Assessment of long-term trends in genetic mean and variance after the introduction of genomic selection in layers: a simulation study
## Pocrnic, Obsteter, Gaynor, Wolc & Grojanc (2023)
## DOI: 10.3389/fgene.2023.1168212 
## [Link to published paper](https://www.frontiersin.org/articles/10.3389/fgene.2023.1168212/full)

- `functions.R` Various support functions needed for the simulation
- `1_Burnin_and_PTS.R` Generating base population genomes; 10 years of  conventional truncation selection on BLUP and random mating with equal contributions (PTS); 20 years of PTS
- `2_GTS.R` 20 years of genomic truncation selection based on ssGBLUP and random mating with equal contributions (GTS)
- `3_MFGTS` 20 years of GTS with optimized mating minimizing future progeny inbreeding with equal contributions (MFGTS)
- `4_GOCS:` 20 years of genomic optimal contribution selection based on ssGBLUP with a constrained number of sires and random pairing of the optimized contributions (GOCS)
- `5_UGOCS.R:` 20 years of genomic optimal contribution selection based on ssGBLUP with an unconstrained number of sires and random pairing of the optimized contributions (UGOCS)
- `6_Random.R` Random selection program (negative control)
- Each script reads 2 input parameters: replicate number and number of sires (in our case 40 or 120). For example, `Rscript 2_GTS.R 3 40`. In addition, GOCS and UGOCS scripts need the trigonometric penalty degrees as well. For example, `Rscript 4_GOCS 3 40 65`.
- Please note that `AlphaMate` (Gorjanc & Hickey (2018), The Roslin Institute , UK) and `AlphaRelate` programs have to be installed on your system. For download and installation details, please consult [AlphaGenes Github](https://github.com/AlphaGenes). 
- Please note that `BLUPF90` and `RENUMF90` programs (Ignacy Misztal et al., University of Georgia, USA) have to be installed on your system. For download and installation details, please consult [BLUPF90 website](http://nce.ads.uga.edu/software/). 
- You might need to change some lines of code in `functions.R` based on your OS and location of the binaries of aforementioned programs. For example:
```
 system(command = "echo renf90.par | $HOME/bin/blupf90 | tee blup.log")
```
or
```
 system(command = "$HOME/bin/AlphaMate | tee AlphaMate.log")
```
- N.B. simulations take a long time to run and we advise use of server/HPC to run them.
