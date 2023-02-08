# Assessment of long-term trends in genetic mean and variance after the introduction of genomic selection in layers: a simulation study
## Pocrnic, Obsteter, Gaynor, Wolc & Grojanc (2023)

- `functions.R` include various support functions needed for the simulation
- `1_Burnin_and_PTS.R` Founders / Burn-in Truncation BLUP / Standard Truncation BLUP
- `2_GTS.R` Truncation ssGBLUP
- `3_MFGTS` Truncation ssGBLUP with `ModeMinInbreeding`
- `4_GOCS:` OCS on both male and female selection candidates, without `MateAllocation`, and `EvolAlgStopTolerance = 0.01`, NRM in OCS
- `5_UGOCS.R:` OCS on both male and female selection candidates, no constrains on the number of males, without `MateAllocation`, and `EvolAlgStopTolerance = 0.01`, NRM in OCS
- `6_Random.R` Scenario with random selection 
- Each script reads 2 parameters: replicate number and number of sires (in our case 40 or 120). For example, `Rscript 2_GTS.R 3 40` 
- Please note that `AlphaMate` program (Gorjanc & Hickey (2018), The Roslin Institute , UK) and `AlphaRelate` programs have to be installed on your system. For download and installation details, please consult [AlphaGenes Github](https://github.com/AlphaGenes). 
- Please note that `BLUPF90` and `RENUMF90` programs (Ignacy Misztal et al., University of Georgia, USA) have to be installed on your system. For download and installation details, please consult [BLUPF90 website](http://nce.ads.uga.edu/software/). 
- You might need to change the following lines of code in `functions.R` based on your OS and location of the binaries of aforementioned programs. For example:
```
 system(command = "echo renf90.par | $HOME/bin/blupf90 | tee blup.log")
```
or
```
 system(command = "$HOME/bin/AlphaMate | tee AlphaMate.log")
```


