#
# Fsquared - estimation of F reference points for the ICES advice rule
# Authors: participants of workshop on the validation of new tool for refpts
#	  Iago Mosqueira (WMR) <iago.mosqueira@wur.nl>
#	  Ernesto Jardim (IPMA) <ernesto.jardim@ipma.pt>
#	  John Trochta (IMR) <john.tyler.trochta@hi.no>
#	  Arni Magnusson (SPC) <arnim@spc.int>
#	  Max Cardinale (SLU) <massimiliano.cardinale@slu.se>
#	  Dorleta Garcia (ICES) <dorleta.garcia@ices.dk>
#	  Colin Millar (ICES) <colin.millar@ices.dk>
#
# Distributed under the terms of the EUPL 1.2
#===============================================================================
setwd("D:\\ARS_18_19_20_github")
rm(list=ls())

#===============================================================================
# # CLUSTER SETTINGS FOR DIPC
#===============================================================================

cluster <-  FALSE #TRUE

# for cluster change library locations:

ref.dir <- getwd()
if (cluster){ 
  
  maindir <- paste0("/", strsplit(ref.dir, "/")[[1]][2])
  usnam   <- Sys.getenv("USER")
  lib.dir <- file.path(maindir, usnam, "rpackages")
  .libPaths( c(lib.dir, .libPaths()) )
  # print(.libPaths())
  
}

#===============================================================================
# PACKAGES
#===============================================================================

# if needed install packages from r universe
#install.packages(c("FLCore", "mse", "msemodules", "mseviz", "FLSRTMB", "ggplotFL", "FLasher", "TAF"), repos=c(FLR="https://flr.r-universe.dev", CRAN="https://cloud.r-project.org")) 
library(devtools)
install_github('flr/Flasher', ref='devel')

# load libraries
library(mse)
library(msemodules)
library(FLSRTMB)  # to estimate SR parameters and add estimation uncertainty
library(TAF)
source("utilities_ARA_ARS.R")
# source("function_tacis2.R") # now included in utilities_ARA_ARS.R

#===============================================================================
# DIRECTORIES
#===============================================================================

if (cluster){
  id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  run.name <- list.dirs("data", recursive = FALSE, full.names = FALSE)[id]
  print(run.name)
  print(id)
}else{
  run.name <- "OM_steepnees_low"
}
run.dir <- file.path("data", run.name)
out.dir <- file.path("model", run.name)

if (!dir.exists(out.dir)){
  dir.create(out.dir, recursive=T)
}

#===============================================================================
# LOAD DATA
#===============================================================================

# load data created previously using the "data.R" script, 'run' is output FLStock
load(file.path(run.dir, "data.rda"))
# ls()

#===============================================================================
# SETUP
#===============================================================================

# Name of the stock
stkname <- "ARS_18-20"
# Recruitment models to be used in the OM conditioning
srmodels <- c("bevholtss3")
# Years to be used to compute stock-recruitment model
spryrs <- ac(2003:2024)
# Initial year of projections
iy <- dims(run)$maxyear
# Data lag
dl <- 2
# Management lag
ml <- 1
# Assessment frequency
af <- 1
# Data year
dy <- iy - dl
# Final year
fy <- iy + 30
# Probability years
pys <- seq(fy - 9, fy)
# How many years from the past to condition the future
conditioning_ny <- 3
# CV for SSB to add uncertainty in the shortcut estimator. Generated in data.R
# bcv_sa
# CV for F to add uncertainty in the shortcut estimator. Generated in data.R
# fcv_sa
# Years for geometric mean in short term forecast
recyrs_mp <- -3
# Btgt. Generated in data.R
# Btgt

# Blim and Btrigger
Blim <- Btgt*0.25
Btrigger <- Btgt*0.5
refpts <- FLPar(c(Blim = Blim, Btrigger = Btrigger))

# F search grid
fg_mp <- seq(0, 1, 0.05)

# Number of iterations (minimum of 50 for testing, 250 for final)
it <- 250

# Random seed
set.seed(987)

#===============================================================================
# SET UP CORES FOR PARALLEL
#===============================================================================

# CHOOSE no. of cores to use in parallel
cores <- 10

future::plan(future::multisession, workers = cores)

options(future.rng.onMisuse = "ignore")

#===============================================================================
# OM conditioning
#===============================================================================

# Stock-recruitment relationship(s)

# The file utilities_fsquared.R has code examples to condition the OM
# using stock recruitment parameters estimated by the stock assessment
# model, like SS3, SAM and a4a.

srpars <- params(srr)
srpars <- propagate(srpars, it)

tmp <- rlnorm(it, logR0, sdlogR0)
srpars["R0"] <- tmp

# GENERATE future deviances: lognormal autocorrelated
srdevs <- rlnormar1(it, sdlog=0.5, rho=0, years=seq(dy, fy), bias.correct=TRUE)

# BUILD FLom, OM FLR object
# om <- FLom(stock=propagate(run, it), refpts=refpts, model="bevholt", params=srpars, deviances=srdevs, name=stkname)

om <- FLom(stock=propagate(run, it), refpts=refpts, model="bevholtss3", params=srpars, deviances=srdevs, name=stkname)

# SETUP om future: average of most recent years set by conditioning_ny
om <- fwdWindow(om, end=fy, nsq=conditioning_ny)

#===============================================================================
# MP
#===============================================================================

# SET intermediate year + start of runs, lags and frequency
mseargs <- list(iy=iy, fy=fy, data_lag=dl, management_lag=ml, frq=af)

# SET shortcut estimator uncertainty: F and SSB deviances and auto-correlation
# Note your SSB deviances and auto-correlation have very little impact on the P(SB<Blim)
sdevs <- shortcut_devs(om, SSBcv=bcv_sa, Fcv=fcv_sa, Fphi=0)

# SETUP standard GFCM HCR (based on ICES advice rule formulation)

arule <- mpCtrl(
  
  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
                args=list(SSBdevs=sdevs$SSB)),
  
  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
                args=list(lim=0, trigger=refpts(om)$Btrigger, target=1, min=0, 
                          metric="ssb", output="fbar")),
  
  # (i)mplementation (sys)tem: tac.is (C ~ F)
  isys = mseCtrl(method=tac.is2, args=list(recyrs=recyrs_mp, Fdevs=sdevs$F, initac=unitSums(catch(om)[,ac(iy)])))
)

# RUN test 

# kk <- mp(om, ctrl=arule, args=mseargs, parallel=F) # OK
# plot(kk)

#===============================================================================
# Run simulations
#===============================================================================

# print message
message(paste("Starting simulations", Sys.time()))

# to run in SEQUENTIAL mode (undo previous settings for running in parallel)
# future::plan(future::sequential)

fgrid <- list() # this is not necessary, we do it for compatibility with the storage of other rules

# run over Ftarget grid
tmp <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp), names=paste0("GFCM_F_", fg_mp))

# compute performance statistics over pys
# performance(tmp) <- performance(tmp, statistics=gfcm_stats, year=pys, type="arule")
# performance(tmp)[,"mp"] <- paste0("GFCM")
# performance(tmp)[,"om"] <- run.name

# add results to list
fgrid[[paste0(run.name,"_GFCM")]] <- tmp

# SAVE
save(fgrid, om, run, file=file.path(out.dir, "fsquared.rdata"))

#===============================================================================
# End simulations
#===============================================================================

