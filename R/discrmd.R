#  Copyright (c) 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
#  All rights reserved.
#
#  This file is part of the psm3mkv program.
#
#  psm3mkv is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==================================================================
# Discretized restricted mean durations
# discrmd.R
# =====================================
#

#' Discretized Restricted Mean Duration calculation for Partitioned Survival Model
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a Partitioned Survival Model structure.
#' @param ptdata Dataset of patient level data. Must be a tibble with columns named:
#' - ptid: patient identifier
#' - pfs.durn: duration of PFS from baseline
#' - pfs.flag: event flag for PFS (=1 if progression or death occurred, 0 for censoring)
#' - os.durn: duration of OS from baseline
#' - os.flag: event flag for OS (=1 if death occurred, 0 for censoring)
#' - ttp.durn: duration of TTP from baseline (usually should be equal to pfs.durn)
#' - ttp.flag: event flag for TTP (=1 if progression occurred, 0 for censoring).
#'
#' Survival data for all other endpoints (time to progression, pre-progression death, post-progression survival) are derived from PFS and OS.
#' @param dpam List of survival regressions for model endpoints. These must include time to progression (TTP) and pre-progression death (PPD).
#' @param psmtype Either "simple" or "complex" PSM formulation
#' @param Ty Time duration over which to calculate (defaults to 10 years). Assumes input is in years, and patient-level data is recorded in weeks.
#' @param discrate Discount rate (%) per year (defaults to zero).
#' @param lifetable Optional. The lifetable must be a dataframe with columns named time and lx. The first entry of the time column must be zero. Data should be sorted in ascending order by time, and all times must be unique.
#' @param timestep Optional, defaults to one (week).
#' @return List containing:
#' - pf: RMD in PF state
#' - pd: RMD in PD state
#' - os: RMD in either alive state
#' @seealso [drmd_stm_cf()] [drmd_stm_cr()]
#' @noRd
# Examples
# # Create dataset and fit survival models (splines)
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_spl(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#   os = find_bestfit_spl(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
# )
# drmd_psm(ptdata=bosonc, dpam=params)
# # Add a lifetable constraint
# ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
# drmd_psm(ptdata=bosonc, dpam=params, lifetable=ltable)
drmd_psm <- function(ptdata, dpam, psmtype="simple", Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- pfprob <- osprob <- adjosprob <- adjfac <- adjprob <- vn <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(0:floor(Tw/timestep)) + timestep/2
  # Obtain all the hazards
  allh <- calc_haz_psm(timevar=tvec, ptdata=ptdata, dpam=dpam, psmtype=psmtype)$adj
  # PFS and OS probabilties from PSM
  pfprob <- prob_pf_psm(tvec, dpam)
  osprob <- prob_os_psm(tvec, dpam)
  # OS and PFS may be constrained already by definitions of PPD and PPS
  maxos <- exp(-cumsum(allh$os))
  maxpfs <- exp(-cumsum(allh$pfs))
  # Further constrain OS by lifetable
  conos <- constrain_survprob(survprob1=maxos, survprob2=osprob, lifetable=lifetable, timevec=tvec)
  conpfs <- constrain_survprob(survprob1=maxpfs, survprob2=pmin(pfprob, osprob), lifetable=lifetable, timevec=tvec)
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(conpfs*vn) * timestep
  os <- sum(conos*vn) * timestep
  # Return values
  return(list(pf=pf, pd=os-pf, os=os))
}

#' Discretized Restricted Mean Duration calculation for State Transition Model Clock Forward structure
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a State Transition Model Clock Forward structure.
#' @inherit drmd_psm params return
#' @seealso [drmd_psm()] [drmd_stm_cr()]
#' @noRd
# Examples
# # Create dataset and fit survival models (splines)
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_spl(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#   os = find_bestfit_spl(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
# )
# drmd_stm_cf(dpam=params)
# # Add a lifetable constraint
# ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
# drmd_stm_cf(dpam=params, lifetable=ltable)
drmd_stm_cf <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- sppd <- sttp <- sos <- NULL
  adjsppd <- adjos <- vn <- pf <- os <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(0:floor(Tw/timestep)) + timestep/2
  # Obtain hazard and survival functions
  hppd <- tvec |> purrr::map_dbl(~calc_haz(.x, survobj=dpam$ppd))
  http <- tvec |> purrr::map_dbl(~calc_haz(.x, survobj=dpam$ttp))
  hpps <- tvec |> purrr::map_dbl(~calc_haz(.x, survobj=dpam$pps_cf))
  sppd <- tvec |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$ppd))
  sttp <- tvec |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$ttp))
  spps <- tvec |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$pps_cf))
  # Derive the constrained S_PPD and S_PFS
  con.sppd <- constrain_survprob(sppd, lifetable=lifetable, timevec=tvec)
  con.spps <- constrain_survprob(spps, lifetable=lifetable, timevec=tvec)
  # Partial prob PD
  con.partprobpd <- sttp*con.sppd*http/con.spps
  con.partprobpd[con.spps==0] <- 0
  con.probpd <- con.spps * cumsum(con.partprobpd)
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(sttp*con.sppd*vn) * timestep
  pd <- sum(con.probpd*vn) * timestep
  # Return values
  return(list(pf=pf, pd=pd, os=pf+pd))
}

#' Discretized Restricted Mean Duration calculation for State Transition Model Clock Reset structure
#' Calculate restricted mean duration (RMD) in PF, PD and OS states under a State Transition Model Clock Reset structure.
#' @inherit drmd_psm params return
#' @seealso [drmd_stm_cf()] [drmd_psm()]
#' @noRd
# Examples
# # Create dataset and fit survival models (splines)
# bosonc <- create_dummydata("flexbosms")
# fits <- fit_ends_mods_spl(bosonc)
# # Pick out best distribution according to min AIC
# params <- list(
#   ppd = find_bestfit_spl(fits$ppd, "aic")$fit,
#   ttp = find_bestfit_spl(fits$ttp, "aic")$fit,
#   pfs = find_bestfit_spl(fits$pfs, "aic")$fit,
#   os = find_bestfit_spl(fits$os, "aic")$fit,
#   pps_cf = find_bestfit_spl(fits$pps_cf, "aic")$fit,
#   pps_cr = find_bestfit_spl(fits$pps_cr, "aic")$fit
# )
# drmd_stm_cr(dpam=params)
# # Add a lifetable constraint
# ltable <- tibble::tibble(lttime=0:20, lx=1-lttime*0.05)
# drmd_stm_cr(dpam=params, lifetable=ltable)
drmd_stm_cr <- function(dpam, Ty=10, discrate=0, lifetable=NA, timestep=1) {
  # Declare local variables
  Tw <- tvec <- sppd <- sttp <- sos <- NULL
  adjsppd <- adjos <- vn <- pf <- os <- NULL
  # Time horizon in weeks (ceiling)
  Tw <- convert_yrs2wks(Ty)
  # Create time vector, with half-cycle addition
  tvec <- timestep*(0:floor(Tw/timestep)) + timestep/2
  # Obtain unconstrained survival functions
  sppd <- tvec |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$ppd))
  sttp <- tvec |> purrr::map_dbl(~calc_surv(.x, survobj=dpam$ttp))
  # Derive the constrained S_PPD
  c.sppd <- constrain_survprob(sppd, lifetable=lifetable, timevec=tvec)
  # Integrand with constraints on S_PPD and S_PPS
  integrand <- function(u, t) {
    i.http <- calc_haz(u, survobj=dpam$ttp)
    i.sttp <- calc_surv(u, survobj=dpam$ttp)
    i.u.sppd <- calc_surv(u, survobj=dpam$ppd)
    i.u.spps <- calc_surv(t-u, survobj=dpam$pps_cr)
    i.slxu <- calc_ltsurv(convert_wks2yrs(u), lifetable=lifetable)
    i.slxt <- calc_ltsurv(convert_wks2yrs(t), lifetable=lifetable)
    i.c.sppd <- pmin(i.u.sppd, i.slxu)
    i.c.spps <- pmin(i.u.spps, i.slxt/i.slxu)
    i.c.spps[i.slxu==0] <- 0
    i.c.sppd * i.sttp * i.http * i.c.spps
  }
  integrand <- Vectorize(integrand, "u")
  # PD membership probability is the integral
  probpd <- function(t) {
    stats::integrate(integrand, lower=0, upper=t, t=t)$value
  }
  probpd <- Vectorize(probpd, "t")
  # Calculate the PD membership probability for each time
  c.probpd <- probpd(tvec)
  # Discount factor
  vn <- (1+discrate)^(-convert_wks2yrs(tvec+timestep/2))
  # Calculate RMDs
  pf <- sum(sttp*c.sppd*vn) * timestep
  pd <- sum(c.probpd*vn) * timestep
  # Return values
  return(list(pf=pf, pd=pd, os=pf+pd))
}

