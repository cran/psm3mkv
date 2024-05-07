## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----packages, message=FALSE--------------------------------------------------
library("boot")
library("ggsci")
library("psm3mkv")
library("purrr")

## ----dataset------------------------------------------------------------------
# Create and review the dummy dataset
bosonc <- create_dummydata("flexbosms")
head(bosonc)
summary(bosonc)

## ----fit_ends_mods_par--------------------------------------------------------
# Create a vector of distributions of interest (flexsurv notation)
alldists <- c("exp", "weibullPH", "llogis", "lnorm", "gamma", "gompertz", "gengamma")

# Fit all distributions to all endpoints (except gengamma to PPD and TTP)
allfits_par <- fit_ends_mods_par(
  bosonc,
  cuttime = 0,
  ppd.dist = alldists[1:6],
  ttp.dist = alldists[1:6],
  pfs.dist = alldists,
  os.dist = alldists,
  pps_cf.dist = alldists,
  pps_cr.dist = alldists
)

# Example 1 - PFS endpoint, distribution 2 (weibullPH)
allfits_par$pfs[[2]]$result

# Example 2 - Parameter values for PPS-CF and PPS-CR endpoints for distribution 3 (llogis)
allfits_par$pps_cf[[3]]$result$res
allfits_par$pps_cr[[3]]$result$res

## ----find_bestfit_par---------------------------------------------------------
# Pick out best distribution according to min AIC
fitpar.ppd <- find_bestfit(allfits_par$ppd, "aic")
fitpar.ttp <- find_bestfit(allfits_par$ttp, "aic")
fitpar.pfs <- find_bestfit(allfits_par$pfs, "aic")
fitpar.os <- find_bestfit(allfits_par$os, "aic")
fitpar.pps_cf <- find_bestfit(allfits_par$pps_cf, "aic")
fitpar.pps_cr <- find_bestfit(allfits_par$pps_cr, "aic")

# Inspect the selection for PFS
fitpar.pfs

## ----fit_ends_mods_spl--------------------------------------------------------
# Fit 1-3 knot splines with all 3 scales (odds, hazard, normal) to each endpoint
allfits_spl <- fit_ends_mods_spl(bosonc)

# Example - PFS endpoint - 1 knot, odds scale
allfits_spl$pfs[[2]]$result
allfits_spl$pfs[[2]]$result$aux$scale # Scale
allfits_spl$pfs[[2]]$result$aux$knots # Knot locations (log time)

## ----find_bestfit-------------------------------------------------------------
# Pick out best distribution according to min AIC
fitspl.ppd <- find_bestfit(allfits_spl$ppd, "aic")
fitspl.ttp <- find_bestfit(allfits_spl$ttp, "aic")
fitspl.pfs <- find_bestfit(allfits_spl$pfs, "aic")
fitspl.os <- find_bestfit(allfits_spl$os, "aic")
fitspl.pps_cf <- find_bestfit(allfits_spl$pps_cf, "aic")
fitspl.pps_cr <- find_bestfit(allfits_spl$pps_cr, "aic")

# Inspect the selection for PFS
fitspl.pfs

## ----bring_together_fits------------------------------------------------------
# Bring together our preferred fits for each endpoint in a list
params <- list(
  ppd = fitpar.ppd$fit,
  ttp = fitpar.ttp$fit,
  pfs = fitspl.pfs$fit,
  os = fitspl.os$fit,
  pps_cf = allfits_par$pps_cf[[2]]$result,
  pps_cr = allfits_spl$pps_cr[[2]]$result
)

## ----count_params-------------------------------------------------------------
# Pull out number of parameters used for each endpoint
count_npar <- map_vec(1:6, ~ params[[.x]]$npars)

# PSM uses PFS (3) and OS (4) endpoints
sum(count_npar[c(3, 4)])

# STM_CF uses PPD (1), TTP (2) and PPS_CF (5) endpoints
sum(count_npar[c(1, 2, 5)])

# STM_CR uses PPD (1), TTP (2) and PPS_CR (6) endpoints
sum(count_npar[c(1, 2, 6)])

## ----calc_likes---------------------------------------------------------------
ll_all <- calc_likes(bosonc, params)
ll_all

## ----calc_allrmds-------------------------------------------------------------
# Call the RMD functions
rmd_all <- calc_allrmds(bosonc, dpam = params)

# Then review the mean duration in PF, PD and total alive (OS)
rmd_all$results

## ----boot---------------------------------------------------------------------
# Bootstrap to calculate SE over 10 bootstrap samples
boot::boot(
  data = bosonc,
  statistic = calc_allrmds,
  R = 10, # Number of samples
  cuttime = 0,
  Ty = 10,
  dpam = params,
  boot = TRUE
)

## ----graph_gen----------------------------------------------------------------
# Generate graphs (can take time)
ptdgraphs <- graph_survs(bosonc, params)

## ----graph_pf-----------------------------------------------------------------
# State membership probabilities for PF state
ptdgraphs$graph$pf + scale_color_npg()

## ----graph_pd-----------------------------------------------------------------
# State membership probabilities for PD state
ptdgraphs$graph$pd + scale_color_npg()

## ----graph_os-----------------------------------------------------------------
# State membership probabilities for OS
ptdgraphs$graph$os + scale_color_npg()

## ----graph_pps----------------------------------------------------------------
# Probabilities of PPS
ptdgraphs$graph$pps + scale_color_npg()

