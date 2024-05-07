## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----packages, message=FALSE--------------------------------------------------
library("dplyr")
library("psm3mkv")
library("tibble")

## ----ltable1, eval=FALSE------------------------------------------------------
#  # HMD HFD package
#  library("HMDHFDplus")
#  
#  # Assumed population age at baseline (time=0)
#  baseage <- 51.0
#  
#  # Mortality data - England & Wales, Female, 2019
#  mort <- readHMDweb(
#    CNTRY = "GBRTENW",
#    item = "fltper_1x1",
#    username = "",
#    password = ""
#  ) |>
#    filter(Year == 2019) |>
#    select(Age, lx, dx) |>
#    mutate(Timey = ceiling(Age - baseage)) |>
#    filter(Timey >= 0)
#  
#  # The table needs to end with lx=0
#  mort <- add_row(mort, Age = 111, lx = 0, Timey = 60)

## ----ltable2------------------------------------------------------------------
mort <- tibble(
  Timey = 0:30,
  lx = 10000 * exp(-0.03 * Timey^1.1),
  dx = lx - lead(lx),
  qx = dx / lx
)

## ----ltable3------------------------------------------------------------------
# Assumed Standardized Mortality Ratio
SMR <- 2

# Recalculate the lifetable with the SMR applied
mort$adjlx <- mort$lx
for (i in 2:length(mort$lx)) {
  mort$adjlx[i] <- mort$adjlx[i - 1] * (1 - SMR * mort$qx[i - 1])
}

# Ensure lx>=0
mort$adjlx[mort$adjlx < 0] <- 0

# Create and view lifetable
ltable <- tibble(lttime = mort$Timey, lx = mort$adjlx)
head(ltable)

## ----datafit------------------------------------------------------------------
# Get some data
bosonc <- create_dummydata("flexbosms")

# We'll make the durations a lot longer
# so that they will be definitely constrained by the lifetable
bosonc <- bosonc |>
  mutate(
    pfs.durn = 20 * pfs.durn,
    os.durn = 20 * os.durn,
    ttp.durn = 20 * ttp.durn
  )

# Fit parametric models to each endpoint
allfits <- fit_ends_mods_par(bosonc)

# Bring together preferred fits for each endpoint
params <- list(
  ppd = find_bestfit(allfits$ppd, "aic")$fit,
  ttp = find_bestfit(allfits$ttp, "aic")$fit,
  pfs = find_bestfit(allfits$pfs, "aic")$fit,
  os = find_bestfit(allfits$os, "aic")$fit,
  pps_cf = find_bestfit(allfits$pps_cf, "aic")$fit,
  pps_cr = find_bestfit(allfits$pps_cr, "aic")$fit
)

## ----proj1--------------------------------------------------------------------
# Set time horizon
thoz <- 10

# Run the calculations
proj1 <- calc_allrmds(bosonc, dpam = params, Ty = thoz)

# Present the results
res1 <- proj1$results |> mutate(method = "int", lxadj = "no", psmmeth = "simple")
res1 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)

## ----proj2--------------------------------------------------------------------
# Run the calculation
proj2 <- calc_allrmds(bosonc, dpam = params, Ty = thoz, rmdmethod = "disc")

# Present the results
res2 <- proj2$results |>
  mutate(method = "disc", lxadj = "no", psmmeth = "simple")
res2 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)

## ----proj3--------------------------------------------------------------------
# Run the calculations
proj3 <- calc_allrmds(bosonc, dpam = params, Ty = thoz, rmdmethod = "disc", lifetable = ltable)

# Present the results
res3 <- proj3$results |>
  mutate(method = "disc", lxadj = "yes", psmmeth = "simple")
res3 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)

# Proportion of time alive spent in PF
proppf3 <- res3$pf / res3$os

## ----proj4--------------------------------------------------------------------
# Run the calculations
proj4 <- calc_allrmds(bosonc, dpam = params, Ty = thoz, psmtype = "complex", rmdmethod = "disc", lifetable = ltable)

# Present the results
res4 <- proj4$results |>
  mutate(method = "disc", lxadj = "yes", psmmeth = "complex")
res4 |> mutate(
  across(c(pf, pd, os), ~ num(.x, digits = 1))
)

# Proportion of time alive spent in PF
proppf4 <- res4$pf / res4$os

