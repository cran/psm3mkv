#  Copyright (c) 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates.
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
# ======================================================================
#
# Datasets.R
# These functions are used to create dummy datasets to illustrate package use
# create_dummydata
# create_dummydata_survcan
# create_dummydata_pharmaonc
# create_dummydata_flexbosms
#
# ======================================

#' Create dummy dataset for illustration
#' @description Create dummy dataset to illustrate [psm3mkv]
#' @param dsname Dataset name, as follows:
#' * `flexbosms` provides a dataset based on [flexsurv::bosms3()]. This contains all the fields necessary for [psm3mkv]. Durations have been converted from months in the original dataset to weeks.
#' * `pharmaonc` provides a dataset based on [pharmaverseadam::adsl] and [pharmaverseadam::adrs_onco] to demonstrate how this package can be used with ADaM ADTTE datasets.
#' * `survcan` provides a dataset based on [survival::cancer()]. This contains the necessary ID and overall survival fields only. Durations have been converted from days in the original dataset to weeks. You will additionally need to supply PFS and TTP data (fields pfs.durn, pfs.flag, ttp.durn and ttp.flag) to use [psm3mkv].
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @export
#' @examples
#' create_dummydata("survcan") |> head()
#' create_dummydata("flexbosms") |> head()
#' create_dummydata("pharmaonc") |> head()
create_dummydata <- function(dsname) {
  dsname <- stringr::str_to_lower(dsname)
  if (dsname=="survcan") {create_dummydata_survcan()}
  else if (dsname=="flexbosms") {create_dummydata_flexbosms()}
  else if (dsname=="pharmaonc") {create_dummydata_pharmaonc()}
  else {stop("Incorrect dataset specified. Must be survcan, flexbosms or pharmaonc.")}
}

#' Create survcan dummy dataset for illustration
#' @description Create 'survcan' dummy dataset to illustrate [psm3mkv]
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @seealso [create_dummydata()]
#' @importFrom rlang .data
#' @noRd
create_dummydata_survcan <- function() {
  # Declare local variables
  ds <- ptid <- os.durn <- os.flag <- NULL
  # Create dataset
  ds <- survival::cancer |>
    dplyr::mutate(
      ptid = dplyr::row_number(),
      os.durn = .data$time/7,
      os.flag = .data$status-1
    ) |>
    dplyr::select(ptid, os.durn, os.flag)
  # Label dataset
  attr(ds$ptid, "label") <- "Patient ID = row number"
  attr(ds$os.durn, "label") <- "Duration of overall survival, =time/7"
  attr(ds$os.flag, "label") <- "Event flag for overall survival (1=event, 0=censor), =status-1"
  return(ds)
}

#' Create flexbosms dataset for illustration
#' @description Create 'flexbosms' dummy dataset to illustrate [psm3mkv]
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @seealso [create_dummydata()]
#' @importFrom rlang .data
#' @noRd
create_dummydata_flexbosms <- function() {
  # Declare local variables
  ds <- id <- pfs.durn <- pfs.flag <- NULL
  os.durn <- os.flag <- ttp.durn <- ttp.flag <- NULL
  # Create dataset
  ds <- flexsurv::bosms3 |>
        tidyr::pivot_wider(id_cols = "id",
                     names_from = "trans",
                     values_from = c("Tstart", "Tstop", "status")
        ) |>
        dplyr::mutate(
          conv = 365.25/(12*7), # Convert from months to weeks
          pfs.durn = .data$conv * .data$Tstop_2,
          pfs.flag = .data$status_1 + .data$status_2,
          os.durn = .data$conv * ifelse(is.na(.data$Tstop_3), .data$Tstop_2, .data$Tstop_3),
          os.flag = ifelse(is.na(.data$Tstop_3), .data$status_2, .data$status_3),
          ttp.durn = .data$pfs.durn,
          ttp.flag = .data$status_1
        ) |>
        dplyr::select(id, pfs.durn, pfs.flag,
                      os.durn, os.flag, ttp.durn, ttp.flag) |>
        dplyr::rename(ptid = "id")
  # Label dataset
  attr(ds$ptid, "label") <- "Patient ID"
  attr(ds$pfs.durn, "label") <- "Duration of PFS"
  attr(ds$pfs.flag, "label") <- "Event flag for PFS (1=event, 0=censor)"
  attr(ds$os.durn, "label") <- "Duration of OS"
  attr(ds$os.flag, "label") <- "Event flag for PFS (1=event, 0=censor)"
  attr(ds$ttp.durn, "label") <- "Duration of TTP"
  attr(ds$ttp.flag, "label") <- "Event flag for TTP (1=event, 0=censor)"
  return(ds)
}

#' Create pharmaonc dataset for illustration
#' @description Create 'pharmaonc' dummy dataset to illustrate [psm3mkv]. This dataset is derived from `pharmaverse::adsl` and `pharmaverse::adrs_onco`. Overall Survival and Time To Progression are derived using `admiral::derive_param_tte()`, then durations are calculated in weeks.
#' @return Tibble dataset, for use with [psm3mkv] functions
#' @seealso [create_dummydata()]
#' @importFrom rlang .data
#' @noRd
create_dummydata_pharmaonc <- function() {
  # Create local variables
  DTHFL <- DTHDT <- LSTALVDT <- AVALC <- ADT <- ASEQ <- RANDDT <- STARTDT <- NULL
  CNSR <- USUBJID <- PARAMCD <- DURN <- EVFLAG <- DURN_OS <- EVFLAG_OS <- DURN_TTP <- EVFLAG_TTP <- NULL
  ttp.durn <- os.durn <- ttp.flag <- os.flag <- ptid <- pfs.durn <- pfs.flag <- NULL
  # Obtain ADSL and ADRS datsets from pharmaverseadam
  adsl <- pharmaverseadam::adsl
  adrs <- pharmaverseadam::adrs_onco
  # Define event: death
  death <- admiral::event_source(
    dataset_name = "adsl",
    filter = DTHFL == "Y",
    date = DTHDT,
    set_values_to = admiral::exprs(
      EVNTDESC = "DEATH",
      SRCDOM = "ADSL",
      SRCVAR = "DTHDT"
    )
  )
  # Define event: last date alive
  last_alive_dt <- admiral::censor_source(
    dataset_name = "adsl",
    date = LSTALVDT,
    set_values_to = admiral::exprs(
      EVNTDESC = "LAST DATE KNOWN ALIVE",
      SRCDOM = "ADSL",
      SRCVAR = "LSTALVDT"
    )
  )
  # Define event: progression
  pd <- admiral::event_source(
    dataset_name = "adrs",
    filter = AVALC == "PD",
    date = ADT,
    set_values_to = admiral::exprs(
      EVENTDESC = "PD",
      SRCDOM = "ADRS",
      SRCVAR = "ADTM",
      SRCSEQ = ASEQ
    )
  )
  # Start creating dataset
  # Derive OS date
  admiral::derive_param_tte(
    dataset_adsl = adsl,
    start_date = RANDDT,
    event_conditions = list(death),
    censor_conditions = list(last_alive_dt),
    source_datasets = list(adsl = adsl, adrs = adrs),
    set_values_to = admiral::exprs(PARAMCD = "OS", PARAM = "Overall Survival")
  ) |>
  # Derive TTP date
    admiral::derive_param_tte(
      dataset_adsl = adsl,
      start_date = RANDDT,
      event_conditions = list(pd),
      censor_conditions = list(last_alive_dt),
      source_datasets = list(adsl = adsl, adrs = adrs),
      set_values_to = admiral::exprs(PARAMCD = "TTP", PARAM = "Time to Progression")
    ) |>
  # Derive durations of TTP and PFS
    dplyr::mutate(
      DURN = admiral::compute_duration(
          start_date = STARTDT,
          end_date = ADT,
          trunc_out = FALSE,
          out_unit = "weeks",
          add_one = FALSE
          ),
      EVFLAG = 1-CNSR
    ) |>
  # Keep only necessary fields
    dplyr::select(USUBJID, PARAMCD, DURN, EVFLAG) |>
  # Pivot wide the duration and event flag fields
    tidyr::pivot_wider(
      id_cols = "USUBJID",
      names_from = "PARAMCD",
      values_from = c("DURN", "EVFLAG")
    ) |>
  # Rename to required field names
    dplyr::rename(
      ptid = USUBJID,
      os.durn = DURN_OS,
      os.flag = EVFLAG_OS,
      ttp.durn = DURN_TTP,
      ttp.flag = EVFLAG_TTP
    ) |>
  # Add a PFS field
    dplyr::mutate(
      pfs.durn = pmin(ttp.durn, os.durn),
      pfs.flag = 1-(1-ttp.flag)*(1-os.flag)
    ) |>
    dplyr::select(ptid, ttp.durn, ttp.flag, pfs.durn, pfs.flag, os.durn, os.flag)
}

#' Check consistency of PFS definition
#' Check that PFS is defined consistently with TTP and OS in a dataset. This convenience function compares `pfs.durn` with the lower of `ttp.durn` and `os.durn`, and checks that the event field `pfs.flag` is consistent with `ttp.flag` and `os.flag` (is 1 when either `ttp.flag` or `os.flag` is one).
#' @param ds Tibble of complete patient-level dataset
#' - `ttp.durn`, `pfs.durn`, and `os.durn` are the durations of TTP (time to progression), PFS (progression-free survival), and OS (overall survival).
#' - `ttp.flag`, `pfs.flag`, and `os.flag`, and `pps.flag` are event flag indicators for TTP, PFS, and OS respectively (1=event, 0=censoring).
#' @export
#' @return List containing:
#' - `durn`: Logical vector comparing expected and actual PFS durations
#' - `flag`: Logical vector comparing expected and actual PFS event flags
#' - `all`: Single logical value of TRUE if all durations and flags match as expected, FALSE otherwise
#' @export
#' @examples
#' ponc <- create_dummydata("pharmaonc")
#' check_consistent_pfs(ponc)
check_consistent_pfs <- function(ds) {
  durn <- flag <- NULL
  durn <- ds$pfs.durn==pmin(ds$ttp.durn, ds$os.durn)
  flag <- ds$pfs.flag==1-(1-ds$ttp.flag)*(1-ds$os.flag)
  list(durn=durn, flag=flag, all=all(c(durn,flag)))
}

