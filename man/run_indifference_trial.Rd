% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_trial_indifference.R
\name{run_indifference_trial}
\alias{run_indifference_trial}
\title{Run an indifference trial with interims}
\usage{
run_indifference_trial(
  id,
  mu,
  delta_indiff = 0.25,
  delta_ctrl = delta_indiff,
  delta_best = 0,
  kappa_sup = 0.95,
  kappa_ind = 0.75,
  kappa_ctr = 0.95,
  ctrl_alloc = 1/13,
  zero_alloc = 0.01,
  N = 10000,
  K = 20
)
}
\description{
An indifference stopping...
}
