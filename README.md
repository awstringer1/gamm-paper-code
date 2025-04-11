# Code for paper Inference for generalized additive mixed models via penalized marginal likelihood

To run the simulations, use the `Makefile`:

`make numsims=100 version=1 date=20250306 simulations`

This will create output files named using `date` and `version` in whatever directory is returned by `getwd()` in `R`. 

The file `summarize-simulations.R` will load all the so-named results in this directory and collate
them. This allows you to perform the simulations in multiple runs on multiple different machines. I created
it this way because running these simulations is extremely onerous. 
To get the `1000` results for `01` reported
in the main text and the `500` for `02-07` reported in the supplement took over a week of running in multiple
batches in parallel on two M4 Mac desktops.
