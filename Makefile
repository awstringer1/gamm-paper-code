# Makefile for gamm simulations
# Arguments:
# 	numsims: number of simulations to run
# 	date: datestamp for output files
# 	version: version number for output files
# Example use: make numsims=100 version=1 date=20250306 simulations


simulations:
	# Rscript 01-main.R $(numsims) "sims-main-$(date)-v$(version)" $(version)
	# Rscript 02-multiple-smooth.R $(numsims) "sims-multiple-$(date)-v$(version)" $(version)
	# Rscript 03-gamm4-poisson.R $(numsims) "sims-poisson-$(date)-v$(version)" $(version)
	# Rscript 04-sims-complexfunction.R $(numsims) "sims-complexfunction-$(date)-v$(version)" $(version)
	# Rscript 05-sims-sigma.R $(numsims) "sims-sigma-$(date)-v$(version)" $(version)
	# Rscript 06-sims-k.R $(numsims) "sims-k-$(date)-v$(version)" $(version)
	# Rscript 07-sims-small-m.R $(numsims) "sims-small-$(date)-v$(version)" $(version)
	Rscript summarize-simulations.R $(date)
	
	
