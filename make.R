#...............................
# Purpose of this script is to
# make this analysis from scratch
# Note, there are extra scripts for reports
# and reports are not made here
# e.g. 00_maps.R is needed for reports
# but not the main analysis
#.................................
source("analyses/00_data_input.R")
source("analyses/01_rehh.R")
source("analyses/02_comparepopehh.R")
source("analyses/03_analyze_crt_dhps_drc.R")
source("analyses/04_XPEHH_stats_crt_dhps.R")
source("analyses/05_haplotype_plots.R")
source("analyses/06_final_figures.R")
