# 2026-Life-History-Code-Archive
Assorted scripts for analyzing pooled genomic data from midwestern European Corn Borer populations.

These aren't particularly clean  - the product of a busy graduate student. I tried to clean up and annotate them as best I could, but there are some limitations. Feel free to message me if you for some reason want to utilize them but can't work out how to run them.

A quick description of each:

D_calculations.R - used to calculate and plot Tajima's D data

baypass_analysis.R - code to run the baypass program itself

betascan_analysis.R - example of filtering and plotting raw betascan data

betascan_to_prep_jtk_v2.R - reformatting popoolation sync data to create input files for JTK cycle

calc_D.sh - script to run Tajima's D estimates in sliding windows

calc_pi.sh - script to run Pi estimations in sliding windows

calc_pi_with_indels.sh - script to run Pi estimations but with indels included

cleanup_rc.R - script to clean and filter the RC file for broader applicability

filter_indels.sh - script for removing indels and creating SNP only files

fst_calculations.R - used to calculate and plot FST

run_poolseq_cbsu.R - example of using JTK cycle to look for oscillations

variance_sliding.sh - single population Pi calculation example

z_scaff_order.txt - referenced by some other scripts, used for filtering and reordering scaffolds on the Z chromosome
