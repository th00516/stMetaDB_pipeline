ls ../2.gatk_pathseq/outputs/*.raw_matrix.rds -d | perl -ne '/.*\/(.*).raw_matrix.rds/; print "$1\n"' | while read p; do echo "Rscript ../.libexec/post_filter.R $p"; done
