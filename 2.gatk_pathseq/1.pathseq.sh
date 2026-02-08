ls ../1.spaceranger/outputs/ | while read p; do p=$(grep $p ../sampleList | sed 's/\t/ /') && echo "../.libexec/pathseq.sh $p"; done
