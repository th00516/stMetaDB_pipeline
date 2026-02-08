cat ../sampleList | sed 's/\t/ /' | while read p; do echo "../.libexec/spaceranger.sh ${p}"; done
