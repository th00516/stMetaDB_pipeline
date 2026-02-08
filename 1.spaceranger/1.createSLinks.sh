fqPath='/XXX/PUB/yuhao/stmetadb_data/fastq'

ls ${fqPath}/*/*_2.fastq.gz -d | perl -ne '/(.*)\/.*/; print"$1\n"' | while read p; do echo "mkdir -p ../../fastq/$(basename $p) && ln -sf ${fqPath}/$(basename $p)/$(basename $p)_1.fastq.gz ../../fastq/$(basename $p)/$(basename $p)_S1_L001_R1_001.fastq.gz && ln -sf ${fqPath}/$(basename $p)/$(basename $p)_2.fastq.gz ../../fastq/$(basename $p)/$(basename $p)_S1_L001_R2_001.fastq.gz &"; done
