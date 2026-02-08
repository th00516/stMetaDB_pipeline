pubDir='/home/dataset-assist-0/PUB/yuhao/stmetadb_tmp/PUB/yuhao'

pathseqdb="${pubDir}/gatk/pathseqdb"

sample_path="../1.spaceranger/outputs"
sample="$1"

case $2 in
    'Homo sapiens')
        spec='hs1'
    ;;
    'Mus musculus')
        spec='mm10'
    ;;
    'Rattus norvegicus')
        spec='rn7'
    ;;
    'Rhesus macaque')
        spec='rheMac10'
    ;;
    *)
        exit 0
esac

export SPARK_LOCAL_IP="127.0.0.1"
mkdir -p tmp/${sample}

${pubDir}/gatk/gatk-4.5.0.0/gatk --java-options "-Xmx1024g" PathSeqPipelineSpark \
    --input ${sample_path}/${sample}/outs/possorted_genome_bam.bam \
    --filter-bwa-image ${pathseqdb}/pathseq_host.${spec}.fa.img \
    --min-clipped-read-length 60 \
    --microbe-dict ${pathseqdb}/pathseq_microbe.dict \
    --microbe-bwa-image ${pathseqdb}/pathseq_microbe.fa.img \
    --taxonomy-file ${pathseqdb}/pathseq_microbe_taxonomy.db \
    --output tmp/${sample}/pathseq.complete.bam \
    --scores-output tmp/${sample}/pathseq.complete.csv \
    --is-host-aligned false \
    --filter-duplicates false \
    --min-score-identity .7
#    --kmer-file ${pathseqdb}/pathseq_host.${spec}.bfi \

../.libexec/UMI_annotator.py \
    ${sample_path}/${sample}/outs/possorted_genome_bam.bam \
    '' \
    ${sample_path}/${sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz \
    tmp/${sample}/pathseq.complete.bam \
    tmp/${sample}/pathseq.complete.csv \
    tmp/${sample}/visium.raw_matrix.readname \
    tmp/${sample}/visium.raw_matrix.unmap_cbub.bam \
    tmp/${sample}/visium.raw_matrix.unmap_cbub.fasta \
    tmp/${sample}/visium.raw_matrix.list \
    tmp/${sample}/visium.raw.raw_matrix.readnamepath \
    tmp/${sample}/visium.raw_matrix.genus.cell \
    tmp/${sample}/visium.raw_matrix.genus.csv \
    tmp/${sample}/visium.raw_matrix.validate.csv
