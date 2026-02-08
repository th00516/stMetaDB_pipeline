pubDir='/home/dataset-assist-0/PUB/yuhao/stmetadb_tmp/PUB/yuhao'


case $2 in
    'Homo sapiens')
        refdata='refdata-gex-GRCh38-2024-A'
    ;;
    'Mus musculus')
        refdata='refdata-gex-mm10-2020-A'
    ;;
    'Rattus norvegicus')
        refdata='refdata-gex-mRatBN7-2-2024-A'
    ;;
    'Rhesus macaque')
        refdata='refdata-gex-rheMac10-2025-A'
    ;;
    *)
        exit 0
esac

${pubDir}/tenx/tool/spaceranger-4.0.1/spaceranger count --id $1 --transcriptome ${pubDir}/tenx/refdata-gex/${refdata}/ --create-bam true --sample $1 --fastqs ../../fastq/$1 --unknown-slide visium-1 --image ../.libexec/misc/fiducial.tif --nosecondary --loupe-alignment ../.libexec/misc/alignment.json && mv $1 outputs/
