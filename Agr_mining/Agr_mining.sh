mkdir -p ./genome_mining/ && cd ./genome_mining/
#BGC mining (need to download Refseq genomes in advance)
mkdir -p ./out_antiSMASH

for i in `ls ./RefSeq/GCF*.fna.gz`;
do
    sample=`basename $i .fna.gz`
    antismash -c 10 --output-dir ./out_antiSMASH/$sample --output-basename $sample --html-title $sample --genefinding-tool prodigal ./RefSeq/${sample}.fna.gz
done
