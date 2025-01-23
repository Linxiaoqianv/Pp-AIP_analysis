mkdir -p ./metagenome/ && cd ./metagenome/

#Quality control (raw data must be downloaded to ./01.RawData in advance)
mkdir -p ./02.fastp/

for i in 01.RawData/HK*
do
    file=`ls $i | sed -n '1p'`
    base=${file%_1.fq.gz}
    mkdir ./02.fastp/$base/
    fastp \
    --in1 $i/${base}_1.fq.gz \
    --in2 $i/${base}_2.fq.gz \
    --out1 ./02.fastp/$base/${base}_1.fq.gz \
    --out2 ./02.fastp/$base/${base}_2.fq.gz \
    --json ./02.fastp/$base/${base}.json \
    --html ./02.fastp/$base/${base}.html \
    --trim_poly_g --poly_g_min_len 10 \
    --trim_poly_x --poly_x_min_len 10 \
    --cut_front --cut_tail \
    --cut_window_size 4 --cut_mean_quality 20 \
    --qualified_quality_phred 15 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --length_required 30 \
    --thread 40
done

#Remove host
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_g
enomic.fna.gz
route_host="../database/host/Arabidopsis_thaliana/TAIR10.1"

mkdir -p ./03.clean/
for i in 02.fastp/HK*
do
    file=`ls $i | sed -n '1p'`
    base=${file%_1.fq.gz}
    mkdir ./03.clean/$base/
    bowtie2 -q -p 32 \
    -1 $i/${base}_1.fq.gz \
    -2 $i/${base}_2.fq.gz \
    -x $route_host \
    -S 03.clean/$base/${base}_clean.sam \
    --un-conc-gz 03.clean/$base/${base}_clean.fq.gz \
    > 03.clean/$base/${base}_clean.log
    rm 03.clean/$base/${base}_clean.sam
done

#Taxonomy annotation (https://github.com/jenniferlu717/KrakenTools must be downloaded to in advance)
standard_db="../database/kraken2_database/standard_db"

mkdir -p ./04.taxo_kraken/result/
for i in ./03.clean/HK*
do
    file=`ls $i | sed -n '1p'`
    base=${file%_clean.fq.1.gz}
    kraken2 --use-names --threads 40 \
            --db $standard_db \
            --report ./04.taxo_kraken/result/${base}.report.txt \
            --gzip-compressed \
            --output ./04.taxo_kraken/result/${base}.result.kraken \
            $i/${base}_clean.fq.1.gz \
            $i/${base}_clean.fq.2.gz
done

mkdir ./04.taxo_kraken/bracken
for i in `ls ./04.taxo_kraken/result/*.report.txt`;
do
    sample=`basename $i .report.txt`
    bracken -d ../database/kraken2_database/standard_db/ -i $i -t 4 -o ./04.taxo_kraken/bracken/${sample}_bracken.out -r 100 -t 20
    kreport2mpa.py -r ./04.taxo_kraken/result/${sample}.report_bracken.txt -o ./04.taxo_kraken/bracken/${sample}_mpa --display-header
    num=`cat $i | head -n 2 | awk '{SUM+=$2}END{print SUM}'`
    echo -e "total_reads\t$num" >> ./04.taxo_kraken/bracken/${sample}_mpa
done

combine_mpa.py -i ./04.taxo_kraken/bracken/*_mpa -o ./04.taxo_kraken/combine_mpa_total.txt
awk '/#/;/d__Bacteria.*s__/;/total_reads/ {print $0}' ./04.taxo_kraken/combine_mpa_total.txt > ./04.taxo_kraken/combine_mpa_bacteria_species.txt
awk '/#/;/d__Bacteria.*g__/;/total_reads/ {print $0}' ./04.taxo_kraken/combine_mpa_total.txt > ./04.taxo_kraken/combine_mpa_bacteria_genus.txt
sed -i '/s__/d' ./04.taxo_kraken/combine_mpa_bacteria_genus.txt
awk '/#/;/d__Bacteria.*p__/;/total_reads/ {print $0}' ./04.taxo_kraken/combine_mpa_total.txt > ./04.taxo_kraken/combine_mpa_bacteria_phylum.txt
sed -i '/[cogs]__/d' ./04.taxo_kraken/combine_mpa_bacteria_phylum.txt
