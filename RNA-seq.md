## RNA-seq
##### PS:Drosophila RNA-seq
脚本位置：
G:\OneDrive - 西湖大学\Ma Lab\Projects\RNA-seq Analysis
F:\RNA_Sequencing_Data


### 1. Linux Part
##### 1.1 Docker Initiation

```
docker run -it -m 32G --cpus 4 --rm -v F:/RNA_Sequencing_Data/rnaseq:/work omicsclass/rnaseq
```

##### 1.2 设置工作路径,注意，每次重启都需要重新设置一次
```
workdir=/work/my_rnaseq  
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts
```

##### 1.3 设置参考基因相关文件变量，方便后续使用,REF_INDEX的文件名为.1.ht2之前的全部文本
```
REF_INDEX=$refdir/Drosophila_melanogaster.BDGP6.32.dna.toplevel
GFF=$refdir/Drosophila_melanogaster.BDGP6.32.103.gff3
GTF=$refdir/Drosophila_melanogaster.BDGP6.32.103.gtf
GENE_BED=$refdir/gene.bed
GENE_LENGTH=$refdir/gene_length.txt
```
##### 1.4 测序数据下机质控(单个质控)
```
cd $workdir
mkdir 1.fastqc

cd $workdir/1.fastqc

fastqc $datadir/*.gz  -o $workdir/1.fastqc
```

##### 1.5 测序数据下机质控(去接头与过滤低质量reads)
```
cd $workdir  
mkdir 2.data_qc
cd $workdir/2.data_qc
```
#Group Information
```
DHR3_1 DHR3_2 DHR3_3 E75+DHR3-1 E75+DHR3-2 E75+DHR3-3 E75-1 E75-2 E75-3 E75-IR-1 E75-IR-2 E75-IR-3 WT-1 WT-2 WT-3
Ras1A Ras2A Ras3A Ras-lgl1A Ras-lgl2A Ras-lgl3A Ras-lgl-GPX1A Ras-lgl-GPX2A Ras-lgl-GPX3A;
R-FB-1A R-FB-2A R-FB-3A RS-FB-1A RS-FB-2A RS-FB-3A
M6-1B M6-2B M6-3B Ras_2B Ras_AA Ras_BA Ras_CA
A1 A2 A3 E1 E2 E3
ptp61f-lar-1 ptp61f-lar-2 Ras-ptp61f-mut-1 Ras-ptp61f-mut-2 Ras-v12-1 Ras-v12-2 W1118-lar-1 W1118-lar-2
4-03A 4-03B 4-03C 119-6DA 119-6DB 119-6DC 119-13DA 119-13Db 119-13DC 141-13DA 141-13DB 141-13DC 189-13DA 189-13DB 189-13DC
685 686 687 688 689 690 
82B-eye-1A 82B-eye-2A 82B-eye-3A 82B-Fat-1A 82B-Fat-2A 82B-Fat-3A scrib-eye-1A scrib-eye-2A scrib-eye-3A scrib-Fat-1A scrib-Fat-2A scrib-Fat-3A TL6-scrib-1A TL6-scrib-2A TL6-scrib-3A 
E75-scrib-wts-1A E75-scrib-wts-2A E75-scrib-wts-3A Case1A Case2A Ras1A Ras2A
119-10DA_1 119-10DA_2 119-10DA_3 61189A_1 61189A_2 61189A_3 Z403A_1 Z403A_2 Z403A_3
QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A
```

```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do 
echo "RUN CMD: fastp --thread 12 --qualified_quality_phred 10 \
--unqualified_percent_limit 50 \
--n_base_limit 10 \
-i $datadir/${i}_1.fq.gz \
-I $datadir/${i}_2.fq.gz \
-o ${i}_1.clean.fq.gz \
-O ${i}_2.clean.fq.gz \
--adapter_fasta $workdir/data/illumina_multiplex.fa \
-h ${i}.html -j ${i}.json"

fastp --thread 14 --qualified_quality_phred 10 \
--unqualified_percent_limit 50 \
--n_base_limit 10 \	
-i $datadir/${i}_1.fq.gz \
-I $datadir/${i}_2.fq.gz \
-o ${i}_1.clean.fq.gz \
-O ${i}_2.clean.fq.gz \
--detect_adapter_for_pe \
-h ${i}.html -j ${i}.json	
done

```
##### 1.6 质控数据统计汇总
```
python $scriptdir/qc_stat.py -d $workdir/2.data_qc/ -o $workdir/2.data_qc/ -p all_sample_qc
```

##### 1.7 Hisat2比对到reference genome
```
cd $workdir/3.map/hisat2
```
Step1-hisat2比对
#hisat2软件，链特异性文库设置
#--rna-strandness RF or FR
#非链特异性文库，不设置此参数
```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do
echo "RUN CMD: hisat2 -p 10 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta --rna-strandness RF \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary"

hisat2 -p 10 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta --rna-strandness RF \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary 
done

```
Step2-sam 转换成bam 格式并排序
```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do
echo "RUN CMD: samtools sort --threads 2 -m 16G -o ${i}.bam ${i}.sam"
samtools sort --threads 2 -m 16G -o ${i}.bam ${i}.sam
done
```
Step3-bam index建立索引
```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do
echo "RUN CMD: samtools index ${i}.bam"
samtools index ${i}.bam
done

```

##### 1.8 对map结果进行QC分析
采用RSeQC 对比对结果文件进行质控分析
```
mkdir -p $workdir/3.map/map_QC

cd $workdir/3.map/map_QC
```
Step1-片段inner size，片段选择是否异常

```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do
echo "RUN CMD: inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size"
inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size
done
```

Step2-基因覆盖情况，RNA是否降解
```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do
echo "RUN CMD: geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody"	
geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody
done
```

##### 1.9 基因表达定量及结果展示
Step1-采用Htseq-count 对已知的基因进行表达定量
```
mkdir 4.expression
cd $workdir/4.expression
```
#--order pos 按reads比对的position进行排序; -name 是按reads比对的基因的name进行排序
#--mode 可以选择union; intersection_strict; intersection_nonempty
#--stranded yes or no   文库类型设置，是否为链特异性文库; 
#--minaqual 为比对质量，默认为10;
#--idattr 若要查看每个转录本的表达情况，则需要讲gene_id改成transcript_id
```
for i in QS8-1A QS8-2A QS8-4A QS11-1A QS11-2A QS11-3A R-FB-1A R-FB-2A R-FB-3A W-FB-1A W-FB-2A W-FB-3A; do
echo "RUN CMD: htseq-count --format bam --order pos --mode union \
--stranded reverse --minaqual 10 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv"

htseq-count --format bam --order pos --mode union \
--stranded no --minaqual 10 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv
done
```

Step2-合并不同样品表达定量结果，方便后续基因表达分析
```
python $scriptdir/merge_gene_count.py -p all_gene_count \
-f  QS8-1A_gene.tsv -l QS8-1A  \
-f  QS8-2A_gene.tsv -l QS8-2A \
-f  QS8-4A_gene.tsv -l QS8-4A \
-f  QS11-1A_gene.tsv -l QS11-1A \
-f  QS11-2A_gene.tsv -l QS11-2A \
-f  QS11-3A_gene.tsv -l QS11-3A \
-f  R-FB-1A_gene.tsv -l R-FB-1A \
-f  R-FB-2A_gene.tsv -l R-FB-2A \
-f  R-FB-3A_gene.tsv -l R-FB-3A \
-f  W-FB-1A_gene.tsv -l W-FB-1A \
-f  W-FB-2A_gene.tsv -l W-FB-2A \
-f  W-FB-3A_gene.tsv -l W-FB-3A

```
##手动删除all_gene_count.tsv中没有数据对应，显示为NA的行





### 2. R-Studio Part








