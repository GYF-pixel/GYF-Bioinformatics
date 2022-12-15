##Drosophila RNA-seq Procedures
#linux 基础：
# https://study.163.com/course/introduction/1006346005.htm?share=1&shareId=1030291076
#更多docker使用及Linux服务器搭建：
# https://study.163.com/course/introduction/1209757831.htm?share=1&shareId=1030291076
######################################################################################
#打开PowerShell
#搜索docker镜像
#docker search omicsclass
#查看本地的docker镜像
#docker images
#下载转录组分析docker镜像
#docker pull omicsclass/rnaseq:v1.0    
#启动docker容器并交互式进入
#docker desktop方法
#docker run -it -m 32G --cpus 4 --rm -v F:/RNA_Sequencing_Data/rnaseq:/work omicsclass/rnaseq
#docker toolbox方法
#docker run --rm -it -m 4G --cpus 1  -v /d/rnaseq:/work omicsclass/rnaseq:v1.0 

#########
###1.参考基因组的准备
#从ensemble下载参考基因组的GTF和GFF3文件--人homo_sapiens; 小鼠mus_musculus; 和果蝇drosophila_melanogaster 
#https://asia.ensembl.org/info/about/species.html
#https://zhuanlan.zhihu.com/p/112047210
#解压ref
gunzip Drosophila_melanogaster.BDGP6.32.103.gff3.gz 
gunzip *.gz
#删除gff3文件当中的gene:和transcript:
sed 's#gene:##' Drosophila_melanogaster.BDGP6.32.103.gff3| sed 's#transcript:##' >Drosophila_melanogaster.BDGP6.32.103.gff31
mv Drosophila_melanogaster.BDGP6.32.103.gff31 Drosophila_melanogaster.BDGP6.32.gff3
# 准备一些路径变量，方便后续调用
#设置工作路径,注意，每次重启都需要重新设置一次
workdir=/work/my_rnaseq  
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts
echo $refdir #检查路径变量是否设置成功
cd $refdir  #进入参考基因组ref目录

###2.Hisat2建立索引
#按顺序录入sh $scriptdir/index.sh      DNA.fa文件      .gff3文件     .gtf文件 ##在没有.gtf文件的情况下也可以通过.fa和.gff3文件去根据index.sh自动生成.gtf文件
sh $scriptdir/index.sh Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa Drosophila_melanogaster.BDGP6.32.103.gff3 Drosophila_melanogaster.BDGP6.32.103.gtf
#若不想自己建立索引，则可以通过官网下载： http://daehwankimlab.github.io/hisat2/download/
#检查.gtf文件当中是否包含"exon" “gene_id”和“transcript_id”
ee Drosophila_melanogaster.BDGP6.32.103.gtf
#检查gene_length.txt文件
ee gene_length.txt
#设置参考基因相关文件变量，方便后续使用,REF_INDEX的文件名为.1.ht2之前的全部文本
REF_INDEX=$refdir/Drosophila_melanogaster.BDGP6.32.dna.toplevel
GFF=$refdir/Drosophila_melanogaster.BDGP6.32.103.gff3
GTF=$refdir/Drosophila_melanogaster.BDGP6.32.103.gtf
GENE_BED=$refdir/gene.bed
GENE_LENGTH=$refdir/gene_length.txt

echo $GENE_LENGTH #检查参考基因相关文件变量是否设置成功


###3.测序质控
#Step1-fastqc所有data
cd $workdir
mkdir 1.fastqc
fastqc $datadir/*.gz  -o $workdir/1.fastqc
##使用echo命令以查看完整的命令 echo“fastqc $datadir/*.gz  -o $workdir/1.fastqc”
#Step2-fastp工具去除adapter
cd $workdir  #回到工作目录
mkdir 2.data_qc
cd  2.data_qc
#大部分的sequencing_data都是已经去过adaptor的，因此把--adapter_fasta $workdir/data/illumina_multiplex.fa \注释掉
#如果含有adaptor但是又未知，则可以用--detect_adapter_for_pe   设置软件自动识别常见接头
for i in KA-1A KA-2A KB-1A KB-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do 
echo "RUN CMD: fastp --thread 12 --qualified_quality_phred 10 \
--unqualified_percent_limit 50 \
--n_base_limit 10 \
-i $datadir/${i}_r1.fq.gz \
-I $datadir/${i}_r2.fq.gz \
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

 #--adapter_fasta $workdir/data/illumina_multiplex.fa \   加在--detect_adapter_for_pe \   位置
 
#Step3-质控数据统计汇总：
python $scriptdir/qc_stat.py -d $workdir/2.data_qc/ -o $workdir/2.data_qc/ -p all_sample_qc

###4.Hisat2比对参考基因组ref
cd $workdir  #回到工作目录
mkdir -p $workdir/3.map/hisat2
cd $workdir/3.map/hisat2
#hisat2软件，链特异性文库设置
#--rna-strandness RF or FR
#非链特异性文库，不设置此参数

for i in KA-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do
echo "RUN CMD: hisat2 -p 10 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta --rna-strandness RF \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary"

hisat2 -p 14 --rg-id=${i} --rg SM:${i} --rg LB:${i} --rg PL:ILLUMINA \
-x $REF_INDEX --dta \
-1 $workdir/2.data_qc/${i}_1.clean.fq.gz \
-2 $workdir/2.data_qc/${i}_2.clean.fq.gz \
-S ${i}.sam 2>${i}.summary
done

# sam 转换成bam 格式并排序
for i in KA-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do
echo "RUN CMD: samtools sort --threads 10 -m 32G -o ${i}.bam ${i}.sam"
samtools sort  --threads 10 -o ${i}.bam ${i}.sam
done
# bam index建立索引
for i in KA-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do
echo "RUN CMD: samtools index ${i}.bam"
samtools index ${i}.bam
done


###5.对map结果进行QC分析
#采用RSeQC 对比对结果文件进行质控分析
cd $workdir  #回到工作目录
mkdir -p $workdir/3.map/map_QC
cd $workdir/3.map/map_QC

#Step1-片段inner size，片段选择是否异常，主要观察其是否满足正态分布，有无杂峰
for i in KA-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do
echo "RUN CMD: inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size"
inner_distance.py -i $workdir/3.map/hisat2/${i}.bam  -r $GENE_BED  -o ${i}_inner_size
done

#Step2-基因覆盖情况，RNA是否降解，若3‘端明显翘起，而5’明显变低，则发生了降解
for i in KA-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do
echo "RUN CMD: geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody"	
geneBody_coverage.py -r $GENE_BED -i $workdir/3.map/hisat2/${i}.bam  -o ${i}.genebody
done


###6.基因表达定量及结果展示
#Step1-采用Htseq-count 对已知的基因进行表达定量
cd $workdir/ #回到工作目录
mkdir 4.expression
cd $workdir/4.expression
#--order pos 按reads比对的position进行排序; -name 是按reads比对的基因的name进行排序
#--mode 可以选择union; intersection_strict; intersection_nonempty
#--stranded yes or no   文库类型设置，是否为链特异性文库; 
#--minaqual 为比对质量，默认为10;
#--idattr 若要查看每个转录本的表达情况，则需要讲gene_id改成transcript_id
for i in KA-2A KC-1A KC-2A KD-2A KE-1A KE-2A; do
echo "RUN CMD: htseq-count --format bam --order pos --mode union \
--stranded reverse --minaqual 10 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv"

htseq-count --format bam --order pos --mode union \
--stranded no --minaqual 10 --type exon \
--idattr gene_id $workdir/3.map/hisat2/${i}.bam $GTF > ${i}_gene.tsv
done

#Step2-合并不同样品表达定量结果，方便后续基因表达分析
#-p后面接输出文件的名字; -f为输入文件的名字; -l为输入文件在输出文件中单独成列的表头
python $scriptdir/merge_gene_count.py -p all_gene_count \
-f  KA-2A_gene.tsv -l KA-2A \
-f  KC-1A_gene.tsv -l KC-1A \
-f  KC-2A_gene.tsv -l KC-2A \
-f  KD-2A_gene.tsv -l KD-2A \
-f  KE-1A_gene.tsv -l KE-1A \
-f  KE-2A_gene.tsv -l KE-2A
      
	  ##手动删除all_gene_count.tsv中没有数据对应，显示为NA的行
	  
#Step3-基因表达定量结果展示
#1.各样本表达量密度 图
#2.各样本表达量box分布图
#3.&4.各样本表达相关性分析热图与聚类图
#注意这里用的是/usr/bin/目录下的Rscript
/usr/bin/Rscript $scriptdir/fpkm_and_plot.R -i all_gene_count.tsv  -l $GENE_LENGTH  -o ./

### 7.基因差异表达分析(DESeq2)，并绘制火山图与MA图 ，以及差异基因表达热图
#到这一步就可以在R上进行操作，也可以继续在linux当中操作
cd $workdir/
mkdir 5.deg
cd 5.deg
#Step1-准备分组文件，在linux里touch一个文件，然后到excel里将建好的group_information复制到其中
touch group_information.txt
#eg：注意ID要和前面count矩阵里文件all_gene_count.tsv的行名一致
#ID	group
#Scrib_Wts-1	Scrib_Wts
#Scrib_Wts-2	Scrib_Wts
#Scrib_Wts-3	Scrib_Wts
#Scrib_Ras-1	Scrib_Ras
#Scrib_Ras-2	Scrib_Ras
#Scrib_Ras-3	Scrib_Ras

#Step2-DESeq2分析差异基因使用的是基因的count值
#-g 分组文件
#-r 指明哪个是对照组，这里指明是Scrib_Wts
#-p 输出文件的前缀，这里指明是Scrib_Ras_vs_Scrib_Wts
/usr/bin/Rscript $scriptdir/deseq_analysis.r \
-i $workdir/4.expression/all_gene_count.tsv \
-g group_information.txt  \
-k $workdir/4.expression/all_gene_fpkm.tsv \
-r Scrib_Wts --fdr 0.01 --fc 2 \
-p Scrib_Ras_vs_Scrib_Wts
#eg：得到的文件如下
#
#
#

