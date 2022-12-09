mkdir 0.fastqc 3.trim 4.hisat2 5.stringtie
cd 2.cleandata
ls | sed /name/d > name
cp name ../3.trim;cp name ../4.hisat2;cp name ../5.stringtie;
for i in `cat name`;do cd $i; md5sum -c MD5.txt >> ../../alllog 2>&1;mv *.fq.gz ..;cd ..;done
~/softwares/FastQC/fastqc -t 15 -o ../0.fastqc *.fq.gz
cd 3.trim
cp ~/softwares/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa .
for i in `cat name`;do 
echo $i;
java -jar ~/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 15 -phred33 ../2.cleandata/${i}_1.clean.fq.gz ../2.cleandata/${i}_2.clean.fq.gz -baseout ${i}.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 >> ../alllog 2>&1;
done
cd 4.hisat2
for i in `cat name`;do
echo $i >> ../alllog;
hisat2 -p 15 --dta --rg-id $i --rg SM:$i --rna-strandness RF --fr -x ~/genome/tair10seq \
-1 ../3.trim/${i}_1P.fq.gz -2 ../3.trim/${i}_2P.fq.gz -U ../3.trim/${i}_1U.fq.gz,../3.trim/${i}_2U.fq.gz -S ${i}.sam >> ../alllog 2>&1;
samtools view -bSF 4 ${i}.sam > ${i}.bam;samtools sort -@ 15 -o ${i}.sorted.bam ${i}.bam;samtools index ${i}.sorted.bam;rm ${i}.sam;
done
for i in `cat name`;do bamCoverage -p 10 -bs 1 -b ${i}.sorted.bam -o ${i}.bw;done
cd 5.stringtie
for i in `cat name`;do echo $i;stringtie ../4.hisat2/${i}.sorted.bam --rf -e -B -A ${i}.gene.abundance -G TAIR10_GFF3_genes_exons.gff -p 18 -o ${i}/${i}.gtf;done
ls */*.gtf >sample_lst.txt
python prepDE.py3 -i sample_lst.txt
###get merged transcriptome
for i in `cat name`;do echo $i;stringtie ../4.hisat2/${i}.sorted.bam --rf -G TAIR10_GFF3_genes_transposons.gff -p 28 -o ${i}.te/${i}.gtf;done
stringtie --merge -G TAIR10_GFF3_genes_transposons.gff -o merge.gtf ddm1-1.te/ddm1-1.gtf ddm1-2.te/ddm1-2.gtf \
ccr4A-ddm1-1.te/ccr4A-ddm1-1.gtf ccr4A-ddm1-2.te/ccr4A-ddm1-2.gtf ccr4B-ddm1-1.te/ccr4B-ddm1-1.gtf ccr4B-ddm1-2.te/ccr4B-ddm1-2.gtf \
cafA1cafB3-ddm1-1.te/cafA1cafB3-ddm1-1.gtf cafA1cafB3-ddm1-2.te/cafA1cafB3-ddm1-2.gtf
gffread -w transcripts.fa -g ~/genome/tair10.fa merge.gtf
samtools faidx transcripts.fa
makeblastdb -in transcripts.fa -dbtype nucl -parse_seqids -out merge
blastn -query transcripts.fa -out merge2mergeall -outfmt 6 -db merge
awk '$1~/MSTR/ && $2!~/MSTR/' merge2mergeall > new
awk '!a[$1]++{print}' new > new2gene
perl ~/bin/merge_total.pl ~/genome/gene_transid_te_type 2 new2gene 2 o
cut -f 1 o |sort  > new_has_type
grep ">" transcripts.fa | sed 's#>##' |grep "MSTRG" |sort> new.id
comm -1 -3 new_has_type new.id > new_not_type
awk '{print $1"\t"$1"\tNA"}' new_not_type > new_has_no_type
awk '{print $13"\t"$1"\t"$15}' o > new_type
cat ~/genome/gene_transid_te_type new_type new_has_no_type > ~/genome/gene_transid_te_assemble_type
