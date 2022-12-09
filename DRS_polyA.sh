~/digging/nanopore_DRS_m6A/minimap2/minimap2 -ax map-ont -L -p 0 -N 10 ~/digging/h101sc2102074/stringtie_assembly/transcripts.fa \
../1.data/ddm/20221104_1045_X2_FAS68456_dfb923be/fastq_pass/pass.fq | samtools view -bh -F 4 |samtools sort -O bam > ddm1_merge_transcriptome.bam
~/digging/nanopore_DRS_m6A/minimap2/minimap2 -ax map-ont -L -p 0 -N 10 ~/digging/h101sc2102074/stringtie_assembly/transcripts.fa \
../1.data/ccr4/20221104_1044_X1_FAV40861_00ced430/fastq_pass/pass.fq | samtools view -bh -F 4 | samtools sort -O bam > ccr4a_merge_transcriptome.bam
samtools index ddm1_merge_transcriptome.bam
samtools index ccr4a_merge_transcriptome.bam

#Collecting ployA length information from raw fast5
nanopolish index -d ../1.data/ddm/20221104_1045_X2_FAS68456_dfb923be/fast5_pass/ -s ../1.data/ddm/20221104_1045_X2_FAS68456_dfb923be/sequencing_summary_FAS68456_4ee67edd.txt \ 
../1.data/ddm/20221104_1045_X2_FAS68456_dfb923be/fastq_pass/pass.fq
nanopolish index -d ../1.data/ccr4/20221104_1044_X1_FAV40861_00ced430/fast5_pass/ -s ../1.data/ccr4/20221104_1044_X1_FAV40861_00ced430/sequencing_summary_FAV40861_b4166b73.txt \
../1.data/ccr4/20221104_1044_X1_FAV40861_00ced430/fastq_pass/pass.fq
nanopolish polya -t 30 -r ../1.data/ddm/20221104_1045_X2_FAS68456_dfb923be/fastq_pass/pass.fq --bam ../2.mapping/ddm1_merge_transcriptome.bam \
-g /public/home/wangling/digging/h101sc2102074/stringtie_assembly/transcripts.fa > ddm1_assembly_reference.tsv
nanopolish polya -t 30 -r ../1.data/ccr4/20221104_1044_X1_FAV40861_00ced430/fastq_pass/pass.fq --bam ../2.mapping/ccr4a_merge_transcriptome.bam \
-g /public/home/wangling/digging/h101sc2102074/stringtie_assembly/transcripts.fa > ccr4a_assembly_reference.tsv

#Check geneBody coverage using RSeQC
~/digging/nanopore_DRS_m6A/minimap2/minimap2 -ax splice -k14 -uf ~/genome/tair10.fa ../1.data/ddm/20221104_1045_X2_FAS68456_dfb923be\
 /fastq_pass/pass.fq |samtools view -bh -F 4 |samtools sort -O bam > ddm1.sorted.bam
~/digging/nanopore_DRS_m6A/minimap2/minimap2 -ax splice -k14 -uf ~/genome/tair10.fa ../1.data/ccr4/20221104_1044_X1_FAV40861_00ced430\
 /fastq_pass/pass.fq |samtools view -bh -F 4 |samtools sort -O bam > ccr4a.sorted.bam
geneBody_coverage.py -r ~/genome/ara_gene.bed -i ddm1.sorted.bam,ccr4a.sorted.bam -o geneBody

#Get expression level
NanoCount --extra_tx_info -i ../2.mapping/ccr4a_merge_transcriptome.bam -o ccr4a_merge.counts.tsv -b ccr4a.merge.bam --max_dist_3_prime -1
NanoCount --extra_tx_info -i ../2.mapping/ddm1_merge_transcriptome.bam -o ddm1_merge.counts.tsv -b ddm1.merge.bam --max_dist_3_prime -1

#Get ployA reads genome mapping
python label.py -b ddm1_merge_transcriptome.bam -o ddm1_merge_transcriptome.polya.bam -p ../3.polya/ddm1_assembly_reference.tsv
samtools view ddm1_merge_transcriptome.polya.bam |awk '$3=="MSTRG.7800.2"' >ddm1_merge_AT2G13160.sam
samtools view -h ddm1_merge_transcriptome.bam | head -86246 > ddm1_merge_transcriptome.header
cat ddm1_merge_transcriptome.header ddm1_merge_AT2G13160.sam | samtools view -bS - > ddm1_merge_AT2G13160.bam
bamToFastq -i ddm1_merge_AT2G13160.bam -fq ddm1_merge_AT2G13160.fq
~/digging/nanopore_DRS_m6A/minimap2/minimap2 -ax splice -k14 -uf ~/genome/tair10.fa ddm1_merge_AT2G13160.fq > ddm1_merge_AT2G13160.bam
python label.py -b ddm1_merge_AT2G13160.bam -o ddm1_merge_AT2G13160.polya.bam -p ../3.polya/ddm1.tsv
samtools sort ddm1_merge_AT2G13160.polya.bam -o ddm1_merge_AT2G13160.polyas.bam
samtools index ddm1_merge_AT2G13160.polyas.bam
 
