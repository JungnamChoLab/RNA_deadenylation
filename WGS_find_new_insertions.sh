###SPLITREADER related are from https://github.com/baduelp/public/tree/master/SPLITREADER
mkdir 0.fastqc 2.mapping 3.splitreader
cd 1.data
ls | sed /name/d > name
for i in `cat name`;do cd $i;md5sum -c MD5.txt;mv *.gz ../;cd ..;done
for i in `cat name`;do ~/softwares/FastQC/fastqc -t 30 ${i}*.fq.gz -o ../0.fastqc/;done
cd 2.mapping
for i in `cat name`;do 
echo $i >> mapping_record;
bowtie2 -p 25 --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive  -x ~/genome/tair10seq -1 ../1.data/${i}_R1.fq.gz -2 ../1.data/${i}_R2.fq.gz -S ${i}.sam >> mapping_record 2>&1;
samtools view -bS ${i}.sam > ${i}_pe.bam;
done
cd ..

for i in `cat name`;do 
sh SPLITREADER-beta2.5_part1.sh $i part1 2.mapping _pe tair10 public/home/wangling/digging/deadenylation_resequencing/3.splitreader TE_all_Athaliana;
done

for i in `cat name`;do sh SPLITREADER-beta2.5_part2.sh $i tair10 public/home/wangling/digging/deadenylation_resequencing/3.splitreader tair10seq te_gene;
done

for i in `cat name`;do 
sh SPLITREADER-sort.sh $i /public/home/wangling/digging/deadenylation_resequencing/3.splitreader /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10 $i tair10;
done
mv /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/*-insertion-sites.*bed /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/ALL

for i in `cat TE_family_name`;do echo $i;
  for j in `cat name`;do sh SRfam_wrapper.sh /public/home/wangling/digging/deadenylation_resequencing/3.splitreader $i $j /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10;
  done;
done

for i in `cat TE_family_name`;do echo $i;
  for j in `cat /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/TE_sequence/${i}.TElist`;
  do sh Intersect_insertions_splitreader.sh ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/${j} ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/${j} tair10 ${j};
  done;
done

for i in `cat TE_family_name`;do echo $i;
  for j in `cat /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/TE_sequence/${i}.TElist`;
  do cd ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/${j};
  sh tampon.intersect.${j}.tair10.sh;
  done;
done

for i in `cat TE_family_name`;do echo $i;
  for j in `cat /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/TE_sequence/${i}.TElist`;
  do cd ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/${j};
  cut -f 1-4 ${j}.tair10-intersect.bed | sortBed -i -|mergeBed -d 5 -i - |sortBed -i - > ${j}.tair10-intersect.mrg.bed;
  done;
done

cd /public/home/wangling/digging/deadenylation_resequencing/

for i in `cat TE_family_name`;do echo $i;
  for j in `cat /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/TE_sequence/${i}.TElist`;
  do cd ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/${j};
  cut -f 1-4 ${j}.tair10-intersect.bed | sortBed -i - > ${j}.tair10-intersect.sort.bed;
  done;
done

cd /public/home/wangling/digging/deadenylation_resequencing/

for i in `cat TE_family_name`;do echo $i;
  for j in `cat /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/TE_sequence/${i}.TElist`;
  do perl Filter_insertions_splitreader.pl tair10 3 ~/digging/deadenylation_resequencing/3.splitreader $j ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/${j} \
  10 50 caf1a1bddm1-L1-1 caf1a1bddm1-L1-2 caf1a1bddm1-L1-3 caf1a1bddm1-L1-4 caf1a1bddm1-L1-5 \
  caf1a1bddm1-L2-1 caf1a1bddm1-L2-2 caf1a1bddm1-L2-3 caf1a1bddm1-L2-4 caf1a1bddm1-L2-5 \
  ccr4a-1ddm1-L1-1 ccr4a-1ddm1-L1-2 ccr4a-1ddm1-L1-3 ccr4a-1ddm1-L1-4 ccr4a-1ddm1-L1-5 \
  ccr4a-1ddm1-L2-1 ccr4a-1ddm1-L2-2 ccr4a-1ddm1-L2-3 ccr4a-1ddm1-L2-4 ccr4a-1ddm1-L2-5 \
  ccr4b-1ddm1-L1-1 ccr4b-1ddm1-L1-2 ccr4b-1ddm1-L1-3 ccr4b-1ddm1-L1-4 ccr4b-1ddm1-L1-5 \
  ccr4b-1ddm1-L2-1 ccr4b-1ddm1-L2-2 ccr4b-1ddm1-L2-3 ccr4b-1ddm1-L2-4 ccr4b-1ddm1-L2-5 \
  ddm1-L1-1 ddm1-L1-2 ddm1-L1-3 ddm1-L1-4 ddm1-L1-5 \
  ddm1-L2-1 ddm1-L2-2 ddm1-L2-3 ddm1-L2-4 ddm1-L2-5 > ${j}.tair10.SRfilt.e;
  done;
done

for i in `cat TE_family_name`;do echo $i;cat /public/home/wangling/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10/${i}/*.tair10-insertions.filt.DP3.bed > 3.splitreader/BEDfiles/SPLITREADER/tair10/${i}.tair10-insertions.filt.DP3.bed;done

cat 3.splitreader/BEDfiles/SPLITREADER/tair10/*.tair10-insertions.filt.DP3.bed > 3.splitreader/BEDfiles/SPLITREADER/tair10/tair10-insertions.filt.DP3.bed

sed -i '/start/d'  3.splitreader/BEDfiles/SPLITREADER/tair10/tair10-insertions.filt.DP3.bed

sed '/start/d' 3.splitreader/BEDfiles/SPLITREADER/tair10/tair10-insertions.filt.DP3.bed |cut -f 1-3 > tair10-insertion-sites.0.bed
sed '/start/d' 3.splitreader/BEDfiles/SPLITREADER/tair10/tair10-insertions.filt.DP3.bed |cut -f 1-3 |awk -v OFS='\t' '{$2 = $2 - 100; $3 = $3 - 100 ; if ( $3 > 0 ) print $1, $2, $3}' > tair10-insertion-sites.100up.bed
sed '/start/d' 3.splitreader/BEDfiles/SPLITREADER/tair10/tair10-insertions.filt.DP3.bed |cut -f 1-3|awk -v OFS='\t' '{$2 = $2 + 100; $3 = $3 + 100 ; if ( $3 > 0 ) print $1, $2, $3}' > tair10-insertion-sites.100down.bed
for i in `cat name`;do sh BAM-readcount_wrapper.sh $i ~/digging/deadenylation_resequencing/2.mapping .sorted $i 3 tair10 ~/digging/deadenylation_resequencing/3.splitreader/tair10/${i}/BAMrc ~/digging/deadenylation_resequencing ~/digging/deadenylation_resequencing/3.splitreader tair10.fa;done
perl Filter_negative_calls_splitreader.pl tair10 3 ~/digging/deadenylation_resequencing/3.splitreader/BEDfiles/SPLITREADER/tair10 ~/digging/deadenylation_resequencing/3.splitreader filt \
caf1a1bddm1-L1-1 caf1a1bddm1-L1-2 caf1a1bddm1-L1-3 caf1a1bddm1-L1-4 caf1a1bddm1-L1-5 \
caf1a1bddm1-L2-1 caf1a1bddm1-L2-2 caf1a1bddm1-L2-3 caf1a1bddm1-L2-4 caf1a1bddm1-L2-5 \
ccr4a-1ddm1-L1-1 ccr4a-1ddm1-L1-2 ccr4a-1ddm1-L1-3 ccr4a-1ddm1-L1-4 ccr4a-1ddm1-L1-5 \
ccr4a-1ddm1-L2-1 ccr4a-1ddm1-L2-2 ccr4a-1ddm1-L2-3 ccr4a-1ddm1-L2-4 ccr4a-1ddm1-L2-5 \
ccr4b-1ddm1-L1-1 ccr4b-1ddm1-L1-2 ccr4b-1ddm1-L1-3 ccr4b-1ddm1-L1-4 ccr4b-1ddm1-L1-5 \
ccr4b-1ddm1-L2-1 ccr4b-1ddm1-L2-2 ccr4b-1ddm1-L2-3 ccr4b-1ddm1-L2-4 ccr4b-1ddm1-L2-5 \
ddm1-L1-1 ddm1-L1-2 ddm1-L1-3 ddm1-L1-4 ddm1-L1-5 \
ddm1-L2-1 ddm1-L2-2 ddm1-L2-3 ddm1-L2-4 ddm1-L2-5 
