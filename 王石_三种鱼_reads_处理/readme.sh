#!/bin/bash
set -e
# wkdir="/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA"
module purge
module load samtools

# 将三次 mapped 结果 软链至此
mkdir Bam_Data
for species in liyu tuantoufang jiyu
do
    wkdir="/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/${species}/reseq_v3/result/02.mapping/picard/"
    for sample in CB-3  F1-1  F1-2  F1-3  F2-1  F2-2  F3-3  NCRC-F1-1
    do

        ln -s ${wkdir}${sample}/${sample}.sorted.mkdup.bam ./Bam_Data/${sample}_${species}.bam
        ln -s ${wkdir}${sample}/${sample}.sorted.mkdup.bam.bai ./Bam_Data/${sample}_${species}.bam.bai

    done
done





# 所有reads 信息
mkdir reads_alignment
for species in liyu tuantoufang jiyu
do
    for sample in CB-3  F1-2  F1-3  F2-1  F2-2  F3-3  F1-1  NCRC-F1-1
    do
        samtools view ./Bam_Data/${sample}_${species}.bam|awk '{if($3!="*" && $6 != "*" && $5 != 0){print $1"\t"$6}}' > ./reads_alignment/${sample}_${species}_mapped_al.txt &
    done
    wait
done

# 统计比对情况 用于后修正
mkdir reads_name
cd reads_alignment
for i in *txt
do
    name=$(echo $i|sed 's/_al//g')
    cut -f1 $i > ../reads_name/$name
done
/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA/script/cmd.py

cd ..

# 统计overlap read name 
mkdir overlap_name
cd mkdir overlap_name
/public/store5/DNA/Project/Reseq/2023/DZOE2023041442-b1/PA/script/stat_reads_mapping.py
cd ..

#nohup ./stat_reads_mapping.py > stat_reads_mapping.log &

#bam_dir="/public/store5/DNA/Project/Reseq/2022/DOE202213215-b1/PA/30_percent/stat/data"
./script/read_stat.py