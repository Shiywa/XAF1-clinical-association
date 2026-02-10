## 转bam为homer可用的tag格式
makeTagDirectory /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-3-HSV-1 \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/alignment/Input-3-HSV-1_1.fastq.gz/Input-3-HSV-1_1.fastq.gz_unique.bam \
    -tbp 1 -fragLength given

makeTagDirectory /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-4-HSV-1 \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/alignment/Input-4-HSV-1_1.fastq.gz/Input-4-HSV-1_1.fastq.gz_unique.bam \
    -tbp 1 -fragLength given

makeTagDirectory /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1 \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/alignment/IP-3-HSV-1_1.fastq.gz/IP-3-HSV-1_1.fastq.gz_unique.bam \
    -tbp 1 -fragLength given

makeTagDirectory /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1 \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/alignment/IP-4-HSV-2_1.fastq.gz/IP-4-HSV-2_1.fastq.gz_unique.bam \
    -tbp 1 -fragLength given

## homer call peak
findPeaks /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1 \
    -i /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-3-HSV-1 \
    -o auto -style clip -size 100 -F 1.1 -L 1.1 -fdr 0.1

findPeaks /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1 -i /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-4-HSV-1 -o auto -style clip -F 2 -L 2 -fdr 0.05

## peak文件转bed
/home/bio/biosoft/homer/bin/pos2bed.pl /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1/peaks.txt > IP-3-HSV-1_peak_homer.bed
/home/bio/biosoft/homer/bin/pos2bed.pl /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1/peaks.txt > IP-4-HSV-1_peak_homer.bed

## homer生成bigwig
makeUCSCfile /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1 \
    -o auto -norm 100000000
gunzip /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1/IP-3-HSV-1.ucsc.bedGraph.gz
bedGraphToBigWig /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1/IP-3-HSV-1.ucsc.bedGraph \
    /home/bio/biosoft/reference/GRCh38_UCSC/star_index/chrNameLength.txt \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-3-HSV-1/IP-3-HSV-1_homer.bw

makeUCSCfile /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1 \
    -o auto -norm 100000000
gunzip /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1/IP-4-HSV-1.ucsc.bedGraph.gz
bedGraphToBigWig /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1/IP-4-HSV-1.ucsc.bedGraph \
    /home/bio/biosoft/reference/GRCh38_UCSC/star_index/chrNameLength.txt \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/IP-4-HSV-1/IP-4-HSV-1_homer.bw


makeUCSCfile /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-3-HSV-1 \
    -o auto -norm 100000000
gunzip /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-3-HSV-1/Input-3-HSV-1.ucsc.bedGraph.gz
bedGraphToBigWig /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-3-HSV-1/Input-3-HSV-1.ucsc.bedGraph \
    /home/bio/biosoft/reference/GRCh38_UCSC/star_index/chrNameLength.txt \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-3-HSV-1/Input-3-HSV-1_homer.bw

makeUCSCfile /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-4-HSV-1 \
    -o auto -norm 100000000
gunzip /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-4-HSV-1/Input-4-HSV-1.ucsc.bedGraph.gz
bedGraphToBigWig /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-4-HSV-1/Input-4-HSV-1.ucsc.bedGraph \
    /home/bio/biosoft/reference/GRCh38_UCSC/star_index/chrNameLength.txt \
    /home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP/all_res/homer_peak/Input-4-HSV-1/Input-4-HSV-1_homer.bw
