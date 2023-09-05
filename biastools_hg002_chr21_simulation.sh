#samtools addreplacerg -r "ID:chr21_hapA" -o chr21_hapA.vg_major.sorted.rg.bam  /home/mlin77/scr4_blangme2/maojan/biastools/vg_test/chr21_hapA.vg_major.sorted.bam
#samtools addreplacerg -r "ID:chr21_hapB" -o chr21_hapB.vg_major.sorted.rg.bam  /home/mlin77/scr4_blangme2/maojan/biastools/vg_test/chr21_hapB.vg_major.sorted.bam
#
#bt2_sorted_bam="vg.major.bam"
#samtools merge -f ${bt2_sorted_bam} chr21_hapA.vg_major.sorted.rg.bam chr21_hapB.vg_major.sorted.rg.bam
#
##echo "[BIASTOOLS] Filter the heterozygous site in vcf file"
##python3 ../../filter_het_VCF.py -v normalized_chr21_HG002.vcf.gz -o chr21_het.vcf.gz
##tabix -p vcf chr21_het.vcf.gz
#
##echo "[BIASTOOLS] Get the wgs alignment bam file of chr21"
##bt2_sorted_chr21_bam="./bt2.refflow.wgs.chr21.sorted.bam"
#bias_report="vg.major.chr21.bias"
##samtools view -@ 16 ${bt2_sorted_bam} "chr21" -o ${bt2_sorted_chr21_bam}
#
#echo "[BIASTOOLS] Intersect the bam file and vcf file"
#bedtools intersect -a ${bt2_sorted_bam} -b chr21_het.vcf.gz | samtools view -bo vg.major.chr21.het.bam
#samtools index vg.major.chr21.het.bam
#
#echo "[BIASTOOLS] WGS testing"
#python3 ../../ref_bi_wgs.py -s vg.major.chr21.het.bam -v chr21_het.vcf.gz -f GRCh38_chr21.fa -p chr21_golden_distribution.rpt.pickle -o ${bias_report}

bias_report="./bias_report/vg.pop5.chr21.bias"
echo "[BIASTOOLS] Plotting"
# making the report of different categories of bias
python3 ../biastools/golden_graph_report.py   -mb ${bias_report}".SNP"
#python3 ../biastools/golden_graph_report.py   -mb ${bias_report}".gap"
## plot the measures with NMB and NAB
#python3 ../biastools/golden_graph.py          -mb ${bias_report}".SNP"
#python3 ../biastools/golden_graph.py          -mb ${bias_report}".gap"
## generate the balance files without golden / show the mapQ and balance difference between bias sites if -SIM is specified
#python3 ../biastools/depth_balance_grapher.py -mb ${bias_report}".SNP" -sim
#python3 ../biastools/depth_balance_grapher.py -mb ${bias_report}".gap" -sim

