


CrypSplice="/mnt/localstorage/michelle/data/Projects/CrypSplice/CrypSplice_Editing/Most_Updated_Version/Load_Filter_Change_3_26_24/CrypSplice2.0_modular.py"
gtf="/mnt/localstorage/michelle/data/References/gencode_v28_hg38/gencode.v28.annotation.gtf"
fasta="/mnt/localstorage/michelle/data/References/gencode_v45_GRCh38.p14/GRCh38.p14.genome.fa"




# Test 1
python3 $CrypSplice CrypticLoad Clust -samples "./Bams/Control_1.CrypSplice.bam,./Bams/Control_2.CrypSplice.bam,./Bams/U2Surp_1.CrypSplice.bam,./Bams/U2Surp_2.CrypSplice.bam" -gtf $gtf -fasta $fasta -s 2 -o ./CrypticLoad_Test1/ 


# Test 2
python3 $CrypSplice CrypticLoad Clust -samples "./Bams/Control_1.CrypSplice.bam,./Bams/Control_2.CrypSplice.bam,./Bams/Rbm17_1.CrypSplice.bam,./Bams/Rbm17_2.CrypSplice.bam" -gtf $gtf -fasta $fasta -s 2 -o ./CrypticLoad_Test2/

