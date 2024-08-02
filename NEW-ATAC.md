[toc]

##ATAC Sequencing new 

all data were accidently removed 

###Conversion of BAM files to fastq.gz 

```
module load samtools
	
samtools fastq -1 103061_stage20_2_1.fastq.gz -2 103061_stage20_2_2.fastq.gz ./raw_data/103061_stage20_2.bam
```
###FastQC on raw data 

```
module load fastqc 

fastqc filename 
```

###Trimming with BBDuk

```
for i in $(ls *fastq.gz |rev|cut -c 12-|rev|uniq)
do
/proj/rpz/slurm_scripts/tmpfile_env --in ${i}_1.fastq.gz --in ${i}_2.fastq.gz --out ${i}_1.bbduk.fastq.gz --out ${i}_2.bbduk.fastq.gz '/scratch/Schmidbaur/bbmap/bbmap/./bbduk.sh in1=%i1 in2=%i2 out1=%o1 out2=%o2 ktrim=r qtrim=15 k=21 mink=8 hdist=0 ref=Adapter.fa' > ${i}.bbduk;
done

for i in *.bbduk
do
/proj/rpz/slurm_scripts/nq -m 25 -c 10 -f $i;
done
```
* all jobs were send with 25 GB except 103061_2 with 40 GB 
* qtrim for quality trim of 15 
	* in fastqc we see reads with bad quality at the end (first quality graph) so. we trim those 
* k mer size is 20 
* added /usr/bin/time 
* this job takes not so long ~15 min 
* memory usage of ~17 GB 

```
/scratch/Schmidbaur/bbmap/bbmap/./bbduk.sh in1=103061_stage20_1_1.fastq.gz in2=103061_stage20_1_2.fastq.gz out1=103061_stage20_1_1.bbduk.fastq.gz out2=103061_stage20_1_2.bbduk.fastq.gz ktrim=r qtrim=15 k=21 mink=8 hdist=0 ref=Adapter.fa
```

###Creating an index file 

```
bowtie2-build /scratch/Schmidbaur/Euprymna_HiC/020318/lachesis_assembly/Eup_chr48.fasta Eupsc_chr48
```
* ~1.5 h 

**Remember the * here, if done as slurm job**

###Alignment with Bowtie2 

* this takes 18h 

~~~b 
for i in $(ls *bbduk.fastq.gz | rev | cut -c 17- |rev | uniq)
do
/proj/rpz/slurm_scripts/tmpfile_env --in ${i}1.bbduk.fastq.gz --in ${i}2.bbduk.fastq.gz --external Eupsc_chr48 --out ${i}bbduk_aln.bam 'bowtie2 --very-sensitive -k 10 -p 8 -x %e1 -1 %i1 -2 %i2 |samtools view -u |samtools sort -n -o %o1' > ${i}bbduk_bwt; 
done 

for i in *bbduk_bwt
do
/proj/rpz/slurm_scripts/nq -m 40 -c 10 -f $i;
done
~~~


###Peak calling Genrich

Stage 20

```
/proj/rpz/slurm_scripts/tmpfile_env --in 97307_Stage20_bbduk_aln.bam --in 103061_stage20_2_bbduk_aln.bam --out Stage20_97307_103061.2.genrich '/scratch/Schmidbaur/Genrich/Genrich-master/Genrich -t %i1,%i2 -o %o1 -j -y -r -v' > stage20_2.genrich

/proj/rpz/slurm_scripts/nq -f stage20.genrich -m 80 -o stage20_2.genrich

```

```
/proj/rpz/slurm_scripts/tmpfile_env --in 97307_Stage20_bbduk_aln.bam --out Stage20_97307.genrich '/scratch/Schmidbaur/Genrich/Genrich-master/Genrich -t %i1 -o %o1 -j -y -r -v' > stage20_only_97307.genrich

/proj/rpz/slurm_scripts/nq -f stage20_97307.genrich -m 80 -o stage20_only_97307.genrich

```

```
/proj/rpz/slurm_scripts/tmpfile_env --in 97307_Stage20_bbduk_aln.bam --in 103061_stage20_1_bbduk_aln.bam --out Stage20_97307_103061.1.genrich '/scratch/Schmidbaur/Genrich/Genrich-master/Genrich -t %i1,%i2 -o %o1 -j -y -r -v' > stage20_1.genrich

/proj/rpz/slurm_scripts/nq -f stage20_1.genrich -m 80 -o stage20_1.genrich
```

Stage 25 

```
/proj/rpz/slurm_scripts/tmpfile_env --in 97308_Stage25_bbduk_aln.bam --in 103062_Stage25_bbduk_aln.bam --out Stage25_97308_103062.genrich '/scratch/Schmidbaur/Genrich/Genrich-master/Genrich -t %i1,%i2 -o %o1 -j -y -r -v' > stage25.genrich 

/proj/rpz/slurm_scripts/nq -f stage25.genrich -m 80 -o stage25.genrich
```

Stage 29

```
/proj/rpz/slurm_scripts/tmpfile_env --in 97309_Stage29_bbduk_aln.bam --in 103063_Stage29_bbduk_aln.bam --out Stage29_97309_103063.genrich '/scratch/Schmidbaur/Genrich/Genrich-master/Genrich -t %i1,%i2 -o %o1 -j -y -r -v' > stage29.genrich 

/proj/rpz/slurm_scripts/nq -f stage29.genrich -m 80 -o stage29.genrich
```
> **IMPORTANT NOTE: THE TWO INPUT FILES SHOULD ONLY BE SEPARATED BY ONE , NO SPACE INBETWEEN OTHERWISE CANNOT ALLOCATE THE TWO INPUT FILES, CONTINUES WITH AN ERROR**

###Check Insertsizes by plotting in R

```
for i in *bbduk_aln.bam;
do
/proj/rpz/slurm_scripts/tmpfile_env --in ${i} --out ${i}_fixed 'samtools fixmate %i1 %o1' > ${i}_fixed;
done

for i in *aln.bam_fixed;
do
/proj/rpz/slurm_scripts/nq -f $i -m 50 -o $i.slurm;
done
```
```
/proj/rpz/slurm_scripts/tmpfile_env --in 103063_Stage29_trim2_aln.bam --out 103063_Stage29_trim2_aln.bam_fixed 'samtools fixmate %i1 %o1' > 103063_Stage29_trim2_aln.bam.fixed

/proj/rpz/slurm_scripts/nq -f 103063_Stage29_trim2_aln.bam.fixed -m 50 -o 103063_Stage29_trim2_aln.bam.fixed.slurm
```

* to get insert sizes: 

```
for i in *aln_fixed.bam
do 
echo "samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) ${i} | awk '{print $9}' > ${i}_insertsizes.txt";
done
```

```
samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) 103061_stage20_1_bbduk_aln_fixed.bam | awk '{print $9}' > 103061_stage20_1_insertsizes.txt

samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) 103062_Stage25_bbduk_aln_fixed.bam | awk '{print $9}' > 103062_Stage25_insertsizes.txt

samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) 103063_Stage29_bbduk_aln_fixed.bam | awk '{print $9}' > 103063_Stage29_insertsizes.txt

samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) 97308_Stage25_bbduk_aln_fixed.bam | awk '{print $9}' > 97308_Stage25_insertsizes.txt

samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) 97307_Stage20_bbduk_aln_fixed.bam | awk '{print $9}' > 97307_Stage20_insertsizes.txt

samtools view -f $(samtools flags PROPER_PAIR,READ1 | cut -f 1) 97309_Stage29_bbduk_aln_fixed.bam | awk '{print $9}' > 97309_Stage29_insertsizes.txt
```


* plotting in R:

~~~r 
data=read.table("97309_stage29_insertsizes.txt")
#change data to positive values: 
data=abs(data) 
#histogramm uses matrix, dataframe conversion to matrix with as.matrix
h<-hist(as.matrix(data),breaks = 200)
#show with h enter
plot(h$mids,h$counts,log="y")
~~~

###GGPlot in R
 
~~~ R 
#simple histogram 

fragLenPlot <- table(data) %>% data.frame %>% rename(InsertSize = data, 
               Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
               Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
              geom_line()
#visualize 
fragLenPlot + theme_bw()

#log scaled 
fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()

#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + 
    geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw()

#and the same for log scaled

fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
~~~

###Fragment size distribution 

* negative fragmentsizes indicate reverse reads 
* rev reads can be taken out 
* to get rid of reverse reads: 

``` 
awk '$1>0 {print $1}' 97309_Stage29_insertsizes.txt > pos_97309_Stage29_insertsizes.txt

awk '$1>0 {print $1}' 103062_Stage25_insertsizes.txt > pos_103062_Stage25_insertsizes.txt

awk '$1>0 {print $1}' 97307_Stage20_insertsizes.txt > pos_97307_Stage20_insertsizes.txt

awk '$1>0 {print $1}' 97308_Stage25_insertsizes.txt > pos_97308_Stage25_insertsizes.txt

awk '$1>0 {print $1}' 103063_Stage29_insertsizes.txt > pos_103063_Stage29_insertsizes.txt
```


###Peak calling via Genrich of single alignments 


```
for i in $(ls *bam |rev |cut -c 14- |rev)
do 
/proj/rpz/slurm_scripts/tmpfile_env --in ${i}bbduk_aln.bam --out ${i}only.genrich '/scratch/Schmidbaur/Genrich/Genrich-master/Genrich -t %i1 -o %o1 -j -y -r -v' > ${i}genrich2.1; 
done

for i in *genrich2.1
do 
/proj/rpz/slurm_scripts/nq -f ${i} -m 80 -o ${i}.slurm;
done
```

> **IMPORTANT: SET MEMORY TO 80**


###Peak calling using MACS2 

* [protocol](https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html)
* M[ACS handbook](https://github.com/taoliu/MACS)
* adjustment of alignment files: 

```
for i in *bbduk_aln.bam;
do
/proj/rpz/slurm_scripts/tmpfile_env --in ${i} --out ${i}_MACS2 'macs2 callpeak  -t %i1  -f BAMPE  -n %o1 --keep-dup 1' > ${i}_MACS;
done

for i in *_MACS;
do
/proj/rpz/slurm_scripts/nq -f $i -m 80 -o $i.slurm;
done
```

* -f BAMPE: MACS will use real fragments from alignment results for reads pileup 
	* otherwise MACS will only keep left mate 5'end tag 
	* it will use insert size pairs of reads to build fragment pileup 
	* otherwise prediction of fragment size by distribution of plus and minus strand reads
* keep-dup 1 will remove all potential duplicates  
* MACS2 will generate several files 
* **DO NOT FORGET TO ADD * TO THE SAVING STEP AT THE END OF THE COMMAND**
* effective genome size: 
	* size of the genome which is mappable 
	* to get the size use khmer
* I counted the N numbers in genome fasta file: Eup_chr48.fasta 
 ``` grep -c 'N' Eup_chr48.fasta ```	
	* effective genome size: 5.1e9 

```
for i in *bbduk_aln.bam;
do
/proj/rpz/slurm_scripts/tmpfile_env --in ${i} --out ${i}_MACS2 'macs2 callpeak  -t %i1  -f BAMPE  -n %o1 -g 5.1e9 --keep-dup 1' > ${i}_MACS2new;
done

for i in *_MACS2new;
do
/proj/rpz/slurm_scripts/nq -f $i -m 80 -o $i.slurm;
done
```

###Peak calling using HOMER 

* [HOMER guide](http://homer.ucsd.edu/homer/ngs/peaks.html)
* make Tag Directory and find Peaks: 

```
in directory /scratch/Schmidbaur/homer
module load samtools
module load homer
makeTagDirectory tagdir20180124 ../Genrich/20180124.noadapt.bam -single 

findPeaks tagdir20180124 -o 20180124_homer.peaks -gsize 5.1e9 -minDist 150 -region
```
* R package requires the NarrowPeak format 
* python script for renaming the Lachesis_Group to chromosome and formating to [narrowPeak file](/Users/pui/Documents/Lab_Vienna/ATAC/Peakfiles/Annotations/HOMER_peakfiles/HOMER_rename.py) 
* new file has the formate: chromosome, start, end, peak_ID, Normalized Tag, strand, Count	size	findPeaks Score	Clonal Fold Change, first column with the initial names given by homer 
* R package requires also 10 columns 
* strand information should be in column 6 

###Motif search using HOMER

findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options] 

* motifsearch with masked genome: ```/scratch/Schmidbaur/genomes_softmasked19/Lachesis_assembly.fasta.masked```
* columns were rearranged to a bed file format in R 
	 * programme did not recognized the file format although correct 
	 * reformat using: ```awk '{printf ("%s\t%s\t%s\t \t%s\n", $1,$2,$3,$4,$6)}```

```
slurm for masked genome 
for i in $(ls Stage2*_peaks_bothRep_*.bed| rev| cut -c 5- |rev) ;do  echo "/proj/rpz/slurm_scripts/tmpfile_env --in ${i}.bed --out /${i)_Motif_Results/ '/apps/homer/4.9/bin/findMotifsGenome.pl %i1 /scratch/Schmidbaur/genomes_softmasked19/Lachesis_assembly.fasta.masked %o1 -size 100 -preparsedDir Lachesis_preparsed' > ${i}.slurm"; done

for i in *slurm;
do
/proj/rpz/slurm_scripts/nq -f $i -m 20 -o $i.slurm;
done

```
* this did not worked properly 
* repeat motifsearch with unmasked genome 
	* unmasked genome: ```/scratch/hoang/PASA/Eup_chr48.fasta ```
	* rearranging the columns to a bed file first 
	* some changes in slurm script had to be done 
	* also because unmasked genome has longer names in the fasta file
	* "Lachesis_group0__68_contigs__length_203855924" 
	* in the bed files the names are short "Lachesis_group0" 
	* names of bed files needed to be changed to the longnames
	* using a python script for that: ```/Users/pui/Documents/Lab_Vienna/ATAC/Cluster/Bedfiles_rename_Chr_longnames.py```
		* with dictionary for longnames 
		* and os library for looping through the files   

```
slurm for unmasked genome: 

for i in $(ls Stage2*_peaks_bothRep_*.bed| rev| cut -c 5- |rev); do awk '{printf ("%s\t%s\t%s\t \t%s\n", $1,$2,$3,$4,$6)}' ${i}.bed > ${i}; done


for i in $(ls Stage2*_peaks_bothRep_*.bed | rev|cut -c 5-|rev); do /proj/rpz/slurm_scripts/tmpfile_env --in $i --out ${i}_Motif_Results_unmasked '/apps/homer/4.9/bin/findMotifsGenome.pl %i1 /scratch/hoang/PASA/Eup_chr48.fasta %o1 -size 100 -preparsedDir Lachesis_preparsed_unmasked' > ${i}_unmasked.slurm; done

for i in *_unmasked.slurm; do /proj/rpz/slurm_scripts/nq -f $i -m 25 -o ${i}.slurm; done
```

* output example:

``` 
Stage25_peaks_bothRep_ceph $TMPDIR; /apps/homer/4.9/bin/findMotifsGenome.pl $TMPDIR/Stage25_peaks_bothRep_ceph /scratch/hoang/PASA/Eup_chr48.fasta $TMPDIR/Stage25_peaks_bothRep_ceph_Motif_Results_unmasked -size 100 -preparsedDir Lachesis_preparsed_unmasked; cp -r $TMPDIR/Stage25_peaks_bothRep_ceph_Motif_Results_unmasked .
```

* Next add synteny ID to the peak/bed files: 
* python script: [longnames_syntenyID_HS.py](/Users/pui/Downloads/longnames_syntenyID_HS.py)
	* first make dictionaries for all synteny IDs per chromosome for start and stop 
	*  also using a list of all chromosomes and defaultdict for dictionary of list 
	*  using range function and dictionaries to check whether start and stop are in a Range and then add synteny ID respectively 

* rerun Homer Motifsearch for new files:

```
for i in *.2txt; do /proj/rpz/slurm_scripts/tmpfile_env --in $i --out ${i}_MotifResults '/apps/homer/4.9/bin/findMotifsGenome.pl %i1 /scratch/hoang/PASA/Eup_chr48.fasta %o1 -size 100 -preparsedDir ../Lachesis_preparsed_unmasked/' > ${i}.slurm; done

for i in *.slurm; do /proj/rpz/slurm_scripts/nq -f $i -m 25 -o ${i}.slurm; done
```
 


###Differential Peak

* create bed files with peaks that are present in both replicates 
* then concatenate all stages into one bed file, sort them after chromosome, merge peaks that overlap : 

```
cat /scratch/hoang/ATACSeq/HOMER_Motifsearch/Stage2*_peaks_bothRep_all.longnames.bed2 > all_stages

bedtools sort -i all_stages > all_stages_sorted

bedtools merge -c 4,5,6,7 -o distinct,distinct,distinct,distinct -i all_stages_sorted > all_stages_merged_allcols.bed

awk '{printf "%s\t%s\t%s\t%s\t-\n", $6,$1,$2,$3}' all_stages_merged_allcols.bed >all_stages_merged.saf

```
* sort: sorts the file after chromosomes, needed for merging 
* merging reads:
	* c: columns that should also be merged 
	* o: must be given when using -c option, here using distinct to keep peak ID, empty column, strand, synteny ID without duplicates 
* awk allows to reorder columns
* run featureCounts, from subreads 

```
module load subread
featureCounts -F 'SAF' -a all_stages_merged.saf -o featureCounts_results /scratch/hoang/ATACSeq/alignment_results/103061_stage20_1_bbduk_aln.bam /scratch/hoang/ATACSeq/alignment_results/103061_stage20_1_bbduk_aln.bam scratch/hoang/ATACSeq/alignment_results/97308_Stage25_bbduk_aln.bam /scratch/hoang/ATACSeq/alignment_results/97307_Stage20_bbduk_aln.bam /scratch/hoang/ATACSeq/alignment_results/103062_Stage25_bbduk_aln.bam /scratch/hoang/ATACSeq/alignment_results/97309_Stage29_bbduk_aln.bam /scratch/hoang/ATACSeq/alignment_results/103063_Stage29_bbduk_aln.bam
```
* featureCounts results looked weird 
* something went wrong with concatenating the files: 

```
cat *.genrich > all_genrich.bed

sort -k1,1 -k2,2n all_genrich.bed > allpeaks.sorted.bed 

bedtools merge -c 4,5,6,7,8,10 -o collapse,collapse,collapse,collapse,collapse,collapse, -i allpeaks.sorted.bed >allpeaks.merged.bed

awk '{printf "%s\t%s\t%s\t%s\t-\n", $4,$1,$2,$3}' allpeaks.merged.bed  > allpeaks.merged.saf
```
* then run again featureCount
* done by Hannah: 
```/scratch/Schmidbaur/Euprymna_atac/atac_sequencing_2019/Genrich_peakcalling_single/featureCounts```

------
Meeting with Hannah:

new_names <- c(rep("Stage20_103061",6)),
                rep("Stage25_97308",6),
                rep("Stage20_97307",6),
            "Stage25_103062","Stage29_97309","Stage29_103063")
colnames(mat_data) <- new_names

condition <- c("Stage20","Stage25","Stage20","Stage25","Stage29","Stage29")


* all genes with peaks, ceph meta, etc, 
* list with genes only in ceph, meta, non-syntenic

* scan("ceph_list", what=character(1), sep = ",")
* ontology biological processes, molecular function

###Handling replicates with IDR 

* until now annotation was done by selecting overlapping peaks in the two biological replicates 
* variability in experiments especially in high throughput methods 
* therefore 2-3 biological replicates
* more infos [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html)
* statistical testing and evaluation of reproduceability between replicates  using the irreproducibility discovery rate framework 
* or bedtools analysis (link to tutorial in linke above) 
* comparison of ranked lists of regionspeaks and assigning values reflecting its reproducibility 
* replicates with same underlying biology with significant peaks = real signals expected to have high consistency between replicates 
* peaks with low significance = noise have low consistency
* consistency between a pair of ranked peaks contain significant and insignificant plotted a transition in consistency is expected 
* consistency transition provides internal indicator of change from signal to noise and suggests how many peaks have been reliably detected 
* WHY? : 
	* avoids choices of initial cutoffs, not compareable between callers 
	* not dependent on arbitrary tresholds all regionspeaks are considered 
	* based on ranks, no input signal calibration required 
* Components: 
	* correspondence curve: matching peaks along the ranked list and selecting, qualitative 
	* inference procedure: summarizes proportion of reproducible and irreproducible signals   
	* irreproducible discovery rate: derives significance value from the inference procedure similar to FDR, can be used to control the level of irreproducibility rate when selecting signals 
	* 0.05 IDR means peak has a 5% chance of being an irreproducible discovery

* first sort genrich narrowPeak files by p-values
* p-values already calculated by genrich 
* sort genrich files by column 8:

```
sort -nk8,8nr ../genrich_peak_files/97309_Stage29_only.narrowPeak >  97309_Stage29_only.sorted.np
```
* then get reproducible peaks with p<= 0.05 : 

```
idr --samples 103063_Stage29_only.sorted.np 97309_Stage29_only.sorted.np  --input-file-type narrowPeak --rank p.value --output-file stage_29_idr --plot --log-output-file stage_29_idr.log  --idr-threshold 0.05
```
* done by Hannah
######Results: 

* upper left: replicate 1 vs replicate 2, red = peaks with higher IDR as 0.05
* upper left: with log10 scores rep1 vs rep2 
* bottom row: Peak rank versus IDR scores are plotted in black. The overlayed boxplots display the distribution of idr values in each 5% quantile. The IDR values are thresholded at the optimization precision - 1e-6 by default.

**Stage 20**

*  Number of peaks passing IDR cutoff of 0.05 - 11454/22894 (50.0%)

**Stage 25**

* Number of peaks passing IDR cutoff of 0.05 - 12873/20835 (61.8%)

**Stage 29**

* Number of peaks passing IDR cutoff of 0.05 - 10807/20171 (53.6%)

* in ATAC Sequencing generally high variability due to: 
	* Tn5 sensitivity, changing activity 
	* differences in Illumina kits itself 
	* depends on single cell dissociation: dead cells, cell clumps 
	* no sorting of cells taken here, whole embryo dissociated 

###GO analysis: 

* carried out using R script 
* do GO analysis independent of the different stages again 
* results only show very general terms 


###Useful stuff: 

* [Command-Line Fundamentals ](https://www.youtube.com/watch?v=HVsySz-h9r4&t=1s&frags=pl%2Cwn)
* [Analysis of ATAC-seq data in R and Bioconductor](https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html)
* [Video](https://vimeo.com/428802218) from spiralia base on ATAC Seq

##LAB MEETING ATAC Seq

* to infer open chromatin regions
* transposase cuts open accessible regions 
* sequence and realignment 
* stage 20 early organogenesis 
* stage 25 late organogenesis 
* stage 29 ready to hatch 

* 30% aligned to several areas but expected because of repetitive 
* excluding or including multimappers 
* genrich not more filtering after alignment 
* removing unmapped reads duplicates 
* fragment length expected peaks 
* sequence artefact because sequence 125 
* internal assembler artefacts 
* MACS2 peak calling needs shifting of reads because of transposase affinity 
* IDR: compares peak files after peak calling 
	* same peaks in both replicates and a peak score of reliability between two differen samples 
* diff binds
* frip score
* samßple correlation 
* FDR 
* differentially open regions: heatmap shows counts 
* ceph 
* in conserved synteny chitin and actin binding 
	* chitin important in cephalopod synteny 
	* differentially open regions to do on GO 
	* protein binding what are the genes  
* tobias: atac seq data to look at TF chromatin, protein bound or nucleosome 
* regions of depleted cutting transposases not cutting 
* estimation of expected transposase cut site 
* correlation of rtesults with TF 