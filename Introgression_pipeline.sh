###############################################################################################
# draft pipeline to test for signals of past introgression between species
# using target sequence data. The code is developed using secapr (Andermann et al., 2019)
# as a base. This code is organized and maintained by Maria Fernanda Torres
# this pipeline will become a snakemake recipe in the future
###############################################################################################

# 0. download SRA data
for i in SRA_ARRAY; do bin/fasterq-dump $il done

#----------------------------------------------------------------------------------------------
# 1. set up the environment and create the index for the reference file
# conda environment = snakemake + secapr (development)
REF='../LOC_sequences.fasta'

bwa index $REF
picard CreateSequenceDictionary R=$REF O=$REF

# renaming files from SRR code to species + SRR code. renaming_files is a tab delimited file with two columns
while IFS='' read -r line; do NR=$(echo "$line" | cut -f 1); NAME=$(echo "$line" | cut -f 2); cp "${NR}" "${NAME}"; done < renaming_files.txt

#----------------------------------------------------------------------------------------------
# 2. Loops to sort, fix, mark PCR duplicates, and extract reads back to fastq files
# 2.1. Map reads to reference
for i in *_1.fastq; do SAMPLE=$(ls ${i} | sed -r 's/_1.fastq//g'); echo ${SAMPLE}; bwa mem -t 2 $REF ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq | picard SortSam I=/dev/stdin O=${SAMPLE}.srtd.bam SORT_ORDER=coordinate; done

# 2.2. Add group information to the BAM files
for i in *.srtd.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.bam//g'); echo ${SAMPLE}; GROUP=$(echo ${SAMPLE}); picard AddOrReplaceReadGroups I=${SAMPLE}.srtd.bam O=${SAMPLE}.srtd.wgrps.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${GROUP} $2> ${SAMPLE}_groupfix.out; done

# 2.3. Sort BAMs by name
for i in *.srtd.wgrps.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.wgrps.bam//g'); echo ${SAMPLE}; samtools sort -@ 2 -n -o "${SAMPLE}".srtdname.wgrps.bam ${SAMPLE}.srtd.wgrps.bam; done

# 2.4. Fix read mate information
for i in *.srtd.wgrps.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.wgrps.bam//g'); echo ${SAMPLE}; samtools fixmate -@ 2 -m ${SAMPLE}.srtdname.wgrps.bam "${SAMPLE}".srtdname.wgrps.fxmt.bam; done

# 2.5. Sort BAMs by coordinate
for i in *.srtd.wgrps.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.wgrps.bam//g'); echo ${SAMPLE}; samtools sort -@ 2 -o "${SAMPLE}".wg.srtcrd.fix.bam ${SAMPLE}.srtdname.wgrps.fxmt.bam; done

# 2.6. Mark duplicates
for i in *.srtd.wgrps.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.wgrps.bam//g'); echo ${SAMPLE}; samtools markdup -@ 2 ${SAMPLE}.wg.srtcrd.fix.bam "${SAMPLE}".wg.srtcrd.fix.mkdp.bam; done

# 2.7 Sort BAMs by name
for i in *.srtd.wgrps.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.wgrps.bam//g'); echo ${SAMPLE}; samtools sort -@ 2 -n -o "${SAMPLE}".wgfixm.srtname.mkdp.bam ${SAMPLE}.wg.srtcrd.fix.mkdp.bam; done

# 2.8. Extract unique reads to fastq
for i in *.srtd.wgrps.bam; do SAMPLE=$(ls ${i} | sed -r 's/.srtd.wgrps.bam//g'); echo ${SAMPLE}; samtools fastq -@ 2 -1 "${SAMPLE}"_nodup_1.fastq -2 "${SAMPLE}"_nodup_2.fastq ${SAMPLE}.wgfixm.srtname.mkdp.bam; done

#----------------------------------------------------------------------------------------------
# 3.1 Quality checks and trimming adapters if necesary
# move nodup reads to the nodupreads folder
secapr quality_check --input nodupreads/ --output raw/

# 3.2. summarize things beyond secapr's summary
multiqc raw/

# 3.3 use index3file.py to create the config.txt file required by secapr. Separate single and double index files in different folders
# https://github.com/mftorres/randomTools/blob/master/index2file.py
# inside the folder with the reads -- nodupreadssingle:
python index2file.py

# trim adapters with secapr and trimmomatic, e.g.
secapr clean_reads --input nodupreadssingle/ --config nodupreadssingle/config.txt --output cleanreadssingle/ --index single [double]

#----------------------------------------------------------------------------------------------
# 4. Make data treated outside secapr compatible with secapr (this won't be needed in the future...)
# 4.1 If reads were cleaned separatedly, make sure they are inside a folder called ${sample}_clean and reads are named ${sample}_clean-READ1.fastq

for i in *_nodup_1.fastq; do
  SAMPLE=$(echo ${i} | sed -r 's/_nodup_1.fastq//g')
  mkdir "${SAMPLE}"_clean
  mv $SAMPLE* ${SAMPLE}_clean
done

# 4.2 rename read files to fit secapr's pattern
# rename 's/nodup_/clean-READ/g' *_clean/*fastq # depending in your rename version
rename nodup_ clean-READ *_clean/*fastq

#----------------------------------------------------------------------------------------------
# 5. Assemble reads using SPAdes in a similar manner to how secapr does it (if samples treated outside secapr, SPAdes must be used separatedly
# secapr requires a stats file. This will be fixed in the future
# run it from the folder with the ${sample}_clean folders

for i in *_clean; do
	SAMPLE=$(echo ${i} | sed -r 's/_clean//g'); echo $SAMPLE;
	spades.py -1 ${SAMPLE}_clean/${SAMPLE}_clean-READ1.fastq -2 ${SAMPLE}_clean/${SAMPLE}_clean-READ2.fastq --only-assembler -o ${SAMPLE}_clean/${SAMPLE}_out;
done

#----------------------------------------------------------------------------------------------
# 6. reformat folders and rename contig files to make sure secapr can extract target contigs

mkdir contigs;
for i in *_clean; do
	SAMPLE=$(echo ${i} | sed -r 's/_clean//g'); echo $SAMPLE;
  cp ${SAMPLE}_clean/*out/contigs.fasta contigs/"${SAMPLE}".fa # replace cp wit mv if you feel adventurous
done

#----------------------------------------------------------------------------------------------
# 7. Extract contigs with find_target_contigs keeping the contig with the best bitscore and keeping paralogs
# to do this, you need to use the hyb_assembly branch of secapr on my repository:
# https://github.com/mftorres/seqcap_processor/tree/hyb_assembly

secapr find_target_contigs --contigs contigs --reference $REF --output extracted_folder --min_identity 90 --keep_bestbitscore --keep_paralogs --blast_threads 14

#----------------------------------------------------------------------------------------------
# 8. visualise the paralogs per sample per locus somehow

# in development
