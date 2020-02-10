#!/bin/sh


#######################################################################################################################
################### prepare genomes ###################

mkdir genome

cd genome

# Download genome

wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.17_FB2017_04/fasta/dmel-all-chromosome-r6.17.fasta.gz

wget ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_r1.07_FB2018_05/fasta/dvir-all-chromosome-r1.07.fasta.gz

# Extract fasta file

gunzip dmel-all-chromosome-r6.17.fasta.gz

gunzip dvir-all-chromosome-r1.07.fasta.gz

# rename fasta to fa

samtools faidx dmel-all-chromosome-r6.17.fasta 2L 2R 3L 3R 4 X Y | sed 's/^>/>chr/' > dmel-all-chromosome-r6.17.fa

mv dvir-all-chromosome-r1.07.fasta dvir-all-chromosome-r1.07.fa

samtools faidx dvir-all-chromosome-r1.07.fa

# build bowtie2 index

bowtie2-build dmel-all-chromosome-r6.17.fa dmel_genome

bowtie2-build dvir-all-chromosome-r1.07.fa dvir_genome


#######################################################################################################################



#######################################################################################################################
################### get data drom SRA ###################

### update when SRR ids are available !!!





#######################################################################################################################



#######################################################################################################################
################### bowtie2 alignment ###################


cd ../FastQ

for file in *.txt.gz; do

	base=`echo ${file} | sed -e "s/.txt.gz//g"`

	# align to dmel

	bowtie_index="../genome/dmel_genome"

	bowtie2 -p 12 -x $bowtie_index -U ${base}.txt.gz  > ${base}.dmel.sam 2> ${base}.dmel.txt

	samtools view -bS -@ 12 ${base}.dmel.sam | samtools sort - -@ 12 | tee ${base}.dmel.bam | samtools index - ${base}.dmel.bam.bai

	rm ${base}.dmel.sam

	# align to dvir

	bowtie_index="../genome/dvir_genome"

	bowtie2 -p 12 -x $bowtie_index -U ${base}.txt.gz  > ${base}.dvir.sam 2> ${base}.dvir.txt

	samtools view -bS -@ 12 ${base}.dvir.sam | samtools sort - -@ 12 | tee ${base}.dvir.bam | samtools index - ${base}.dvir.bam.bai

	rm ${base}.dvir.sam


done


#######################################################################################################################


#######################################################################################################################
################### makeTagDirectory ###################


mkdir ../Homer

cd ../Homer

mv ../FastQ/*.bam .
mv ../FastQ/*.bai .

for bam in *.dmel.bam; do

    base=`echo $bam | sed 's/.dmel.bam//'`

	dir_base=`echo $base | sed 's/_[G,A,T,C][G,A,T,C][G,A,T,C][G,A,T,C].*//'`

	# calculate total number of mapped reads

    mapped=`samtools view -c -F 4 ${base}.dvir.bam | bc`

	# use totalReads for normalization

    makeTagDirectory ${dir_base}.dvirNorm.dir -totalReads $mapped ${base}.dmel.bam

	makeTagDirectory ${dir_base}.dmelNorm.dir ${base}.dmel.bam

done


#######################################################################################################################




#######################################################################################################################
################### makeUCSCfile ###################


### update when SRR ids are available !!!

# iterate through IP directories

for chip in *_IP*.dir; do

	# replace IP by Input

	input=`echo $chip | sed 's/_IP/_INP/'`

	echo $chip $input

	chip_base=`echo $chip | sed 's/.dir//'`

	makeUCSCfile ${chip} -fragLength 150 -i ${input}  -o ${chip_base}.INPnorm.bedgraph

done

#######################################################################################################################





#######################################################################################################################
