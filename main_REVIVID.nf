#! usr/bin/env nextflow
// script parameters

workDir = '/staging/leuven/stg_00086/Laurens/FNRCP/work'
params.genome ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.fasta'
params.genomefai ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai'
params.genomedict ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.dict'
params.indexes = '/home/laurens/data/resources/reference_genome/broadhg38/Homo_sapiens_assembly38.fasta.*'
params.ped ='/staging/leuven/stg_00086/Laurens/REVIVID/F1.ped'
params.mask ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/mask.bed'
params.home ='/staging/leuven/stg_00086/Laurens/FNRCP'
params.fastqgz = '/staging/leuven/stg_00086/Laurens/FNRCP/raw_data/FASTQGZ/*/*/*.fastq.gz'
fastqgz_ch = channel.fromPath(params.fastqgz)

//script

myFile = file(params.ped)
myReader = myFile.newReader()
String line
familymap = [:]
ids = []
while( line = myReader.readLine() ) {
(empty, family, id, father, mother, sex, phenotype) = (line =~ /(^F\d{1,2})\t(GC\d+)\t(\w+)\t(\w+)\t(\d+)\t(\d+)/)[0]
        familymap[id]=family
        ids << id
}
myReader.close()

Channel.fromList(ids).map { it -> [it, familymap[it]] }into{ idfamily_ch1; idfamily_ch2; idfamily_ch3; idfamily_ch4; idfamily_ch5 }

//processes



fastqgz_ch.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().set{gzipped_ch}
gzipped_ch.into{temp_ch1;temp_ch2}
temp_ch2.view()


process pear {

        tag "$lane"
        time '30m'
		memory '2 GB'
        errorStrategy 'retry'
        maxRetries 3

 //       publishDir '/staging/Leuven/stg_00086/Laurens/FNRCP/paired'

        input:
        tuple val(id), val(lane),file(R1), file(R2) from temp_ch1.flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}

        output:
        tuple val(id), val(lane), file('*.assembled.fastq'),file('*.forward.fastq'), file('*.reverse.fastq') into paired_ch

        """
        pear -j ${task.cpus} -f $R1 -r $R2 -o $lane
        """
}
 process alignment {

		time '2h'
        maxForks 3
        tag "$id"

        indexes_ch = Channel.fromPath(params.indexes)

//        publishDir '/staging/Leuven/stg_00086/Laurens/FNRCP/mapped'

        input:
        tuple val(id),val(lane), file(assembled), file(forward), file(reverse) from paired_ch
        path genome from params.genome
        path indexes from indexes_ch.toList()
        path home from params.home

        output:
        tuple val(id), file('*.indexed.bam') into mapped_ch

        """
        bwa mem -r 0.5 -t ${task.cpus} $genome $assembled | samtools view -@ ${task.cpus} -bS > ${lane}.assembled.unsorted.bam
        bwa mem -r 0.5 -t ${task.cpus} $genome $forward | samtools view -@ ${task.cpus} -bS > ${lane}.forward.unsorted.bam
        bwa mem -r 0.5 -t ${task.cpus} $genome $reverse | samtools view -@ ${task.cpus} -bS > ${lane}.reverse.unsorted.bam
        samtools merge -@ ${task.cpus} ${lane}.indexed.unsorted.bam  ${lane}.assembled.unsorted.bam ${lane}.forward.unsorted.bam ${lane}.reverse.unsorted.bam
        samtools sort -@ ${task.cpus} -o ${lane}.indexed.bam ${lane}.indexed.unsorted.bam
        """
}

mapped_ch.into{mapped_ch1;mapped_ch2}

mapped_ch2.groupTuple().set{mappedgrouped_ch}

process mergebams {

	tag "$id"
        

	input:
	tuple val(id),file(bams) from mappedgrouped_ch
		
	output:
	tuple val(id),file("${id}.bam") into mergedbam_ch

	"""
	samtools merge -@ 8 ${id}.bam $bams
	"""

}

//process duplicates { }
mergedbam_ch.into{mergedbam1_ch;mergedbam2_ch}
mergedbam2_ch.view()

process genotype {

        tag "$id"
		cpus 4

        input:
        tuple val(id), path(bams) from mergedbam1_ch
        path genome from params.genome
        path dict from params.genomedict
//        path index from params.indexes
        path faidx from params.genomefai
		path mask from params.mask

        output:
        tuple val(id), file("${id}.vcf") into vcf_uncallibrated_ch

        """
	gatk AddOrReplaceReadGroups -I $bams -O ${id}.RG.bam -LB REVIVID -PL ILLUMINA -PU $id -SM $id --MAX_RECORDS_IN_RAM 200000000
	samtools index -@ 8 ${id}.RG.bam
        gatk HaplotypeCaller --verbosity INFO -XL $mask -R $genome -I ${id}.RG.bam -O ${id}.vcf --sequence-dictionary ${dict}
        """

}