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
        cpus 2
        errorStrategy 'retry'
        maxRetries 3

        publishDir '/staging/Leuven/stg_00086/Laurens/FNRCP/paired'

        input:
        tuple val(id), val(lane),file(R1), file(R2) from temp_ch1.flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}

        output:
        tuple val(id), val(lane), file('*.assembled.fastq'),file('*.forward.fastq'), file('*.reverse.fastq') into paired_ch

        """
        pear -v 0 -f $R1 -r $R2 -o $lane
        """
}
 process alignment {

        maxForks 3
        tag "$x"

        indexes_ch = Channel.fromPath(params.indexes)

        publishDir '/staging/Leuven/stg_00086/Laurens/FNRCP/mapped'

        input:
        tuple val(id),val(lane), file(assembled), file(forward), file(reverse) from paired_ch
        path genome from params.genome
        path indexes from indexes_ch.toList()
        path home from params.home

        output:
        tuple val(x), file('*.indexed.bam') into mapped_ch

        """
        bwa mem -r 0.5 -t 4 $genome $assembled | samtools view -@ 4 -bS > ${x}.assembled.unsorted.bam
        bwa mem -r 0.5 -t 4 $genome $forward | samtools view -@ 4 -bS > ${x}.forward.unsorted.bam
        bwa mem -r 0.5 -t 4 $genome $reverse | samtools view -@ 4 -bS > ${x}.reverse.unsorted.bam
        samtools merge -@ 4 ${x}.indexed.unsorted.bam  ${x}.assembled.unsorted.bam ${x}.forward.unsorted.bam ${x}.reverse.unsorted.bam
        samtools sort -@ 4 -o ${x}.indexed.bam ${x}.indexed.unsorted.bam
        """
}

process genotype {

        memory '16 GB'
        tag "$bams.baseName"
        maxForks 4

        publishDir '/staging/Leuven/stg_00086/Laurens/FNRCP/genotyped'

        input:
        tuple val(x), file(bams) from mapped_ch
        path genome from params.genome
        path dict from params.genomedict
        path index from params.indexes
        path faidx from params.genomefai


        output:
        tuple val(x), file('*.g.vcf') into gvcfs_ch

        """
        /usr/roadrunner/programs/gatk-4.1.8.1/gatk HaplotypeCaller --verbosity WARNING  -R $genome -I ${bams} -O ${bams}.g.vcf --sequence-dictionary ${dict} -ERC GVCF
        """

}
