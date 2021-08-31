#! usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\

 R E V I V I D - W G S - P I P E L I N E
 =======================================
 #######################################
 ############### DSL 2 #################
 ############### W I P #################
 ############### alpha #################
 #######################################
 =======================================
 """

include { importfastq } from './importfastq.nf'
include { fastQC } from './fastQC.nf'
include { pear } from './pear.nf'
include { alignment } from './alignment.nf'
include { readgroups } from './readgroups.nf'
include { duplicates } from './duplicates.nf'
include { mergebams } from './mergebams.nf'
include { generateCRAM } from './generateCRAM.nf'
include { baserecalibrator } from './baserecalibrator.nf'
include { delete_file } from './delete_file.nf'
include { checkbam } from './checkbam.nf'
include { applyBQSR } from './applyBQSR.nf'
include { genotype } from './genotype.nf'
include { variantrecalibration } from './variantrecalibration.nf'
include { compressandindex } from './compressandindex.nf'
include { mergevcf } from './mergevcf.nf'
include { test } from './test.nf'

// script parameters
params.fastqgz = '/mnt/hdd2/data/lau/phd/testyard/FNRCP/FASTQ/*/*/*.fastq.gz'
fastqgz_ch = channel.fromPath(params.fastqgz)


indexes_ch = Channel.fromPath(params.indexes).toList()

donebams_ch = channel.fromPath('./results/bams/*.bam*').toSortedList().flatten().collate( 2 ).map{bam,bai -> tuple(bam.simpleName,bam,bai)}.flatten().collate( 3 )

garbage_ch = Channel.empty()

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

Channel.fromList(ids).map { it -> [it, familymap[it]] }set{ idfamily_ch }


//workflows

workflow download_fastq_to_bam_and_cram {
take: idfam
main:

//importfastq(idfam, params.home) 
//importfastq.out.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}

fastqgz_ch.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3),file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
//.map{file -> tuple(file.getBaseName(3), file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
//gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(R1,R2)}.flatten().toList().size.view()



fastQC(gzipped_ch)
pear(gzipped_ch, params.home)
//if ( pear.out[0].flatten().toList().size.view() > 0 ) {
//	delete_file(gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(R1,R2)}.flatten())
//}
alignment(pear.out, params.genome,indexes_ch, params.home)
readgroups(alignment.out,params.home)
duplicates(readgroups.out,params.home)

mergebams(duplicates.out[0].groupTuple(),params.home)
generateCRAM(mergebams.out[0],params.genome,indexes_ch)

garbage_ch.concat(gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.flatten().toList(),pear.out.flatten().collate( 5 ).map{id,lane,paired,forward,reverse -> tuple(lane,paired,forward,reverse)},alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane,bam)},readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(lane,bam,bai)}).groupTuple().flatten().unique().toList().dump(tag:garbage).view().set{workflow1garbage}

//delete_file(workflow1garbage,mergebams.out[1])

emit:
bams = mergebams.out[0]
crams = generateCRAM.out[0]
garbage = workflow1garbage

}
workflow testwf {
take: 
data
bams
main:

test(data)

}

workflow createvcfs {
take: bam 

main:
baserecalibrator(bam,params.genome, indexes_ch, params.genomedict, params.snps, params.snpsindex)
applyBQSR(baserecalibrator.out,params.genome,indexes_ch,params.genomedict)
genotype(applyBQSR.out,params.genome,indexes_ch,params.genomedict,params.mask)
variantrecalibration(genotype.out,params.genome,params.genomedict,indexes_ch,params.snps, params.snpsindex,params.indels,params.indelsindex,params.mask)
compressandindex(variantrecalibration.out)

mergevcf(idfamily_ch.join(compressandindex.out).map{ id, family, vcf, index -> tuple(family,vcf,index)}.groupTuple())


emit:
mergevcf.out
}

workflow { 
main:
checkbam(idfamily_ch)
checkbam.out.test_ch.filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{done_ch}
done_ch.toSortedList().flatten().collate(1).combine(donebams_ch, by:0).map{id,bam,bai -> tuple(id,bam,bai)}.set{alldone_ch}
download_fastq_to_bam_and_cram(checkbam.out.test_ch.filter( ~/.*todo.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id,family)})
download_fastq_to_bam_and_cram.out.bams.concat(alldone_ch).set{mixed}
testwf(download_fastq_to_bam_and_cram.out.garbage,download_fastq_to_bam_and_cram.out.bams)
createvcfs(mixed)

}