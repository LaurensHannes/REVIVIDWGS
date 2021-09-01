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
include { SelectVariantsdenovo } from './SelectVariantsdenovo.nf'
include { SelectVariantsAR } from './SelectVariantsAR.nf'
include { annotate as annotatedenovo; annotate as annotateAR } from './annotate.nf'


// script parameters


indexes_ch = Channel.fromPath(params.indexes).toList()
donebams_ch = channel.fromPath('./results/bams/*.bam*').toSortedList().flatten().collate( 2 ).map{bam,bai -> tuple(bam.simpleName,bam,bai)}.flatten().collate( 3 )
garbage_ch = Channel.empty()
testcollection = Channel.empty()
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
importfastq(idfam, params.home) 
importfastq.out.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
fastQC(gzipped_ch)
pear(gzipped_ch, params.home)
gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.join(fastQC.out.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane)}).join(pear.out.flatten().collate( 5 ).map{id,lane,paired,forward,reverse -> tuple(lane)}).set{testgarbage_ch}

alignment(pear.out, params.genome,indexes_ch, params.home)
readgroups(alignment.out,params.home)
duplicates(readgroups.out,params.home)

mergebams(duplicates.out[0].groupTuple(),params.home)
generateCRAM(mergebams.out[0],params.genome,indexes_ch)

garbage_ch.concat(gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.flatten().toList(),pear.out.flatten().collate( 5 ).map{id,lane,paired,forward,reverse -> tuple(lane,paired,forward,reverse)},alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane,bam)},readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(lane,bam,bai)}).groupTuple().dump(tag:"garbage").set{workflow1garbage}
duplicates.out[0].flatten().collate ( 2 ).map{lane,bam -> tuple(bam.getBaseName(2))}.join(workflow1garbage).flatten().dump(tag:"merged").set{garbagemerge}
testcollection.concat(testgarbage_ch)

emit:
bams = mergebams.out[0]
crams = generateCRAM.out[0]
garbage = garbagemerge
testgarbage = testcollection
}
workflow testwf {
take: 
data
test
main:
test(test)
delete_file(data)

}

workflow createvcfs {
take: bam 

main:
baserecalibrator(bam,params.genome, indexes_ch, params.genomedict, params.snps, params.snpsindex)
applyBQSR(baserecalibrator.out,params.genome,indexes_ch,params.genomedict)
genotype(applyBQSR.out,params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)
variantrecalibration(genotype.out,params.genome,params.genomedict,indexes_ch,params.snps, params.snpsindex,params.indels,params.indelsindex,params.mask)
compressandindex(variantrecalibration.out)

mergevcf(idfamily_ch.join(compressandindex.out).map{ id, family, vcf, index -> tuple(family,vcf,index)}.groupTuple())


emit:
triovcf = mergevcf.out
}

workflow trioVCFanalysis {
take:vcf

main:
SelectVariantsdenovo(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(SelectVariantsAR.out[0],params.programs,params.humandb,params.annovardbs)
}



workflow { 
main:
checkbam(idfamily_ch)
checkbam.out.test_ch.filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{done_ch}
done_ch.toSortedList().flatten().collate(1).combine(donebams_ch, by:0).map{id,bam,bai -> tuple(id,bam,bai)}.set{alldone_ch}
download_fastq_to_bam_and_cram(checkbam.out.test_ch.filter( ~/.*todo.*/ ).dump(tag:"todo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id,family)})
download_fastq_to_bam_and_cram.out.bams.concat(alldone_ch).set{mixed}
testwf(download_fastq_to_bam_and_cram.out.garbage,download_fastq_to_bam_and_cram.out.testgarbage)
createvcfs(mixed)
trioVCFanalysis(createvcfs.out.triovcf)

}