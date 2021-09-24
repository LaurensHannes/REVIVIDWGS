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

include { importfastq } from './modules/importfastq.nf'
include { fastQC } from './modules/fastQC.nf'
include { pear } from './modules/pear.nf'
include { alignment } from './modules/alignment.nf'
include { readgroups } from './modules/readgroups.nf'
include { duplicates } from './modules/duplicates.nf'
include { mergebams } from './modules/mergebams.nf'
include { generateCRAM } from './modules/generateCRAM.nf'
include { baserecalibrator } from './modules/baserecalibrator.nf'
include { delete_file } from './modules/delete_file.nf'
include { checkbam } from './modules/checkbam.nf'
include { applyBQSR } from './modules/applyBQSR.nf'
include { genotype } from './modules/genotype.nf'
include { variantrecalibration } from './modules/variantrecalibration.nf'
include { compressandindex } from './modules/compressandindex.nf'
include { mergevcf } from './modules/mergevcf.nf'
include { test } from './testmodules/test.nf'
include { SelectVariantsdenovo } from './modules/SelectVariantsdenovo.nf'
include { SelectVariantsAR } from './modules/SelectVariantsAR.nf'
include { annotate as annotatedenovo; annotate as annotateAR } from './modules/annotate.nf'


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
importfastq(idfam, params.home
importfastq.out.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.dump(tag:"test").groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
fastQC(gzipped_ch)
pear(gzipped_ch, params.home)
gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.dump(tag:"gzippedG").join(fastQC.out.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane)}.dump(tag:"fastQCG")).join(pear.out.flatten().collate( 5 ).map{id,lane,paired,forward,reverse -> tuple(lane)}.dump(tag:"pearG")).dump(tag:"G1").set{testgarbage_ch}

alignment(pear.out, params.genome,indexes_ch, params.home)
pear.out.flatten().collate( 5 ).map{id,lane,paired,forward,reverse -> tuple(lane,paired,forward,reverse)}.join(alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane)}).dump(tag:"garbage2").set{testgarbage_ch2}

readgroups(alignment.out,params.home)
alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane,bam)}.join(readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(lane)}).dump(tag:"garbage3").set{testgarbage_ch3}

duplicates(readgroups.out,params.home)

readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(id,bam,bai)}.dump(tag:"garbage4part1").join(duplicates.out[0].flatten().collate ( 2 ).map{id,bam -> tuple(bam.getBaseName(2)).dump(tag:"garbage4part2")}.dump(tag:"garbage4part3")).dump(tag:"garbage4").set{testgarbage_ch4}



mergebams(duplicates.out[0].groupTuple(),params.home)
generateCRAM(mergebams.out[0],params.genome,indexes_ch)
duplicates.out[0].flatten().collate ( 2 ).map{id,bam -> tuple(id,bam)}.join(mergebams.out[0].flatten().collate ( 3 ).map{id,bam,bai -> tuple(id)}).join(generateCRAM.out[0].flatten().collate ( 3 ).map{id,cram,crai -> tuple(id)}).dump(tag:"garbage5").set{testgarbage_ch5}
garbage_ch.concat(gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.flatten().toList(),pear.out.flatten().collate( 5 ).map{id,lane,paired,forward,reverse -> tuple(lane,paired,forward,reverse)},alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane,bam)},readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(lane,bam,bai)}).groupTuple().dump(tag:"garbage").set{workflow1garbage}
duplicates.out[0].flatten().collate ( 2 ).map{lane,bam -> tuple(bam.getBaseName(2))}.join(workflow1garbage).flatten().dump(tag:"merged").set{garbagemerge}
testcollection.concat(testgarbage_ch,testgarbage_ch2,testgarbage_ch3).dump(tag:"G12345")

emit:
bams = mergebams.out[0]
crams = generateCRAM.out[0]
garbage = garbagemerge
testgarbage = testgarbage_ch
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
testwf(download_fastq_to_bam_and_cram.out.garbage,download_fastq_to_bam_and_cram.out.testgarbage.flatten())
createvcfs(mixed)
trioVCFanalysis(createvcfs.out.triovcf)

}