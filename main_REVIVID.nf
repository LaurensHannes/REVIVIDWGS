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
include { alignment } from './modules/alignment.nf'
include { readgroups } from './modules/readgroups.nf'
include { duplicates } from './modules/duplicates.nf'
include { mergebams } from './modules/mergebams.nf'
include { generateCRAM } from './modules/generateCRAM.nf'
include { CollectWgsMetrics } from './modules/CollectWgsMetrics.nf'
include { baserecalibrator } from './modules/baserecalibrator.nf'
include { delete_file } from './modules/delete_file.nf'
include { checkbam } from './modules/checkbam.nf'
include { checkvcf } from './modules/checkvcf.nf'
include { checkfamilyvcf } from './modules/checkfamilyvcfnew.nf'
include { applyBQSR } from './modules/applyBQSR.nf'
include { genotype } from './modules/genotype.nf'
include { leftalignandtrim } from './modules/leftalignandtrim.nf'
include { variantrecalibration } from './modules/variantrecalibration.nf'
include { compressandindex } from './modules/compressandindex.nf'
include { mergevcf } from './modules/mergevcf.nf'
include { combineGVCFs } from './modules/combineGVCFs.nf'
include { genotypeGVCFs } from './modules/genotypeGVCFs.nf'
include { test } from './testmodules/test.nf'
include { SelectVariantsdenovo } from './modules/SelectVariantsdenovo.nf'
include { SelectVariantsAR } from './modules/SelectVariantsAR.nf'
include { SelectVariantsX } from './modules/SelectVariantsX.nf'
include { annotate as annotatedenovo; annotate as annotateAR; annotate as annotateX } from './modules/annotate.nf'
include { parliament2 } from './modules/parliament2.nf'
include { vcftoolshardfilter } from './modules/vcftoolshardfilter.nf'

// script parameters


indexes_ch = Channel.fromPath(params.indexes).toList()
donebams_ch = channel.fromPath('./results/bams/*.bam*').toSortedList().flatten().collate( 2 ).map{bam,bai -> tuple(bam.simpleName,bam,bai)}.flatten().collate( 3 )
donevcfs_ch = channel.fromPath('./results/vcfs/*.vcf*').toSortedList().flatten().collate( 1 ).map{vcf -> tuple(vcf.simpleName,vcf)}.flatten().collate( 2 )
donefamilyvcfs_ch = channel.fromPath('./results/familyvcfs/*.alignedandtrimmed.vcf.gz').toSortedList().flatten().collate( 1 ).map{vcf -> tuple(vcf.simpleName,vcf)}.flatten().collate( 2 )
donefamilyvcfs_ch.view()


garbage_ch = Channel.empty()
testcollection = Channel.empty()

vcfcollection = Channel.empty()

//script

myFile = file(params.ped)
myReader = myFile.newReader()
String line
familymap = [:]
ids = []
while( line = myReader.readLine() ) {
(empty, family, id, father, mother, sex, phenotype) = (line =~ /(^.*F\d{1,2})\t(GC\d+)\t(\w+)\t(\w+)\t(\d+)\t(\d+)/)[0]
        familymap[id]=family
        ids << id
}
myReader.close()

Channel.fromList(ids).map { it -> [it, familymap[it]] }set{ idfamily_ch }
Channel.fromList(ids).map { it -> familymap[it] }.unique().collate( 1 ).dump(tag:"family").set{ family_ch }
family_ch.view()

//workflows

workflow download_fastq_to_bam_and_cram {
take: idfam
main:
importfastq(idfam, params.home,params.arch)
importfastq.out.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
fastQC(gzipped_ch)

alignment(gzipped_ch, params.genome,indexes_ch, params.home)
gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.dump(tag:"gzippedG").join(fastQC.out.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane)}.dump(tag:"fastQCG")).join(alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane)}).dump(tag:"G1").set{testgarbage_ch}

readgroups(alignment.out,params.home)
alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane,bam)}.join(readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(lane)}).dump(tag:"garbage3").set{testgarbage_ch3}

duplicates(readgroups.out,params.home)

//readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(id,bam,bai)}.dump(tag:"garbage4part1").join(duplicates.out[0].flatten().collate ( 2 ).map{id,bam -> tuple(id)}.dump(tag:"garbage4part2")).dump(tag:"garbage4").set{testgarbage_ch4}





mergebams(duplicates.out[0].groupTuple(),params.home)
generateCRAM(mergebams.out[0],params.genome,indexes_ch)
CollectWgsMetrics(mergebams.out[0],params.genome)
//duplicates.out[0].flatten().collate ( 2 ).map{id,bam -> tuple(id,bam)}.join(mergebams.out[0].flatten().collate ( 3 ).map{id,bam,bai -> tuple(id)}).join(generateCRAM.out[0].flatten().collate ( 3 ).map{id,cram,crai -> tuple(id)}).dump(tag:"garbage5").set{testgarbage_ch5}
//garbage_ch.concat(gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(lane,R1,R2)}.flatten().toList(),alignment.out.flatten().collate( 3 ).map{id,lane,bam -> tuple(lane,bam)},readgroups.out.flatten().collate( 4 ).map{id,lane,bam,bai -> tuple(lane,bam,bai)}).groupTuple().dump(tag:"garbage").set{workflow1garbage}
//duplicates.out[0].flatten().collate ( 2 ).map{lane,bam -> tuple(bam.getBaseName(2))}.join(workflow1garbage).flatten().dump(tag:"merged").set{garbagemerge}
//testcollection.concat(testgarbage_ch,testgarbage_ch3,testgarbage_ch4,testgarbage_ch5).dump(tag:"G12345").set{concatedtestcollection}

emit:
bams = mergebams.out[0]
crams = generateCRAM.out[0]
//garbage = garbagemerge
//testgarbage = concatedtestcollection
}
workflow testwf {
take: 
data
test
main:

delete_file(data)
test(test)

}

workflow createindividualvcfs {
take: bam 

main:
baserecalibrator(bam,params.genome, indexes_ch, params.genomedict, params.snps, params.snpsindex)
applyBQSR(baserecalibrator.out,params.genome,indexes_ch,params.genomedict)
baserecalibrator.out[0].flatten().collate ( 4 ).map{id,bam,bai,recaltable -> tuple(id,bam,bai)}.join(applyBQSR.out[0].flatten().collate ( 3 ).map{id,bam,bai -> tuple(id)}).set{testgarbage_ch6}

genotype(applyBQSR.out,params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)
applyBQSR.out[0].flatten().collate ( 3 ).map{id,bam,bai -> tuple(id,bam,bai)}.join(genotype.out[0].flatten().collate ( 2 ).map{id,vcf -> tuple(id)}).set{testgarbage_ch7}



emit:
individualvcf = genotype.out[0]

}

workflow createfamilyvcfs {
take: vcf

main: 
combineGVCFs((idfamily_ch.join(vcf).map{ id, family, vcf -> tuple(family, vcf)}.groupTuple().flatten().collate ( 4 ).dump(tag:"combine")),params.genome,indexes_ch,params.genomedict)
genotypeGVCFs(combineGVCFs.out[0],params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)

variantrecalibration(genotypeGVCFs.out[0],params.genome,params.genomedict,indexes_ch,params.snps, params.snpsindex,params.indels,params.indelsindex,params.mask)
//genotype.out[0].flatten().collate ( 2 ).join(variantrecalibration.out[0].flatten().collate ( 2 ).map{id,vcf -> tuple(id)}).set{testgarbage_ch8}

//compressandindex(variantrecalibration.out[0])

//variantrecalibration.out[0].flatten().collate ( 2 ).join(compressandindex.out[0].flatten().collate ( 3 ).map{id,vcfgz,vcfgztbi -> tuple(id)}).set{testgarbage_ch9}



//mergevcf(idfamily_ch.join(compressandindex.out).map{ id, family, vcf, index -> tuple(family,vcf,index)}.groupTuple())
//compressandindex.out[0].flatten().collate ( 3 ).map{id,vcfgz,vcfgztbi -> tuple(familymap[id]),id,vcfgz,vcfgztbi}.join(mergevcf.out[0].flatten().collate ( 3 ).map{family,vcfgz,vcfgztbi -> tuple(family)}).set{testgarbage_ch10}
vcftoolshardfilter(variantrecalibration.out[0])
//vcftoolshardfilter(compressandindex.out[0])
leftalignandtrim(vcftoolshardfilter.out[0],params.genome,indexes_ch,params.genomedict)
//mergevcf.out[0].flatten().collate ( 3 ).join(leftalignandtrim.out[0].flatten().collate ( 3 ).map{family,vcfgz,vcfgztbi -> tuple(family)}).set{testgarbage_ch9}

//vcfcollection.concat(testgarbage_ch9,testgarbage_ch10).set{concatedvcfcollection}


emit:
triovcf = leftalignandtrim.out
//vcfgarbage = concatedvcfcollection
}

workflow trioVCFanalysis {
take:vcf

main:
SelectVariantsdenovo(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsX(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(SelectVariantsAR.out[0],params.programs,params.humandb,params.annovardbs)
annotateX(SelectVariantsX.out[0],params.programs,params.humandb,params.annovardbs)
}



workflow { 
main:

checkbam(idfamily_ch)
checkbam.out.bamcheck_ch.dump(tag:"done").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{bamdone_ch}
bamdone_ch.toSortedList().flatten().collate(1).combine(donebams_ch, by:0).map{id,bam,bai -> tuple(id,bam,bai)}.set{bamalldone_ch}
download_fastq_to_bam_and_cram(checkbam.out.bamcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"todo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id,family)})
download_fastq_to_bam_and_cram.out.bams.concat(bamalldone_ch).set{bammixed}

parliament2(bammixed,params.genome,indexes_ch)

checkvcf(idfamily_ch)
checkvcf.out.vcfcheck_ch.dump(tag:"vcfdone").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{vcfdone_ch}
vcfdone_ch.toSortedList().flatten().collate(1).combine(donevcfs_ch, by:0).map{id,vcf -> tuple(id,vcf)}.set{vcfalldone_ch}
createindividualvcfs(checkvcf.out.vcfcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"vcftodo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id)}.join(bammixed))
createindividualvcfs.out.individualvcf.concat(vcfalldone_ch).set{vcfmixed}

checkfamilyvcf(family_ch.flatten()).view()
checkfamilyvcf.out.familyvcfcheck_ch.dump(tag:"done").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 2 ).map{family,status -> family}.set{familyvcfdone_ch}
familyvcfdone_ch.toSortedList().flatten().collate(1).combine(donefamilyvcfs_ch, by:0).map{family,vcf -> tuple(family,vcf)}.set{familyvcfalldone_ch}
familyvcfalldone_ch.view()
createfamilyvcfs(checkfamilyvcf.out.familyvcfcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"todo").groupTuple().flatten().collate( 2 ).map{family,status -> tuple(family)}.join(vcfmixed))
createfamilyvcfs.out.triovcf.concat(familyvcfalldone_ch).set{familyvcfmixed} 
familyvcfmixed.view()

//right-one testwf(download_fastq_to_bam_and_cram.out.testgarbage.flatten(),download_fastq_to_bam_and_cram.out.testgarbage.flatten())
//testwf(download_fastq_to_bam_and_cram.out.testgarbage.flatten(),createfamilyvcfs.out.vcfgarbage.flatten())


trioVCFanalysis(familyvcfmixed)


}