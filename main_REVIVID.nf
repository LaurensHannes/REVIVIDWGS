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
include { deeptrio; glnexusdt; glnexusprocessing} from './modules/deeptrio.nf'
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
include { leftalignandtrim as  leftalignandtrimgatk;  leftalignandtrim as  leftalignandtrimdeepvariant} from './modules/leftalignandtrim.nf'
include { variantrecalibration } from './modules/variantrecalibration.nf'
include { compressandindex } from './modules/compressandindex.nf'
include { concatvcf; mergevcf;intersectvcf; normalizeindels as normalizeindelsdeepvariant;normalizeindels as normalizeindelsgatk } from './modules/bcftools.nf'
include { combineindividualGVCFs } from './modules/combineindividualGVCFs.nf'
include { combineGVCFs } from './modules/combineGVCFs.nf'
include { genotypeGVCFs } from './modules/genotypeGVCFs.nf'
include { test } from './testmodules/test.nf'
include { SelectVariantsdenovo } from './modules/SelectVariantsdenovo.nf'
include { SelectVariantsAR } from './modules/SelectVariantsAR.nf'
include { SelectVariantsX } from './modules/SelectVariantsX.nf'
include { annotate as annotatedenovo; annotate as annotateAR; annotate as annotateX } from './modules/annotate.nf'
include { parliament2 } from './modules/parliament2.nf'
include { AnnotSV } from './modules/AnnotSV.nf'
include { mergeCNV } from './modules/mergeCNV.nf'
include { vcftoolshardfilter } from './modules/vcftoolshardfilter.nf'
include { splitbamlanes } from './modules/splitbamlanes.nf'
include { splitbamindividuals } from './modules/splitbamindividuals.nf'
include { createfilterbedfileCNV } from './modules/createfilterbedfileCNV.nf'

// script parameters

chromosomes_ch = Channel.of('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM')



indexes_ch = Channel.fromPath(params.indexes).toList()
donebams_ch = channel.fromPath('./results/bams/*.bam*').toSortedList().flatten().collate( 2 ).map{bam,bai -> tuple(bam.simpleName,bam,bai)}.flatten().collate( 3 )
donevcfs_ch = channel.fromPath('./results/vcfs/*.vcf*').toSortedList().flatten().collate( 1 ).map{vcf -> tuple(vcf.simpleName,vcf)}.flatten().collate( 2 )
donefamilyvcfs_ch = channel.fromPath('./results/familyvcfs/*.vcf.gz').toSortedList().flatten().collate( 1 ).map{vcf -> tuple(vcf.simpleName,vcf)}.flatten().collate( 2 )
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
fathers = []
mothers = []
while( line = myReader.readLine() ) {
(empty, family, id, father, mother, sex, phenotype) = (line =~ /(^.*F\d{1,2})\t(GC\d+)\t(\w+)\t(\w+)\t(\d+)\t(\d+)/)[0]
        familymap[id]=family
        ids << id
		fathers << father
		mothers << mother
}
myReader.close()
a =  Channel.fromList(ids)
b = Channel.fromList(fathers)
c = Channel.fromList(mothers)

a.first().concat(b.first(),c.first()).flatten().collate( 3 ).set{ shortped_ch }

Channel.fromList(ids).map { it -> [it, familymap[it]] }.set{ idfamily_ch }
Channel.fromList(ids).map { it -> familymap[it] }.unique().collate( 1 ).dump(tag:"family").set{ family_ch }
family_ch.view()

//workflows

workflow download_fastq_to_bam_and_cram {
take: idfam
main:

importfastq(idfam, params.home,params.arch,params.download)
importfastq.out.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}

fastQC(gzipped_ch)
alignment(gzipped_ch, params.genome,indexes_ch, params.home)
mergebams(alignment.out[0].map{id,lane,bam,bai -> tuple(id,bam)}.groupTuple(),params.home)
generateCRAM(mergebams.out[0],params.genome,indexes_ch)
CollectWgsMetrics(mergebams.out[0],params.genome)

emit:
mergedbams = mergebams.out[0]
crams = generateCRAM.out[0]

}



workflow createindividualbams {
take: bam

main: 

splitbamindividuals(bam,chromosomes_ch)

emit:
individualsplitbam = splitbamindividuals.out[0]

}

workflow deepvariant {
take: bam 


main:

deeptrio(bam.map{id,chr,bam,bai -> tuple(chr,tuple(bam,bai))}.groupTuple().flatten().collate( 7 ).combine(shortped_ch).flatten().collate( 10 ),params.genome,indexes_ch)
concatvcf(deeptrio.out[0].concat( deeptrio.out[1], deeptrio.out[2]).groupTuple(sort:true).flatten().collate( 51 ))
glnexusdt(idfamily_ch.join(concatvcf.out[0]).map{ id, family, vcf ,vcftbi -> tuple(family,vcf,vcftbi)}.groupTuple().flatten().collate( 7 ).combine(shortped_ch).flatten().collate( 10 ))
glnexusprocessing(glnexusdt.out[0])
normalizeindelsdeepvariant(glnexusprocessing.out[0],params.genome)

emit:
deepvariantvcf = normalizeindelsdeepvariant.out[0]

}



workflow createindividualvcfs {
take: bamperchr



main:



baserecalibrator(bamperchr,params.genome, indexes_ch, params.genomedict, params.snps, params.snpsindex)
applyBQSR(baserecalibrator.out,params.genome,indexes_ch,params.genomedict)
genotype(applyBQSR.out,params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)
combineindividualGVCFs(genotype.out[0].groupTuple(size: 25).flatten().collate ( 26 ),params.genome,indexes_ch,params.genomedict)

emit:
individualvcf = combineindividualGVCFs.out[0]

}

workflow createfamilyvcfs {
take: vcf

main: 
combineGVCFs(vcf,params.genome,indexes_ch,params.genomedict)
genotypeGVCFs(combineGVCFs.out[0],params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)

variantrecalibration(genotypeGVCFs.out[0],params.genome,params.genomedict,indexes_ch,params.snps, params.snpsindex,params.indels,params.indelsindex,params.mask)
vcftoolshardfilter(variantrecalibration.out[0])
normalizeindelsgatk(vcftoolshardfilter.out[0],params.genome)

emit:
triovcf = normalizeindelsgatk.out
}

workflow CNVanalysis {
take:bam

main:
parliament2(bam,params.genome,indexes_ch)
createfilterbedfileCNV(bam)
idfamily_ch.join(parliament2.out[0].join(createfilterbedfileCNV.out[0])).map{ id, family, vcf ,lowmq -> tuple(family,vcf,lowmq)}.groupTuple().flatten().collate( 7 ).view()
mergeCNV(idfamily_ch.join(parliament2.out[0].join(createfilterbedfileCNV.out[0])).map{ id, family, vcf ,lowmq -> tuple(family,vcf,lowmq)}.groupTuple().flatten().collate( 7 ))

emit:
mergeCNV.out[0]
}

workflow consensus {

take:
vcf1
vcf2

main:
intersectvcf(vcf1,vcf2)

emit:
intersectvcf.out[0].concat(intersectvcf.out[1],intersectvcf.out[2])
}
workflow triovcfanalysis {
take:
vcf
cnvvcf

main:
SelectVariantsdenovo(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsX(vcf,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(SelectVariantsAR.out[0],params.programs,params.humandb,params.annovardbs)
annotateX(SelectVariantsX.out[0],params.programs,params.humandb,params.annovardbs)
AnnotSV(vcf.join(cnvvcf).groupTuple().flatten().collate( 4 ))
}


workflow { 
main:

checkbam(idfamily_ch)
checkbam.out.bamcheck_ch.dump(tag:"done").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{bamdone_ch}
bamdone_ch.toSortedList().flatten().collate(1).combine(donebams_ch, by:0).map{id,bam,bai -> tuple(id,bam,bai)}.set{bamalldone_ch}
download_fastq_to_bam_and_cram(checkbam.out.bamcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"todo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id,family)})
download_fastq_to_bam_and_cram.out.mergedbams.concat(bamalldone_ch).set{bammixed}

CNVanalysis(bammixed)

checkvcf(idfamily_ch)
checkvcf.out.vcfcheck_ch.dump(tag:"vcfdone").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{vcfdone_ch}
vcfdone_ch.toSortedList().flatten().collate(1).combine(donevcfs_ch, by:0).map{id,vcf -> tuple(id,vcf)}.set{vcfalldone_ch}
createindividualbams(checkvcf.out.vcfcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"vcftodo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id)}.join(bammixed))
deepvariant(createindividualbams.out)
createindividualvcfs(createindividualbams.out)
createindividualvcfs.out.individualvcf.concat(vcfalldone_ch).map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}

checkfamilyvcf(family_ch.flatten())
checkfamilyvcf.out.familyvcfcheck_ch.dump(tag:"done").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 2 ).map{family,status -> family}.set{familyvcfdone_ch}
familyvcfdone_ch.toSortedList().flatten().collate(1).combine(donefamilyvcfs_ch, by:0).map{family,vcf -> tuple(family,vcf)}.set{familyvcfalldone_ch}
createfamilyvcfs(checkfamilyvcf.out.familyvcfcheck_ch.filter( ~/.*todo.*/ ).groupTuple().flatten().collate( 2 ).map{family,status -> tuple(family)}.join(vcfmixed))
createfamilyvcfs.out.triovcf.concat(familyvcfalldone_ch).set{familyvcfmixed} 


//right-one testwf(download_fastq_to_bam_and_cram.out.testgarbage.flatten(),download_fastq_to_bam_and_cram.out.testgarbage.flatten())
//testwf(download_fastq_to_bam_and_cram.out.testgarbage.flatten(),createfamilyvcfs.out.vcfgarbage.flatten())

consensus(deepvariant.out.deepvariantvcf,familyvcfmixed)

triovcfanalysis(consensus.out,CNVanalysis.out[0])


}

workflow consensusentry {

tbi_ch = Channel.fromPath(params.tbi).map{tbi -> tuple(tbi.simpleName,tbi.getBaseName(),tbi)}
short_ch = Channel.fromPath(params.tbi).map{tbi -> tuple(tbi.simpleName)}
callers = Channel.from('deepvariant','GATK')
combined = short_ch.combine(callers)
tbi_ch.combine(combined, by:0).map{fam,vcf,vcftbi,caller -> tuple(fam,caller,vcf,vcftbi)}.view().set{finishedconsensusentry_ch}

main:
intersectvcf(finishedconsensusentry_ch.first(),finishedconsensusentry_ch.last())
intersectvcf.out[0].concat(intersectvcf.out[1],intersectvcf.out[2]).set{isec_ch}
SelectVariantsdenovo(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsX(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(SelectVariantsAR.out[0],params.programs,params.humandb,params.annovardbs)
annotateX(SelectVariantsX.out[0],params.programs,params.humandb,params.annovardbs)

}