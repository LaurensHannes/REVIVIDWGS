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

include { importfastq;importexomefastq} from './modules/importfastq.nf'
include { deeptrio; glnexusdt; glnexusprocessing} from './modules/deeptrio.nf'
include { fastQC } from './modules/fastQC.nf'
include { alignment } from './modules/alignment.nf'
include { readgroups } from './modules/readgroups.nf'
include { duplicates } from './modules/duplicates.nf'
include { stmergebams;mergebams } from './modules/mergebams.nf'
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
include { variantrecalibration;variantcohortrecalibration;variantchrrecalibration;variantcohortchrrecalibration } from './modules/variantrecalibration.nf'
include { compressandindex } from './modules/compressandindex.nf'
include { concatvcf; mergevcf;intersectvcf; normalizeindels as normalizeindelsdeepvariant;normalizeindels as normalizeindelsgatk } from './modules/bcftools.nf'
include { combineindividualGVCFs;combinechrGVCFs } from './modules/combineindividualGVCFs.nf'
include { combineGVCFs;combinecohortGVCFs;combinechrVCFs } from './modules/combineGVCFs.nf'
include { genotypeGVCFs;genotypechrGVCFs } from './modules/genotypeGVCFs.nf'
include { test } from './testmodules/test.nf'
include { SelectVariantsdenovo } from './modules/SelectVariantsdenovo.nf'
include { SelectVariantsAR } from './modules/SelectVariantsAR.nf'
include { SelectVariantsX } from './modules/SelectVariantsX.nf'
include { annovar as annotatedenovo; annovar as annotateAR; annovar as annotateX } from './modules/annotate.nf'
include { VEP as VEPdenovo; VEP as VEPAR; VEP as VEPX } from './modules/annotate.nf'
include { parliament2;indelible } from './modules/CNV.nf'
include { peddy } from './modules/peddy.nf'
include { AnnotSV } from './modules/AnnotSV.nf'
include { mergeCNV } from './modules/mergeCNV.nf'
include { vcftoolshardfilter } from './modules/vcftoolshardfilter.nf'
include { splitbamlanes } from './modules/splitbamlanes.nf'
include { splitbamindividuals } from './modules/splitbamindividuals.nf'
include { createfilterbedfileCNV } from './modules/createfilterbedfileCNV.nf'

// script parameters

		chromosomes_ch = Channel.of('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')


switch (params.flow) {
	case 'wgs':
		chromosomes_ch = Channel.of('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM')
		break
	case 'wes':
		chromosomes_ch = Channel.of('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
		break
	case 'sim':
		chromosomes_ch = Channel.of('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')
		break


} 




// bams_ch = channel.fromPath('')

ped_ch = Channel.fromPath(params.ped).map{ ped -> ped.getSimpleName()}.flatten()
indexes_ch = Channel.fromPath(params.indexes).toList()
donebams_ch = channel.fromPath('./results/bams/*.bam*').toSortedList().flatten().collate( 2 ).map{bam,bai -> tuple(bam.simpleName,bam,bai)}.flatten().collate( 3 )
donevcfs_ch = channel.fromPath('./results/vcfs/*.vcf*').toSortedList().flatten().collate( 1 ).map{vcf -> tuple(vcf.simpleName,vcf)}.flatten().collate( 2 )
donefamilyvcfs_ch = channel.fromPath('./results/familyvcfs/*.vcf.gz').toSortedList().flatten().collate( 1 ).map{vcf -> tuple(vcf.simpleName,vcf)}.flatten().collate( 2 )


garbage_ch = Channel.empty()
testcollection = Channel.empty()

vcfcollection = Channel.empty() 

//script

myFile = file(params.ped)
myReader = myFile.newReader()
String line
familymap = [:]
trios = []
families = []
individuals = []

while( line = myReader.readLine() ) {
(empty, family, id, father, mother, sex, phenotype) = (line =~ /(^.*F\d{1,5})\t(\w+)\t(\w+)\t(\w+)\t(\d+)\t(\d+)/)[0]
        familymap[id]=family
		individuals << id
        if(father!="0") {
		trios << tuple(id,father,mother)
		families << tuple(family,id,father,mother)
		}
}
myReader.close()
//shortped_ch = Channel.fromList(trios)
//shortped_ch.flatten().set{ ids }
ids = Channel.fromList(individuals)
familytrio_ch = Channel.fromList(families)
triofamilywithchr_ch = chromosomes_ch.combine(familytrio_ch).map{chr,fam,index,father,mother -> tuple(fam,chr,index,father,mother)}

ids.map { it -> [it, familymap[it]] }.set{ idfamily_ch }
ids.map { it -> familymap[it] }.unique().collate( 1 ).dump(tag:"family").set{ family_ch }
family_ch.view()

//workflows

workflow download_fastq_to_bam_and_cram {
take: idfam
main:

importfastq(idfam, params.home,params.arch,params.download,params.bucket)
importfastq.out.flatten().filter(~/.*R[12]\.*.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}

fastQC(gzipped_ch,params.arch)
alignment(gzipped_ch, params.genome,indexes_ch, params.home,params.arch)
mergebams(alignment.out[0].map{id,lane,bam,bai -> tuple(id,bam)}.groupTuple(),params.home)
generateCRAM(mergebams.out[0],params.genome,indexes_ch,params.arch)
CollectWgsMetrics(mergebams.out[0],params.genome,params.arch)

emit:
mergedbams = mergebams.out[0]
crams = generateCRAM.out[0]

}



workflow createindividualbams {
take: bam

main: 

splitbamindividuals(bam,chromosomes_ch,params.genome,indexes_ch)

emit:
individualsplitbam = splitbamindividuals.out[0]

}

workflow deepvariant {
take: bam 


main:
bam.map{id,chr,bam,bai -> tuple(familymap[id],chr,bam,bai)}.groupTuple(by:[0,1]).join(triofamilywithchr_ch,by:[0,1]).flatten().collate(11).set{ deeptrioinput_ch }

deeptrio(deeptrioinput_ch,params.genome,indexes_ch)
deeptrio.out[0].concat( deeptrio.out[1], deeptrio.out[2]).groupTuple(by:[0,1],sort:true).flatten().collate( 52 ).view()
concatvcf(deeptrio.out[0].concat( deeptrio.out[1], deeptrio.out[2]).groupTuple(by:[0,1],sort:true).flatten().collate( 52 ))
glnexusdt(concatvcf.out[0].map{ family, id, vcf ,vcftbi -> tuple(family,vcf,vcftbi)}.groupTuple().flatten().collate( 7 ).join(familytrio_ch).flatten().collate( 10 ))
glnexusprocessing(glnexusdt.out[0])
normalizeindelsdeepvariant(glnexusprocessing.out[0],params.genome)

emit:
deepvariantvcf = normalizeindelsdeepvariant.out[0]

}



workflow createindividualvcfs {
take: bamperchr



main:

baserecalibrator(bamperchr,params.genome, indexes_ch, params.genomedict, params.snps, params.snpsindex,params.broadinterval)
applyBQSR(baserecalibrator.out,params.genome,indexes_ch,params.genomedict,params.broadinterval)
genotype(applyBQSR.out,params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)
if( params.cohort ) {
	genotype.out[0].flatMap{ id, vcf, idx -> tuple(vcf, idx) }.toList().set{cohortgenotypeoutput}
	combinechrGVCFs(cohortgenotypeoutput,chromosomes_ch,params.genome,indexes_ch,params.genomedict,params.broadinterval,ped_ch)
	combinechrGVCFs.out[0].set{individualgvcfsoutput_ch}
}
else {
genotype.out[0].groupTuple(size: 24).flatten().collate ( 25 ).view().set{collatedgenotypes_ch}
combineindividualGVCFs(collatedgenotypes_ch,params.genome,indexes_ch,params.genomedict)
combineindividualGVCFs.out[0].set{individualgvcfsoutput_ch}
}

emit:
individualvcf = individualgvcfsoutput_ch

}

workflow createfamilyvcfs {
take: vcf

main: 
if( params.cohort) {
genotypechrGVCFs(vcf.flatten().toList(),chromosomes_ch,params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)
variantcohortchrrecalibration(genotypechrGVCFs.out[0].toSortedList(),genotypechrGVCFs.out[1].toSortedList(),params.genome,params.genomedict,indexes_ch,params.snps,params.snpsindex,params.snpstruth,params.snpstruthindex,params.indels,params.indelsindex,chromosomes_ch)
combinechrVCFs(variantcohortchrrecalibration.out[0].toSortedList(),params.genome,indexes_ch,params.genomedict)
combinechrVCFs.out[0].set{vrecal_ch}

}
else {
combineGVCFs(vcf,params.genome,indexes_ch,params.genomedict)
genotypeGVCFs(combineGVCFs.out[0],params.genome,indexes_ch,params.broadinterval,params.genomedict,params.mask)
variantrecalibration(genotypeGVCFs.out[0],params.genome,params.genomedict,indexes_ch,params.snps,params.snpsindex,params.indels,params.indelsindex)
variantrecalibration.out[0].set{vrecal_ch}
}
vcftoolshardfilter(vrecal_ch)
normalizeindelsgatk(vcftoolshardfilter.out[0],params.genome)

emit:
triovcf = normalizeindelsgatk.out
}

workflow CNVanalysis {
take:bam

main:
if ( params.exome == 'false' ) {
parliament2(bam,params.genome,indexes_ch)
createfilterbedfileCNV(bam)
idfamily_ch.join(parliament2.out[0].join(createfilterbedfileCNV.out[0])).map{ id, family, vcf ,lowmq -> tuple(family,vcf,lowmq)}.groupTuple().flatten().collate( 7 ).view()
mergeCNV(idfamily_ch.join(parliament2.out[0].join(createfilterbedfileCNV.out[0])).map{ id, family, vcf ,lowmq -> tuple(family,vcf,lowmq)}.groupTuple().flatten().collate( 7 ))
mergeCNV.out[0].set{CNV_ch}
}
else {
indelible(familytrio_ch,bam.map{ id, bam, bai -> tuple(bam,bai)}.flatten().toList())
indelible.out[0].set{CNV_ch}
}
emit:
CNV_ch
}

workflow consensus {

take:
vcf1
vcf2

main:
intersectvcf(vcf1,vcf2)
peddy(vcf1.concat(vcf2),params.ped)
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
SelectVariantsAR.out[0].map{fam,analysis,mode,vcfgz,tbi -> tuple(fam,analysis,vcfgz,tbi)}.groupTuple(by:1).flatten().unique().view()
mergevcf(SelectVariantsAR.out[0].map{fam,analysis,mode,vcfgz,tbi -> tuple(fam,analysis,vcfgz,tbi)}.groupTuple(by:1).flatten().unique().collate( 8 ))
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(mergevcf.out[0],params.programs,params.humandb,params.annovardbs)
annotateX(SelectVariantsX.out[0],params.programs,params.humandb,params.annovardbs)
AnnotSV(mergevcf.out[0].join(cnvvcf).groupTuple().flatten().collate( 6 ).map{fam,analysis,mode,snvvcf,snvvcftbi,cnvvcf -> tuple(fam,snvvcf,mode,cnvvcf)})
if ( params.VEPCALL == 'true' ) {
VEPdenovo(SelectVariantsdenovo.out[0],params.VEP)
VEPAR(mergevcf.out[0],params.VEP)
VEPX(SelectVariantsX.out[0],params.VEP)
}

}


workflow { 
main:

checkbam(idfamily_ch)
checkbam.out.bamcheck_ch.dump(tag:"done").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{bamdone_ch}
bamdone_ch.toSortedList().flatten().collate(1).combine(donebams_ch, by:0).map{id,bam,bai -> tuple(id,bam,bai)}.set{bamalldone_ch}
download_fastq_to_bam_and_cram(checkbam.out.bamcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"todo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id,family)})
download_fastq_to_bam_and_cram.out.mergedbams.concat(bamalldone_ch).set{bammixed}



checkvcf(idfamily_ch)
checkvcf.out.vcfcheck_ch.dump(tag:"vcfdone").filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.set{vcfdone_ch}
vcfdone_ch.toSortedList().flatten().collate(1).combine(donevcfs_ch, by:0).map{id,vcf -> tuple(id,vcf)}.set{vcfalldone_ch}
createindividualbams(checkvcf.out.vcfcheck_ch.filter( ~/.*todo.*/ ).dump(tag:"vcftodo").groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id)}.join(bammixed))
deepvariant(createindividualbams.out)
CNVanalysis(bammixed)
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









// CUSTOM WORKFLOWS for RANDOM ENTRY 

workflow consensusentry {

dv_ch = Channel.fromPath(params.tbi).filter( ~/.*deeptrio.*/ ).toSortedList().flatten().collate( 2 )
gatk_ch = Channel.fromPath(params.tbi).filter( ~/.*GATK.*/ ).toSortedList().flatten().collate( 2 )
dv_complete = dv_ch.map{vcf,tbi -> tuple(tbi.simpleName,"deepvariant",vcf,tbi)}
gatk_complete = gatk_ch.map{vcf,tbi -> tuple(tbi.simpleName,"GATK",vcf,tbi)}

main:
//peddy(dv_complete.concat(gatk_complete),params.ped)
intersectvcf(dv_complete,gatk_complete)
intersectvcf.out[0].concat(intersectvcf.out[1],intersectvcf.out[2]).set{isec_ch}
SelectVariantsdenovo(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsX(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR.out[0].map{fam,analysis,mode,vcfgz,tbi -> tuple(fam,analysis,vcfgz,tbi)}.groupTuple(by:1).flatten().unique().collate( 8 ).view()
mergevcf(SelectVariantsAR.out[0].map{fam,analysis,mode,vcfgz,tbi -> tuple(fam,analysis,vcfgz,tbi)}.groupTuple(by:1).flatten().unique().collate( 8 ))
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(mergevcf.out[0],params.programs,params.humandb,params.annovardbs)
annotateX(SelectVariantsX.out[0],params.programs,params.humandb,params.annovardbs)
VEPdenovo(SelectVariantsdenovo.out[0],params.VEP)
VEPAR(mergevcf.out[0],params.VEP)
VEPX(SelectVariantsX.out[0],params.VEP)
}

workflow download {


main:
importfastq(idfamily_ch, params.home,params.arch,params.download,params.bucket)

}


workflow createbams {


main:

importfastq(idfamily_ch, params.home,params.arch,params.download,params.bucket,params.exome,params.lsbucket)
importfastq.out.flatten().filter(~/.*R\d+.*.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
//importfastq.out.flatten().filter(~/.*R\d+.*.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").flatten().collate( 3 ).splitFastq(by: 10_000_000, pe: true, file:true).map{lane,R1,R2 -> tuple(R1.simpleName,R1.getBaseName(1),R1,R2)}.view().set{gzipped_ch}

fastQC(gzipped_ch,params.arch)
alignment(gzipped_ch, params.genome,indexes_ch, params.home,params.arch)
stmergebams(alignment.out[0].map{id,lane,bam,bai -> tuple(id,bam)}.groupTuple(),params.home,params.arch)

stmergebams.out[0].map{file -> tuple(file.simpleName,file)}.groupTuple().view().set{mergedbamstemp1_ch}
stmergebams.out[1].map{file -> tuple(file.simpleName,file)}.groupTuple().view().set{mergedbamstemp2_ch}
mergedbamstemp1_ch.join(mergedbamstemp2_ch).flatten().collate( 3 ).view().set{mergedbamstemp_ch}
CollectWgsMetrics(mergedbamstemp_ch,params.genome,params.arch)
generateCRAM(mergedbamstemp_ch,params.genome,indexes_ch,params.arch)
}

workflow annotatewithcnventry {


dv_ch = Channel.fromPath(params.tbi).filter( ~/.*deeptrio.*/ ).toSortedList().flatten().collate( 2 )
gatk_ch = Channel.fromPath(params.tbi).filter( ~/.*GATK.*/ ).toSortedList().flatten().collate( 2 )
cnv_ch = Channel.fromPath(params.cnv).map{cnv -> tuple(cnv.getSimpleName(),cnv)}
dv_complete = dv_ch.map{vcf,tbi -> tuple(tbi.simpleName,"deepvariant",vcf,tbi)}
gatk_complete = gatk_ch.map{vcf,tbi -> tuple(tbi.simpleName,"GATK",vcf,tbi)}

main:	
peddy(dv_complete.concat(gatk_complete),params.ped)
intersectvcf(dv_complete,gatk_complete)
intersectvcf.out[0].concat(intersectvcf.out[1],intersectvcf.out[2]).set{isec_ch}
SelectVariantsdenovo(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsX(isec_ch,params.genome,params.genomedict,indexes_ch,params.ped,params.mask)
SelectVariantsAR.out[0].map{fam,analysis,mode,vcfgz,tbi -> tuple(fam,analysis,vcfgz,tbi)}.groupTuple(by:1).flatten().unique().collate( 8 )
mergevcf(SelectVariantsAR.out[0].map{fam,analysis,mode,vcfgz,tbi -> tuple(fam,analysis,vcfgz,tbi)}.groupTuple(by:1).flatten().unique().collate( 8 ))
annotatedenovo(SelectVariantsdenovo.out[0],params.programs,params.humandb,params.annovardbs)
annotateAR(mergevcf.out[0],params.programs,params.humandb,params.annovardbs)
annotateX(SelectVariantsX.out[0],params.programs,params.humandb,params.annovardbs)
AnnotSV(mergevcf.out[0].join(cnv_ch).groupTuple().flatten().collate( 6 ).map{fam,analysis,mode,snvvcf,snvvcftbi,cnvvcf -> tuple(fam,snvvcf,mode,cnvvcf)})
VEPdenovo(SelectVariantsdenovo.out[0],params.VEP)
VEPAR(mergevcf.out[0],params.VEP)
VEPX(SelectVariantsX.out[0],params.VEP)

}

workflow createbamsandcallvariants {


main:


importfastq(idfamily_ch, params.home,params.arch,params.download,params.bucket,params.exome,params.lsbucket)
importfastq.out.flatten().filter(~/.*R\d+.*.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").view().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
//importfastq.out.flatten().filter(~/.*R\d+.*.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().dump(tag:"test").flatten().collate( 3 ).splitFastq(by: 10_000_000, pe: true, file:true).map{lane,R1,R2 -> tuple(R1.simpleName,R1.getBaseName(1),R1,R2)}.view().set{gzipped_ch}

fastQC(gzipped_ch,params.arch)
alignment(gzipped_ch, params.genome,indexes_ch, params.home,params.arch)
stmergebams(alignment.out[0].map{id,lane,bam,bai -> tuple(id,bam)}.groupTuple(),params.home,params.arch)

stmergebams.out[0].map{file -> tuple(file.simpleName,file)}.groupTuple().view().set{mergedbamstemp1_ch}
stmergebams.out[1].map{file -> tuple(file.simpleName,file)}.groupTuple().view().set{mergedbamstemp2_ch}
mergedbamstemp1_ch.join(mergedbamstemp2_ch).flatten().collate( 3 ).view().set{mergedbamstemp_ch}
CollectWgsMetrics(mergedbamstemp_ch,params.genome,params.arch)
generateCRAM(mergedbamstemp_ch,params.genome,indexes_ch,params.arch)

createindividualbams(mergedbamstemp_ch)
if ( params.CNV == 'true' ) {
CNVanalysis(bammixed)
}
if (params.caller == 'both' ) {
deepvariant(createindividualbams.out)
createindividualvcfs(createindividualbams.out)
createindividualvcfs.out.map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}
createfamilyvcfs(vcfmixed)
}
else if (params.caller == 'deepvariant' ) {
deepvariant(createindividualbams.out)
}
else if (params.caller == 'gatk' ) {
createindividualvcfs(createindividualbams.out)

if (params.cohort) {
	createindividualvcfs.out.set{vcfmixed}
}
		else {
createindividualvcfs.out.map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}
		}
createfamilyvcfs(vcfmixed)
}
}

workflow callgvariantsexome {

take:
bam
main:


createindividualbams(bam)

if (params.caller == 'both' ) {
deepvariant(createindividualbams.out)
createindividualvcfs(createindividualbams.out)
createindividualvcfs.out.map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}
//createfamilyvcfs(vcfmixed)
}
else if (params.caller == 'deepvariant' ) {
deepvariant(createindividualbams.out)
}
else if (params.caller == 'gatk' ) {
createindividualvcfs(createindividualbams.out)
//createindividualvcfs.out.map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}




}
emit:
createindividualvcfs.out[0]

}



workflow cohortexomecalling {

main:

if (params.individualchrgvcfs) {

createfamilyvcfs(Channel.fromPath(params.individualchrgvcfs))
}

else {
channel.fromPath(params.exomemapped).map{file -> tuple(file.simpleName,file)}.groupTuple().set{mergedbamstemp1_ch}
ids.join(mergedbamstemp1_ch).flatten().collate( 2 ).groupTuple().view().set{mergedbamstemp2_ch}
channel.fromPath(params.exomeindexes).map{file -> tuple(file.simpleName,file)}.groupTuple().set{mergedbamstemp3_ch}
mergedbamstemp2_ch.join(mergedbamstemp3_ch).flatten().collate( 3 ).view().set{mergedbamstemp_ch}
callgvariantsexome(mergedbamstemp_ch)
if ( params.CNV == 'true' ) {
CNVanalysis(mergedbamstemp_ch)
}
createfamilyvcfs(callgvariantsexome.out[0])

}

}


workflow callvariantsfrombams {


main:



createindividualbams(mergedbamstemp_ch)
if ( params.CNV == 'true' ) {
CNVanalysis(bammixed)
}
if (params.caller == 'both' ) {
deepvariant(createindividualbams.out)
createindividualvcfs(createindividualbams.out)
createindividualvcfs.out.map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}
createfamilyvcfs(vcfmixed)
}
else if (params.caller == 'deepvariant' ) {
deepvariant(createindividualbams.out)
}
else if (params.caller == 'gatk' ) {
createindividualvcfs(createindividualbams.out)
createindividualvcfs.out.map{id,vcf -> tuple(familymap[id], vcf)}.groupTuple().flatten().collate( 4 ).set{vcfmixed}
createfamilyvcfs(vcfmixed)
}
}