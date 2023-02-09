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
include { annovar as annotatedenovo; annovar as annotateAR; annovar as annotateX } from './modules/annotate.nf'
include { VEP as VEPdenovo; VEP as VEPAR; VEP as VEPX } from './modules/annotate.nf'
include { parliament2 } from './modules/parliament2.nf'
include { peddy } from './modules/peddy.nf'
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
(empty1, family1, id1, father1, mother1, sex1, phenotype1) = (line =~ /(^.*F\d{1,2})\t(GC\d+)\t(GC\d+)\t(GC\d+)\t(\d+)\t(\d+)/)[0]
		ids1 << id1
		fathers1 << father1
		mothers1 << mother1
}
myReader.close()
a =  Channel.fromList(ids)
b = Channel.fromList(fathers)
c = Channel.fromList(mothers)

d =  Channel.fromList(ids1)
e = Channel.fromList(fathers1)
f = Channel.fromList(mothers1)


a.first().concat(b.first(),c.first()).flatten().collate( 3 ).set{ shortped_ch }

d.merge(e).merge(f).view()

