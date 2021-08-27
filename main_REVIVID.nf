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
include { pear } from './pear.nf'
include { alignment } from './alignment.nf'
include { readgroups } from './readgroups.nf'
include { duplicates } from './duplicates.nf'
include { mergebams } from './mergebams.nf'
include { generateCRAM } from './generateCRAM.nf'
include { baserecalibrator } from './baserecalibrator.nf'
include { delete_file } from './delete_file.nf'
include { checkbam } from './checkbam.nf'


// script parameters



indexes_ch = Channel.fromPath(params.indexes).toList()



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
//gzipped_ch.view()

pear(gzipped_ch, params.home)
if (pear.out)
	delete_file(gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(R1,R2)}.flatten())
gzipped_ch.flatten().collate( 4 ).map{id,lane,R1,R2 -> tuple(R1,R2)}.flatten().view()
alignment(pear.out, params.genome,indexes_ch, params.home)
readgroups(alignment.out,params.home)
duplicates(readgroups.out,params.home)

mergebams(duplicates.out[0].groupTuple(),params.home)
generateCRAM(mergebams.out,params.genome,indexes_ch)
emit:
mergebams.out
generateCRAM.out[0]


}

workflow createvcfs {
take: bam

main:
baserecalibrator(bam,params.genome, indexes_ch, params.genomedict, params.snps, params.snpsindex)
emit:
baserecalibrator.out
}

workflow { 
Channel.empty().set{ createvcfsinput_ch }
checkbam(idfamily_ch)
checkbam.out.test_ch.filter( ~/.*done.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> id}.view()
download_fastq_to_bam_and_cram(checkbam.out.test_ch.filter( ~/.*todo.*/ ).groupTuple().flatten().collate( 3 ).map{id,family,status -> tuple(id,family)})
Channel.empty().set{ mixed }
createvcfsinput_ch.view()
mixed.mix(createvcfsinput_ch,download_fastq_to_bam_and_cram.out[0])
createvcfs(mixed)

}