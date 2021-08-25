#! usr/bin/env nextflow
nextflow.enable.dsl=2

include { importfastq } from './importfastq.nf'
include { pear } from './pear.nf'
include { alignment } from './alignment.nf'
include { readgroups } from './readgroups.nf'
include { duplicates } from './duplicates.nf'
include { mergebams } from './mergebams.nf'
include { generateCRAM } from './generateCRAM.nf'
include { baserecalibrator } from './baserecalibrator.nf'


// script parameters



indexes_ch = Channel.fromPath(params.indexes).toList()

//script
workflow {
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

importfastq(idfamily_ch, params.home) 

importfastq.out.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
//gzipped_ch.view()

pear(gzipped_ch, params.home)
alignment(pear.out, params.genome,indexes_ch, params.home)
alignment.out.view()
readgroups(alignment.out,params.home)
duplicates(readgroups.out,params.home)
mergebams(duplicates.out[0].groupTuple(),params.home)
generateCRAM(mergebams.out,params.genome,indexes_ch)
baserecalibrator(mergebams.out,params.genome, params.genomedict, params.snps, params.snpsindex)
baserecalibrator.out.view()

}
