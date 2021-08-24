#! usr/bin/env nextflow
nextflow.enable.dsl=2

include { importfastq } from './importfastq.nf'
// script parameters


params.ped ='/mnt/hdd2/data/lau/phd/testyard/FNRCP/FNRCP.ped'
params.home ='/mnt/hdd2/data/lau/phd/testyard/FNRCP'


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

importfastq (idfamily_ch, params.home) 

fastqgz_ch.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
gzipped_ch.into{temp_ch1;temp_ch2}
//temp_ch2.view()
} 
