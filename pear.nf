nextflow.enable.dsl=2
include { delete_file } from './delete_file.nf'
process pear {



        tag "$lane"
        time { 1.hour * task.attempt }
		memory '2 GB'
		cpus 15
        errorStrategy 'retry'
        maxRetries 3

        input:
        tuple val(id), val(lane),file(R1), file(R2)
		path home
		
        output:
        tuple val(id), val(lane), file("${lane}.assembled.fastq"), file("${lane}.unassembled.forward.fastq"), file("${lane}.unassembled.reverse.fastq")

		script:
		"""
		pear -j ${task.cpus} -f $R1 -r $R2 -o $lane
		"""
		delete_file(R1)
		}