process splitbamindividuals {

		cpus { 1 * task.attempt }

		memory { 2.GB * task.attempt }
		tag "$id $chr"
			 time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		
	input:
	tuple val(id),file(bam),file(bai) 
	each chr
	
	output:
	tuple val(id), file("${id}.${chr}.bam"), file("${id}.${chr}.bam.bai")

	
"""
	samtools view -@ ${task.cpus} -b $bam $chr > ${id}.${chr}.bam
	samtools index -@ ${task.cpus} ${id}.${chr}.bam 
"""
}
