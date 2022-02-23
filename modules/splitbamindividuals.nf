process splitbamindividuals {

		cpus { 2 * task.attempt }

		memory { 8.GB * task.attempt }
		tag "$id $chr"
			 time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		
	input:
	tuple val(id),file(bam),file(bai) 
	each val(chr)
	
	output:
	tuple val(id), val(id), file("${id}.${chr}.bam"), file("${id}.${chr}.bam.bai")

	
"""
	samtools view -@ ${task.cpus} -b $bam $chr > ${id}.${chr}.bam
	samtools index -@ ${task.cpus} ${id}.${chr}.bam 
"""
}
