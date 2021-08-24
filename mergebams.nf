process mergebams {

	tag "$id"
    cpus 6
	 time { 2.hour * task.attempt }
	 errorStrategy 'retry' 
	maxRetries 3
	
	input:
	tuple val(id),file(bams)
	path home

	
	output:
	tuple val(id),file("${id}.bam"),file("${id}.bam.bai")

	"""
	samtools merge -@ ${task.cpus} ${id}.bam $bams
	samtools index -@ ${task.cpus} ${id}.bam
	"""

}
