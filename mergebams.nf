process mergebams {

	tag "$id"
    cpus 16
	 time { 2.hour * task.attempt }
	 errorStrategy 'retry' 
	maxRetries 3
	publishDir '/results/bams', mode: 'copy', overwrite: false
	
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
