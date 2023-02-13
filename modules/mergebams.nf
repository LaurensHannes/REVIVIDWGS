process mergebams {

	tag "$id"
    cpus 18
	 time { 2.hour * task.attempt }
	 errorStrategy 'retry' 
	maxRetries 3
	myDir = file('./results/bams')
	myDir.mkdirs()
	publishDir './results/bams', mode: 'copy', overwrite: true
	input:
	tuple val(id),file(bams)
	path home

	
	output:
	tuple val(id),file("${id}.bam"),file("${id}.bam.bai")
	tuple val(id),val("true")

	"""
	samtools merge -@ ${task.cpus} ${id}.bam $bams
	samtools index -@ ${task.cpus} ${id}.bam
	"""

}
