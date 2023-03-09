process mergebams {

	tag "$id"
    cpus 18
	 time { 2.hour * task.attempt }
	 errorStrategy 'retry' 
	maxRetries 3
	myDir = file('./results/bams')
	myDir.mkdirs()
	publishDir 'vsc-hard-mounts/leuven-archive/arc_00086/results/bams', mode: 'copy', overwrite: true
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

process stmergebams {

	tag "$id"
    cpus 18
	 time { 2.hour * task.attempt }
	 errorStrategy 'retry' 
	maxRetries 3
	myDir = file('./results/bams')
	myDir.mkdirs()
	storeDir './results/bams'
	input:
	tuple val(id),file(bams)
	path home

	
	output:
	file("${id}.bam")
	file("${id}.bam.bai")


	"""
	samtools merge -@ ${task.cpus} ${id}.bam $bams
	samtools index -@ ${task.cpus} ${id}.bam
	"""

}