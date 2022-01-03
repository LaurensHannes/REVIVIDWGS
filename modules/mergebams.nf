process mergebams {

	tag "$id"
    cpus 16
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
	samtools merge -@ ${task.cpus} ${id}.temp.bam $bams
	samtools sort -@ ${task.cpus} ${id}.temp.bam > ${id}.bam
	samtools index -@ ${task.cpus} ${id}.bam
	"""

}
