process mergevcf {

	tag "$id"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3

	input:
	tuple val(id), file(vcf)

	output:
	tuple val(id), path("${id}.vcf.gz"), path("${id}.vcf.gz.tbi")
	
"""
	bcftools merge -o $id.vcf.gz -O z --threads ${task.cpus} $vcf 
	tabix $id.vcf.gz
"""

}