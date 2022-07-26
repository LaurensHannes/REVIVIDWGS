process mergevcf {

	tag "$id"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3

	input:
	tuple val(id), tuple(vcf)

	output:
	tuple val(family), path("${family}.vcf.gz"), path("${family}.vcf.gz.tbi")
	
"""
	bcftools merge -o ${family}.vcf.gz -O z --threads ${task.cpus} $vcf 
	tabix ${family}.vcf.gz
"""

}