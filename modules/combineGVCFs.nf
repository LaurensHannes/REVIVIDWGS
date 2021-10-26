process combineGVCFs {

	tag "$family"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
		
	input:
	tuple val(family), path(vcf1), path(vcf2), path(vcf3) 
	path genome 
	
	output:
	tuple val(family), path("${family}.g.vcf.gz")
	
"""
	gatk CombineGVCFs -R $genome -V $vcf1 -V $vcf2 -V $vcf3 -O ${family}.g.vcf.gz
"""

}