process combineGVCFs {

	tag "$family"
	cpus 4
	time { 2.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
		scratch '$VSC_SCRATCH_NODE'
	stageInMode	'copy'
		
	input:
	tuple val(family), path(vcf1), path(vcf2), path(vcf3)
	path genome 
	path indexes
	path dict
	
	output:
	tuple val(family), path("${family}.g.vcf.gz"), path("${family}.g.vcf.gz.tbi")
	
"""
	gatk IndexFeatureFile -I $vcf1
	gatk IndexFeatureFile -I $vcf2
	gatk IndexFeatureFile -I $vcf3
	gatk CombineGVCFs -R $genome -V $vcf1 -V $vcf2 -V $vcf3 -O ${family}.g.vcf.gz
"""

}