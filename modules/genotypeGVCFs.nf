process genotypeGVCFs {

	tag "$family"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	cpus 4
			scratch '$VSC_SCRATCH_NODE'
	stageInMode	'copy'
	
	input:
	tuple val(family), path(vcf), path(vcftbi)
        path genome 
        path indexes
		path broadinterval
		path dict		
		path mask 
	
	output:
	tuple val(family), path("${family}.vcf.gz"), path("${family}.vcf.gz.tbi")
	
"""
	gatk GenotypeGVCFs -R $genome -V $vcf -O ${family}.vcf.gz -L $broadinterval --sequence-dictionary $dict -XL $mask
"""

}