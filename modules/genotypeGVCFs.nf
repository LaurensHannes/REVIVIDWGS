process genotypeGVCFs {

	tag "$family"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
		
	input:
	tuple val(family), path(vcf) 
        path genome 
        path indexes
		path broadinterval
		path dict		
		path mask 
	
	output:
	tuple val(family), path("${family}.vcf.gz"), path("${family}.vcf.gz.tbi")
	
"""
	gatk genotypeGVCFs -R $genome -V -$vcf -O ${family}.vcf.gz -L $broadinterval --sequence-dictionary $dict
"""

}