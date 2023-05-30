process genotypeGVCFs {

	tag "$family"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	cpus 4
		cache 'lenient'
	
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
	gatk GenotypeGVCFs -R $genome -V $vcf -O ${family}.vcf.gz -L $broadinterval --sequence-dictionary $dict
"""

}

process genotypechrGVCFs {

	tag "$chr"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
		storeDir './results/vcfs'
	memory { 8.GB * task.attempt }
	cpus 4
		cache 'lenient'
	
	input:
		path vcf
		val chr
        path genome 
        path indexes
		path broadinterval
		path dict		
		path mask 
	
	output:
	path("${chr}.vcf.gz")
	path("${chr}.vcf.gz.tbi")
	
"""
gatk IndexFeatureFile -I ${chr}.g.vcf.gz
egrep -i -w "^${chr}" ${broadinterval} > ${chr}.bed
	gatk GenotypeGVCFs -R $genome -V *${chr}.g.vcf.gz -O ${chr}.vcf.gz --sequence-dictionary $dict -L ${chr}.bed
"""

}