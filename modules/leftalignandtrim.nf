process leftalignandtrim {

        tag "$family"
		 time { 1.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	cpus 1



        input:
        tuple val(family), file(vcfgz), file(vcfgztbi)
        path genome 
		path indexes
		path dict
		
        output:
        tuple val(family), file("${family}.alignedandtrimmed.vcf.gz"), file("${family}.alignedandtrimmed.vcf.gz.tbi") 
        
		
		"""
		gatk LeftAlignAndTrimVariants -R $genome -V $vcfgz -O ${family}.alignedandtrimmed.vcf.gz
 
		"""		
		
}