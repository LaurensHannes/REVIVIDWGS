process leftalignandtrim {

        tag "$id"
		 time { 1.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	cpus 1



        input:
        tuple val(id), file(vcf)
        path genome 
	
		
        output:
        tuple val(id), file("${id}.alignedandtrimmed.vcf") 
        
		
		"""
		gatk LeftAlignAndTrimVariants -R $genome -V $vcf -O ${id}.alignedandtrimmed.vcf
 
		"""		
		
}