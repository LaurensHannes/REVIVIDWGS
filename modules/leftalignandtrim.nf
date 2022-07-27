process leftalignandtrim {

        tag "$family"
		 time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 4.GB * task.attempt }
	cpus 2
publishDir "./results/familyvcfs", mode: 'copy', overwrite: true


        input:
        tuple val(family),val(caller), file(vcfgz), file(vcfgztbi)
        path genome 
		path indexes
		path dict
		path broadinterval
		
        output:
        tuple val(family), file("${family}.${caller}.alignedandtrimmed.vcf.gz")
        
		
		"""
		gatk LeftAlignAndTrimVariants -R $genome -V $vcfgz -O ${family}.${caller}.alignedandtrimmed.vcf.gz -L $broadinterval
 
		"""		
		
}