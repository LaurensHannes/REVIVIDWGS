process duplicates { 

        tag "$lane"
		
		container "docker://broadinstitute/gatk"
        errorStrategy 'retry'
         maxRetries 9
		memory { 8.GB * task.attempt }
		cpus { 2 * task.attempt }
			 time { 15.minute * task.attempt }

		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 
	path home

	output:
	tuple val(id),file("${lane}.dups.bam")
	tuple val(id),file("${lane}.metrics.txt")

	
	"""

	gatk MarkDuplicates -I $bam -O ${lane}.dups.bam -M ${lane}.metrics.txt --TAGGING_POLICY OpticalOnly 


	"""
}