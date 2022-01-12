process duplicates { 

        tag "$lane"
		
		container "docker://broadinstitute/gatk"
        errorStrategy 'retry'
         maxRetries 3
		memory { 20.GB * task.attempt }
		cpus { 2 * task.attempt }
			 time { 30.minute * task.attempt }

		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 
	path home

	output:
	tuple val(id),file("${lane}.dups.bam")
	tuple val(id),file("${lane}.metrics.txt")

	
	"""
	gatk MarkDuplicatesSpark -I $bam -O ${lane}.dups.bam -M ${lane}.metrics.txt --duplicate-tagging-policy OpticalOnly --tmp-dir $home/temp

	"""
}