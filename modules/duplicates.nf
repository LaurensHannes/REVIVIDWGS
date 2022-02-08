process duplicates { 

        tag "$lane"
		
		container "docker://broadinstitute/gatk"
        errorStrategy 'retry'
         maxRetries 9
		memory { 32.GB * task.attempt }
		cpus 2
			 time { 30.minute * task.attempt }
		scratch '/lustre1/project/stg_00086/scratch'
			stageInMode	'copy'

		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 


	output:
	tuple val(id),file("${lane}.dups.bam")
	tuple val(id),file("${lane}.metrics.txt")

	
	"""

	gatk MarkDuplicates -I $bam -O ${lane}.dups.bam -M ${lane}.metrics.txt --TAGGING_POLICY OpticalOnly 


	"""
}