process duplicates { 

        tag "$lane"
        errorStrategy 'retry'
         maxErrors 3
		memory { 2.GB * task.attempt }
		cpus { 18 * task.attempt }
			 time { 4.hour * task.attempt }
		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 
	path home

	output:
	tuple val(id),file("${lane}.dups.bam")
	tuple val(id),file("${lane}.metrics.txt")

	
	"""
	gatk MarkDuplicates -I $bam -O ${lane}.dups.bam -M ${lane}.metrics.txt

	"""
}