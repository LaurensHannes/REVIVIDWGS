process duplicates { 

        tag "$lane"
        errorStrategy 'retry'
         maxErrors 3
		memory { 60.GB * task.attempt }
		cpus { 18 * task.attempt }
			 time { 4.hour * task.attempt }
		scratch './scratch'
		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 
	path home

	output:
	tuple val(id),file("${lane}.dups.bam")
	tuple val(id),file("${lane}.metrics.txt")

	
	"""
	gatk MarkDuplicates -I $bam -O ${lane}.dups.bam -M ${lane}.metrics.txt --MAX_RECORDS_IN_RAM 2000000 --COMPRESSION_LEVEL 9

	"""
}