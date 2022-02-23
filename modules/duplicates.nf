process duplicates { 

        tag "$lane"
		
		container "docker://broadinstitute/gatk"
        errorStrategy 'retry'
         maxRetries 9
		memory { 8.GB * task.attempt }
		cpus 2
			 time { 15.minute * task.attempt }


		
	input:
	tuple val(id),val(lane),val(chr),file(bam),file(bai) 


	output:
	tuple val(id),file("${lane}.${chr}.dups.bam")
	tuple val(id),file("${lane}.metrics.txt")

	
	"""

	gatk MarkDuplicates -I $bam -O ${lane}.${chr}.dups.bam -M ${lane}.metrics.txt --TAGGING_POLICY OpticalOnly -MAX_FILE_HANDLES 2000 --SORTING_COLLECTION_SIZE_RATIO 0.025 --MAX_OPTICAL_DUPLICATE_SET_SIZE -1 --MAX_RECORDS_IN_RAM 50000 --ASSUME_SORT_ORDER coordinate


	"""
}