process readgroups {

	tag "$lane"
        errorStrategy 'retry'
         maxRetries 3
		cpus { 4 * task.attempt }

			 time { 4.hour * task.attempt }
	input:
	tuple val(id), val(lane), file(bam) 
	path home
	
	output:
	tuple val(id),val(lane), file("${lane}.RG.bam"),file("${lane}.RG.bam.bai")
		


	"""
	gatk AddOrReplaceReadGroups -I $bam -O ${lane}.RG.bam -LB REVIVID -PL ILLUMINA -PU $lane -SM $id 
	samtools index -@ ${task.cpus} ${lane}.RG.bam
	"""

}