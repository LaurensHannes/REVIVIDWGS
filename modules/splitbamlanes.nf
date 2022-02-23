process splitbamlanes {

		cpus { 2 * task.attempt }

		memory { 8.GB * task.attempt }
		tag "$lane $chr"
			 time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 
	each val(chr)
	
	output:
	tuple val(id), val(lane), file("${lane}.${chr}.bam"), file("${lane}.${chr}.bam.bai")

	
"""
	samtools view -@ ${task.cpus} -b $bam $chr > ${lane}.${chr}.bam
	samtools index -@ ${task.cpus} ${lane}.${chr}.bam 
"""
}
