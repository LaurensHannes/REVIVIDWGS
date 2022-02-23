process splitbamlanes {

		cpus { 1 * task.attempt }

		memory { 2.GB * task.attempt }
		tag "$lane $chr"
			 time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		
	input:
	tuple val(id),val(lane),file(bam),file(bai) 
	each chr
	
	output:
	tuple val(id), val(lane),val(chr), file("${lane}.${chr}.bam"), file("${lane}.${chr}.bam.bai")

	
"""
	samtools view -@ ${task.cpus} -b $bam $chr > ${lane}.${chr}.bam
	samtools index -@ ${task.cpus} ${lane}.${chr}.bam 
"""
}
