process splitbamindividuals {

		cpus { 1 * task.attempt }

		memory { 2.GB * task.attempt }
		tag "$id $chr"
			 time { 30.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		
	input:
	tuple val(id),file(bam),file(bai) 
	each chr
	path genome 
	path indexes
	
	output:
	tuple val(id), val(chr), file("${id}.${chr}.bam"), file("${id}.${chr}.bam.bai")

	
"""
	samtools view -@ ${task.cpus} -T $genome -b $bam $chr > ${id}.${chr}.bam
	samtools index -@ ${task.cpus} ${id}.${chr}.bam 
"""
}
