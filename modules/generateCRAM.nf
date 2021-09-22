process generateCRAM {

	tag "$id"
		myDir = file('./results/crams')
		myDir.mkdirs()
			
		input:
		tuple val(id),file(bam),file(bai) 
		path genome
        path indexes
				
		output:
		tuple val(id),file("${id}.cram"),file("${id}.cram.crai")
		tuple val(id),file(bam),file(bai)
		
		"""
		samtools view -@ ${task.cpus} -C -o ${id}.cram -T $genome $bam 
		samtools index -@ ${task.cpus} -c ${id}.cram
		"""
		
}