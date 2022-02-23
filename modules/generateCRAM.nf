process generateCRAM {

	tag "$id"
	cpus 4
		myDir = file('./results/crams')
		myDir.mkdirs()
		publishDir './results/crams', mode: 'copy', overwrite: true

			
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