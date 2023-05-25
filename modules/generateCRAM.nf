process generateCRAM {

	tag "$id"
	cpus 4
		myDir = file('./results/crams')
		myDir.mkdirs()
		publishDir "$arch/results/crams/$id", mode: 'copy', overwrite: false
        errorStrategy 'retry'
         maxRetries 3
			 time { 2.hours * task.attempt }
			
		input:
		tuple val(id),file(bam),file(bai) 
		path genome
        path indexes
		val arch		

		output:
		tuple val(id),file("${id}.cram"),file("${id}.cram.crai")

		script:
        outputcram = file("${arch}/results/crams/${id}/${id}.bam")
        if (!outputcram.exists())
		"""
		samtools view -@ ${task.cpus} -C -o ${id}.cram -T $genome $bam 
		samtools index -@ ${task.cpus} -c ${id}.cram
		"""
		else
		"""
		touch ${id}.cram
        touch ${id}.cram.crai
		"""
}