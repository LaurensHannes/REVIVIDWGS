 process alignment {


		tag "$lane"
			 time { 12.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		


		input:
        tuple val(id), val(lane),file(R1), file(R2)
        path genome
        path indexes
        path home
        val arch

        output:
        tuple val(id), val(lane), file("${lane}.indexed.bam"),file("${lane}.indexed.bam.bai")
        
        
        
        script:
        outputbam = "${arch}/results/bams/${id}/${id}.bam"
        if (!outputbam.exists())
        """
        bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${id}\\tSM:${id}\\tLB:REVIVID\\tPL:ILLUMINA\\tPU:${lane}" $genome $R1 $R2 | samtools sort -@ ${task.cpus} -o ${lane}.indexed.bam
        samtools index -@ ${task.cpus} ${lane}.indexed.bam
		"""
        else 
        """
        touch ${lane}.indexed.bam
        touch ${lane}.indexed.bam.bai
        """
}