 process alignment {


		cpus { 72 }

		memory { 128.GB * task.attempt }
		tag "$lane"
			 time { 8.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		


		input:
        tuple val(id), val(lane),file(R1), file(R2)
        path genome
        path indexes
        path home

        output:
        tuple val(id), val(lane), file("${lane}.indexed.bam"),file("${lane}.indexed.bam.bai")

        """
        bwa-mem2 mem -t 60 -R "@RG\\tID:${id}\\tSM:${id}\\tLB:REVIVID\\tPL:ILLUMINA\\tPU:${lane}" $genome $R1 $R2 | samtools sort -@ 12 -o ${lane}.indexed.bam
        samtools index -@ ${task.cpus} ${lane}.indexed.bam
		"""
}