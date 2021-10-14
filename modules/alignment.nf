 process alignment {


		cpus { 16 * task.attempt }

		memory { 8.GB * task.attempt }
		tag "$lane"
			 time { 4.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3


		input:
        tuple val(id),val(lane), file(assembled), file(forward), file(reverse)
        path genome
        path indexes
        path home

        output:
        tuple val(id), val(lane), file("${lane}.indexed.bam")

        """
        bwa mem -t ${task.cpus} $genome $assembled | samtools view -@ ${task.cpus} -bS > ${lane}.assembled.bam
        bwa mem -t ${task.cpus} $genome $forward | samtools view -@ ${task.cpus} -bS > ${lane}.forward.bam
        bwa mem -t ${task.cpus} $genome $reverse | samtools view -@ ${task.cpus} -bS > ${lane}.reverse.bam
        samtools merge -@ ${task.cpus} ${lane}.indexed.unsorted.bam  ${lane}.assembled.bam ${lane}.forward.bam ${lane}.reverse.bam
        samtools sort -@ ${task.cpus} -o ${lane}.indexed.bam ${lane}.indexed.unsorted.bam
		rm ${lane}.assembled.bam ${lane}.forward.bam ${lane}.reverse.bam ${lane}.indexed.unsorted.bam
        """
}