 process alignment {

		
		memory { 8.GB * task.attempt }
		tag "$lane"
			 time { 1.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3


		storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
        indexes_ch = Channel.fromPath(params.indexes)


        input:
        tuple val(id),val(lane), file(assembled), file(forward), file(reverse) from paired_ch
        path genome from params.genome
        path indexes from indexes_ch.toList()
        path home from params.home

        output:
        tuple val(id), val(lane), file("${lane}.indexed.bam")into mapped_ch

        """
		if [ -f $home/tempstorage/${id}/${lane}.indexed.bam ] || [ -f $home/tempstorage/${id}/${lane}.RG.bam ] || [ -f $home/tempstorage/${id}/${lane}.dups.bam ] || [ -f $home/tempstorage/${id}/${id}.bam ] || [ -f $home/tempstorage/${id}/${id}.recallibrated.bam ] || [ -f $home/tempstorage/${id}/${id}.vcf ] || [ -f $home/tempstorage/${id}/${id}.filtered.vcf ]  || [ -f $home/tempstorage/${id}/${id}.filtered.vcf.gz ]
		then
		echo "done" >  /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${lane}.indexed.bam
		else
        bwa mem -t ${task.cpus} $genome $assembled | samtools view -@ ${task.cpus} -bS > ${lane}.assembled.bam
        bwa mem -t ${task.cpus} $genome $forward | samtools view -@ ${task.cpus} -bS > ${lane}.forward.bam
        bwa mem -t ${task.cpus} $genome $reverse | samtools view -@ ${task.cpus} -bS > ${lane}.reverse.bam
        samtools merge -@ ${task.cpus} ${lane}.indexed.unsorted.bam  ${lane}.assembled.bam ${lane}.forward.bam ${lane}.reverse.bam
        samtools sort -@ ${task.cpus} -o ${lane}.indexed.bam ${lane}.indexed.unsorted.bam
		echo "done" >   ${lane}.assembled.bam  ${lane}.forward.bam ${lane}.reverse.bam  ${lane}.indexed.unsorted.bam
		echo "done" >  /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${lane}*.fastq
		fi
        """
}