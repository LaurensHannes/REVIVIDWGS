process applyBQSR {

        tag "$id"
		 time { 5.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		cpus 4
		scratch '$VSC_SCRATCH_NODE'
	stageInMode	'copy'

        input:

		tuple val(id), file(bam), file(bai), file(table)
        path genome
		path indexes 
        path dict 


        output:
        tuple val(id),file("${id}.recallibrated.bam"),file("${id}.recallibrated.bam.bai") 

        """
        gatk ApplyBQSR -R $genome -I $bam -bqsr-recal-file $table -O ${id}.recallibrated.bam
		samtools index -@ ${task.cpus} ${id}.recallibrated.bam
		
		
        """

}