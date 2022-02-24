process applyBQSR {

        tag "$id $chr"
		 time { 30.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		cpus 2

        input:

		tuple val(id), val(chr) , file(bam), file(bai), file(table)
        path genome
		path indexes 
        path dict 


        output:
        tuple val(id), val(chr) ,file("${id}.${chr}.recallibrated.bam"),file("${id}.${chr}.recallibrated.bam.bai") 

        """
        gatk ApplyBQSR -R $genome -I $bam -bqsr-recal-file $table -O ${id}.${chr}.recallibrated.bam
		samtools index -@ ${task.cpus} ${id}.${chr}.recallibrated.bam
		
		
        """

}