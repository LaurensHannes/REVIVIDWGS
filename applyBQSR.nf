process applyBQSR {

        tag "$id"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
		 time { 5.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		cpus 4


        input:
        //tuple val(id), file(bam), file(bai) from mergedbam3_ch
		tuple val(id), file(bam), file(bai), file(table) from recal_data_ch
        path genome from params.genome
     //   path table from recal_data_ch
	path faidx from params.genomefai
        path dict from params.genomedict


        output:
        tuple val(id),file("${id}.recallibrated.bam"),file("${id}.recallibrated.bam.bai") into BQSR_applied_ch

        """
        gatk ApplyBQSR -R $genome -I $bam -bqsr-recal-file $table -O ${id}.recallibrated.bam
		samtools index -@ ${task.cpus} ${id}.recallibrated.bam
		echo "done" >  /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${id}.bam /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${id}.bam.bai
		
        """

}