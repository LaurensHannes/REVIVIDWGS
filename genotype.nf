process genotype {

        tag "$id"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"


        input:
        tuple val(id), file(bam),file(bai) from BQSR_applied_ch
        path genome from params.genome
        path dict from params.genomedict
//        path index from params.indexes
        path faidx from params.genomefai
	path mask from params.maskrepeats

        output:
        tuple val(id), file("${id}.vcf") into vcf_uncallibrated_ch

        """
        gatk HaplotypeCaller --verbosity INFO -XL $mask -R $genome -I $bam -O ${id}.vcf --sequence-dictionary ${dict} --native-pair-hmm-threads ${task.cpus}
        """

}