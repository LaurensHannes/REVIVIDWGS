process genotype {

        tag "$id"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	cpus 8



        input:
        tuple val(id), file(bam),file(bai) 
        path genome 
        path indexes
		path broadinterval
		path dict		
		path mask 
		
		
        output:
        tuple val(id), file("${id}.vcf") 
        """
        gatk HaplotypeCaller --verbosity INFO -L $broadinterval -XL $mask -R $genome -I $bam -O ${id}.vcf --sequence-dictionary ${dict} -pcr_indel_model NONE --native-pair-hmm-threads ${task.cpus}
        """

}