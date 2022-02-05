process baserecalibrator {

	tag "$id"
	cpus 4	
	time { 5.hour * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
	container "docker://broadinstitute/gatk"
	memory { 16.GB * task.attempt }
	scratch '$VSC_SCRATCH'
	
        input:
        tuple val(id), file(merged), file(bai) 
        path genome 
		path indexes
        path dict 
		path snps
        path snpsindex 

        output:
        tuple val(id), file(merged), file(bai), file("${id}.recal_data.table") 

        """
        gatk BaseRecalibrator -I $merged -R $genome -O ${id}.recal_data.table --known-sites $snps --verbosity WARNING
        """
}