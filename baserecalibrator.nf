process baserecalibrator {

	tag "$id"
	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
	cpus 4	
	time { 5.hour * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
	container "docker://broadinstitute/gatk"
	memory { 16.GB * task.attempt }

        input:
        tuple val(id), file(merged), file(bai) 
        path genome 
        path dict 
		path snps
        path snpsindex 

        output:
        tuple val(id), file(merged), file(bai), file("${id}.recal_data.table") 

        """
        gatk BaseRecalibrator -I $merged -R $genome -O ${id}.recal_data.table --known-sites $snps --verbosity WARNING
        """
}