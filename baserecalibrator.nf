process baserecalibrator {

	tag "$id"
	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
	cpus 4	
	time { 5.hour * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
	container "docker://broadinstitute/gatk"
	memory { 32.GB * task.attempt }

        input:
        tuple val(id), file(merged), file(bai) from mergedbam2_ch
        path genome from params.genome
       path faidx from params.genomefai
        path dict from params.genomedict
	 path snps from params.snps
        path snpsindex from params.snpsindex

        output:
        tuple val(id), file(merged), file(bai), file("${id}.recal_data.table") into recal_data_ch

        """
        gatk BaseRecalibrator -I $merged -R $genome -O ${id}.recal_data.table --known-sites $snps --verbosity WARNING
        """
}