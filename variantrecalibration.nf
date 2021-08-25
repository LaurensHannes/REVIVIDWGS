process variantrecalibration {

	tag "$id"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
	time { 30.minutes * task.attempt }
	cpus 6
		 errorStrategy 'retry' 
		maxRetries 3
	
	container "docker://broadinstitute/gatk"

	input:
	tuple val(id), path(vcf) from vcf_uncallibrated_ch
	path genome from params.genome
        path dict from params.genomedict
//      path index from params.indexes
        path faidx from params.genomefai
	path snps from params.snps
	path snpsindex from params.snpsindex
	path indels from params.indels
	path indelsindex from params.indelsindex
	path mask from params.maskrepeats


	output:
	tuple val(id), file("${id}.filtered.vcf")  into individual_vcf_ch

	"""
	gatk CNNScoreVariants -V $vcf -R $genome -O ${id}.pretranched.vcf
	gatk FilterVariantTranches -V ${id}.pretranched.vcf --resource $snps --resource $indels -O ${id}.filtered.vcf --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 
	"""
}