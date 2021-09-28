process variantrecalibration {

	tag "$id"
	time { 30.minutes * task.attempt }
	cpus 6
		 errorStrategy 'retry' 
		maxRetries 3
	
	container "docker://broadinstitute/gatk"

	input:
	tuple val(id), path(vcf)
	path genome
    path dict 
	path indexes
	path snps 
	path snpsindex 
	path indels 
	path indelsindex 
	path mask


	output:
	tuple val(id), file("${id}.filtered.vcf")  

	"""
	gatk CNNScoreVariants -V $vcf -R $genome -O ${id}.pretranched.vcf
	gatk FilterVariantTranches -V ${id}.pretranched.vcf --resource $snps --resource $indels -O ${id}.filtered.vcf --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 
	"""
}