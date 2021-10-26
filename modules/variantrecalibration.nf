process variantrecalibration {

	tag "$family"
	time { 30.minutes * task.attempt }
	cpus 6
		 errorStrategy 'retry' 
		maxRetries 3
	
	container "docker://broadinstitute/gatk"

	input:
	tuple val(family), path(vcf)
	path genome
    path dict 
	path indexes
	path snps 
	path snpsindex 
	path indels 
	path indelsindex 
	path mask


	output:
	tuple val(family), file("${family}.filtered.vcf.gz")  

	"""
	gatk CNNScoreVariants -V $vcf -R $genome -O ${family}.pretranched.vcf.gz
	gatk FilterVariantTranches -V ${family}.pretranched.vcf.gz --resource $snps --resource $indels -O ${family}.filtered.vcf.gz --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 
	"""
}