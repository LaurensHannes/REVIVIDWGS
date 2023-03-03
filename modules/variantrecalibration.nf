process variantrecalibration {

	tag "$family"
	time '1h'
	cpus 6
	cache 'lenient'
	
	container "docker://broadinstitute/gatk"

	input:
	tuple val(family), path(vcf), path(vcftbi)
	path genome
    path dict 
	path indexes
	path snps 
	path snpsindex 
	path indels 
	path indelsindex 


	output:
	tuple val(family), file("${family}.filtered.vcf.gz"), file("${family}.filtered.vcf.gz.tbi")

	script:
	"""
	gatk CNNScoreVariants -V $vcf -R $genome -O ${family}.pretranched.vcf.gz
	gatk FilterVariantTranches -V ${family}.pretranched.vcf.gz --resource $snps --resource $indels -O ${family}.filtered.vcf.gz --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 
	"""
	

}

process variantcohortrecalibration {

	tag "$family"
	time '1h'
	cpus 6
	cache 'lenient'
	
	container "docker://broadinstitute/gatk"

	input:
	tuple val(family), path(vcf), path(vcftbi)
	path genome
    path dict 
	path indexes
	path snps 
	path snpsindex 
	path indels 
	path indelsindex 


	output:
	tuple val(family), file("${family}.filtered.vcf.gz"), file("${family}.filtered.vcf.gz.tbi")

	"""
	gatk VariantRecalibrator -V $vcf -R $genome --resource $snps --resource $indels -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode BOTH -O output.recal --tranches-file output.tranches 
	gatk ApplyVQSR  -V $vcf -R $genome -O ${family}.filtered.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file output.tranches --recal-file output.recal -mode BOTH
	"""
}