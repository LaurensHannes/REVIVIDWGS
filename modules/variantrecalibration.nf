process variantrecalibration {

	tag "$family"
	time { 30.minute * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
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
	time { 30.minute * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
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
	gatk VariantRecalibrator -V $vcf -R $genome -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 --resource:1000G,known=false,training=true,truth=false,prior=10.0 $snps -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O output_snps.recal --tranches-file output_snp.tranches 
	gatk VariantRecalibrator -V $vcf -R $genome -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 --resource:mills,known=false,training=true,truth=true,prior=12 $indels -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL -O output_indels.recal --tranches-file output_indel.tranches 
	gatk ApplyVQSR  -V $vcf -R $genome -O snp.filtered.vcf.gz --truth-sensitivity-filter-level 99.7 --tranches-file output_snp.tranches --recal-file output_snps.recal -mode SNP
	gatk ApplyVQSR  -V snp.filtered.vcf.gz -R $genome -O ${family}.filtered.vcf.gz --truth-sensitivity-filter-level 99.7 --tranches-file output_indels.tranches --recal-file output_indels.recal -mode INDEL
	"""
}