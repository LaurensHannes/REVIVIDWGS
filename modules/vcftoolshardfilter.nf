process vcftoolshardfilter {

		tag "$family"
			cpus 2
		
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi)
		
		output:
        tuple val(family),val("GATK"), file("${family}.filtered.vcf.gz"), file("${family}.filtered.vcf.gz.tbi")

		"""
		vcftools --gzvcf $vcfgz --out ${family}.filtered --remove-filtered-all --minGQ 20 --minDP 10 --max-missing 1 --recode
		bcftools convert -O z -o ${family}.filtered.vcf.gz ${family}.filtered.recode.vcf
		tabix ${family}.filtered.vcf.gz
		"""
		
		
}