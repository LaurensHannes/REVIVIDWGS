process vcftoolshardfilter {

		tag "$family"
			cpus 2
		cache 'deep'
		publishDir "./results/familyvcfs", mode: 'copy', overwrite: true
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi)
		
		output:
        tuple val(family),val("GATK"), file("${family}.filtered.vcf.gz"), file("${family}.filtered.vcf.gz.tbi")

		"""
		vcftools --gzvcf $vcfgz --out ${family}.hardfiltered --remove-filtered-all --minGQ 20 --minDP 10 --max-missing 1 --recode
		bcftools convert -O z -o ${family}.hardfiltered.vcf.gz ${family}.hardfiltered.recode.vcf
		tabix ${family}.hardfiltered.vcf.gz
		"""
		
		
}