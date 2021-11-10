process SelectVariantsX {


        tag "$family"
		cpus 2
		
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi)
		path genome
		path dict 
		path indexes
        path ped
        path mask 
		
        output:
        tuple val(family), val("X"), file("${family}.X.vcf.gz"), file("${family}.X.vcf.gz.tbi")

        """
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask -L chrX --exclude-filtered true --remove-unused-alternates true  --restrict-alleles-to BIALLELIC -O ${family}.X.vcf.gz 
        """
}
