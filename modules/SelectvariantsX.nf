process SelectVariantsX {

		
        tag "$family"
		cpus 2
		time { 15.minute * task.attempt }
		errorStrategy 'retry'
         maxRetries 9
		container "docker://broadinstitute/gatk"
		
        input:
        tuple val(family), val(mode), file(vcfgz)
		path genome
		path dict 
		path indexes
        path ped
        path mask 
		
        output:
        tuple val(family), val("X"), val(mode), file("${family}.${mode}.X.vcf.gz"), file("${family}.${mode}.X.vcf.gz.tbi")

        """
		gatk IndexFeatureFile -I $vcfgz
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask -L chrX --exclude-filtered true --remove-unused-alternates true  --restrict-alleles-to BIALLELIC -O ${family}.${mode}.X.vcf.gz 
        """
}
