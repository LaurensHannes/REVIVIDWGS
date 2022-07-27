process SelectVariantsAR {


        tag "$family"
		cpus 2
		time { 30.minute * task.attempt }
		errorStrategy 'retry'
         maxRetries 9
		 memory { 16.GB * task.attempt }
		container "docker://broadinstitute/gatk"
		
        input:
        tuple val(family), file(vcfgz)
		path genome
		path dict 
		path indexes
        path ped
        path mask 
		
        output:
        tuple val(family), val("AR"), val(mode), file("${family}.${mode}.recessive.vcf.gz"), file("${family}.${mode}.recessive.vcf.gz.tbi")

        """
		gatk IndexFeatureFile -I $vcfgz
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask --mendelian-violation true --exclude-filtered true --invert-mendelian-violation true --mendelian-violation-qual-threshold 30  --remove-unused-alternates true  --restrict-alleles-to BIALLELIC -O ${family}.${mode}.recessive.vcf.gz 
        """
}
