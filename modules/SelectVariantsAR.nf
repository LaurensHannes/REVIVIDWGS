process SelectVariantsAR {


        tag "$family"
		cpus 2
		time { 15.minute * task.attempt }
		errorStrategy 'retry'
         maxRetries 9
		 memory { 4.GB * task.attempt }
		container "docker://broadinstitute/gatk"
		
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi)
		path genome
		path dict 
		path indexes
        path ped
        path mask 
		
        output:
        tuple val(family), val("AR"), file("${family}.recessive.vcf.gz"), file("${family}.recessive.vcf.gz.tbi")

        """
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask --mendelian-violation true --exclude-filtered true --invert-mendelian-violation true --mendelian-violation-qual-threshold 30  --remove-unused-alternates true  --restrict-alleles-to BIALLELIC -O ${family}.recessive.vcf.gz 
        """
}
