process SelectVariantsdenovo {

        tag "$family"
		cpus 2
		time { 30.minute * task.attempt }
		errorStrategy 'retry'
         maxRetries 9
		 memory { 16.GB * task.attempt }
		container "docker://broadinstitute/gatk"
	
        input:
        tuple val(family), val(mode),file(vcfgz)
		path genome
		path dict 
		path indexes
        path ped
        path mask 


        output:
        tuple val(family), val("denovo"), val(mode), file("${family}.${mode}.denovo.vcf.gz"), file("${family}.${mode}.denovo.vcf.gz.tbi") 
     
		"""
		gatk IndexFeatureFile -I $vcfgz
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask --exclude-filtered false --mendelian-violation true --mendelian-violation-qual-threshold 30 -O ${family}.${mode}.denovo.vcf.gz --remove-unused-alternates true
        """
//exclude-filtered changed to false --remove-unused-alternates to false  --restrict-alleles-to BIALLELIC (removed)
}