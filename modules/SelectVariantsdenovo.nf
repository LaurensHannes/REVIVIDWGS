process SelectVariantsdenovo {

        tag "$family"
		cpus 2
		time { 15.minute * task.attempt }
		errorStrategy 'retry'
         maxRetries 9
		container "docker://broadinstitute/gatk"
	
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi) 
		path genome
		path dict 
		path indexes
        path ped
        path mask 


        output:
        tuple val(family), val("denovo"), file("${family}.denovo.vcf.gz"), file("${family}.denovo.vcf.gz.tbi") 
     
		"""
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask --exclude-filtered false --mendelian-violation true --mendelian-violation-qual-threshold 30 -O ${family}.denovo.vcf.gz --remove-unused-alternates true
        """
//exclude-filtered changed to false --remove-unused-alternates to false  --restrict-alleles-to BIALLELIC (removed)
}