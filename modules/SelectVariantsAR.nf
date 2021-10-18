process SelectVariantsAR {


        tag "$family"
		cpus 4
		
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
		gatk FilterVcf -I $vcfgz -O filtered${vcfgz} --MIN_DP 10 --MIN_GQ 30
        gatk SelectVariants -V filtered${vcfgz} -ped $ped -R $genome -XL $mask --mendelian-violation true --exclude-filtered true --invert-mendelian-violation true --mendelian-violation-qual-threshold 30  --remove-unused-alternates true  --restrict-alleles-to BIALLELIC -O ${family}.recessive.vcf.gz 
        """
}
