process compressandindex {

	tag "$id"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
input:
	tuple val(id), file(vcf) from individual_vcf_ch 

output:
        tuple val(id), file("${id}.filtered.vcf.gz"), file("${id}.filtered.vcf.gz.tbi") into individual_vcfgz_ch

"""
        bcftools view --threads ${task.cpus} -o ${id}.filtered.vcf.gz -O z ${id}.filtered.vcf
        tabix ${id}.filtered.vcf.gz
"""
}