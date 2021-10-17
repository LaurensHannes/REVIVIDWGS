process compressandindex {

	tag "$id"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
input:
	tuple val(id), file(vcf) 

output:
        tuple val(id), file("${id}.filtered.vcf.gz"), file("${id}.filtered.vcf.gz.tbi") 

"""
        bcftools view --threads ${task.cpus} -o ${id}.filtered.vcf.gz -O z $vcf
        tabix ${id}.filtered.vcf.gz
"""
}
