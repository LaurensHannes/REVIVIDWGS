process mergevcf {

	tag "$id"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3

	input:
	tuple val(id), file(vcf1), file(vcf2), file(vcf3), file(vcf4), file(vcf5), file(vcf6), file(vcf7), file(vcf8), file(vcf9), file(vcf10), file(vcf11), file(vcf12), file(vcf13), file(vcf14), file(vcf15), file(vcf16), file(vcf17), file(vcf18), file(vcf19), file(vcf20), file(vcf21), file(vcf22), file(vcf23), file(vcf24), file(vcf25)

	output:
	tuple val(id), file("${id}.vcf.gz"), file("${id}.vcf.gz.tbi")
	
"""
	bcftools merge -o ${id}.vcf.gz -O z --threads ${task.cpus} $vcf1 $vcf2 $vcf3 $vcf4 $vcf5 $vcf6 $vcf7 $vcf8 $vcf9 $vcf10 $vcf11 $vcf12 $vcf13 $vcf14 $vcf15 $vcf16 $vcf17 $vcf18 $vcf19 $vcf20 $vcf21 $vcf22 $vcf23 $vcf24 $vcf25
	tabix ${id}.vcf.gz
"""

}