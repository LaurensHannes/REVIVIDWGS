process concatvcf {

	tag "$id"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3

	input:
	tuple val(id), file(vcf1), file(vcf2), file(vcf3), file(vcf4), file(vcf5), file(vcf6), file(vcf7), file(vcf8), file(vcf9), file(vcf10), file(vcf11), file(vcf12), file(vcf13), file(vcf14), file(vcf15), file(vcf16), file(vcf17), file(vcf18), file(vcf19), file(vcf20), file(vcf21), file(vcf22), file(vcf23), file(vcf24), file(vcf25),file(vcftbi1), file(vcftbi2), file(vcftbi3), file(vcftbi4), file(vcftbi5), file(vcftbi6), file(vcftbi7), file(vcftbi8), file(vcftbi9), file(vcftbi10), file(vcftbi11), file(vcftbi12), file(vcftbi13), file(vcftbi14), file(vcftbi15), file(vcftbi16), file(vcftbi17), file(vcftbi18), file(vcftbi19), file(vcftbi20), file(vcftbi21), file(vcftbi22), file(vcftbi23), file(vcftbi24), file(vcftbi25)

	output:
	tuple val(id),file("${id}.g.vcf.gz"), file("${id}.g.vcf.gz.tbi")
	
	
"""
	bcftools concat -o ${id}.g.vcf.gz -O z --threads ${task.cpus} $vcf1 $vcf2 $vcf3 $vcf4 $vcf5 $vcf6 $vcf7 $vcf8 $vcf9 $vcf10 $vcf11 $vcf12 $vcf13 $vcf14 $vcf15 $vcf16 $vcf17 $vcf18 $vcf19 $vcf20 $vcf21 $vcf22 $vcf23 $vcf24 $vcf25
	tabix ${id}.g.vcf.gz
"""

}

process mergevcf {

	tag "${fam}_deeptrio"
	cpus 4
	time { 30.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3

	input:
	tuple val(fam), file(vcf1), file(vcf2), file(vcf3), file(vcftbi1), file(vcftbi2), file(vcftbi3),val(index),val(father),val(mother)

	output:
	tuple val(fam),val("deeptrio"),file("${fam}.vcf.gz"), file("${fam}.vcf.gz.tbi")
	
"""
	bcftools merge -o ${fam}.vcf.gz -O z --threads ${task.cpus} ${index}.vcf.gz ${father}.vcf.gz ${mother}.vcf.gz
	tabix ${fam}.vcf.gz
"""

}

process intersectvcf {

	tag "intersecting vcf for family:${fam}"
	cpus 1
	time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
	
	input: 
	tuple val(fam),val(caller1),file(vcf1),file(vcftbi1)
	tuple val(fam),val(caller2),file(vcf2),file(vcftbi2)
	
	output:
	tuple val(fam),val("consensus"),file("output/0003.vcf")
	tuple val(fam),val(caller1),file("output/0000.vcf")
	tuple val(fam),val(caller2),file("output/0001.vcf")
	
	"""
	bcftools isec -p "output"  $vcf1 $vcf2
	"""
	
}

process normalizeindels {

	tag "normalizing vcf for caller:${caller} from family:${fam}"
	cpus 4
	time { 30.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		publishDir "./results/familyvcfs", mode: 'copy', overwrite: true
	
	input: 
	tuple val(fam),val(caller), file(vcfgz), file(vcfgztbi)
	path genome 
	
	output:
	tuple val(fam),val(caller),file("${fam}.${caller}.normalized.vcf.gz"),file("${fam}.${caller}.normalized.vcf.gz.tbi")
	
	"""
	bcftools norm -m- -f $genome $vcfgz -O z -o ${fam}.${caller}.normalized.vcf.gz
	tabix ${fam}.${caller}.normalized.vcf.gz
	"""
	
	
}

	