process combineindividualGVCFs {

	tag "$id"
	cpus 4
	time { 2.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	publishDir './results/vcfs', mode: 'copy', overwrite: false
	cache 'lenient'
		
	input:
	tuple val(id), path(vcf1), path(vcf2), path(vcf3),path(vcf4) ,path(vcf5) ,path(vcf6) ,path(vcf7) ,path(vcf8) ,path(vcf9) ,path(vcf10) ,path(vcf11) ,path(vcf12) ,path(vcf13) ,path(vcf14) ,path(vcf15) ,path(vcf16) ,path(vcf17) ,path(vcf18) ,path(vcf19) ,path(vcf20) ,path(vcf21) ,path(vcf22) ,path(vcf23) ,path(vcf24) ,path(vcf25) 

	path genome 
	path indexes
	path dict
	
	output:
	tuple val(id), path("${id}.g.vcf.gz")
	
"""
find \$PWD -name "*.vcf.gz" > input.list
lines=\$(cat input.list)
N=${task.cpus}
for line in \$lines
do
	((i=i%N)); ((i++==0)) && wait
	gatk IndexFeatureFile -I \$line & 
done
	gatk CombineGVCFs -R $genome -V input.list -O ${id}.g.vcf.gz
"""

}