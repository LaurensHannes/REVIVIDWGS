process combineindividualGVCFs {

	tag "$id"
	cpus 4
	time { 2.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"

		
	input:
	tuple val(id), path(vcf1), path(vcf2), path(vcf3),path(vcf4) ,path(vcf5) ,path(vcf6) ,path(vcf7) ,path(vcf8) ,path(vcf9) ,path(vcf10) ,path(vcf11) ,path(vcf12) ,path(vcf13) ,path(vcf14) ,path(vcf15) ,path(vcf16) ,path(vcf17) ,path(vcf18) ,path(vcf19) ,path(vcf20) ,path(vcf21) ,path(vcf22) ,path(vcf23) ,path(vcf24) ,path(vcf25) 

	path genome 
	path indexes
	path dict
	
	output:
	tuple val(id), path("${id}.g.vcf.gz"), path("${id}.g.vcf.gz.tbi")
	
"""
	gatk IndexFeatureFile -I $vcf1
	gatk IndexFeatureFile -I $vcf2
	gatk IndexFeatureFile -I $vcf3
	gatk IndexFeatureFile -I $vcf4
	gatk IndexFeatureFile -I $vcf5
	gatk IndexFeatureFile -I $vcf6
	gatk IndexFeatureFile -I $vcf7
	gatk IndexFeatureFile -I $vcf8
	gatk IndexFeatureFile -I $vcf9
	gatk IndexFeatureFile -I $vcf10
	gatk IndexFeatureFile -I $vcf11
	gatk IndexFeatureFile -I $vcf12
	gatk IndexFeatureFile -I $vcf13
	gatk IndexFeatureFile -I $vcf14
	gatk IndexFeatureFile -I $vcf15
	gatk IndexFeatureFile -I $vcf16
	gatk IndexFeatureFile -I $vcf17
	gatk IndexFeatureFile -I $vcf18
	gatk IndexFeatureFile -I $vcf19
	gatk IndexFeatureFile -I $vcf20
	gatk IndexFeatureFile -I $vcf21
	gatk IndexFeatureFile -I $vcf22
	gatk IndexFeatureFile -I $vcf23
	gatk IndexFeatureFile -I $vcf24
	gatk IndexFeatureFile -I $vcf25
	gatk CombineGVCFs -R $genome -V $vcf1 -V $vcf2 -V $vcf3 -V vcf4 -V vcf5 -V vcf6 -V vcf7 -V vcf8 -V vcf9 -V vcf10 -V vcf11 -V vcf12 -V vcf13 -V vcf14 -V vcf15 -V vcf16 -V vcf17 -V vcf18 -V vcf19 -V vcf20 -V vcf21 -V vcf22 -V vcf23 -V vcf24 -V vcf25
 -O ${id}.g.vcf.gz
"""

}