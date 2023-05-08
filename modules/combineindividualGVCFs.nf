process combineindividualGVCFs {

	tag "$id"
	cpus 36
	time { 4.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	publishDir './results/vcfs', mode: 'copy', overwrite: false
		cache 'lenient'
	memory { 32.GB * task.attempt }
	
	input:
	tuple val(id), path(vcf1), path(vcf2), path(vcf3),path(vcf4) ,path(vcf5) ,path(vcf6) ,path(vcf7) ,path(vcf8) ,path(vcf9) ,path(vcf10) ,path(vcf11) ,path(vcf12) ,path(vcf13) ,path(vcf14) ,path(vcf15) ,path(vcf16) ,path(vcf17) ,path(vcf18) ,path(vcf19) ,path(vcf20) ,path(vcf21) ,path(vcf22) ,path(vcf23) ,path(vcf24) 

	path genome 
	path indexes
	path dict
	
	output:
	tuple val(id), path("${id}.g.vcf.gz")
	
"""
find \$PWD -name "*.vcf.gz" > input.list
lines=\$(cat input.list)

for line in \$lines
do

	gatk IndexFeatureFile -I \$line & 
done
sleep 300
	gatk CombineGVCFs -R $genome -V input.list -O ${id}.g.vcf.gz
"""

}

process combinechrGVCFs {

	tag "$chr"
	cpus 2
	time { 16.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	storeDir './results/gvcfs'
	memory { 32.GB * task.attempt }
	scratch true
	input:
	
	path vcf
	val chr
	path genome 
	path indexes
	path dict
	path broadinterval

	output:
	path("${chr}.g.vcf.gz")
	
script:


"""
mkdir temp
egrep -i -w "^${chr}" ${broadinterval} > ${chr}.bed

find \$PWD -name "*.${chr}.*.vcf.gz" > initial.list
initials=\$(cat initial.list)

for ini in \$initials
do
	cp \$ini temp/  
done
find \$PWD/temp -name "*.${chr}.*.vcf.gz" > input.list
lines=\$(cat input.list)
for line in \$lines
do
	gatk IndexFeatureFile -I \$line & 
done
sleep 180
	gatk CombineGVCFs -R $genome -V input.list -O ${chr}.g.vcf.gz -L ${chr}.bed
	rm -r temp
"""

}