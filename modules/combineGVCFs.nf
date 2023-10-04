process combineGVCFs {

	
	
	cpus 4
	time { 4.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
		memory { 48.GB * task.attempt }
			cache 'lenient'
	publishDir "./results/familyvcfs", mode: 'copy', overwrite: true
	
	input:
	tuple val(family), path(vcf1), path(vcf2), path(vcf3)
	path genome 
	path indexes
	path dict
	
	output:
	tuple val(family), path("${family}.g.vcf.gz"), path("${family}.g.vcf.gz.tbi")
	
"""
find \$PWD -name "*.vcf.gz" > input.list
lines=\$(cat input.list)

for line in \$lines
do

	gatk IndexFeatureFile -I \$line & 
done
wait
	gatk CombineGVCFs -R $genome -V $vcf1 -V $vcf2 -V $vcf3 -O ${family}.g.vcf.gz
"""

}

process combinecohortGVCFs {

	tag "cohort"
	cpus { 4 * task.attempt }
	time { 6.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
		memory { 46.GB * task.attempt }
			cache 'lenient'
		publishDir "./results/familyvcfs", mode: 'copy', overwrite: true
	input:
	path vcf
	path genome 
	path indexes
	path dict
	
	output:
	path("cohort.g.vcf.gz")
	path("cohort.g.vcf.gz.tbi")
	
"""
find \$PWD -name "*.vcf.gz" > input.list
lines=\$(cat input.list)

for line in \$lines
do

	gatk IndexFeatureFile -I \$line & 
done
wait
	gatk CombineGVCFs -R $genome  -O cohort.g.vcf.gz -V input.list --sequence-dictionary $dict
"""

}

process combinechrVCFs {

	tag "${fam}"
	cpus { 4 * task.attempt }
	time { 6.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
		memory { 46.GB * task.attempt }
			cache 'lenient'
		publishDir "./results/familyvcfs", mode: 'copy', overwrite: true
		
	input:
	path vcf
	path genome 
	path indexes
	path dict
	file ped 

	output:
	
	path("${fam}.vcf.gz")
	path("${fam}.vcf.gz.tbi")
	
	script:

	fam = ped.baseName

"""

find \$PWD -name "*.vcf.gz" > input.list
lines=\$(cat input.list)

for line in \$lines
do

	gatk IndexFeatureFile -I \$line & 
done
wait
	gatk MergeVcfs -R $genome  -O ${fam}.vcf.gz -I input.list -D $dict
"""

}