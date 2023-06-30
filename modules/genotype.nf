process genotype {

        tag "$id $chr"
		 time { 6.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	cpus 4
		myDir = file('./results/vcfs')
	myDir.mkdirs()





        input:
        tuple val(id), val(chr),file(bam),file(bai) 
        path genome 
        path indexes
		path broadinterval
		path dict		
		path mask 
		
		
        output:
        tuple val(id), file("${id}.${chr}.g.vcf.gz"), file("${id}.${chr}.g.vcf.gz.tbi")

		script:
		if(chr == "chrY")
        """
		egrep -i -w "^${chr}" ${broadinterval} > ${chr}.bed
        gatk HaplotypeCaller --verbosity INFO -ERC GVCF -ploidy 1 -L ${chr}.bed -R $genome -I $bam -O ${id}.${chr}.g.vcf.gz --sequence-dictionary ${dict} --pcr-indel-model NONE -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads ${task.cpus}
        """
		else  
		"""
		egrep -i -w "^${chr}" ${broadinterval} > ${chr}.bed
        gatk HaplotypeCaller --verbosity INFO -ERC GVCF -L ${chr}.bed -R $genome -I $bam -O ${id}.${chr}.g.vcf.gz --sequence-dictionary ${dict} --pcr-indel-model NONE -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads ${task.cpus}
        """
}