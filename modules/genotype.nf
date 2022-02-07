process genotype {

        tag "$id"
		 time { 12.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 16.GB * task.attempt }
	cpus 8
		myDir = file('./results/vcfs')
	myDir.mkdirs()
	publishDir './results/vcfs', mode: 'copy', overwrite: false
	scratch '$VSC_SCRATCH_NODE'



        input:
        tuple val(id), file(bam),file(bai) 
        path genome 
        path indexes
		path broadinterval
		path dict		
		path mask 
		
		
        output:
        tuple val(id), file("${id}.g.vcf.gz")
        """
        gatk HaplotypeCaller  --java-options "-Xmx16g" --verbosity INFO -ERC GVCF -L $broadinterval -XL $mask -R $genome -I $bam -O ${id}.g.vcf.gz --sequence-dictionary ${dict} --pcr-indel-model NONE -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads ${task.cpus}
        """

}