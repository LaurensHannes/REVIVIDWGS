process deepmosaic {
	tag "$id"
	cpus 16
	time { 4.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://laurenshannes/deepmosaic:latest"
	publishDir "./results/deepmosaic/$id", mode: 'copy', overwrite: true
		cache 'lenient'
	memory { 32.GB * task.attempt }
	
	input:
	tuple val(id), val(chr),path(vcf),path(bam),path(bai)
	path genome 
	path indexes
	path dict
	val coverage 
	val(programpath)
    val(humandbpath)


	output:
	tuple val(id), path("${id}.${chr}.deepmosaic.txt")
	
"""
echo -e "#sample_name\tbam\tvcf\tdepth\tsex" >> input.txt
echo -e "${id}\t${bam}\t${vcf}\t${coverage}\tM" >> input.txt

/DeepMosaic/deepmosaic/deepmosaic-draw -i input.txt -o ./ -a ${programpath}/annovar/ -b hg38 -db gnomad30 -t ${task.cpus}
/DeepMosaic/deepmosaic/deepmosaic-predict -i features.txt  -o ${id}.${chr}.deepmosaic.txt -m efficientnet-b4_epoch_6.pt -b 10 -gb hg38
"""

}

