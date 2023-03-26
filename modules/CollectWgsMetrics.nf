process CollectWgsMetrics {

	tag "$id"
	cpus 1	
	time { 6.hour * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
	container "docker://broadinstitute/gatk"
	memory { 32.GB * task.attempt }
	publishDir "$arch/QC/$id", mode: 'copy', overwrite: false

        input:
        tuple val(id), file(merged), file(bai) 
        path genome 
        val arch

        output:
        tuple val(id), file("${id}_wgs_metrics.txt") 

        """
        gatk CollectWgsMetrics -I $merged -R $genome -O ${id}_wgs_metrics.txt
        """
}