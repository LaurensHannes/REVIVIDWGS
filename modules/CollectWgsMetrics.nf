process CollectWgsMetrics {

	tag "$id"
	cpus 1	
	time { 3.hour * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
	container "docker://broadinstitute/gatk"
	memory { 4.GB * task.attempt }
	publishDir "./QC/$id", mode: 'copy', overwrite: false

        input:
        tuple val(id), file(merged), file(bai) 
        path genome 

        output:
        tuple val(id), file("${id}_wgs_metrics.txt") 

        """
        gatk CollectWgsMetrics -I $merged -R $genome -O ${id}_wgs_metrics.txt
        """
}