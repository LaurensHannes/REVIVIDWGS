process fastQC {

        tag "$id"
		 time { 3.hour * task.attempt }
		 errorStrategy 'retry' 
		 	cpus 1
		maxRetries 3
		container "biocontainers/fastqc:v0.11.9_cv8"
	     memory { 2.GB * task.attempt }
		publishDir "./QC/$id", mode: 'copy', overwrite: false
	
	input: 
	        tuple val(id), val(lane),file(R1), file(R2)
				
	output:
			tuple val(), val(lane), file("*R1_fastqc.html"), file("*R2_fastqc.html")
	"""
	fastqc $R1 $R2
	"""
	} 
	
	