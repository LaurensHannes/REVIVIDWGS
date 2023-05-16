process fastQC {

        tag "$id"
		 time { 3.hour * task.attempt }
		 errorStrategy 'retry' 
		 	cpus 1
		maxRetries 3
		container "biocontainers/fastqc:v0.11.9_cv8"
	     memory { 32.GB * task.attempt }
		publishDir "$arch/results/QC/$id", mode: 'copy', overwrite: false
	
	input: 
	        tuple val(id), val(lane),file(R1), file(R2)
	val arch

	output:
			tuple val(), val(lane), file("*R1*fastqc.html"), file("*R2*fastqc.html")
	
	script:

	if (!"${arch}/results/QC/${id}/${lane}.R1_fastqc.html".exists() && !"${arch}/results/QC/${id}/${lane}.R2_fastqc.html".exists()) 
	
	"""
	fastqc $R1 $R2
	"""
	else 
	"""
	ln -s ${arch}/results/QC/${id}/${lane}.R1_fastqc.html
	ln -s ${arch}/results/QC/${id}/${lane}.R2_fastqc.html
	} 
	
	