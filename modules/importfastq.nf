process importfastq {

        tag "$id"

      container "docker://laurenshannes/gsutil"
	  
   
		cpus 4
        errorStrategy 'ignore'
		time { 2.hour * task.attempt }
         maxErrors 3
		disk { 50.GB * task.attempt }
		memory { 8.GB * task.attempt }
		storeDir "/vsc-hard-mounts/leuven-archive/arc_00086/FASTQ/$family/$id"
		
        input:
        tuple val(id),val(family) 
        path home
		path arch
		val download
		val bucket
	
        output:
         path "$id*.fastq.gz"
		
		script:
				
		if( download == "true")
        """
 
        gsutil cp -prn gs://${bucket}/**/$id*R[12]*.fastq.gz ./

        """
		
		else if( download == "false")
		
		"""
		ln -s /vsc-hard-mounts/leuven-archive/arc_00086/FASTQ/$family/$id/$id.*R[12].fastq.gz ./
		"""
}