process importfastq {

        tag "$id"

      container "docker://laurenshannes/gsutil"
	  
   
		cpus 4
        errorStrategy 'ignore'
		time { 4.hour * task.attempt }
         maxErrors 3
		disk { 50.GB * task.attempt }
		memory { 8.GB * task.attempt }
		storeDir "${arch}/FASTQ/$family/$id"
		maxForks 25
		
        input:
        tuple val(id),val(family) 
        path home
		val arch
		val download
		val bucket
		val exome
		path lsbucket


        output:
         path "$id*.fastq.gz"
		
		script:
		if( download == "true" && exome == "true") 
		
		"""
		for line in $lsbucket ; do
		gsutil cp -prn $(grep "$id" \$line) ./
		done
		"""
		



		if( download == "true" && exome == "false") 

        """

        gsutil cp -prn gs://${bucket}/**/$id*R[12]*.fastq.gz ./

        """

		else if( download == "false") 
		
		"""
		ln -s /vsc-hard-mounts/leuven-archive/arc_00086/FASTQ/$family/$id/$id.*R[12].fastq.gz ./
		"""


}