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
		if( exome == "true") 
		
		"""
		for line in ${lsbucket} ; do
		grep $id \$line > temp.txt
		done
		for line in temp.txt ; do
		gsutil cp -prn \$line ./
		done
		"""
		

		if( params.download == "true" && params.exome == "false") 

        """

        gsutil cp -prn gs://${bucket}/**/$id*R[12]*.fastq.gz ./

        """

		else if( params.download == "false") 
		
		"""
		ln -s /vsc-hard-mounts/leuven-archive/arc_00086/FASTQ/$family/$id/$id.*R[12].fastq.gz ./
		"""


}

process importexomefastq {

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

		
		"""
		for line in $lsbucket ; do
		grep $id \$line > temp.txt
		done
		for fastq in $(cat temp.txt) ; do
		gsutil cp -prn \$fastq ./
		done
		"""


}