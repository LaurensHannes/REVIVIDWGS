process parliament2 { 

        tag "$id"
        errorStrategy 'retry'
         maxErrors 3
		       container "docker://dnanexus/parliament2:latest"
			   runOptions =  "bin/bash"  -B "./:/home/dnanexus/in" -B "./:/home/dnanexus/out"
		memory { 28.GB * task.attempt }
		cpus { 4 * task.attempt }
			 time { 4.hour * task.attempt }
		


		input:

		tuple val(id), file(bam), file(bai)
		path genome
		path indexes

		"""
		python /home/dnanexus/parliament2.py --bam $bam --bai $bai -r $genome --fai $genome.fai
		"""
		
		}