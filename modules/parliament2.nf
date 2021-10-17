process parliament2 { 

        tag "$id"
		label 'parliament'
        errorStrategy 'ignore'
         maxRetries 3
		       container "docker://sameerdcosta/parliament2:latest"
			   containerOptions '-B `pwd`:/home/dnanexus/in:rw -B `pwd`:/home/dnanexus/out:rw'
		memory { 64.GB * task.attempt }
		cpus { 36 * task.attempt }
			 time { 4.hour * task.attempt }
		


		input:

		tuple val(id), file(bam), file(bai)
		path genome
		path indexes

		"""
		python /home/dnanexus/parliament2.py --bam $bam --bai $bai -r $genome --fai "$genome".fai --genotype --svviz --filter_short_contigs
		"""
		
		}