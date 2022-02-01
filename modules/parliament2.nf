process parliament2 { 

        tag "$id"
		label 'parliament'
        errorStrategy 'retry'
         maxRetries 3
		       container "docker://sameerdcosta/parliament2:latest"
			   containerOptions '-B `pwd`:/home/dnanexus/in:rw -B `pwd`:/home/dnanexus/out:rw'
		memory { 64.GB * task.attempt }
		cpus 36
			 time { 16.hour * task.attempt }
			publishDir "./results/CNV/$id/parliament", mode: 'copy', overwrite: true


		input:

		tuple val(id), file(bam), file(bai)
		path genome
		path indexes

		output:
		
		tuple val(id), file("${id}.survivor_sorted.vcf")

		"""
		python /home/dnanexus/parliament2.py --bam $bam --bai $bai -r $genome --fai "$genome".fai --genotype --filter_short_contigs
		
		"""
		
		}