process parliament2 { 

        tag "$id"
		errorStrategy 'retry'
         maxRetries 3
		       container "docker://dnanexus/parliament2:v0.1.11-0-gb492db6d"
			   containerOptions '-B `pwd`:/home/dnanexus/in:rw -B `pwd`:/home/dnanexus/out:rw'
		memory { 180.GB * task.attempt }
		cpus 26
		clusterOptions '-A lp_revivid'
			 time { 8.hour * task.attempt }

			storeDir "./results/CNV/$id/parliament"


		input:

		tuple val(id), file(bam), file(bai)
		path genome
		path indexes

		output:
		
		tuple val(id), file("${id}.combined.genotyped.vcf")

		"""
		python /home/dnanexus/parliament2.py --bam $bam --bai $bai -r $genome --fai "$genome".fai --genotype --filter_short_contigs
		
		"""
		
		}

process indelible { 

        tag "$fam"
		errorStrategy 'retry'
         maxRetries 3
		       container "docker://mercury/indelible:1.1.3"
			   containerOptions '--pwd /usr/src/app/indelible'
		memory { 180.GB * task.attempt }
		cpus 16
		
			 time { 8.hour * task.attempt }

			storeDir "./results/CNV/$fam/indelbile"


		input:

		tuple val(fam), val(index), val(father), val(mother)
		path bam
		path genome
		path indexes 
		
		output:
		
		tuple val(id), file("${id}.bam.indelible.denovo.tsv")

		"""
		ls > input
		export index=\$(egrep '${index}' < input | egrep 'am\$')
		/usr/src/app/indelible/indelible.py complete --config /usr/src/app/indelible/config.hg38.yml --i  \$(echo \$index) --o ./ --r $genome --priors /usr/src/app/indelible/data/Indelible_db_10k.hg38.bed   --tb ${task.cpus}
		
		"""
		
		}