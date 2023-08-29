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
		       container "docker://laurenshannes/indelible:1.1.3"
			   containerOptions '-B `pwd`:/usr/src/app/indelible:rw'
		memory { 64.GB * task.attempt }
		cpus 4
		clusterOptions '--partition=batch_long'

			 time 7.day

			storeDir "./results/CNV/$fam/indelbile"


		input:

		tuple val(fam), val(index), val(father), val(mother)
		path bam
		path genome
		path indexes 
		
		output:
		
		file("*.tsv")

		"""
		ls > input
		export index=\$(egrep '${index}' < input | egrep 'am\$')
		export father=\$(egrep '${father}' < input | egrep 'am\$')
		export mother=\$(egrep '${mother}' < input | egrep 'am\$')
		ln -s /lustre1/project/stg_00086/resources/programs/indelible/data
		/lustre1/project/stg_00086/resources/programs/indelible/indelible.py complete --config /lustre1/project/stg_00086/resources/programs/indelible/config.hg38.yml --i  \$(echo \$index) --o ./ --r $genome --priors /lustre1/project/stg_00086/resources/programs/indelible/data/Indelible_db_10k.hg38.bed   --tb ${task.cpus} --p \$(echo \$father) --m \$(echo \$mother)
		
		"""
		
		}