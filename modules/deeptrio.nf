process deeptrio { 

        tag "$chr from index $index"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://google/deepvariant:deeptrio-1.3.0"
			   containerOptions '--cleanenv -H $PWD -B /usr/lib/locale/:/usr/lib/locale/,/usr/bin/parallel 		 -B `pwd`:/input:rw -B `pwd`:/output:rw -B `pwd`:/reference:rw  -B ${VSC_SCRATCH},${TMPDIR},${VSC_SCRATCH}/tmp:/tmp'
			memory { 180.GB * task.attempt }
		cpus 36
		executor 'PBS'
		clusterOptions '-A lp_revivid'
			 time { 12.hour * task.attempt }
			publishDir "./results/deeptrio/$id/", mode: 'copy', overwrite: true



		input:

		tuple val(chr), file(bam1), file(bai1),file(bam2),file(bai2),file(bam3),file(bai3),val(index),val(father),val(mother)
		path genome
		path indexes
		

		output:
		
		file("${index}.vcf.gz")
		file("${index}.g.vcf.gz")

		"""
		/opt/deepvariant/bin/deeptrio/run_deeptrio  --regions $chr  --model_type WGS   --ref $genome   --reads_child ${index}.${chr}.bam   --reads_parent1 ${father}.${chr}.bam    --reads_parent2 ${mother}.${chr}.bam    --sample_name_child '${index}'   --sample_name_parent1 '${father}'   --sample_name_parent2 '${mother}' --output_vcf_child ${index}.vcf.gz   --output_vcf_parent1 ${father}.vcf.gz   --output_vcf_parent2 ${mother}.vcf.gz --output_gvcf_child ${index}.g.vcf.gz   --output_gvcf_parent1 ${father}.g.vcf.gz   --output_gvcf_parent2 ${mother}.g.vcf.gz --num_shards ${task.cpus}  --intermediate_results_dir /output/intermediate_results_dir 
		"""
		
		}