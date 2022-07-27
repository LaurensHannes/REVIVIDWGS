process deeptrio { 

        tag "$chr from index $index"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://google/deepvariant:deeptrio-1.3.0"
			   containerOptions '--cleanenv -H $PWD -B /usr/lib/locale/:/usr/lib/locale/,/usr/bin/parallel 		 -B `pwd`:/input:rw -B `pwd`:/output:rw -B `pwd`:/reference:rw  -B ${VSC_SCRATCH},${TMPDIR},${VSC_SCRATCH}/tmp:/tmp'
			memory 180.GB
		cpus 36
		executor 'PBS'
		clusterOptions '-A lp_revivid'
			 time { 12.hour * task.attempt }
			publishDir "./results/deeptrio/", mode: 'copy', overwrite: true



		input:

		tuple val(chr), file(bam1), file(bai1),file(bam2),file(bai2),file(bam3),file(bai3),val(index),val(father),val(mother)
		path genome
		path indexes
		

		output:
		
		tuple val("${index}"),file("${index}.${chr}.vcf.gz"),file("${index}.${chr}.vcf.gz.tbi")
		tuple val("${father}"),file("${father}.${chr}.vcf.gz"),file("${father}.${chr}.vcf.gz.tbi")
		tuple val("${mother}"),file("${mother}.${chr}.vcf.gz"),file("${mother}.${chr}.vcf.gz.tbi")
		tuple val("${index}"),file("${index}.${chr}.g.vcf.gz"),file("${index}.${chr}.g.vcf.gz.tbi")
		tuple val("${father}"),file("${father}.${chr}.g.vcf.gz"),file("${father}.${chr}.g.vcf.gz.tbi")
		tuple val("${mother}"),file("${mother}.${chr}.g.vcf.gz"),file("${mother}.${chr}.g.vcf.gz.tbi")
		tuple val("QC"),file("${index}.${chr}.visual_report.html"),file("${father}.${chr}.visual_report.html"),file("${mother}.${chr}.visual_report.html")

		"""
		/opt/deepvariant/bin/deeptrio/run_deeptrio  --regions $chr  --model_type WGS   --ref $genome   --reads_child ${index}.${chr}.bam   --reads_parent1 ${father}.${chr}.bam    --reads_parent2 ${mother}.${chr}.bam    --sample_name_child '${index}_deeptrio'   --sample_name_parent1 '${father}_deeptrio'   --sample_name_parent2 '${mother}_deeptrio' --output_vcf_child ${index}.${chr}.vcf.gz   --output_vcf_parent1 ${father}.${chr}.vcf.gz   --output_vcf_parent2 ${mother}.${chr}.vcf.gz --output_gvcf_child ${index}.${chr}.g.vcf.gz   --output_gvcf_parent1 ${father}.${chr}.g.vcf.gz   --output_gvcf_parent2 ${mother}.${chr}.g.vcf.gz --num_shards ${task.cpus}  --intermediate_results_dir /output/intermediate_results_dir 
		"""
		
		}
		
		process glnexusdt { 

        tag "$chr for family $family"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://quay.io/mlin/glnexus:v1.3.1"
			   containerOptions '--cleanenv -H $PWD -B /usr/lib/locale/:/usr/lib/locale/,/usr/bin/parallel 		 -B `pwd`:/input:rw -B `pwd`:/output:rw -B `pwd`:/reference:rw  -B ${VSC_SCRATCH},${TMPDIR},${VSC_SCRATCH}/tmp:/tmp'
			memory 8.GB
		cpus 1
		executor 'PBS'
		clusterOptions '-A lp_revivid'
			 time { 30.minute * task.attempt }
			publishDir "./results/deeptrio/$id/", mode: 'copy', overwrite: true



		input:

		tuple val(chr), file(bam1), file(bai1),file(bam2),file(bai2),file(bam3),file(bai3),val(index),val(father),val(mother)
		path genome
		path indexes
		

		output:
		
		tuple val("${index}"),file("${index}.${chr}.vcf.gz"),file("${index}.${chr}.vcf.gz.tbi")
		tuple val("${father}"),file("${father}.${chr}.vcf.gz"),file("${father}.${chr}.vcf.gz.tbi")
		tuple val("${mother}"),file("${mother}.${chr}.vcf.gz"),file("${mother}.${chr}.vcf.gz.tbi")
		tuple val("${index}"),file("${index}.${chr}.g.vcf.gz")
		tuple val("${father}"),file("${father}.${chr}.g.vcf.gz")
		tuple val("${mother}"),file("${mother}.${chr}.g.vcf.gz")
		tuple val("QC"),file("${index}.${chr}.visual_report.html"),file("${father}.${chr}.visual_report.html"),file("${mother}.${chr}.visual_report.html")

		"""
		/opt/deepvariant/bin/deeptrio/run_deeptrio  --regions $chr  --model_type WGS   --ref $genome   --reads_child ${index}.${chr}.bam   --reads_parent1 ${father}.${chr}.bam    --reads_parent2 ${mother}.${chr}.bam    --sample_name_child '${index}_deeptrio'   --sample_name_parent1 '${father}_deeptrio'   --sample_name_parent2 '${mother}_deeptrio' --output_vcf_child ${index}.${chr}.vcf.gz   --output_vcf_parent1 ${father}.${chr}.vcf.gz   --output_vcf_parent2 ${mother}.${chr}.vcf.gz --output_gvcf_child ${index}.${chr}.g.vcf.gz   --output_gvcf_parent1 ${father}.${chr}.g.vcf.gz   --output_gvcf_parent2 ${mother}.${chr}.g.vcf.gz --num_shards ${task.cpus}  --intermediate_results_dir /output/intermediate_results_dir 
		"""
		
		}