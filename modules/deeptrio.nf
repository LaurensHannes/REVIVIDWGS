process deeptrio { 

        tag "${bam1}"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://google/deepvariant:deeptrio-1.3.0"
			   containerOptions '-B `pwd`:/input:rw -B `pwd`:/output:rw -B `pwd`:/reference:rw'
		memory { 180.GB }
		cpus 36
			 time { 5.hour * task.attempt }
			publishDir "./results/deeptrio/$id/", mode: 'copy', overwrite: true



		input:

		tuple file(bam1), file(bai1),file(bam2),file(bai2),file(bam3),file(bai3)
		path genome
		

		output:
		
		file("${bam1}.vcf.gz")

		"""
		/opt/deepvariant/bin/deeptrio/run_deeptrio   --model_type WGS   --ref $genome   --reads_child $bam1   --reads_parent1 $bam2   --reads_parent2 $bam3   --output_vcf_child ${bam1}.vcf.gz   --output_vcf_parent1 ${bam2}.vcf.gz   --output_vcf_parent2 ${bam3}.vcf.gz  --num_shards ${task.cpus}  --intermediate_results_dir /output/intermediate_results_dir 
		"""
		
		}