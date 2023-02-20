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
			 time { 3.hour * task.attempt }
			storeDir "./results/deeptrio/"



		input:

		tuple val(fam), val(chr), file(bam1), file(bai1),file(bam2),file(bai2),file(bam3),file(bai3),val(index),val(father),val(mother)
		path genome
		path indexes
		

		output:
		
		tuple val(fam),val(index),file("${index}.${chr}.vcf.gz"),file("${index}.${chr}.vcf.gz.tbi")
		tuple val(fam),val(father),file("${father}.${chr}.vcf.gz"),file("${father}.${chr}.vcf.gz.tbi")
		tuple val(fam),val(mother),file("${mother}.${chr}.vcf.gz"),file("${mother}.${chr}.vcf.gz.tbi")
		tuple val(fam),val(index),file("${index}.${chr}.g.vcf.gz"),file("${index}.${chr}.g.vcf.gz.tbi")
		tuple val(fam),val(father),file("${father}.${chr}.g.vcf.gz"),file("${father}.${chr}.g.vcf.gz.tbi")
		tuple val(fam),val(mother),file("${mother}.${chr}.g.vcf.gz"),file("${mother}.${chr}.g.vcf.gz.tbi")
		tuple val("QC"),file("${index}.${chr}.visual_report.html"),file("${father}.${chr}.visual_report.html"),file("${mother}.${chr}.visual_report.html")

		"""
		/opt/deepvariant/bin/deeptrio/run_deeptrio  --regions $chr  --model_type WGS   --ref $genome   --reads_child ${index}.${chr}.bam   --reads_parent1 ${father}.${chr}.bam    --reads_parent2 ${mother}.${chr}.bam    --sample_name_child '${index}'   --sample_name_parent1 '${father}'   --sample_name_parent2 '${mother}' --output_vcf_child ${index}.${chr}.vcf.gz   --output_vcf_parent1 ${father}.${chr}.vcf.gz   --output_vcf_parent2 ${mother}.${chr}.vcf.gz --output_gvcf_child ${index}.${chr}.g.vcf.gz   --output_gvcf_parent1 ${father}.${chr}.g.vcf.gz   --output_gvcf_parent2 ${mother}.${chr}.g.vcf.gz --num_shards ${task.cpus}  --intermediate_results_dir /output/intermediate_results_dir 
		"""
		
		}
		
		process glnexusdt { 

        tag "family $fam"

		errorStrategy 'retry'
         maxRetries 3
		       container "docker://quay.io/mlin/glnexus:v1.3.1"
			   containerOptions '--cleanenv -H $PWD -B /usr/lib/locale/:/usr/lib/locale/,/usr/bin/parallel 		 -B `pwd`:/data  -B ${VSC_SCRATCH},${TMPDIR},${VSC_SCRATCH}/tmp:/tmp'
			memory 4.GB
		cpus 18
		executor 'PBS'
		clusterOptions '-A lp_revivid'
			 time { 30.minute * task.attempt }
			



		input:

		tuple val(fam), file(vcf1), file(vcf2), file(vcf3), file(vcftbi1), file(vcftbi2), file(vcftbi3),val(index),val(father),val(mother)

		output:
		
		tuple val(fam),file("${fam}.unprocessed.vcf.gz")

		"""
		
		/usr/local/bin/glnexus_cli --config DeepVariantWGS  /data/${index}.g.vcf.gz /data/${father}.g.vcf.gz /data/${mother}.g.vcf.gz > /data/${fam}.unprocessed.vcf.gz
		"""
		
		}
		

process glnexusprocessing {

	tag "$fam"
	cpus 4
	time { 30.minutes * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
	publishDir "./results/deeptrio/$fam/", mode: 'copy', overwrite: true

	input:
	tuple val(fam),file(vcfgz)
	
	output:
	tuple val(fam),val("deeptrio"),file("${fam}.vcf.gz"), file("${fam}.vcf.gz.tbi")
	
	
"""
	cat $vcfgz | bcftools view - | bgzip -c > ${fam}.vcf.gz
	tabix ${fam}.vcf.gz
"""

}