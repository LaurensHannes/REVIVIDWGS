process AnnotSV{ 

        tag "$fam"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://laurenshannes/annotsv:3.1"

		memory { 4.GB * task.attempt }
		cpus 1
			 time { 30.minute * task.attempt }
			publishDir "./results/annotated/CNV/$fam" , mode: 'copy', overwrite: true



		input:

		tuple(val(fam), tuple(  file(snvvcf),file(cnvvcf)))

		output:
		
		tuple val(fam), file("${fam}.merge.vcf")

		"""
         ./AnnotSV/bin/AnnotSV -bedtools /usr/bin/bedtools -bcftools /usr/bin/bcftools -SVinputfile $cnvvcf -candidateSnvIndelFiles $snvvcf
		"""
		
		}