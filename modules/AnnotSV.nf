process AnnotSV{ 

        tag "$fam"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://laurenshannes/annotsv:3.1"

		memory { 4.GB * task.attempt }
		cpus 2
			 time { 30.minute * task.attempt }
			publishDir "./results/annotated/CNV/$fam" , mode: 'copy', overwrite: true



		input:

		tuple val(fam),  file(snvvcf),file(cnvvcf)

		output:
		
		tuple val(fam), file("**/${fam}.annotated.CNV.tsv")

		"""
         /AnnotSV/bin/AnnotSV -bedtools /usr/bin/bedtools -bcftools /usr/bin/bcftools -SVinputfile $cnvvcf -candidateSnvIndelFiles $snvvcf -outputFile ${fam}.annotated.CNV.tsv
		"""
		
		}