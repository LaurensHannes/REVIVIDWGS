process mergeCNV{ 

        tag "$fam"

        errorStrategy 'retry'
         maxRetries 3
		       container "docker://sameerdcosta/parliament2:latest"

		memory { 4.GB * task.attempt }
		cpus 1
			 time { 5.minute * task.attempt }
			publishDir "./results/CNV/$fam(", mode: 'copy', overwrite: true



		input:

		tuple val(fam), file(vcf1), file(vcf2), file(vcf3)

		output:
		
		tuple val(fam), file("${fam}.merge.vcf")

		"""
		ls ${vcf1} >> ${fam}
		ls ${vcf2} >> ${fam}
		ls ${vcf3} >> ${fam}
		survivor merge ${fam} 100 1 1 0 1 100 ${fam}.merge.vcf
		"""
		
		}