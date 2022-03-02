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

		tuple val(fam), file(vcf1), file(vcf2),file(vcf3),file(lowmq1),file(lowmq2),file(lowmq3)

		output:
		
		tuple val(fam), file("${fam}.merge.vcf")

		"""
		survivor bincov $lowmq1 10 2 > ${lowmq1}.bed
		survivor bincov $lowmq2 10 2 > ${lowmq2}.bed
		survivor bincov $lowmq3 10 2 > ${lowmq3}.bed
		survivor filter ${vcf1} ${lowmq1}.bed 50 -1 0.01 10 ${vcf1}.filterd.vcf
		survivor filter ${vcf2} ${lowmq2}.bed 50 -1 0.01 10 ${vcf2}.filterd.vcf
		survivor filter ${vcf3} ${lowmq3}.bed 50 -1 0.01 10 ${vcf3}.filterd.vcf
		ls ${vcf1}.filterd.vcf >> ${fam}
		ls ${vcf2}.filterd.vcf >> ${fam}
		ls ${vcf3}.filterd.vcf >> ${fam}
		survivor merge ${fam} 100 1 1 0 1 100 ${fam}.merge.vcf
		"""
		
		}