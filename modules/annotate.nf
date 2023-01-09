process annovar {
        tag "${analysis} analysis in mode:${mode} for family:${family}"
		cpus 16
			publishDir "./results/annotated/$analysis/$family/$mode", mode: 'copy', overwrite: true
			containerOptions "-B $programpath -B $humandbpath"
		time { 1.hour * task.attempt }
		errorStrategy 'ignore'
         maxRetries 9
		 
        input:
        tuple val(family), val(analysis), val(mode), file(vcfgz),file(vcfgztbi)
		val(programpath)
		val(humandbpath)
		val(annovardbs)
		
        output:
        tuple val(family), val(analysis), val(mode), file("${family}.${analysis}.${mode}.hg38_multianno.vcf"),file("${family}.${analysis}.${mode}.hg38_multianno.txt") 

        """
        perl ${programpath}/annovar/table_annovar.pl $vcfgz ${humandbpath} -intronhgvs 25 --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.${analysis}.${mode} -remove -polish ${annovardbs} -nastring . -polish -vcfinput
        """
}

process VEP {
        tag "${analysis} analysis in mode:${mode} for family:${family}"
		cpus 16
		container "docker://ensemblorg/ensembl-vep:release_108.2"
			publishDir "./results/annotated/$analysis/$family/$mode", mode: 'copy', overwrite: true
			containerOptions "-B $VEPpath"
		time { 1.hour * task.attempt }
		errorStrategy 'ignore'
         maxRetries 9
		 
        input:
        tuple val(family), val(analysis), val(mode), file(vcfgz),file(vcfgztbi)
		val(VEPpath)
				
        output:
        tuple val(family), val(analysis), val(mode), file("${family}.${analysis}.${mode}.hg38_VEP.vcf"),file("${family}.${analysis}.${mode}.hg38_VEP.txt") 

        """
		vep -i $vcfgz --cache --dir ${VEPpath} --fork ${task.cpus} -o ${family}.${analysis}.${mode}.hg38_VEP.vcf --af_gnomadg --hgvs --symbol --canonical --check_frequency --freq_pop 1KG_EUR  --freq_freq 0.01 --freq_filter include --plugin NMD,pLI --vcf
		"""
}