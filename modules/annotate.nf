process annotate {
        tag "${analysis} analysis in mode:${mode} for family:${family}"
		cpus 16
			publishDir "./results/annotated/$analysis/$family/$mode", mode: 'copy', overwrite: true
			containerOptions "-B $programpath -B $humandbpath"
		time { 1.hour * task.attempt }
		errorStrategy 'retry'
         maxRetries 9
		 
        input:
        tuple val(family), val(analysis), val(mode), file(vcfgz),file(vcfgztbi)
		val(programpath)
		val(humandbpath)
		val(annovardbs)
		
        output:
        tuple val(family), val(analysis, val(mode)), file("${family}.${analysis}.hg38_multianno.vcf"),file("${family}.${analysis}.hg38_multianno.txt") 

        """
        perl ${programpath}/annovar/table_annovar.pl $vcfgz ${humandbpath} -intronhgvs 25 --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.${analysis}.${mode} -remove -polish ${annovardbs} -nastring . -polish -vcfinput
        """
}