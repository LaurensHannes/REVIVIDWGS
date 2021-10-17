process annotate {
        tag "$family"
		cpus 16
			publishDir "./results/annotated/$analysis", mode: 'copy', overwrite: false
			containerOptions "-B $programpath $humandbpath"

        input:
        tuple val(family), val(analysis), file(vcfgz),file(vcfgztbi)
		val(programpath)
		val(humandbpath)
		val(annovardbs)
		
        output:
        tuple val(family), val(analysis), file("${family}.${analysis}.hg38_multianno.vcf"),file("${family}.${analysis}.hg38_multianno.txt") 

        """
        perl ${programpath}/annovar/table_annovar.pl $vcfgz ${humandbpath} -intronhgvs 25 --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.${analysis} -remove -polish ${annovardbs} -nastring . -polish -vcfinput
        """
}