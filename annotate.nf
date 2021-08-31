process annotate {
        tag "$family"
		cpus 16
			publishDir "./results/annotated/$analysis", mode: 'copy', overwrite: false
// change input names for gz and gztbi
        input:
        tuple val(family), val(analysis), file(vcfgz),file(vcfgztbi)
		val(programpath)
		val(humandbpath)
        output:
        tuple val(family), val(analysis), file("${family}.${analysis}.hg38_multianno.vcf"),file("${family}.${analysis}.hg38_multianno.txt") 

        """
        perl ${programpath}/annovar/table_annovar.pl $vcfgz ${humandbpath} -intronhgvs 25 --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.${analysis} -remove -polish -protocol refgene,avsnp150,gnomad30_genome,clinvar_20210123,regsnpintron,dbnsfp41a -operation g,f,f,f,f,f -nastring . -polish -intronhgvs 50 -vcfinput
        """
}