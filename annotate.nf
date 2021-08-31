process annotate {
        tag "$family"
		cpus 16
			publishDir "./results/annotated/$analysis", mode: 'copy', overwrite: false
// change input names for gz and gztbi
        input:
        tuple val(family), val(analysis), file(vcfgz),file(vcfgztbi)
        output:
        tuple val(family), val(analysis), file("${family}.${analysis}.hg38_multianno.vcf"),file("${family}.${analysis}.hg38_multianno.txt") 

        """
        perl /staging/leuven/stg_00086/resources/programs/annovar/table_annovar.pl $vcfgz /staging/leuven/stg_00086/resources/humandb/ -intronhgvs 25 --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.${analysis} -remove -polish -protocol refgene,avsnp150,gnomad30_genome,clinvar_20210123,regsnpintron,dbnsfp41a -operation g,f,f,f,f,f -nastring . -polish -intronhgvs 50 -vcfinput
        """
}