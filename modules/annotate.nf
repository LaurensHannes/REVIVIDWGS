process annovar {
        tag "${analysis} analysis in mode:${mode} for family:${family}"
		cpus 16
			publishDir "./results/annotated/$analysis/$family/$mode", mode: 'copy', overwrite: true
			containerOptions "-B $programpath -B $humandbpath"
		time { 1.hour * task.attempt }
		errorStrategy 'ignore'
         maxRetries 9
		memory { 32.GB * task.attempt }

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
		cpus 18
		container "docker://ensemblorg/ensembl-vep:release_108.2"
			publishDir "./results/annotated/$analysis/$family/$mode", mode: 'copy', overwrite: true
			containerOptions "-B $VEPpath"
		time { 8.hour * task.attempt }
		errorStrategy 'ignore'
         maxRetries 9
		 memory { 96.GB * task.attempt }
		 
        input:
        tuple val(family), val(analysis), val(mode), file(vcfgz),file(vcfgztbi)
		val(VEPpath)
				
        output:
        tuple val(family), val(analysis), val(mode) ,file("${family}.${analysis}.${mode}.hg38_VEP.txt") 
		tuple val(family), val(analysis), val(mode) ,file("${family}.${analysis}.${mode}.hg38_VEP.txt_summary.html")
        """
		vep -i $vcfgz --cache --dir ${VEPpath} --fork ${task.cpus}  --check_frequency --freq_pop 1KG_EUR --freq_gt_lt gt --freq_freq 0.01 --freq_filter exclude -o ${family}.${analysis}.${mode}.hg38_VEP.txt --af_gnomadg --hgvs --show_ref_allele --symbol --pick --regulatory --no_intergenic --canonical --gene_phenotype --pubmed --plugin NMD --plugin pLI,${VEPpath}/pLI_values.txt --plugin REVEL,${VEPpath}/new_tabbed_revel_grch38.tsv.gz --tab
		vep -i $vcfgz --cache --dir ${VEPpath} --fork ${task.cpus}  --check_frequency --freq_pop 1KG_EUR --freq_gt_lt gt --freq_freq 0.01 --freq_filter exclude -o ${family}.${analysis}.${mode}.hg38_VEP.txt --af_gnomadg --hgvs --show_ref_allele --symbol --pick --regulatory --no_intergenic --canonical --gene_phenotype --pubmed --plugin NMD --plugin pLI,${VEPpath}/pLI_values.txt --plugin REVEL,${VEPpath}/new_tabbed_revel_grch38.tsv.gz --vcf --no_stats
		"""
}