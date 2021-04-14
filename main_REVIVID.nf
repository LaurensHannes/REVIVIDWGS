#! usr/bin/env nextflow
// script parameters

workDir = '/staging/leuven/stg_00086/Laurens/FNRCP/work'
params.genome ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.fasta'
params.genomefai ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai'
params.genomedict ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.dict'
params.indexes = '/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Homo_sapiens_assembly38.fasta.*'
params.ped ='/staging/leuven/stg_00086/Laurens/REVIVID/FNRCP.ped'
params.mask ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/mask.bed'
params.home ='/staging/leuven/stg_00086/Laurens/FNRCP'
params.maskrepeats ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/repeats.bed'
params.fastqgz = '/staging/leuven/stg_00086/Laurens/FNRCP/raw_data/FASTQGZ/*/*/*.fastq.gz'
params.snps ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
params.snpsindex  ='/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
params.indels = '/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
params.indelsindex = '/staging/leuven/stg_00086/resources/reference_genomes/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
//fastqgz_ch = channel.fromPath(params.fastqgz)

//script

myFile = file(params.ped)
myReader = myFile.newReader()
String line
familymap = [:]
ids = []
while( line = myReader.readLine() ) {
(empty, family, id, father, mother, sex, phenotype) = (line =~ /(^F\d{1,2})\t(GC\d+)\t(\w+)\t(\w+)\t(\d+)\t(\d+)/)[0]
        familymap[id]=family
        ids << id
}
myReader.close()

Channel.fromList(ids).map { it -> [it, familymap[it]] }into{ idfamily_ch1; idfamily_ch2; idfamily_ch3; idfamily_ch4; idfamily_ch5 }

//processes





process importfastq {
        tag "$id"
      container "docker://laurenshannes/gsutil"
   
		cpus 4
        errorStrategy 'retry'
         maxErrors 3
		disk { 50.GB * task.attempt }
		memory { 16.GB * task.attempt }
        input:
        tuple val(id),val(family) from idfamily_ch1
        path home from params.home

        output:
         path "$home/FASTQ/$family/$id/*.fastq.gz" into fastqgz_ch

        """

        [ ! -d "$home/FASTQ" ] && mkdir "$home/FASTQ"
		[ ! -d "$home/tempstorage" ] && mkdir "$home/tempstorage"
        [ ! -d "$home/FASTQ/$family" ] && mkdir "$home/FASTQ/$family"
		[ ! -d "$home/tempstorage/$id" ] && mkdir "$home/tempstorage/$id"
		[ ! -d "$home/tempstorage/$family" ] && mkdir "$home/tempstorage/$family"
        [ ! -d "$home/FASTQ/$family/$id" ] && mkdir "$home/FASTQ/$family/$id"
        gsutil cp -prn gs://gcpi-rkjbr/GC085/$id/uploads/$id.*.fastq.gz $home/FASTQ/$family/$id/
        """
}

fastqgz_ch.flatten().filter(~/.*R\d+.fastq.gz/).map{file -> tuple(file.getBaseName(3), file)}.groupTuple().flatten().collate( 3 ).map{lane,R1,R2 -> tuple(R1.simpleName,lane,R1,R2)}.set{gzipped_ch}
gzipped_ch.into{temp_ch1;temp_ch2}
//temp_ch2.view()

process pear {

        tag "$lane"
        time '1h'
		memory '2 GB'
		cpus 9
        errorStrategy 'retry'
        maxRetries 3
		storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"

        input:
        tuple val(id), val(lane),file(R1), file(R2) optional true from temp_ch1
		path home from params.home
		
        output:
        tuple val(id), val(lane), file("${lane}.assembled.fastq"), file("${lane}.unassembled.forward.fastq"), file("${lane}.unassembled.reverse.fastq") into paired_ch

        """
		if [ -f $home/tempstorage/${id}/${lane}*.fastq ] || [ -f $home/tempstorage/${id}/${lane}.indexed.bam ] || [ -f $home/tempstorage/${id}/${id}.bam ] || [ -f $home/tempstorage/${id}/${id}.recallibrated.bam ] || [ -f $home/tempstorage/${id}/${id}.vcf ] || [ -f $home/tempstorage/${id}/${id.filtered.vcf ] || [ -f $home/tempstorage/${id}/${id}.filtered.vcf.gz ]
		then 
		rm ${lane}.assembled.fastq
		rm ${lane}.unassembled.forward.fastq
		rm ${lane}.unassembled.reverse.fastq
		else
		pear -j ${task.cpus} -f $R1 -r $R2 -o $lane
        fi
		"""
}

 process alignment {

		
		memory { 8.GB * task.attempt }
		tag "$lane"
			 time { 1.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3


		storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
        indexes_ch = Channel.fromPath(params.indexes)


        input:
        tuple val(id),val(lane), file(assembled), file(forward), file(reverse) from paired_ch
        path genome from params.genome
        path indexes from indexes_ch.toList()
        path home from params.home

        output:
        tuple val(id), val(lane), file("${lane}.indexed.bam") optional true into mapped_ch

        """
		if [ -f $home/tempstorage/${id}/${lane}.indexed.bam ] || [ -f $home/tempstorage/${id}/${lane}.RG.bam ] || [ -f $home/tempstorage/${id}/${lane}.dups.bam ] || [ -f $home/tempstorage/${id}/${id}.bam ] || [ -f $home/tempstorage/${id}/${id}.recallibrated.bam ] || [ -f $home/tempstorage/${id}/${id}.vcf ] || [ -f $home/tempstorage/${id}/${id.filtered.vcf ]  || [ -f $home/tempstorage/${id}/${id}.filtered.vcf.gz ]
		then
		rm ${lane}.indexed.bam
		else
        bwa mem -t ${task.cpus} $genome $assembled | samtools view -@ ${task.cpus} -bS > ${lane}.assembled.bam
        bwa mem -t ${task.cpus} $genome $forward | samtools view -@ ${task.cpus} -bS > ${lane}.forward.bam
        bwa mem -t ${task.cpus} $genome $reverse | samtools view -@ ${task.cpus} -bS > ${lane}.reverse.bam
        samtools merge -@ ${task.cpus} ${lane}.indexed.unsorted.bam  ${lane}.assembled.bam ${lane}.forward.bam ${lane}.reverse.bam
        samtools sort -@ ${task.cpus} -o ${lane}.indexed.bam ${lane}.indexed.unsorted.bam
		rm  ${lane}.assembled.bam  ${lane}.forward.bam ${lane}.reverse.bam  ${lane}.indexed.unsorted.bam
		rm /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${lane}*.fastq
		fi
        """
}



process readgroups {

	tag "$lane"
	cpus 2
	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"

	input:
	tuple val(id), val(lane), file(bam) from mapped_ch
	
	output:
	tuple val(id),val(lane), file("${lane}.RG.bam"),file("${lane}.RG.bam.bai") optional true into mapped_RG_ch	
		


	"""
	gatk AddOrReplaceReadGroups -I $bam -O ${lane}.RG.bam -LB REVIVID -PL ILLUMINA -PU $lane -SM $id 
	samtools index -@ ${task.cpus} ${lane}.RG.bam
	rm $bam
	"""

}

process duplicates { 

        tag "$lane"
		storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
        errorStrategy 'retry'
         maxErrors 3
		memory { 2.GB * task.attempt }
		cpus 2

		
	input:
	 tuple val(id),val(lane),file(bam),file(bai) from mapped_RG_ch

	output:
	tuple val(id),file("${lane}.dups.bam") optional true into dups_ch
	tuple val(id),file("${lane}.metrics.txt") optional true into dup_metrics_ch

	
	"""
	gatk MarkDuplicates -I $bam -O ${lane}.dups.bam -M ${lane}.metrics.txt
	rm $bam
	rm $bai
	"""
}

dups_ch.groupTuple().set{mappedgrouped_ch}

process mergebams {

	tag "$id"
    cpus 6
	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
			 time { 2.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
	
	input:
	tuple val(id),file(bams) from mappedgrouped_ch
	path home from params.home

	
	output:
	tuple val(id),file("${id}.bam"),file("${id}.bam.bai") optional true into mergedbam_ch

	"""
	if [ -f $home/tempstorage/${id}/${id}.bam ] || [ -f $home/tempstorage/${id}/${id}.recallibrated.bam ]  || [ -f $home/tempstorage/${id}/${id}.vcf ] || [ -f $home/tempstorage/${id}/${id.filtered.vcf ]  || [ -f $home/tempstorage/${id}/${id}.filtered.vcf.gz ]
		then
		echo "done" > ${id}.bam
	    echo "done" > ${id}.bam.bai
		else
	samtools merge -@ ${task.cpus} ${id}.bam $bams
	samtools index -@ ${task.cpus} ${id}.bam
	fi
	"""

}


mergedbam_ch.into{mergedbam1_ch;mergedbam2_ch}


process generateCRAM {

	tag "$id"
		storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
		
		input:
		tuple val(id),file(bam),file(bai) from mergedbam1_ch
		path genome from params.genome
        path faidx from params.genomefai
				
		output:
		tuple val(id),file("${id}.cram"),file("${id}.cram.crai") optional true into mergedcram_ch
		tuple val(id),file(bam),file(bai) into mergedbam3_ch
		
		"""
		samtools view -@ ${task.cpus} -C -o ${id}.cram -T $genome $bam 
		samtools index -@ ${task.cpus} -c ${id}.cram
		"""
		
}



//merged_dups_ch.into{mergedbam1_ch;mergedbam2_ch}


process baserecalibrator {

	tag "$id"
	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
	cpus 4	
	time { 5.hour * task.attempt }
	errorStrategy 'retry' 
	maxRetries 3
	container "docker://broadinstitute/gatk"
	memory { 32.GB * task.attempt }

        input:
        tuple val(id), file(merged), file(bai) from mergedbam2_ch
        path genome from params.genome
       path faidx from params.genomefai
        path dict from params.genomedict
	 path snps from params.snps
        path snpsindex from params.snpsindex

        output:
        tuple val(id), file(merged), file(bai), file("${id}.recal_data.table") into recal_data_ch

        """
        gatk BaseRecalibrator -I $merged -R $genome -O ${id}.recal_data.table --known-sites $snps --verbosity WARNING
        """
}

process applyBQSR {

        tag "$id"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"
		 time { 5.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		cpus 4


        input:
        //tuple val(id), file(bam), file(bai) from mergedbam3_ch
		tuple val(id), file(bam), file(bai), file(table) from recal_data_ch
        path genome from params.genome
     //   path table from recal_data_ch
	path faidx from params.genomefai
        path dict from params.genomedict


        output:
        tuple val(id),file("${id}.recallibrated.bam"),file("${id}.recallibrated.bam.bai") into BQSR_applied_ch

        """
        gatk ApplyBQSR -R $genome -I $bam -bqsr-recal-file $table -O ${id}.recallibrated.bam
		samtools index -@ ${task.cpus} ${id}.recallibrated.bam
		rm /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${id}.bam /staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}/${id}.bam.bai
		
        """

}



		


process genotype {

        tag "$id"
		 time { 10.hour * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://broadinstitute/gatk"
	memory { 8.GB * task.attempt }
	       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"


        input:
        tuple val(id), file(bam),file(bai) from BQSR_applied_ch
        path genome from params.genome
        path dict from params.genomedict
//        path index from params.indexes
        path faidx from params.genomefai
	path mask from params.maskrepeats

        output:
        tuple val(id), file("${id}.vcf") into vcf_uncallibrated_ch

        """
        gatk HaplotypeCaller --verbosity INFO -XL $mask -R $genome -I $bam -O ${id}.vcf --sequence-dictionary ${dict} --native-pair-hmm-threads ${task.cpus}
        """

}

process variantrecalibration {

	tag "$id"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"

	container "docker://broadinstitute/gatk"

	input:
	tuple val(id), path(vcf) from vcf_uncallibrated_ch
	path genome from params.genome
        path dict from params.genomedict
//      path index from params.indexes
        path faidx from params.genomefai
	path snps from params.snps
	path snpsindex from params.snpsindex
	path indels from params.indels
	path indelsindex from params.indelsindex
	path mask from params.maskrepeats


	output:
	tuple val(id), file("${id}.filtered.vcf")  into individual_vcf_ch

	"""
	gatk CNNScoreVariants -V $vcf -R $genome -O ${id}.pretranched.vcf
	gatk FilterVariantTranches -V ${id}.pretranched.vcf --resource $snps --resource $indels -O ${id}.filtered.vcf --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 
	"""
}

process compressandindex {

	tag "$id"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/${id}"


input:
	tuple val(id), file(vcf) from individual_vcf_ch 

output:
        tuple val(id), file("${id}.filtered.vcf.gz"), file("${id}.filtered.vcf.gz.tbi") into individual_vcfgz_ch

"""
        bcftools view --threads ${task.cpus} -o ${id}.filtered.vcf.gz -O z ${id}.filtered.vcf
        tabix ${id}.filtered.vcf.gz
"""
}

idfamily_ch2.join(individual_vcfgz_ch).map{ id, family, vcf, index -> tuple(family,vcf,index)}.groupTuple().set{individual_vcf_for_merge_ch}
individual_vcf_for_merge_ch.into{individual_vcf_for_merge_ch1;individual_vcf_for_merge_ch2}
individual_vcf_for_merge_ch2.view()


process mergevcf {

	tag "$family"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/$family"



	input:
	tuple val(family), path(vcf), path(index) from individual_vcf_for_merge_ch1

	output:
	tuple val(family), path("${family}.vcf.gz"), path("${family}.vcf.gz.tbi") into merged_vcf_ch
	
"""
	bcftools merge -o ${family}.vcf.gz -O z --threads ${task.cpus} $vcf 
	tabix ${family}.vcf.gz
"""

}

//process gatkfilter {
//
//      tag "$family"
//
//        storeDir "/mnt/hdd2/data/lau/phd/FNRCP/tempstorage/${family}"
//	container "docker://broadinstitute/gatk"
//		
//      input:
//    tuple val(family), file(vcfgz), file(vcfgztbi) from merged_vcf_ch

//
//      path genome from params.genome
//    path index from params.genomefai
//        path dict from params.genomedict
//
//      output:
//        tuple val(family), file("${family}.gatk.filtered.vcf.gz"), file("${family}.gatk.filtered.vcf.gz.tbi") into combinedfilteredvcf_ch
//
//
//      """
//  gatk VariantFiltration -R $genome -V $vcfgz -O ${family}.gatk.filtered.vcf.gz --filter-name "qualfilter" --filter-expression "QUAL < 50" --genotype-filter-name "gqfilter" --genotype-filter-expression "GQ < 30"
//    """
//
//}

//combinedfilteredvcf_ch.into{vcftodenovo;vcftorecessive}
merged_vcf_ch.into{vcftodenovo;vcftorecessive}
process SelectVariantsdenovo {

        tag "$family"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/$family"

        analysis_ch = channel.value("denovo")

	
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi) from vcftodenovo
        path genome from params.genome
        path index from params.genomefai
        path dict from params.genomedict
        path ped from params.ped
        path mask from params.mask
        val y from analysis_ch

        output:
        tuple val(family), val(y), file("${family}.denovo.vcf.gz"), file("${family}.denovo.vcf.gz.tbi") into denovovcf_ch

        """
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask --exclude-filtered false --mendelian-violation true --mendelian-violation-qual-threshold 30 --remove-unused-alternates false -O ${family}.denovo.vcf.gz
        """
//exclude-filtered changed to false --remove-unused-alternates to false  --restrict-alleles-to BIALLELIC (removed)
}

process SelectVariantsAR {


        tag "$family"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/$family"

        analysis_ch = channel.value("AR")

		
        input:
        tuple val(family), file(vcfgz), file(vcfgztbi) from vcftorecessive
        path genome from params.genome
        path index from params.genomefai
        path dict from params.genomedict
        path ped from params.ped
        path mask from params.mask
        val z from analysis_ch

        output:
        tuple val(family), val(z), file("${family}.recessive.vcf.gz"), file("${family}.recessive.vcf.gz.tbi") into recessivevcf_ch

        """
        gatk SelectVariants -V $vcfgz -ped $ped -R $genome -XL $mask --mendelian-violation true --exclude-filtered true --invert-mendelian-violation true --mendelian-violation-qual-threshold 30  --remove-unused-alternates true  --restrict-alleles-to BIALLELIC -O ${family}.recessive.vcf.gz
        """
}

//process gatkannotate {
//
//        tag "$family"
//		container "docker://broadinstitute/gatk"
//		
  //      storeDir "/mnt/hdd2/data/lau/phd/FNRCP/tempstorage/${family}"
//
//
  //      input:
    //    tuple val(family), val(y), file(denovovcfgz), file(denovovcfgztbi) from denovovcf_ch
//        path genome from params.genome
//        path index from params.genomefai
//        path dict from params.genomedict
//        path ped from params.ped
//
//        output:
//        tuple val(family), val(y), file("${family}.denovo.gatk.vcf.gz"), file ("${family}.denovo.gatk.vcf.gz.tbi") into denovogatkvcf_ch
//
//        """
//        gatk VariantAnnotator -A TandemRepeat -A PossibleDeNovo -R $genome -V $denovovcfgz -O ${family}.denovo.temp.gatk.vcf.gz -ped $ped
//        gatk SelectVariants -V  ${family}.denovo.temp.gatk.vcf.gz -O  ${family}.denovo.gatk.vcf.gz -select 'vc.hasAttribute("hiConfDeNovo")'
//        """
//
//}

denovovcf_ch.join(recessivevcf_ch).set{joined_ch}

process annotate {
        tag "$family"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/$family"
		cpus 36
		
// change input names for gz and gztbi
        input:
        tuple val(family), val(denovo), file(denovogatkfvcfgz),file(denovovcfgatkfgztbi),val(AR), file(ARgatkfvcfgz),file(ARvcfgatkfgztbi) from joined_ch
        output:
        tuple val(family), val(denovo), file("${family}.denovo.hg38_multianno.vcf"),file("${family}.denovo.hg38_multianno.txt") into annotated_denovo_ch
        tuple val(family), val(AR), file("${family}.AR.hg38_multianno.vcf") ,file("${family}.AR.hg38_multianno.txt") into annotated_AR_ch
        """
        perl /staging/leuven/stg_00086/resources/programs/annovar/table_annovar.pl $denovogatkfvcfgz /staging/leuven/stg_00086/resources/humandb/ -intronhgvs --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.denovo -remove -polish -protocol refgene,avsnp150,gnomad30_genome,clinvar_20210123,regsnpintron,dbnsfp41a -operation g,f,f,f,f,f -nastring . -polish -intronhgvs 50 -vcfinput
        perl /staging/leuven/stg_00086/resources/programs/annovar/table_annovar.pl $ARgatkfvcfgz /staging/leuven/stg_00086/resources/humandb/ -intronhgvs --maxgenethread ${task.cpus} -thread ${task.cpus} -buildver hg38 -out ${family}.AR -remove -polish -protocol refgene,avsnp150,gnomad30_genome,clinvar_20210123,regsnpintron,dbnsfp41a -operation g,f,f,f,f,f -nastring . -polish -intronhgvs 50 -vcfinput
        """
}

process subset {

        tag "$family"
       	storeDir "/staging/leuven/stg_00086/Laurens/FNRCP/tempstorage/$family"


        input:
        tuple val(family), val(AR), file(ARannotatedVCF), file(ARannotatedTXT) from annotated_AR_ch

        """
        cat $ARannotatedTXT | awk '\$12 < 0.00002 || NR==1' >  ${family}.${AR}.subset.txt
        """

}