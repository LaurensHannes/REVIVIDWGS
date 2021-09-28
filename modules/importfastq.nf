process importfastq {

        tag "$id"

      container "docker://laurenshannes/gsutil"
   
		cpus 4
        errorStrategy 'retry'
         maxErrors 3
		disk { 50.GB * task.attempt }
		memory { 8.GB * task.attempt }
		publishDir "./FASTQ/$family/$id", mode: 'copy', overwrite: false
		
        input:
        tuple val(id),val(family) 
        path home

        output:
         path "$id*.fastq.gz"

        """
		gsutil cp -prn gs://gcpi-rkjbr/**/$id.200724.NovaSeq2.FCA.lane1.s_1_[12].gcap_20_01.R[12].fastq.gz ./
        """
}