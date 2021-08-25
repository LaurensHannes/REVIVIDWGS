process importfastq {

        tag "$id"
      container "docker://laurenshannes/gsutil"
   
		cpus 4
        errorStrategy 'retry'
         maxErrors 3
		disk { 50.GB * task.attempt }
		memory { 8.GB * task.attempt }
		
        input:
        tuple val(id),val(family) 
        path home

        output:
         path "$home/FASTQ/$family/$id/*.fastq.gz"

        """
        [ ! -d "$home/FASTQ" ] && mkdir "$home/FASTQ"
        [ ! -d "$home/FASTQ/$family" ] && mkdir "$home/FASTQ/$family"
        [ ! -d "$home/FASTQ/$family/$id" ] && mkdir "$home/FASTQ/$family/$id"
        gsutil cp -prn gs://gcpi-rkjbr/GC085/$id/uploads/$id.200724.NovaSeq2.FCA.lane2.s_2_2.gcap_20_01.R1.fastq.gz $home/FASTQ/$family/$id/
        """
}