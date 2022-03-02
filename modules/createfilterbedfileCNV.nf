process createfilterbedfileCNV{ 

        tag "$id"

        errorStrategy 'retry'
         maxRetries 3


		memory { 4.GB * task.attempt }
		cpus 9
			 time { 15.minute * task.attempt }




		input:

		tuple val(id), file(bam),file(bai)
	

		output:
		
		tuple val(id), file("${id}.lowMQ.cov")

		"""
	samtools view -@ ${task.cpus} -H $bam > lowMQ.sam
	samtools view -@ ${task.cpus} $bam | awk '$5<5 {print $0}' >>  lowMQ.sam
	samtools view -@ ${task.cpus} -S -b -h lowMQ.sam > lowMQ.bam
	samtools depth lowMQ.bam >  ${id}.lowMQ.cov
		"""
		
		}