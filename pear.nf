process pear {

        tag "$lane"
        time { 1.hour * task.attempt }
		memory '2 GB'
		cpus 9
        errorStrategy 'retry'
        maxRetries 3

        input:
        tuple val(id), val(lane),file(R1), file(R2)
		path home
		
        output:
        tuple val(id), val(lane), file("${lane}.assembled.fastq"), file("${lane}.unassembled.forward.fastq"), file("${lane}.unassembled.reverse.fastq")

        """
		pear -j ${task.cpus} -f $R1 -r $R2 -o $lane
		"""
}