process checkfamilyvcf {


	tag "$id"
	cpus 1
	executor 'local'
     runOptions = "-B $launchDir"


	
	input:
    tuple val(id),val(family) 
	
	output:
	tuple val(id), val(family), env(status), emit: vcfcheck_ch
	
	"""
	if [ -f $launchDir/results/vcfs/"$family".vcf ]
    then 
	status="done"
	else
	status="todo"
	fi
	"""
	
	}