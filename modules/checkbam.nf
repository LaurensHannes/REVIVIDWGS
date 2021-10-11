process checkbam {


	tag "$id"
	cpus 1
	executor 'local'

	
	input:
    tuple val(id),val(family) 
	
	output:
	tuple val(id), val(family), env(status), emit: test_ch
	
	"""
	if [ -f $launchDir/results/bams/"$id".bam ]
    then 
	status="done"
	else
	status="todo"
	fi
	"""
	
	}