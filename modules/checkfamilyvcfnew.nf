process checkfamilyvcf {


	tag "$family"
	cpus 1
	executor 'local'
     runOptions = "-B $launchDir"


	
	input:
    val(family) 
	
	output:
	tuple val(family), env(status), emit: familyvcfcheck_ch
	
	"""
	if [ -f $launchDir/results/"$family"/"$family".alignedandtrimmed.vcf.gz ]
    then 
	status="done"
	else
	status="todo"
	fi
	"""
	
	}