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
	number_of_files=\$(ls $launchDir/results/familyvcfs/"$family"*.vcf.gz | wc -l)
	if [  \$number_of_files > 0 ]
	then 
	status="done"
	else
	status="todo"
	fi
	"""
	
	}