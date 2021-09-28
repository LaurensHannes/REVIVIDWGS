nextflow.enable.dsl=2
process checkbam {


	tag "$id"
	cpus 1
	input:
    tuple val(id),val(family) 
	
	output:
	tuple val(id), val(family), env(status), emit: test_ch
	
	"""
	if [ -f $launchDir/results/bams/"$id".bam ]
    then gatk ValidateSamFile -I $launchDir/results/bams/"$id".bam -O $launchDir/results/bams/"$id"_validate_report.txt
    V1=\$( cat $launchDir/results/bams/"$id"_validate_report.txt)
    if [ "\$V1" == "No errors found" ]; then
    echo "bamtastic"
	status="done"
	fi
	else
	status="todo"
	fi
	"""
	
	}