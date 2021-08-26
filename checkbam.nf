nextflow.enable.dsl=2
process checkbam {


	tag "$id"
	input:
    tuple val(id),val(family) 
	
	output:
	tuple val("$todoID"), val("$todoF"), optional: true
	tuple env(doneID), path('$launchDir/results/bams/"$id".bam'),path('$launchDir/results/bams/"$id".bam.bai'), optional: true
	"""
	[ -f $launchDir/results/bams/"$id".bam ] && gatk ValidateSamFile -I $launchDir/results/bams/"$id".bam -O $launchDir/results/bams/"$id"_validate_report.txt
    V1=\$( cat $launchDir/results/bams/"$id"_validate_report.txt)
    if [ "\$V1" == "No errors found" ]; then
    echo "bamtastic"
	doneID="\$id"
	else
	todoID="\$id"
	todoF="\$family"
	
    fi
	"""
	
	}