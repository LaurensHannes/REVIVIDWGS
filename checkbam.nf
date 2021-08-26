process checkbam {


	tag "$id"
	input:
    tuple val(id),val(family) 
	output:
	tuple val(id),env(outcome)
	
	"""
	gatk ValidateSamFile -I $workDir/results/bams/"$id".bam -O $workDir/results/bams/"$id"_validate_report.txt
    V1=\$( cat $workDir/results/bams/"$id"_validate_report.txt)
    if [ "\$V1" == "No errors found" ]; then
    outcome="true"
	else
	outcome="false"
    fi
	"""
	
	}