process checkbam {


	tag "$id"
	input:
    tuple val(id),val(family) 
	output:
	tuple val(id),env(outcome)
	
	"""
	gatk ValidateSamFile -I ./results/bams/$id.bam -O /results/bams/$id_validate_report.txt
    V1=$( cat /results/bams/$id_validate_report.txt)
    if [ "$V1" == "No errors found" ]; then
    outcome="true"
	else
	outcome="false"
    fi
	"""