process peddy {

	tag "$fam"
	cpus 4
	time { 15.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		container "docker://xmzhuo/peddy:0.4.8"
		publishDir "./QC/$fam", mode: 'copy', overwrite: true
		
	input:
	tuple val(fam), val(caller), file(vcf), file(tbi)
	path ped
	
	output:
	tuple val(family), val(caller), file("${fam}.${caller}.normalized.html")
	
"""
	peddy -p ${task.cpus} --plot --prefix ${fam}.${caller}.normalized ${fam}.${caller}.normalized.vcf.gz $ped
"""

}