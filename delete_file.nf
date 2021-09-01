process 'delete_file' {

	tag "$inputFile"

	time { 1.minute * task.attempt }
		 errorStrategy 'retry' 
		maxRetries 3
		cpus 1
  input:
	file(inputFile)
	tuple val(identifier), val(permission)
	output:
    // There is no output
	
	when: 
	$permission = "true"
	
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}