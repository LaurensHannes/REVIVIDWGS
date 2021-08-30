process 'delete_file' {

	tag "$inputFile"


  input:
	file(inputFile)
	tuple val(identifier), val(permission)
	output:
    // There is no output
	
	when: 
	$permission = "true"
	
  shell:
    '''
    echo $identifier
	echo $inputFile
    '''
}