process 'delete_file' {

	tag "$inputFile"


  input:
	file(inputFile)

	output:
    // There is no output
	
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}