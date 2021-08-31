process 'delete_file' {

	tag "$inputFile"
	cpus 1

  input:
	file(inputFile)

	output:
    // There is no output
	
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}