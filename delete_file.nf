process 'delete_file' {
  input:
    tuple val(id), file(inputFile)
    tuple val(idnextstep), file(files)
	output:
    // There is no output
	
	when:
	idnextstep =~ /id/
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}