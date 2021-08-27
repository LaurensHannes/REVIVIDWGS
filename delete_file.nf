process 'delete_file' {
  input:
    tuple val(id), val(lane), file(inputFile)
    tuple val(idnextstep), val(lanenextstep), file(files)
	output:
    // There is no output
	
	when:
	lanenextstep =~ /lane/
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}