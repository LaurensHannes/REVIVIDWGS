process 'test' {

	tag "$inputFile"
	cpus 1
executor 'local'
  input:
	file(inputFile)

	output:
    // There is no output

	
  shell:
    """
	echo ${inputFile}
    """
}