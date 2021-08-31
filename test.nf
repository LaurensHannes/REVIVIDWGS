process 'test' {

	tag "$inputFile"


  input:
	file(inputFile)

	output:
    // There is no output

	
  shell:
    """
	echo ${inputFile}
    """
}