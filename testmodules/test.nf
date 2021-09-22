process 'test' {

	tag "$inputFile"
	cpcus 1

  input:
	file(inputFile)

	output:
    // There is no output

	
  shell:
    """
	echo ${inputFile}
    """
}