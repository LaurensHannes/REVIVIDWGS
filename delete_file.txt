process 'delete_file' {
  input:
    file inputFile
    val sampleName
    output:
    // There is no output
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}