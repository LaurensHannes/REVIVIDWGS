process 'delete_file' {
  input:
    tuple val(id), file(inputFile)
    output:
    // There is no output
  shell:
    '''
    rm "$(readlink '!{inputFile}')"
    '''
}