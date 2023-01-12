



def _test_input_validation():
    ins = validate_inputs(p["fasta"], p["fastq1"])
    assert p["fasta"] == ins.fasta
    assert ins.csv == ""
    assert ins.is_paired() == False
    assert ins.supplied_csv() == False

    # check to make sure we get the proper errors for supplying file that does
    # not exist
    with pytest.raises(DREEMInputException) as exc_info:
        validate_inputs(p["fasta"], "")
    assert exc_info.value.args[0] == "fastq1 file: does not exist !"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_inputs("fake_path", p["fastq1"])
    assert exc_info.value.args[0] == "fasta file: does not exist fake_path!"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_inputs(p["fasta"], p["fastq1"], "fake_path")
    assert exc_info.value.args[0] == "fastq2 file: does not exist fake_path!"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_inputs(p["fasta"], p["fastq1"], csv="fake_path")
    assert exc_info.value.args[0] == "csv file: does not exist !"


def _test_fasta_checks():
    fasta_test_path = TEST_DIR + "/resources/test_fastas/"
    path = fasta_test_path + "blank_line.fasta"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_fasta_file(path)
    assert (
        exc_info.value.args[0]
        == "blank line found on ln: 1. These are not allowed in fastas."
    )
    path = fasta_test_path + "incorrect_format.fasta"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_fasta_file(path)
    assert (
        exc_info.value.args[0]
        == "reference sequence names are on line zero and even numbers. line 0 "
        "has value which is not correct format in the fasta"
    )
    path = fasta_test_path + "incorrect_sequence.fasta"
    with pytest.raises(DREEMInputException) as exc_info:
        validate_fasta_file(path)
    print(exc_info)
