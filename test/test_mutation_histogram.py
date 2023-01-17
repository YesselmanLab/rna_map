"""
test mutational histogram
"""
import os
import pickle
import json
import pytest

from dreem.mutation_histogram import (
    MutationHistogram,
    get_mut_histos_from_json_file,
    plot_read_coverage,
    plot_modified_bases,
    plot_mutation_histogram,
    plot_population_avg_old,
    plot_population_avg,
    merge_mut_histo_dicts
)

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def get_example_mut_histo() -> MutationHistogram:
    pickle_path = os.path.join(
        TEST_DIR,
        "resources",
        "case_1",
        "output",
        "BitVector_Files",
        "mutation_histos.p",
    )
    with open(pickle_path, "rb") as f:
        mhs = pickle.load(f)
    return mhs["mttr-6-alt-h3"]


def test_mutation_histogram():
    mh = MutationHistogram("construct_1", "ACGTACGT", "DMS")
    assert mh.start == 1
    assert mh.end == 8


def test_merge():
    """
    test merge
    """
    mh = get_example_mut_histo()
    mh2 = get_example_mut_histo()
    mh.merge(mh2)
    assert mh.num_reads == mh2.num_reads * 2
    # picked 10 at random to check
    for i in range(1, 134):
        assert mh.num_of_mutations[i] == mh2.num_of_mutations[i] * 2
        assert mh.mod_bases["A"][i] == mh2.mod_bases["A"][i] * 2
    assert mh.skips["low_mapq"] == mh2.skips["low_mapq"] * 2


def test_read_mutation_histogram():
    mh = get_example_mut_histo()
    assert mh.name == "mttr-6-alt-h3"


def test_mutation_histogram_to_json():
    mh = get_example_mut_histo()
    d = {mh.name: mh.get_dict()}
    with open("test.json", "w") as f:
        json.dump(d, f)
    mhs = get_mut_histos_from_json_file("test.json")
    assert type(mhs[mh.name]) == MutationHistogram
    os.remove("test.json")


def test_plot_read_coverage():
    """
    test plot_read_coverage
    """
    mh = get_example_mut_histo()
    fname = f"{mh.name}_{mh.start}_{mh.end}_read_coverage.html"
    plot_read_coverage(mh.get_nuc_coords(), mh.get_read_coverage(), fname)
    assert os.path.isfile("mttr-6-alt-h3_1_134_read_coverage.html")
    os.remove("mttr-6-alt-h3_1_134_read_coverage.html")


def test_plot_modified_bases():
    """
    test plot_modified_bases
    """
    mh = get_example_mut_histo()
    fname = f"{mh.name}_{mh.start}_{mh.end}_mutations.html"
    plot_modified_bases(mh.get_nuc_coords(), mh.mod_bases, fname)
    assert os.path.isfile("mttr-6-alt-h3_1_134_mutations.html")
    os.remove("mttr-6-alt-h3_1_134_mutations.html")


def test_plot_mutation_histogram():
    """
    test plot_mutation_histogram
    """
    mh = get_example_mut_histo()
    fname = f"{mh.name}_{mh.start}_{mh.end}_mutation_histogram.html"
    plot_mutation_histogram(mh.get_nuc_coords(), mh.num_of_mutations, fname)
    assert os.path.isfile("mttr-6-alt-h3_1_134_mutation_histogram.html")
    os.remove("mttr-6-alt-h3_1_134_mutation_histogram.html")


def test_get_pop_avg_dataframe():
    """
    test get_pop_avg_dataframe
    """
    mh = get_example_mut_histo()
    df = mh.get_pop_avg_dataframe()
    assert len(df) == 134
    assert "position" in df.columns
    assert "mismatches" in df.columns
    df.drop("nuc", inplace=True, axis=1)
    df = df.rename(
        columns={
            "position": "Position",
            "mismatches": "Mismatches",
            "mismatch_del": "Mismatches + Deletions",
        }
    )
    # df.to_csv("test.tsv", sep="\t", index=False)
    # popavg_filename = file_base_name + "popavg_reacts.txt"


def test_plot_population_avg_old():
    """
    test plot_population_avg_old
    """
    mh = get_example_mut_histo()
    df = mh.get_pop_avg_dataframe()
    fname = f"{mh.name}_{mh.start}_{mh.end}_pop_avg.html"
    plot_population_avg_old(df, mh.name, fname)
    assert os.path.isfile("mttr-6-alt-h3_1_134_pop_avg.html")
    os.remove("mttr-6-alt-h3_1_134_pop_avg.html")


def _test_plot_population_avg():
    """
    test plot_population_avg
    """
    mh = get_example_mut_histo()
    df = mh.get_pop_avg_dataframe()
    fname = f"{mh.name}_{mh.start}_{mh.end}_pop_avg.html"
    plot_population_avg(df, mh.name, fname)
    assert os.path.isfile("mttr-6-alt-h3_1_134_pop_avg.html")
    os.remove("mttr-6-alt-h3_1_134_pop_avg.html")
