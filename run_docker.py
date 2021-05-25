import click
import subprocess
import os

@click.command()
# required options
@click.option(
        "-fa",
        "--fasta",
        type=click.Path(exists=True),
        required=True,
        help="reference sequences in fasta format",
)
@click.option(
        "-fq1",
        "--fastq1",
        type=click.Path(exists=True),
        required=True,
        help="fastq sequencing file of mate 1",
)
@click.option(
        "-fq2",
        "--fastq2",
        type=click.Path(exists=True),
        required=False,
        help="fastq sequencing file of mate 1",
)
def main(fasta, fastq1, fastq2):
    # generate docker container and run code
    cmd = "docker run -v "
    files = [fasta,fastq1,fastq2]
    abs_files = []
    for f in files:
        spl = f.split("/")
        abs_files.append(os.path.abspath(f) + ":/data/"+spl[-1])
    cmd += " -v ".join(abs_files)
    cmd += " --name dreem_cont -t dreem dreem -fa test.fasta -fq1 test_mate1.fastq -fq2 test_mate2.fastq"
    subprocess.call(cmd, shell=True)
    subprocess.call("docker cp dreem_cont:/data/output output", shell=True)
    subprocess.call("docker rm dreem_cont", shell=True)


if __name__ == "__main__":
    main()
