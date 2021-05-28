import shutil
import re

from dreem import parameters, settings

static = """
@optgroup.group("main arguments")
@optgroup.option("-fa", "--fasta", type=click.Path(exists=True), required=True,
                 help="reference sequences in fasta format")
@optgroup.option("-fq1", "--fastq1", type=click.Path(exists=True), required=True,
                 help="fastq sequencing file of mate 1")
@optgroup.option("-fq2", "--fastq2", type=click.Path(exists=True),
                 help="fastq sequencing file of mate 2")
@optgroup.group("common options")
@optgroup.option("--db", type=click.Path(exists=True),
                help="A csv formatted file that contains dot bracket info for each sequence")
@optgroup.option("-pf", "--param-file", type=click.Path(exists=True),
                help="A yml formatted file to specify parameters")
@optgroup.option("-ow", "--overwrite", is_flag=True,
                help="overwrites previous results, if not set will keep previous " 
                     "calculation checkpoints")
@optgroup.option("-ll", "--log-level", help="set log level (INFO|WARN|DEBUG|ERROR|FATAL)", default="INFO")
"""

def get_parameters(obj, name, skip):
    args = obj.__dict__
    s = f"@optgroup.group(\"{name} options\")\n"
    for k, v in args.items():
        if k in skip:
            continue
        if type(v) == bool:
            s += f"@optgroup.option(\"--{k}\", is_flag=True,\nhelp=\"{obj.description[k]}\")\n"
        else:
            s += f"@optgroup.option(\"--{k}\", default=None,\nhelp=\"{obj.description[k]}\")\n"
    return s

def main():
    pf = parameters.ParametersFactory()
    run_py_file_path = settings.get_py_path() + "/run.py"
    f = open(run_py_file_path)
    lines = f.readlines()
    f.close()
    s = "".join(lines)
    spl = re.split('@click.command\(\)[\S\s]+def main\(\*\*args\)\:',s)
    f = open(run_py_file_path, "w")
    f.write(spl[0])
    f.write("@click.command()\n")
    f.write(static)
    f.write(get_parameters(pf._Map(), "map", "description".split(",")))
    f.write(get_parameters(pf._BitVector(), "bit vector",
                           "description,miss_info,ambig_info,nomut_bit,del_bit".split(",")))
    f.write("def main(**args):")
    f.write(spl[1])
    f.close()



if __name__ == "__main__":
    main()
