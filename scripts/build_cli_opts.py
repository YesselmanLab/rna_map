import yaml
import cloup
from pathlib import Path

@cloup.command()
@cloup.argument("param_file", type=cloup.Path(exists=True))
def main(param_file):
    """
    main function for script
    """
    with open(param_file) as f:
        cli_opts = yaml.safe_load(f)
    f = open(Path("rna_map/cli_opts.py"), "w")
    f.write("import yaml\n")
    f.write("import cloup\n")
    f.write("from cloup import option_group, option\n")
    f.write("from rna_map.logger import get_logger\n\n")
    f.write("log = get_logger('CLI_OPTS')\n\n")
    t = "    "
    for group, v in cli_opts.items():
        f.write(f"def {group}_options():\n")
        f.write(f"{t}return option_group(\n")
        f.write(f"{t}{t}\"{v['name']}\",\n")
        f.write(f"{t}{t}\"{v['description']}\",\n")
        for opt, opt_vals in v["opts"].items():
            f.write(f'{t}{t}option(\n{t}{t}{t}"{opt}",\n')
            if "alt_name" in opt_vals:
                f.write(f'{t}{t}{t}"{opt_vals["alt_name"]}",\n')
            if "type" in opt_vals:
                f.write(f'{t}{t}{t}type={opt_vals["type"]},\n')
            if "is_flag" in opt_vals:
                f.write(f"{t}{t}{t}is_flag={opt_vals['is_flag']},\n")
            if "default" in opt_vals:
                if opt_vals["default"] == "":
                    f.write(f'{t}{t}{t}default="",\n')
                elif opt_vals["default"] == "None":
                    f.write(f"{t}{t}{t}default=None,\n")
                elif type(opt_vals["default"]) == str:
                    f.write(f'{t}{t}{t}default="{opt_vals["default"]}",\n')
                else:
                    f.write(f"{t}{t}{t}default={opt_vals['default']},\n")
            if "required" in opt_vals:
                f.write(f"{t}{t}{t}required={opt_vals['required']},\n")
            if "help" in opt_vals:
                f.write(f"{t}{t}{t}help=\"{opt_vals['help']}\",\n")
            f.write(f"{t}{t}),\n")
        f.write(f"{t})\n\n")

    f.write(f"\n\n")
    f.write(f"def parse_cli_args(params, args):\n")
    for group, v in cli_opts.items():
        f.write(f"{t}# {group} options\n")
        for opt, opt_vals in v["opts"].items():
            if "param" not in opt_vals:
                continue
            spl = opt_vals["param"].split(":")
            param_statement = "params"
            opt_name = opt
            if opt_name.startswith("--"):
                opt_name = opt_name[2:]
            if opt_name.startswith("-"):
                opt_name = opt_name[1:]
            opt_name = opt_name.replace("-", "_")
            for e in spl:
                param_statement += f'["{e}"]'
            if "is_flag" in opt_vals:
                f.write(f"{t}if args['{opt_name}']:\n")
            elif "default" in opt_vals:
                if opt_vals["default"] == "None":
                    f.write(f"{t}if args['{opt_name}'] is not None:\n")
                else:
                    f.write(f"{t}if args['{opt_name}'] != {opt_vals['default']}:\n")
            else:
                f.write(f"{t}if args['{opt_name}'] is not None:\n")
            if "log_msg" in opt_vals:
                if "{value}" in opt_vals["log_msg"]:
                    f.write(
                        f"{t}{t}log.info(\"{opt_vals['log_msg']}\".format(value=args['{opt_name}']))\n"
                    )
                else:
                    f.write(f"{t}{t}log.info(\"{opt_vals['log_msg']}\")\n")
            f.write(f"{t}{t}{param_statement} = args['{opt_name}']\n")
    f.close()


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
