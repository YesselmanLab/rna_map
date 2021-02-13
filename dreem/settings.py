import os


def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return base_dir
