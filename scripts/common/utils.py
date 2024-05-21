import os
import subprocess


def utils_directory():
    return os.path.dirname(os.path.abspath(__file__))


def artefact_directory():
    return os.path.dirname(os.path.abspath(__file__)) + "/../artefacts"


def strip_ansi_chars(stream):
    cmd = [f"{utils_directory()}/strip_ansi_chars.sh"]
    result = subprocess.run(cmd, input=stream, capture_output=True, text=True)
    return result.stdout
