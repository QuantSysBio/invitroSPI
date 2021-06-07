shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake import utils
from snakemake.utils import min_version
from snakemake import logging

min_version("6.0")
shell.prefix("set -euo pipefail;")

config = yaml.load(open("INPUT/config.yaml", "r+"), Loader=yaml.FullLoader)
dependencies = yaml.load(open("SOURCE/dependencies.yaml", "r+"), Loader=yaml.FullLoader)

snakefiles = "SOURCE/"
include: snakefiles + "rules.py"

rule all:
    input:
        ProteasomeDB = "OUTPUT/ProteasomeDB.csv"
