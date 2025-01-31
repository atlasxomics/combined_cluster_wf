import logging
import subprocess
import re
from latch import large_task
from latch.types import LatchDir

from dataclasses import dataclass
from enum import Enum
from typing import List

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


@dataclass
class Groupings:
    clusterA: str
    conditionA: str
    multipleA: bool
    clusterB: str
    conditionB: str
    multipleB: bool


class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'


@dataclass
class CompareOutput:
    visual_output_dir: LatchDir


def strip_string(input):
    return re.sub(r'\s+', '', input)


def expand_string(input_string):
    """Split the input string into parts
    """
    stripped_string = strip_string(input_string)
    split_comma = stripped_string.split(',')
    final_string = []
    for i in split_comma:
        mult_digit = re.findall(r'\d+', i)
        if len(mult_digit) > 1:

            # Extract the letter part and numeric part
            start_num = int(mult_digit[0])
            end_num = int(mult_digit[1])

            # Generate the list of values
            result_list = [f"C{num}" for num in range(start_num, end_num + 1)]

            # Join the list into a comma-separated string
            result_string = ','.join(result_list)

            final_string += result_list
        else:
            final_string += [i]

    result_string = ','.join(final_string)
    return result_string


def resolve_bool(value: bool):
    if value:
        return 't'
    else:
        return 'f'


@large_task
def compare_task(
    project_name: str,
    groupings: List[Groupings],
    archrproject: LatchDir,
    genome: Genome,
) -> CompareOutput:

    out_dir = f"/root/{project_name}/"
    remote_dir = f"latch:///compare_outs/{project_name}"

    subprocess.run(["mkdir", out_dir])
    _r_cmd = [
        "Rscript",
        "wf/compare_clusters.R",
        project_name,
        expand_string(groupings[0].clusterA),
        strip_string(groupings[0].conditionA),
        resolve_bool(groupings[0].multipleA),
        expand_string(groupings[0].clusterB),
        strip_string(groupings[0].conditionB),
        resolve_bool(groupings[0].multipleB),
        archrproject.local_path,
        genome.value,
        out_dir
    ]

    logging.info("Initiating R script.")
    subprocess.run(_r_cmd, check=True)

    logging.info("Rscript complete; uploading results.")

    # Get rid of unnecessary files
    subprocess.run(
        ["rm", "-r", f"{out_dir}Rplots.pdf", f"{out_dir}ArchRLogs"]
    )

    return CompareOutput(
        visual_output_dir=LatchDir(out_dir, remote_dir)
    )


if __name__ == "__main__":
    compare_task(
        project_name="D1234_default",
        groupings=[Groupings(
            clusterA="C4, C8, C3",
            conditionA="Tumor",
            multipleA=True,
            clusterB="C1",
            conditionB="Pdx",
            multipleB=True
        )],
        archrproject="latch://13502.account/ArchRProjects/Gaykalova/Gaykalova_ArchRProject",
        genome=Genome.mm10
    )
