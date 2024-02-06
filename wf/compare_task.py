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
    clusterB: str
    conditionB: str


class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'


@dataclass
class CompareOutput:
    visual_output_dir: LatchDir


def expand_string(input_string):
    """Split the input string into parts
    """
    stripped_string = re.sub(r'\s+', '', input_string)
    split_comma = stripped_string.split(',')
    final_string = []
    for i in split_comma:
        parts = i.split('-')

        if len(parts) == 2:
            # Extract the start and end values
            start, end = parts

            # Extract the letter part and numeric part
            start_num = int(start[1:])
            end_num = int(end[1:])

            # Generate the list of values
            result_list = [f"C{num}" for num in range(start_num, end_num + 1)]

            # Join the list into a comma-separated string
            result_string = ','.join(result_list)

            final_string += result_list
        else:
            final_string += [i]

    result_string = ','.join(final_string)
    return result_string


@large_task
def compare_task(
    project_name: str,
    groupings: List[Groupings],
    archrproject: LatchDir,
    genome: Genome,
    output_directory: str
) -> CompareOutput:

    project_dir = f"/root/{output_directory}/{project_name}_compare_clusters"
    local_dir = f"/root/{output_directory}"
    remote_dir = f"latch:///compare_wf/{output_directory}"

    subprocess.run(["mkdir", local_dir])
    subprocess.run(["mkdir", project_dir])

    _r_cmd = [
        "Rscript",
        "wf/compare_clusters.R",
        project_name,
        expand_string(groupings[0].clusterA),
        groupings[0].conditionA,
        expand_string(groupings[0].clusterB),
        groupings[0].conditionB,
        archrproject.local_path,
        genome.value,
        project_dir
    ]

    logging.info("Initiating R script.")
    subprocess.run(_r_cmd)
    logging.info("Rscript complete; uploading results.")

    return CompareOutput(
        visual_output_dir=LatchDir(local_dir, remote_dir)
    )


if __name__ == "__main__":
    compare_task(
        project_name="D1234_default",
        groupings=[Groupings(
            clusterA="",
            conditionA="Young",
            clusterB="C5-C6",
            conditionB="Old"
        )],
        archrproject="latch://13502.account/ArchRProjects/Babayev_2/Babayev_2_ArchRProject",
        genome=Genome.mm10,
        output_directory="dev_test"
    )
