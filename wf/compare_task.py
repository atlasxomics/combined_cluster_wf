import logging
import subprocess
import re
from latch import custom_task
from latch.types import LatchDir

from dataclasses import dataclass
from enum import Enum

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


def strip_string(input: str):
    return re.sub(r'\s+', '', input)


def expand_string(input_string: str):
    """Split the input string into parts
    """

    # Remove whitespace
    stripped_string = strip_string(input_string)

    # 'C1,C2,C5-C7' -> ['C1', 'C2', 'C5-C7']
    split_comma = stripped_string.split(',')
    final_string = []

    for i in split_comma:

        # Handle 'C5-C7', return [list] + ['C5', 'C6', C7']
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

    # Convert list back into string: 'C1,C2,C3,C5,C6,C7'
    result_string = ','.join(final_string)
    return result_string


def resolve_bool(value: bool):
    if value:
        return 't'
    else:
        return 'f'


@custom_task(cpu=8, memory=64, storage_gib=1000)
def compare_task(
    project_name: str,
    groupings: Groupings,
    archrproject: LatchDir,
    genome: Genome,
) -> LatchDir:

    out_dir = f"/root/{project_name}/"
    remote_dir = f"latch:///compare_outs/{project_name}"

    # Have to ensure the project path is the same as when save in ArchR1.0.3
    archrproject_name = archrproject.remote_path.split("/")[-1]
    archrproj_dest = f"/root/{archrproject_name}"
    subprocess.run(
        ["mv", archrproject.local_path, archrproj_dest]
    )

    subprocess.run(["mkdir", out_dir])

    _r_cmd = [
        "Rscript",
        "wf/compare_clusters.R",
        project_name,
        expand_string(groupings.clusterA),
        strip_string(groupings.conditionA),
        resolve_bool(groupings.multipleA),
        expand_string(groupings.clusterB),
        strip_string(groupings.conditionB),
        resolve_bool(groupings.multipleB),
        archrproj_dest,
        genome.value,
        out_dir
    ]

    logging.info("Initiating R script.")
    subprocess.run(_r_cmd, check=True)

    logging.info("Rscript complete; uploading results.")

    subprocess.run(["mv", archrproj_dest, out_dir])

    # Get rid of unnecessary files
    subprocess.run(
        ["rm", "-r", f"{out_dir}Rplots.pdf", f"{out_dir}ArchRLogs"]
    )

    return LatchDir(out_dir, remote_dir)


if __name__ == "__main__":
    compare_task(
        project_name="Pieper_154_brain_TBI-Sham_TBI-24_Cluster2_develop",
        groupings=Groupings(
            clusterA="2",
            conditionA="Sham_24hr",
            multipleA=False,
            clusterB="2",
            conditionB="TBI_24hr",
            multipleB=False
        ),
        archrproject="latch://13502.account/snap_outs/Pieper_154_brain_ArchRFull_10Core/Pieper_154_brain_ArchRFull_10Core_ArchRProject",
        genome=Genome.mm10
    )
