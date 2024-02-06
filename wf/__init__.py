"""Latch workflow for comparing cluster/condition groupings in an ArchRProject
"""

from wf.compare_task import compare_task, Groupings, Genome, CompareOutput

from typing import List

from latch.resources.workflow import workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchMetadata,
    LatchParameter,
    LatchRule,
)


metadata = LatchMetadata(
    display_name="compare clusters",
    author=LatchAuthor(
        name="AtlasXomics Inc.",
        email="joshuab@atlasxomics.com",
        github="https://github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/combined_cluster_wf",
    license="MIT",
    parameters={
        "project_name": LatchParameter(
            display_name="project name",
            description="Name or identifier for project.",
            batch_table_column=True,
        ),
        "archrproject": LatchParameter(
            display_name="ArchRProject",
            description="Path to ArchRProject folder containing a \
                        Save-ArchR-Project.rds file.",
            batch_table_column=True
        ),
        "groupings": LatchParameter(
            display_name="Specifications of groupings",
            description="Cluster and condition specifications for the two \
                        cell groupings to be compared.",
            samplesheet=True
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome',
            batch_table_column=True,
        ),
        "output_directory": LatchParameter(
            display_name="output directory",
            batch_table_column=True,
            description="Name of Latch directory for output files; files \
                        will be saved to /compare_outs/{output directory}.",
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="output directory name cannot start with a '/'"
                ),
                LatchRule(
                    regex="^\S+$",
                    message="directory name cannot contain whitespace"
                )
            ]
        ),
    },
    tags=[],
)


@workflow(metadata)
def compare_workflow(
    project_name: str,
    groupings: List[Groupings],
    archrproject: LatchDir,
    genome: Genome,
    output_directory: str
) -> CompareOutput:

    return compare_task(
        project_name=project_name,
        groupings=groupings,
        archrproject=archrproject,
        genome=genome,
        output_directory=output_directory
    )


LaunchPlan(
    compare_workflow,
    "default",
    {
        "project_name": "default",
        "groupings": [
            Groupings(
                clusterA="C2-4",
                conditionA="young",
                clusterB="C5",
                conditionB="old"
            )
        ],
        "archrproject": LatchDir(
            "latch://13502.account/ArchRProjects/Babayev_2/Babayev_2_ArchRProject",
        ),
        "genome": Genome.mm10,
        "output_directory": "default"
    },
)
