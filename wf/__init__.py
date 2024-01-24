"""Latch workflow for normalizing hot rows and columns in spatial ATAC-seq data
"""

from wf.compare_task import compare_task, Groupings, Genome, CompareOutput

from typing import List

from latch import workflow, map_task
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchFile,
    LatchDir,
    LatchMetadata,
    LatchParameter,
    LatchRule,
)

    
    
metadata = LatchMetadata(
    display_name="compare cluster wf",
    author=LatchAuthor(
        name="AtlasXomics Inc.",
        email="joshuab@atlasxomics.com",
        github="https://github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/combined_cluster_wf",
    parameters={
        "project_name": LatchParameter(
            display_name="run id",
            description="ATX Run ID with optional prefix, default to \
                        Dxxxxx_NGxxxxx format.",
            batch_table_column=True,
        ),
        "archrproject": LatchParameter(
            display_name="ArchRProject",
            description="",
            batch_table_column=True
        ),
        'groupings': LatchParameter(
            display_name='Specifications of groupings',
            description='Comparisons between clusters and conditions within subset',
            samplesheet=True
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome to be used for geneAnnotation and \
                        genomeAnnotation',
            batch_table_column=True,
        ),
        "output_directory": LatchParameter(
            display_name="output directory",
            batch_table_column=True,
            description="Name of Latch directory for merge fastq files; files \
                        will be saved to /impute/{output directory}.",
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
                clusterA="C1",
                conditionA="young",
                clusterB="C2-4",
                conditionB="old"
            )
        ],
        "archrproject": LatchDir(
            "latch:///compare_wf/D1234/demo_ArchRProject"
        ),
        "genome": Genome.mm10,
        "output_directory": "demo"
    },
)
