"""Latch workflow for comparing cluster/condition groupings in an ArchRProject
"""
from typing import Union

from wf.compare_task import compare_task, Barcodes, Groupings, Genome

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
        email="jamesm@atlasxomics.com",
        github="https://github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/combined_cluster_wf",
    license="MIT",
    parameters={
        "project_name": LatchParameter(
            display_name="project name",
            description="Name or identifier for project; outputs will be \
                        saved to /compare_outs/[project_name].",
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                ),
                LatchRule(
                    regex="^\S+$",
                    message="project name cannot contain whitespace"
                )
            ],
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
            batch_table_column=True
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome',
            batch_table_column=True,
        ),
    },
    tags=[],
)


@workflow(metadata)
def compare_workflow(
    project_name: str,
    groupings: Union[Groupings, Barcodes],
    archrproject: LatchDir,
    genome: Genome,
) -> LatchDir:

    '''Explore differences in genes, peaks, and motifs within an ArchRProject.

    # compare clusters

    **compare clusters** is a [latch.bio](https://latch.bio/) workflow for
    comparing differences in genes, peaks, and motifs between user-defined
    cluster/condition groupings within an ArchRProject.  Provided an
    ArchRProject and grouping specifications, **compare clusters** generates,
    * for genes
        * volcano plot
        * all_genes.csv: all genes with test results from ArchR::getMarkerFeatures
        * marker_genes.csv: all_genes filtered and scored with significance
        thresholds; data for the volcano plot.
    * for [peaks](https://www.archrproject.com/bookdown/pairwise-testing-between-groups.html)
        * MA Plot
        * all_peaks: all peaks with test results from ArchR::getMarkerFeatures
        * marker_peaks.csv: all_genes filtered and scored with significance
        thresholds; data for the volcano plot.
    * for [motifs](https://www.archrproject.com/bookdown/motif-enrichment-in-differential-peaks.html​)
        * [up/down]Regulated_motifs.csv: Up or down regulated motifs ranked by 
        significance values; see ArchR [docs](https://www.archrproject.com/bookdown/motif-enrichment.html#motif-enrichment-in-differential-peaks)
        * [up/down] enrichment plot: scatter plot of motifs ranked by -log10(FDR)
        * all_motifs.csv: all motifs with test results from ArchR::getMarkerFeatures
        * marker_motifs.csv: all_motifs filtered and scored with significance
        thresholds; data for the volcano plot.
        * volcano plot

    ## Inputs
    All input files for **compare clusters** must be on the latch.bio
    [file system](https://wiki.latch.bio/wiki/data/overview). Each run in the
    workflow takes the following parameters,
    * project name: A name for the output folder
    * ArchRProject: A file path on the latch.bio file system pointing to a
    directory generated via [ArchR::saveArchRProject()](https://www.archrproject.com/reference/saveArchRProject.html)
    * genome: A reference genome for the ArchRProject
    * Specifications of groupings: Cluster and condition labels defining the
    subsets of cells to be compared.
        * ClusterA: Cluster specifications for first subset of cells.
        * ConditionA: Condition specifications for first subset of cells.
        * MultipleA: If an ArchRProject contains multiple conditions, returns
        all cells associated with CondationA (see below).
        * ClusterB: Cluster specifications for second subset of cells.
        * ConditionB: Condition specifications for first subset of cells.
        * MultipleB: If an ArchRProject contains multiple conditions, returns
        all cells associated with CondationB (see below).

    #### Rules for groupings
    * Clusters can be specified as a common separated list (C2,C3,C5), a range
    (C2-C4), or a combination of the two (C2-C4,C6).
    * The Condition **must** match that provided when generating the
    ArchRProject (i.e., in the **create ArchRProject** Latch Workflow).
    * If no Condition is specified, groupings **cannot** share clusters.
    For example, C1-C3 versus C3,C4 is not allowed; however, it is allow if the
    groupings have different clusters.
    * If multiple conditions are supplied when creating the ArchRProject, use
    the Multiple button to select all samples with that condition label.  As
    an example, take an ArchRProject with the following samples:
        * Sample1: old-control
        * Sample2: young-control
        * Sample3: old-experiment
        * Sample4: young-experiment

    To select cells from Sample1 and Sample2, enter 'young' into the Condition
    field and set the Multiple button to True (on).  Otherwise, Sample1 and
    Sample2 must be specified by entering 'old-control,old-experiment' into
    the Condition field.

    ## Outputs
    Outputs from **create ArchRProject** are loaded into latch.bio
    [Data module](https://wiki.latch.bio/wiki/data/overview) in the
    `compare_outs` directory.

    ## Next Steps
    Downsteam analysis can be performed locally or in a latch.bio
    [Pod](https://wiki.latch.bio/wiki/pods/overview).  For access to
    ATX-specific Pods, please contact your AtlasXomics Support Scientist.

    ## Support
    Questions? Comments?  Contact support@atlasxomics.com.

    '''

    return compare_task(
        project_name=project_name,
        groupings=groupings,
        archrproject=archrproject,
        genome=genome,
    )


LaunchPlan(
    compare_workflow,
    "default",
    {
        "project_name": "default",
        "groupings": Groupings(
            clusterA="C1-C3",
            conditionA="WT",
            multipleA=False,
            clusterB="C4,C6,C7",
            conditionB="Lupus",
            multipleB=False
        ),
        "archrproject": LatchDir(
            "s3://latch-public/test-data/13502/compare_ArchRProject",
        ),
        "genome": Genome.mm10,
    },
)
