"""Latch workflow for comparing cluster/condition groupings in an ArchRProject
"""
from typing import Annotated, Union
from flytekit.core.annotation import FlyteAnnotation

from wf.compare_task import compare_task, Barcodes, Groupings

from latch.resources.workflow import workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule,
)
from latch.types.metadata import Fork, ForkBranch, Params, Spoiler, Text


flow = [
    Params("project_name"),
    Params("archrproject"),
    Spoiler(
        "Specify Sample Size",
        Text(
            "Larger sample sizes will take longer "
            "to run; default 500 cells."
        ),
        Fork(
            "Specify sample size",
            "Specify sample size",
            max_cells=ForkBranch("Number of cells", Params("max_cells")),
            use_max_possible_cells=ForkBranch(
                "Maximumize sample size", Params("use_max_possible_cells")
            ),
        ),
    ),
    Params("groupings"),
]


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
        "max_cells": LatchParameter(
            display_name="max cells per comparison group",
            description="Upper bound on the number of cells passed to "
                        "ArchR::getMarkerFeatures for each group.",
            batch_table_column=True
        ),
        "use_max_possible_cells": LatchParameter(
            display_name="use maximum possible cells",
            description="Ignore the explicit cap and use the largest "
                        "possible matched group size for the comparison.",
            batch_table_column=True
        ),
        "groupings": LatchParameter(
            display_name="Specifications of groupings",
            description="Cluster, condition, and sample specifications for \
                        the two cell groupings to be compared.",
            batch_table_column=True
        ),
    },
    tags=[],
    flow=flow
)


@workflow(metadata)
def compare_workflow(
    project_name: str,
    groupings: Annotated[
        Union[Groupings, Barcodes, LatchFile],
        FlyteAnnotation({
            "union_variant_names": ["Manual", "Barcodes", "File"]
        })
    ],
    archrproject: LatchDir,
    use_max_possible_cells: bool,
    max_cells: int = 500,
) -> LatchDir:

    '''Explore differences in genes, peaks, and motifs within an ArchRProject.

    # compare clusters

    **compare clusters** is a [latch.bio](https://latch.bio/) workflow for
    comparing differences in genes, peaks, and motifs between user-defined
    cluster/condition groupings within an ArchRProject.  Provided an
    ArchRProject and grouping specifications, **compare clusters** generates,
    * `UpdateClustName_by_barcode.csv`: selected cell barcodes and their
    assigned `UpdateClustName` values (`GroupA` or `GroupB`).
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
    * ArchRProject: A file path on the latch.bio file system pointing to an
    ArchR project directory containing a `Save-ArchR-Project.rds` file.
    * max cells: Upper bound on the number of cells used per group in
    ArchR::getMarkerFeatures comparisons.
    * use maximum possible cells: If true, use the largest possible matched
    group size from the selected cells and ignore the explicit cap.
    * Specifications of groupings:
        * **Manual**: Cluster, condition, and sample labels defining the
        groups of cells to be compared.
            * ClusterA: Cluster name(s) for first group of cells.
            * ConditionA: Condition name(s) for first group of cells.
            * SampleA: Sample name(s) for the first group of cells.
            * MultipleA: If an ArchRProject contains multiple conditions,
            returns all cells associated with CondationA (see below).
            * ClusterB: Cluster name(s) for second group of cells.
            * ConditionB: Condition names(s) for second group of cells.
            * SampleB: Sample name(s) for the second group of cells.
            * MultipleB: If an ArchRProject contains multiple conditions,
            returns all cells associated with CondationB (see below).
        * **Barcodes**: Comma-separated lists of barcodes (i.e.,
        'D00000#AAAAAAAAAAAAAAAA-1,D00000#TTTTTTTTAAAAAAAA-1...')
        * **File**: JSON file saved in Latch Data specifying barcodes for
        each group:
            ```
            {
                "groupA": [
                "D00000#GCAAGAAGGTCGTAGA-1",
                "D00000#ACTAGGTCGTGTTCTA-1",
                "D00000#TCCATTGGACACGACC-1",
                ],
                "groupB": [
                "D00001#ATCGCAGTGTGTTCTA-1",
                "D00001#TCACTGCATAGGATGA-1",
                "D00001#TCACTGCAAGCACCTC-1",
                ]
            }
            ```

    #### Rules for  manual  groupings
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
        max_cells=max_cells,
        use_max_possible_cells=use_max_possible_cells,
        groupings=groupings,
        archrproject=archrproject,
    )


LaunchPlan(
    compare_workflow,
    "default",
    {
        "project_name": "default",
        "max_cells": 500,
        "use_max_possible_cells": False,
        "groupings": Groupings(
            clusterA="C1-C3",
            conditionA="WT",
            sampleA="",
            multipleA=False,
            clusterB="C4,C6,C7",
            conditionB="Lupus",
            sampleB="",
            multipleB=False
        ),
        "archrproject": LatchDir(
            "s3://latch-public/test-data/13502/compare_ArchRProject",
        ),
    },
)
