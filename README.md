,# compare clusters

    **compare clusters** is a [latch.bio](https://latch.bio/) workflow for
    comparing differences in genes, peaks, and motifs between user-defined
    cluster/condition groupings within an ArchRProject.  Provided an
    ArchRProject and grouping specifications, **compare clusters** generates,
    * for genes
        * volcano plot
        * gene_markers.csv
    * for [peaks](https://www.archrproject.com/bookdown/pairwise-testing-between-groups.html)
        * volcano plot
        * peak_markers.csv
    * for motifs (https://www.archrproject.com/bookdown/motif-enrichment-in-differential-peaks.html​
        * enrichment plot
        * motif_enrichment.csv

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
    * Groupings **cannot** share clusters.  For example, C1-C3 versus C3,C4 is
    not allowed.
    * The Condition **must** match that provided when generating the
    ArchRProject (i.e., in the **create ArchRProject** Latch Workflow).
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