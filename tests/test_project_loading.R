# Run from the repository root: Rscript --vanilla tests/test_project_loading.R
expressions <- parse("wf/compare_clusters.R")
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("load_comparison_project"))) eval(expr)
}

setClass("ComparisonProjectFixture", slots = c(projectMetadata = "list"))
path <- tempfile("comparison-project-")
dir.create(path)
project_file <- file.path(path, "Save-ArchR-Project.rds")
project <- new("ComparisonProjectFixture", projectMetadata = list(
  GroupCoverages = list(Clusters = list(coverageMetadata = data.frame(
    File = "/old/machine/GroupCoverages/missing.h5"
  ))),
  outputDirectory = "/old/machine"
))

# Simulate the failing ArchR cache check, while retaining other validation.
loadArchRProject <- function(path) {
  loaded <- readRDS(file.path(path, "Save-ArchR-Project.rds"))
  for (coverage in loaded@projectMetadata$GroupCoverages) {
    zfiles <- coverage$coverageMetadata$File
    stopifnot(all(file.exists(zfiles)))
  }
  if (!file.exists(file.path(path, "required.arrow"))) stop("Missing Arrow")
  loaded
}

saveRDS(project, project_file)
original_checksum <- tools::md5sum(project_file)
file.create(file.path(path, "required.arrow"))
loaded <- load_comparison_project(path)
stopifnot(
  is.null(loaded@projectMetadata$GroupCoverages),
  identical(loaded@projectMetadata$outputDirectory, "/old/machine"),
  identical(tools::md5sum(project_file), original_checksum)
)

# Required data errors must propagate, and the original RDS must be restored.
unlink(file.path(path, "required.arrow"))
error <- tryCatch(load_comparison_project(path), error = identity)
stopifnot(
  inherits(error, "error"),
  identical(conditionMessage(error), "Missing Arrow"),
  identical(tools::md5sum(project_file), original_checksum)
)

# Projects with no cache should load normally without rewriting their RDS.
project@projectMetadata$GroupCoverages <- NULL
saveRDS(project, project_file)
original_checksum <- tools::md5sum(project_file)
file.create(file.path(path, "required.arrow"))
stopifnot(
  identical(load_comparison_project(path), project),
  identical(tools::md5sum(project_file), original_checksum),
  length(list.files(path)) == 2L
)
unlink(path, recursive = TRUE)
cat("Project loading regression checks passed.\n")
