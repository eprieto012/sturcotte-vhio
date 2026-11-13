# Run inside the OOD RStudio image (R 4.5.2, ~/R/library + site-library): records in renv.lock the
# versions of the packages the sturcotte scripts use, as installed in the environment the analyses were run in.
project <- "/mnt/bioinfnas/immuno/eprieto/analyses/sturcotte/renv-sturcotte"
libs <- .libPaths()[1:2]
ip <- installed.packages(lib.loc = .libPaths())
ip <- ip[!duplicated(ip[, "Package"]), ]
scripts <- list.files(file.path(project, "..", c("scripts", "setup")), pattern = "\\.(R|Rmd)$", full.names = TRUE, recursive = TRUE)
scripts <- scripts[!grepl("/setup/(MAGIC|env-lock)/", scripts)]
txt <- unlist(lapply(scripts, readLines, warn = FALSE))
direct <- unique(c(
  unlist(regmatches(txt, gregexpr("(?<=library\\()[A-Za-z0-9.]+(?=\\))", txt, perl = TRUE))),
  unlist(regmatches(txt, gregexpr("\\b[A-Za-z][A-Za-z0-9.]*(?=::)", txt, perl = TRUE))),
  "knitr", "rmarkdown", "UCell", "STACAS", "SingleCellExperiment", "BiocManager", "BiocVersion", "renv"
))
direct <- intersect(direct, rownames(ip))
# Optional (Suggests) packages of the analysis packages that were installed: some change behaviour when
# present (presto, glmGamPoi, ...). Development/testing/documentation tools and spatial stacks are left out.
analysis <- intersect(c("Seurat", "SeuratObject", "ProjecTILs", "scGate", "UCell", "STACAS", "scDblFinder", "SoupX",
                        "djvdj", "Rmagic", "pdfCluster", "ggplot2", "patchwork", "reticulate"), rownames(ip))
sug <- unique(unlist(tools::package_dependencies(analysis, db = ip, which = "Suggests")))
skip <- c("sf", "lwgeom", "rstan", "covr", "testthat", "vdiffr", "spelling", "roxygen2", "devtools", "pak", "remotes",
          "profvis", "shiny", "miniUI", "downlit", "tinytex", "sessioninfo", "gitcreds", "pkgload", "brio", "callr",
          "waldo", "rtracklayer", "GenomicAlignments", "Rsamtools", "here", "bslib", "sass", "pkgdown", "usethis",
          "rstudioapi", "plotly", "htmlwidgets", "DBI", "RSQLite", "dbplyr", "blob", "maps", "mapproj", "multcomp",
          "quantreg", "ragg", "xts", "profmem", "bench", "lobstr", "BiocStyle")
sug <- setdiff(intersect(sug, rownames(ip)), skip)
sug <- sug[is.na(ip[sug, "Priority"])]
pkgs <- sort(unique(c(direct, sug)))
writeLines(direct, "lock-direct-packages.txt"); writeLines(setdiff(sug, direct), "lock-suggested-packages.txt")
renv::snapshot(project = project, library = libs, lockfile = file.path(project, "renv.lock"),
               packages = pkgs, prompt = FALSE, force = TRUE)
