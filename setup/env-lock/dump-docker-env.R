# Dump the package set of the environment the analyses were run in (OOD RStudio image + ~/R/library)
ip <- as.data.frame(installed.packages(fields = c("Repository", "RemoteType", "RemoteSha")), stringsAsFactors = FALSE)
write.table(ip[, c("Package", "LibPath", "Version", "Priority", "Built", "Repository", "RemoteType", "RemoteSha")],
            "docker-installed-packages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
sink("docker-session.txt")
print(R.version.string); print(.libPaths()); print(extSoftVersion()); print(La_version()); print(La_library())
print(sessionInfo()); print(getOption("repos")); print(Sys.getenv(c("R_LIBS_USER","R_LIBS_SITE")))
sink()
