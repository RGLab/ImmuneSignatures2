args = commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Argumant must be supplied GITHUB_PAT")
}

# re-install preprocessCore to avoid issue with openblas >= 0.3.5
# https://support.bioconductor.org/p/122925/
BiocManager::install("preprocessCore", configure.args = "--disable-threading", force = TRUE)

devtools::install_github("RGLab/ImmuneSignatures2", auth_token = args[1], dependencies = TRUE)
