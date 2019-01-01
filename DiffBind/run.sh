# 1. Installation
enter R:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffBind", version = "3.8")

# 2. Obtaining differentially bound sites
# 2.1 Reading in the peaksets: 
# The first step is to read in a set of peaksets and associated metadata.
library(DiffBind)
tamoxifen <- dba(sampleSheet="tamoxifen.csv",dir=system.file("extra", package="DiffBind"))

# 2.2 Counting reads: 
# The next step is to calculate a binding matrix with scores based on read counts for every sample (affinity scores), rather than confidence scores for only those peaks called in a specific sample (occupancy scores). 
tamoxifen <- dba.count(tamoxifen,summits=250)

# 2.3 Establishing a contrast: 
# Before running the differential analysis, we need to tell DiffBind which cell lines fall in which groups. This is done using the dba.contrast function, as follows:
tamoxifen <- dba.contrast(tamoxifen, categories=DBA_CONDITION)

# 2.4 Performing the differential analysis
# The main differential analysis function is invoked as follows:
tamoxifen <- dba.analyze(tamoxifen)

# 2.5 Retrieving the differentially bound sites
# The final step is to retrieve the differentially bound sites as follows:
tamoxifen.DB <- dba.report(tamoxifen)


