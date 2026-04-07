# Base Image
FROM bioconductor/bioconductor_docker:devel

# Install R Packages
RUN R -e "pak::pak('rikenbit/symTensor', dependencies=TRUE);\
    tools::testInstalledPackage('symTensor')"
