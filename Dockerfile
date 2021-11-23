FROM condaforge/mambaforge:4.10.3-7

COPY . .

RUN mamba create -n signac

RUN mamba install --yes -c conda-forge -c bioconda \
    bats==0.4.0 \
    r-optparse==1.7.1 \
    libpng==1.6.37 \
    r-cairo==1.5_12.2 \
    r-hdf5r==1.3.5 \
    r-workflowscriptscommon==0.0.7 \
    r-signac==1.4.0