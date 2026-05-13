FROM continuumio/miniconda3:24.1.2-0

RUN conda install -y -c bioconda -c conda-forge \
        sra-tools=3.1.1 \
        fastqc=0.12.1 \
        trimmomatic=0.39 \
    && conda clean -afy

# trimmomatic adapter files accessible by filename
RUN ln -s /opt/conda/share/trimmomatic/adapters/TruSeq3-PE.fa /opt/conda/bin/TruSeq3-PE.fa

ENV PATH=/opt/conda/bin:$PATH
