process download_reads {
    input:
        val accession
    output:
        tuple val(accession), path("${accession}_1.fastq"), path("${accession}_2.fastq")
    script:
    """
    fasterq-dump --split-files $accession
    """
}

process run_qc {
    input:
        tuple val(reads_label), path(r1), path(r2)
    output:
        path "${reads_label}_qc_report/", type: 'folder'
    script:
    """
    mkdir ${reads_label}_qc_report
    fastqc -o ${reads_label}_qc_report/ $r1 $r2
    """
}

process trimm {
    input:
        tuple val(reads_label), path(r1), path(r2)
    output:
        tuple val(reads_label), path("${reads_label}_R1_p.fq.gz"), path("${reads_label}_R2_p.fq.gz")
    script:
    """
    trimmomatic PE $r1 $r2 \
        ${reads_label}_R1_p.fq.gz ${reads_label}_R1_u.fq.gz \
        ${reads_label}_R2_p.fq.gz ${reads_label}_R2_u.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

process assemble {
    input:
        tuple val(reads_label), path(r1), path(r2)
    output:
        path "spades_out/contigs.fasta"
    script:
    """
    spades.py -1 $r1 -2 $r2 -o spades_out
    """
}

process index_reference {
    input:
        path ref
    output:
        tuple path(ref), path("${ref}.*")
    script:
    """
    bwa index $ref
    """
}

process samtools_faidx {
    input:
        path ref
    output:
        tuple path(ref), path("${ref}.fai")
    script:
    """
    samtools faidx $ref
    """
}

process map_reads {
    input:
        tuple val(reads_label), path(r1), path(r2)
        tuple path(ref), path(index_files)
    output:
        tuple val(reads_label), path("${reads_label}.bam"), path("${reads_label}.bam.bai")
    script:
    """
    bwa mem $ref $r1 $r2 \
        | samtools sort -o ${reads_label}.bam
    samtools index ${reads_label}.bam
    """
}

process coverage_plot {
    input:
        tuple val(reads_label), path(bam), path(bai)
    output:
        path "${reads_label}_coverage.png"
    script:
    """
    samtools depth -a $bam > ${reads_label}_depth.txt
    plot_coverage.py ${reads_label}_depth.txt ${reads_label}_coverage.png
    """
}
