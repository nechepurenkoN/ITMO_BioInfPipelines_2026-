params.input_reads_folder = ''
params.ncbi_id = ''
params.reference = ''

process download_reads {
    input:
        val accession
    output:
        tuple val(accession), path("${accession}_{1,2}.fastq")
    script:
    """
    fasterq-dump --split-files $accession
    """
}

process run_initial_qc {
    input:
        tuple val(reads_label), path(reads)
    output:
        path "initial_qc_report/", type: 'folder'
    script:
    """
    mkdir initial_qc_report
    fastqc -o initial_qc_report/ $reads
    """
}

process run_trimmed_qc {
    input:
        tuple val(reads_label), path(reads)
    output:
        path "trimmed_qc_report/", type: 'folder'
    script:
    """
    mkdir trimmed_qc_report
    fastqc -o trimmed_qc_report/ $reads
    """
}

process trimm {
    input:
        tuple val(reads_label), path(reads)
    output:
        tuple val(reads_label), path("out_R?_p.fq.gz")
    script:
    """
    trimmomatic PE $reads \
        out_R1_p.fq.gz out_R1_u.fq.gz \
        out_R2_p.fq.gz out_R2_u.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

process assemble {
    input:
        tuple val(reads_label), path(reads)
    output:
        path "spades_out/contigs.fasta"
    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out
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

process map_reads {
    input:
        tuple val(reads_label), path(reads)
        tuple path(ref), path(index_files)
    output:
        tuple val(reads_label), path("${reads_label}.bam"), path("${reads_label}.bam.bai")
    script:
    """
    bwa mem $ref ${reads[0]} ${reads[1]} \
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

workflow {
    main:
        if (params.ncbi_id) {
            reads_ch = download_reads(params.ncbi_id)
        } else {
            reads_ch = channel.fromFilePairs("${params.input_reads_folder}/*_{1,2}.f*q*")
        }

        initial_qc_result = run_initial_qc(reads_ch)
        trimmed_ch = trimm(reads_ch)
        trimmed_qc_result = run_trimmed_qc(trimmed_ch)

        if (params.reference) {
            ref_ch = channel.fromPath(params.reference)
        } else {
            ref_ch = assemble(trimmed_ch)
        }

        indexed_ref_ch = index_reference(ref_ch).first()
        bam_ch = map_reads(trimmed_ch, indexed_ref_ch)
        coverage_result = coverage_plot(bam_ch)

    publish:
        initial_qc = initial_qc_result
        trimmed_reads = trimmed_ch
        trimmed_qc = trimmed_qc_result
        mapping = bam_ch
        coverage = coverage_result
}

output {
    initial_qc {
        path 'initial_qc'
    }
    trimmed_reads {
        path 'trimmed_reads'
    }
    trimmed_qc {
        path 'trimmed_qc'
    }
    mapping {
        path 'mapping'
    }
    coverage {
        path 'coverage'
    }
}
