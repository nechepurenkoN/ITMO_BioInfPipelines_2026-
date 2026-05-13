include { download_reads; run_qc; trimm; assemble; index_reference; samtools_faidx; map_reads; coverage_plot } from './modules/hw2_processes'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CALL } from './modules/nf-core/bcftools/call/main'

workflow {
    main:
        if (params.ncbi_id) {
            reads_ch = download_reads(params.ncbi_id)
        } else {
            reads_ch = Channel.fromFilePairs("${params.input_reads_folder}/*_{1,2}.f*q*", flat: true)
        }

        trimmed_ch = trimm(reads_ch)

        reads_tagged = reads_ch.map   { label, r1, r2 -> tuple("initial_${label}", r1, r2) }
        trimmed_tagged = trimmed_ch.map { label, r1, r2 -> tuple("trimmed_${label}", r1, r2) }
        qc_result = run_qc(reads_tagged.mix(trimmed_tagged))

        if (params.reference) {
            ref_ch = Channel.value(file(params.reference))
        } else {
            ref_ch = assemble(trimmed_ch)
        }

        indexed_ref_ch = index_reference(ref_ch)
        fai_result = samtools_faidx(ref_ch)
        bam_ch = map_reads(trimmed_ch, indexed_ref_ch)
        coverage_result = coverage_plot(bam_ch)

        bam_meta_ch = bam_ch.map { label, bam, bai -> tuple([id: label], bam, [], []) }
        fai_meta_ch = fai_result.map { ref, fai -> tuple([id: 'ref'], ref, fai) }.first()
        mpileup_result = BCFTOOLS_MPILEUP(bam_meta_ch, fai_meta_ch, params.save_mpileup)
        vcf_input_ch = mpileup_result.vcf.join(mpileup_result.tbi)
        vcf_result = BCFTOOLS_CALL(vcf_input_ch, [], [], [])

    publish:
        qc = qc_result
        trimmed_reads = trimmed_ch
        mapping = bam_ch
        coverage = coverage_result
        variants = vcf_result.vcf
}

output {
    qc {
        path 'qc'
    }
    trimmed_reads {
        path 'trimmed_reads'
    }
    mapping {
        path 'mapping'
    }
    coverage {
        path 'coverage'
    }
    variants {
        path 'variants'
    }
}
