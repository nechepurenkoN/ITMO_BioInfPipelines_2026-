include { run_qc; trimm; index_reference; samtools_faidx; map_reads; coverage_plot } from './modules/hw2_processes'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CALL } from './modules/nf-core/bcftools/call/main'
include { filter_variants } from './modules/filter_variants'

workflow {
    main:
        samples_ch = Channel.fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, row.group, file(row.fastq_1), file(row.fastq_2)) }

        reads_ch = samples_ch.map { id, group, r1, r2 -> tuple(id, r1, r2) }
        trimmed_ch = trimm(reads_ch)

        reads_tagged = reads_ch.map   { id, r1, r2 -> tuple("initial_${id}", r1, r2) }
        trimmed_tagged = trimmed_ch.map { id, r1, r2 -> tuple("trimmed_${id}", r1, r2) }
        qc_result = run_qc(reads_tagged.mix(trimmed_tagged))

        ref_ch = Channel.value(file(params.reference))
        indexed_ref_ch = index_reference(ref_ch).first()
        fai_result = samtools_faidx(ref_ch)
        fai_meta_ch = fai_result.map { ref, fai -> tuple([id: 'ref'], ref, fai) }.first()

        group_ch = samples_ch.map { id, group, r1, r2 -> tuple(id, group) }
        trimmed_with_group = trimmed_ch.join(group_ch)

        branched = trimmed_with_group.branch {
            virus_a: it[3] == 'virus_a'
            virus_b: it[3] == 'virus_b'
        }
        reads_a = branched.virus_a.map { id, r1, r2, group -> tuple(id, r1, r2) }
        reads_b = branched.virus_b.map { id, r1, r2, group -> tuple(id, r1, r2) }

        bam_a_ch = map_reads(reads_a, indexed_ref_ch)
        bam_b_ch = map_reads(reads_b, indexed_ref_ch)
        coverage_result = coverage_plot(bam_a_ch.mix(bam_b_ch))

        vcf_a_ch = call_variants(bam_a_ch, fai_meta_ch).vcf
        vcf_b_ch = call_variants(bam_b_ch, fai_meta_ch).vcf

        all_vcfs_ch = vcf_a_ch.mix(vcf_b_ch)

        filtered_ch = filter_variants(all_vcfs_ch)

    publish:
        qc  = qc_result
        trimmed_reads = trimmed_ch
        mapping = bam_a_ch.mix(bam_b_ch)
        coverage = coverage_result
        variants = all_vcfs_ch.map { meta, vcf -> vcf }
        filtered = filtered_ch.map { meta, vcf -> vcf }
}

workflow call_variants {
    take:
        bam_ch
        fai_meta_ch
    main:
        bam_meta_ch    = bam_ch.map { label, bam, bai -> tuple([id: label], bam, [], []) }
        mpileup_result = BCFTOOLS_MPILEUP(bam_meta_ch, fai_meta_ch, params.save_mpileup)
        vcf_input_ch   = mpileup_result.vcf.join(mpileup_result.tbi)
        vcf_result     = BCFTOOLS_CALL(vcf_input_ch, [], [], [])
    emit:
        vcf = vcf_result.vcf
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
    filtered {
        path 'filtered_variants'
    }
}
