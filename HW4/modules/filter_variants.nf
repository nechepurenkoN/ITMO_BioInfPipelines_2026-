process filter_variants {
    input:
        tuple val(meta), path(vcf)
    output:
        tuple val(meta), path("${meta.id}_filtered.vcf.gz")
    stub:
    """
    touch ${meta.id}_filtered.vcf.gz
    """
    script:
    """
    bcftools filter -e 'QUAL<20 || DP<10' -O z -o ${meta.id}_filtered.vcf.gz $vcf
    """
}
