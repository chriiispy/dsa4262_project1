// parameters
params.refFa = '/path/to/ref.fa'
params.refGtf = '/path/to/ref.gtf'
params.reads = '/path/to/reads.fq'
params.outdir = 'results'

// minimap2
process MINIMAP2_ALIGN {
    input:
        path refFa
        path reads
    output:
        path "aligned_reads.sam"

    """
        minimap2 -ax splice -uf -k14 "${refFa}" "${reads}" > aligned_reads.sam
    """;
}

// sam to bam
process SAM_TO_BAM {
    publishDir params.outdir
    
    input:
        path reads_sam
    output:
        path "aligned_reads.bam"

    """
        samtools view -b $reads_sam > aligned_reads.bam
    """;
}

// bambu
process BAMBU {
    publishDir params.outdir

    input:
        path refFa
        path refGtf
        path reads_bam
    output:
        path "counts_transcript.txt"
        path "counts_gene.txt"
        path "extended_annotations.gtf"

    """
        #!/usr/bin/env Rscript --vanilla
        library(bambu)
        show("$refGtf")
        annotations <- prepareAnnotations("$refGtf")
        se <- bambu(reads = "$reads_bam", annotations = annotations, genome = "$refFa", ncore = 1)
        writeBambuOutput(se, path = "./")
    """;
}


// WORKFLOW
workflow {
    MINIMAP2_ALIGN(params.refFa, params.reads)
    SAM_TO_BAM(MINIMAP2_ALIGN.out)
    BAMBU(params.refFa, params.refGtf, SAM_TO_BAM.out)
}
