/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_isoxpressor_pipeline'
include { getGenomeAttribute      } from '../subworkflows/local/utils_nfcore_isoxpressor_pipeline'



include { FASTQ_ALIGN_STAR         } from '../subworkflows/nf-core/fastq_align_star'
include { FASTQ_TRIM_FASTP_FASTQC  } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'

include { STAR_GENOMEGENERATE       } from '../modules/nf-core/star/genomegenerate/main'
include { HTSEQ_COUNT               } from '../modules/nf-core/htseq/count/main'
include { ISOXPRESSOR_ISOSEGMENTER  } from '../modules/local/isoxpressor/isosegmenter/main'
include { ISO2GTF                   } from '../modules/local/iso2gtf/main'
include { ISOSORT                   } from '../modules/local/isosort/main.nf'
include { SEQKIT_SPLIT              } from '../modules/local/seqkit/split/main'



params.fasta            = getGenomeAttribute('fasta')
params.star_index       = getGenomeAttribute('star')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOXPRESSOR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input


    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_fasta=Channel.fromPath(params.fasta)



    if(params.skip_isosegmenter){
        ch_gtf=Channel.value(file(params.gtf))
    }else {

        SEQKIT_SPLIT(ch_fasta)
        ch_chromosomes=SEQKIT_SPLIT.out.chromosomes.flatten().map{ n ->[n.baseName,n ]}

        ISOXPRESSOR_ISOSEGMENTER(ch_chromosomes)
        ch_isochores=ISOXPRESSOR_ISOSEGMENTER.out.isochores

        header = Channel.of("Chromosome,Isochore Class,GC Level,Isochore Start,Isochore End,Isochore Size,\n")

        ch_isochores
        .map { file ->
            def name = file.baseName
            def content = file.splitCsv(skip: 1)
            content.collect { row -> "${name},${row[3]},${row[4]},${row[0]},${row[1]},${row[2]},\n" }
        }
        .flatten()
        .set{ch_read_head_temp}

        header.concat(ch_read_head_temp)
            .collectFile(name: "unsorted.csv", sort:{})
            .set { ch_read_head_unsorted }


        ISOSORT(ch_read_head_unsorted)
        ISO2GTF(ISOSORT.out.csv)


        ch_gtf=ISO2GTF.out.gtf

    }

    if(params.skip_genome_generate){

        ch_star_index=Channel.fromPath(params.star_index)
        ch_star_index=ch_star_index.map{ [ [:], it ] }

    }else{

        STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } )
        ch_star_index = STAR_GENOMEGENERATE.out.index
    }



    //
    // MODULE: FASTQ_TRIM_FASTP_FASTQC
    //

    FASTQ_TRIM_FASTP_FASTQC (
        ch_samplesheet,
        params.adapter_fasta ?: [],
        params.save_trimmed_fail,
        params.discard_trimmed_pass,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)
    ch_filtered_reads = FASTQ_TRIM_FASTP_FASTQC.out.reads


    //
    // SUBWORKFLOW: align with STAR
    // alignments and run BAM_SORT_STATS_SAMTOOLS for each
    //

    FASTQ_ALIGN_STAR(
        ch_filtered_reads,
        ch_star_index,
        ch_gtf.map { [ [:], it ] },
        params.star_ignore_sjdbgtf,
        '',
        '',
        ch_fasta.map { [ [:], it ] },
        ''
    )

    ch_genome_bam              = FASTQ_ALIGN_STAR.out.bam
    ch_genome_bam_index        = FASTQ_ALIGN_STAR.out.bai
    ch_versions                = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

    //
    // Create counts with htseq
    //

    ch_bam_bai = ch_genome_bam.join(ch_genome_bam_index)
    ch_grouped = ch_bam_bai.groupTuple().collect()
    HTSEQ_COUNT(ch_grouped,ch_gtf.map{[[:],it]})

    //
    // Create reads.csv file
    //

    ch_counts=HTSEQ_COUNT.out.txt

    ch_counts
    .tap{count_header}
    .map{it[0]}
    .toList()
    .tap{data}
    .map{it[1]}
    .splitText()
    .groupTuple()
    .toSortedList()
    .map{it[1]}

    count_header.concat(data)
    .collectFile(name:"all_counts.csv")

    // TODO: Join with reads_head while omiting the last 5 lines




    //
    // Create condition file from the samplesheet
    //
    Channel.fromPath(params.input)
    .splitCsv(header:true,skip:1)
    .map { row -> "${row[0]},${row[1]},${row[2]}\n" }
    .collectFile(name:"conditions.csv")

    //
    // Run analysis.py
    //
    // TODO: Add the analysis module


    //
    // Create the profiles
    //
    // TODO: Add the graphs module
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'isoxpressor_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
