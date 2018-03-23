#!/usr/bin/env nextflow

if( params.help ) {
    return help()
}

if( !nextflow.version.matches('0.25+') ) {
    return nextflow_version_error()
}

if( params.reference ) {
    reference = file(params.reference)
    if( !reference.exists() ) return reference_error(reference)
}

if( params.adapters ) {
    adapters = file(params.adapters)
    if( !adapters.exists() ) return adapter_error(adapters)
}

if( params.fqc_adapters ) {
    fqc_adapters = file(params.fqc_adapters)
    if( !fqc_adapters.exists() ) return fastqc_error(fqc_adapters)
}

threads = params.threads

// Trimmomatic options
leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { return fastq_error(params.reads) }
    .into { read_pairs; fastqc_pairs }

process FastQC {
    tag { dataset_id }
    
    publishDir "${params.output}/QualityMetrics", mode: "symlink"
    
    input:
        set dataset_id, file(forward), file(reverse) from fastqc_pairs
    
    output:
        set dataset_id, file("*_fastqc.zip") into (fastqc_logs)
    
    """
    mkdir output
    fastqc -f fastq ${forward} ${reverse} -t ${threads} -o output
    mv output/*.zip .
    """
}

process QualityControl {
    tag { dataset_id }
    
    publishDir "${params.output}/QualityControlOutput", mode: "symlink",
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
            else if(filename.indexOf(".log") > 0) "Log/$filename"
            else {}
        }
    
    input:
        set dataset_id, file(forward), file(reverse) from read_pairs
    
    output:
        set dataset_id, file("${dataset_id}.1P.fastq"), file("${dataset_id}.2P.fastq") into (paired_fastq_alignment, paired_fastq_assembly)
        set dataset_id, file("${dataset_id}.1U.fastq"), file("${dataset_id}.2U.fastq") into (unpaired_fastq)
        set dataset_id, file("${dataset_id}.trimmomatic.stats.log") into (trimmomatic_logs)
    
    """
    /usr/lib/jvm/java-7-openjdk-amd64/bin/java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar \
        PE \
        -threads ${threads} \
        $forward $reverse -baseout ${dataset_id} \
        ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
        LEADING:${leading} \
        TRAILING:${trailing} \
        SLIDINGWINDOW:${slidingwindow} \
        MINLEN:${minlen} \
        2> ${dataset_id}.trimmomatic.stats.log
    
    mv ${dataset_id}_1P ${dataset_id}.1P.fastq
    mv ${dataset_id}_2P ${dataset_id}.2P.fastq
    mv ${dataset_id}_1U ${dataset_id}.1U.fastq
    mv ${dataset_id}_2U ${dataset_id}.2U.fastq
    """
}

process BuildReferenceIndex {
    tag { reference.baseName }
    
    publishDir "${params.output}/BuildReferenceIndex", mode: "symlink"
    
    input:
        file(reference)
    
    output:
        file 'genome.index*' into index
    
    """
    bwa index -p genome.index ${reference}
    """
}

process AlignReadsToGenome {
    tag { dataset_id }
    
    publishDir "${params.output}/AlignedFiles", mode: "symlink"
    
    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq_alignment
        file ref_index from index.first()
    
    output:
        set dataset_id, file("${dataset_id}.aligned.sam") into genome_sam
    
    """
    bwa mem genome.index ${forward} ${reverse} -t ${threads} > ${dataset_id}.aligned.sam
    """
}

process AssembleReads {
    tag { dataset_id }
    
    publishDir "${params.output}/AssembledFiles", mode: "symlink"
    
    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq_assembly
    
    output:
        set dataset_id, file("${dataset_id}.contigs.fa") into (spades_contigs)
        set dataset_id, file("${dataset_id}.graph.fastg") into (spades_graphs)
    
    script:
    """
    spades.py \
        -t ${threads} \
        --cov-cutoff auto \
        -1 ${forward} \
        -2 ${reverse} \
        -o output
    mv output/contigs.fasta ./${dataset_id}.contigs.fa
    mv output/assembly_graph.fastg ./${dataset_id}.graph.fastg
    """
}

def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastqc_error(def input) {
    println ""
    println "[params.fqc_adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def reference_error(def input) {
    println ""
    println "[params.reference] fail to open: '" + reference + "' : No such file or directory"
    println ""
    return 1
}

def help() {
    println ""
    println "Program: bear-wgs"
    println "Version: $workflow.repository - $workflow.revision [$workflow.commitId]"
    println "Documentation: https://github.com/lakinsm/bear-wgs"
    println "Contact: Steven Lakin <steven.m.lakin@gmail.com>"
    println ""
    println "Usage:    nextflow bear-wgs.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --adapters      STR      path to FASTA formatted adapter sequences"
    println "    --reference     STR      path to reference genome FASTA file"
    println "    --output        STR      directory to write process outputs to"
    println ""
    println "Trimming options:"
    println ""
    println "    --leading       INT      cut bases off the start of a read, if below a threshold quality"
    println "    --minlen        INT      drop the read if it is below a specified length"
    println "    --slidingwindow INT      perform sw trimming, cutting once the average quality within the window falls below a threshold"
    println "    --trailing      INT      cut bases off the end of a read, if below a threshold quality"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of threads to use for each process"
    println "    --ploidy        INT      genome copy number"
    println "    --min-alt-count INT      requires this number of observations supporting an alternate allele"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}





















