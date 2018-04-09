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

/*
if( params.db ) {
    db = params.db
    if(
}
*/

threads = params.threads
db = params.db

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
        set dataset_id, file("${dataset_id}.1P.fastq"), file("${dataset_id}.2P.fastq") into (paired_fastq_alignment, paired_fastq_assembly, paired_fastq_ariba)
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
        set dataset_id, file("${dataset_id}.aligned.sam") into prokka_sam
        set dataset_id, file("${dataset_id}.aligned.sam") into freebayes_sam
    
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
        set dataset_id, file("${dataset_id}.contigs.fa") into (spades_contigs_ariba, spades_contigs_phaster)
        set dataset_id, file("${dataset_id}.graph.fastg") into (spades_graphs)
        set dataset_id, file("${dataset_id}.contigs.paths") into (spades_paths)
    
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
    mv output/contigs.paths ./${dataset_id}.contigs.paths
    """
}

process FindAribaMarkers {
    tag { dataset_id }
    
    publishDir "${params.output}/AribaFiles", mode: "symlink"
    
    input:
        set dataset_id, file(forward), file(reverse) from paired_fastq_ariba
    
    output:
        set dataset_id, file("${dataset_id}.card.assemblies.fa.gz") into (ariba_card_assemblies)
        set dataset_id, file("${dataset_id}.megares.assemblies.fa.gz") into (ariba_megares_assemblies)
        set dataset_id, file("${dataset_id}.plasmidfinder.assemblies.fa.gz") into (ariba_plasmidfinder_assemblies)
        set dataset_id, file("${dataset_id}.virulencefinder.assemblies.fa.gz") into (ariba_virulencefinder_assemblies)
        file("${dataset_id}.card.report.tsv") into (ariba_card_reports)
        file("${dataset_id}.megares.report.tsv") into (ariba_megares_reports)
        file("${dataset_id}.plasmidfinder.report.tsv") into (ariba_plasmidfinder_reports)
        file("${dataset_id}.virulencefinder.report.tsv") into (ariba_virulencefinder_reports)
        
        """
        ariba run \
        --threads ${threads} \
        /opt/out.card.prepareref \
        ${forward} \
        ${reverse} ariba_out_card
        mv ariba_out_card/assemblies.fa.gz ./${dataset_id}.card.assemblies.fa.gz
        mv ariba_out_card/report.tsv ./${dataset_id}.card.report.tsv
        
        ariba run \
        --threads ${threads} \
        /opt/out.megares.prepareref \
        ${forward} \
        ${reverse} ariba_out_megares
        mv ariba_out_megares/assemblies.fa.gz ./${dataset_id}.megares.assemblies.fa.gz
        mv ariba_out_megares/report.tsv ./${dataset_id}.megares.report.tsv
        
        ariba run \
        --threads ${threads} \
        /opt/out.plasmidfinder.prepareref \
        ${forward} \
        ${reverse} ariba_out_plasmidfinder
        mv ariba_out_plasmidfinder/assemblies.fa.gz ./${dataset_id}.plasmidfinder.assemblies.fa.gz
        mv ariba_out_plasmidfinder/report.tsv ./${dataset_id}.plasmidfinder.report.tsv
        
        ariba run \
        --threads ${threads} \
        /opt/out.virulencefinder.prepareref \
        ${forward} \
        ${reverse} ariba_out_virulencefinder
        mv ariba_out_virulencefinder/assemblies.fa.gz ./${dataset_id}.virulencefinder.assemblies.fa.gz
        mv ariba_out_virulencefinder/report.tsv ./${dataset_id}.virulencefinder.report.tsv
        """
}

process SummarizeAribaReports {
    tag { "ARIBA Reports" }
    
    publishDir "${params.output}/AribaFiles", mode: "symlink"
    
    input:
        file card_reports from ariba_card_reports.toList()
        file megares_reports from ariba_megares_reports.toList()
        file plasmidfinder_reports from ariba_plasmidfinder_reports.toList()
        file virulencefinder_reports from ariba_virulencefinder_reports.toList()
    output:
        set file("ariba.card.summary.csv"), file("ariba.megares.summary.csv"), file("ariba.plasmidfinder.summary.csv"), file("ariba.virulencefinder.summary.csv") into (ariba_summary_files)
    
    """
    #!/bin/bash
    ariba summary ariba.card.summary $card_reports
    ariba summary ariba.megares.summary $megares_reports
    ariba summary ariba.plasmidfinder.summary $plasmidfinder_reports
    ariba summary ariba.virulencefinder.summary $virulencefinder_reports
    """
}

process SeparatePlasmidContigs {
    tag { "SPAdes Assemblies" }
    
    publishDir "${params.output}/AssembledFiles", mode: "symlink"
    
    input:
        set dataset_id, file(spades_contigs) from spades_contigs_ariba
        set dataset_id, file(ariba_plasmid_assemblies) from ariba_plasmidfinder_assemblies
    output:
        set dataset_id, file("${dataset_id}.genome.contigs.fa") into (genome_contigs)
        set dataset_id, file("${dataset_id}.plasmid.contigs.fa") into (plasmid_contigs)
    
    """
    #!/bin/bash
    if [[ ! -s ${ariba_plasmid_assemblies} ]]; then
        gzip -d -c $ariba_plasmid_assemblies > plasmid_unzipped.fa
        mummer -b $spades_contigs plasmid_unzipped.fa > ${dataset_id}.plasmid.alignment.out
        separate_plasmid_contigs.py $spades_contigs ${dataset_id}.plasmid.alignment.out ${dataset_id}.plasmid.contigs.fa ${dataset_id}.genome.contigs.fa
    else
        touch ${dataset_id}.plasmid.contigs.fa
        cp $spades_contigs ${dataset_id}.genome.contigs.fa
    fi
    """
}

if( params.phage ) {
    process FindPhages {
        tag { dataset_id }
        
        input:
            set dataset_id, file(spades_contigs) from spades_contigs_phaster
        output:
            set dataset_id, file("${dataset_id}.phaster.results.zip") into (phaster_results)
        
        """
        mkdir temp_in
        mkdir temp_out
        
        cp $spades_contigs temp_in/query.fa
        phaster_query.py temp_in temp_out
        mv temp_out/zipfiles/*.zip ./${dataset_id}.phaster.results.zip
        
        rm temp_in/*
        rmdir temp_in
        rmdir temp_out/zipfiles
        rm temp_out/*
        rmdir temp_out
        """
    }
}

\*
// Need to generate consensus here for prokka
process AnnotateGenomeAlignments {
    tag { dataset_id }
    
    publishDir "${params.output}/AnnotatedGenomeAlignments", mode: "symlink"
    
    input:
        set dataset_id, file(sam_file) from prokka_sam
    output:
        set dataset_id, file("${dataset_id}.*") into (annotated_genome_alignments)
    
    """
    #!/bin/bash
    prokka --outdir annotations --usegenus --genus $db --cpus $threads --prefix ${dataset_id}.alignment $sam_file
    mv annotations/* .
    """
}
\*

process AnnotateGenomeAssemblies {
    tag { dataset_id }
    
    publishDir "${params.output}/AnnotatedGenomeAssemblies", mode: "symlink"
    
    input:
        set dataset_id, file(genome_contig_file) from genome_contigs
    output:
        set dataset_id, file("${dataset_id}.*") into (annotated_genome_assemblies)
    
    """
    #!/bin/bash
    prokka --outdir annotations --usegenus --genus $db --cpus $threads --prefix ${dataset_id}.genome $genome_contig_file
    mv annotations/* .
    """
}

process AnnotatePlasmidAssemblies {
    tag { "Plasmid Contigs" }
    
    publishDir "${params.output}/AnnotatedPlasmidAssemblies", mode: "symlink"
    
    input:
        set dataset_id, file(plasmid_contig_file) from plasmid_contigs
    output:
        set dataset_id, file("${dataset_id}.*") into (annotated_plasmid_assemblies)
    
    """
    #!/bin/bash
    if [[ -s ${plasmid_contig_file} ]]; then
        cat $plasmid_contig_file
        prokka --outdir annotations --usegenus --genus Plasmid --cpus $threads --prefix ${dataset_id}.plasmid $plasmid_contig_file
        mv annotations/* .
    else
        touch ${dataset_id}.empty.no.annotations
    fi
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

/*
def db_error(def input) {
    println ""
    println "[params.db] invalid argument: '" + db "' : Choose one of [Salmonella, Efaecalis]"
    println ""
    return 1
}
*/

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
    println "    --db            STR      database to use for annotation [Salmonella, Efaecalis]"
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
    println "    --phage         FLAG     query assemblies to find phages, false by default"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}





















