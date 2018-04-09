Usage
-----

### Quick Start Example

To run against a Salmonella database with approximately 16 threads:
```
nextflow bear-wgs.nf -w "/path/to/working/directory" --reference "./containers/data/references/CP016573.1_Salmonella_Heidelberg.fasta" --output "/path/to/output/folder" --threads 4 --ploidy 1 --db Salmonella -resume
```

### Display Help Message

The `help` parameter displays the available options and commands.
```
$ nextflow run auir.nf --help
```

### File Inputs

#### Set custom sequence data

The `reads` parameter accepts sequence files in standard fastq and gz format.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq"
```

### File Outputs

#### Set output directory

The `output` parameter writes output files to the desired directory.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --output "test"
```

### Other Parameters

#### Set custom trimming parameters

```
$ nextflow run auir.nf \
    --reads "data/raw/*_R{1,2}.fastq" \
    --leading 3 \
    --trailing 3 \
    --minlen 36 \
    --slidingwindow 4 \
    --adapters "data/adapters/nextera.fa" \
    --output "test"
```

#### Set number of threads

The `threads` parameter specifies the number of threads to use for each process.
```
$ nextflow run auir.nf --reads "data/raw/*_R{1,2}.fastq" --threads 16 --output "test"
```
