version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem


struct dragmapResources {
    String modules
    String hashTable
}

workflow dragmap {

    input {
        File fastqR1
        File? fastqR2
        String outputFileNamePrefix
        Int numChunk = 1
        Boolean doUMIextract = false
        Boolean doTrim = false
        String reference
        Int? numReads
    }

    parameter_meta {
        fastqR1: "Fastq file for read 1"
        fastqR2: "Fastq file for read 2"
        outputFileNamePrefix: "Prefix for output files"
        numChunk: "Number of chunks to split fastq file [1, no splitting]"
        doUMIextract: "If true, UMI will be extracted before alignment [false]"
        doTrim: "If true, adapters will be trimmed before alignment [false]"
        reference: "The genome reference build. For example: hg19, hg38, mm10"
        numReads: "Number of reads"
    }

    Map[String,dragmapResources] resourceByGenome = { 
        "hg19": {
            "modules": "samtools/1.9 dragmap/1.2.1 dragmap-hash-table-hg19/p13", 
            "hashTable": "$DRAGMAP_HASH_TABLE_HG19_ROOT/hg19"
        },
        "hg38": {
            "modules": "samtools/1.9 dragmap/1.2.1 dragmap-hash-table-hg38/p12", 
            "hashTable": "$DRAGMAP_HASH_TABLE_HG38_ROOT/hg38"
        },
        "mm10": {
            "modules": "samtools/1.9 dragmap/1.2.1 dragmap-hash-table-mm10/p6", 
            "hashTable": "$DRAGMAP_HASH_TABLE_MM10_ROOT/mm10"
        }
    }

    String dragmapModules = resourceByGenome[reference].modules
    String dragmapHashTable = resourceByGenome[reference].hashTable

    if (numChunk > 1) {
        call bwaMem.countChunkSize as countChunkSize {
            input:
            fastqR1 = fastqR1,
            numChunk = numChunk,
            numReads = numReads
        }
    
        call bwaMem.slicer as slicerR1 { 
            input: 
            fastqR = fastqR1,
            chunkSize = countChunkSize.chunkSize
        }
        if (defined(fastqR2)) {
            # workaround for converting File? to File
            File fastqR2File = select_all([fastqR2])[0]
            call bwaMem.slicer as slicerR2 {
                input:
                fastqR = fastqR2File,
                chunkSize = countChunkSize.chunkSize
            }
        }
    }

    Array[File] fastq1 = select_first([slicerR1.chunkFastq, [fastqR1]])

    if(defined(fastqR2)) {
        Array[File?] fastq2 = select_first([slicerR2.chunkFastq, [fastqR2]])
        Array[Pair[File,File?]] pairedFastqs = zip(fastq1,fastq2)
    }

    if(!defined(fastqR2)) {
        Array[Pair[File,File?]] singleFastqs = cross(fastq1,[fastqR2])
    }

    Array[Pair[File,File?]] outputs = select_first([pairedFastqs, singleFastqs])

    scatter (p in outputs) {
        if (doUMIextract) {
            call bwaMem.extractUMIs as extractUMIs { 
                input:
                fastq1 = p.left,
                fastq2 = p.right,
            }
        }

        if (doTrim) {
            call bwaMem.adapterTrimming as adapterTrimming { 
                input:
                fastqR1 = select_first([extractUMIs.fastqR1, p.left]),
                fastqR2 = if (defined(fastqR2)) then select_first([extractUMIs.fastqR2, p.right]) else fastqR2,
            }
        }

        call runDragmap  { 
            input: 
            read1s = select_first([adapterTrimming.resultR1, extractUMIs.fastqR1, p.left]),
            read2s = if (defined(fastqR2)) then select_first([adapterTrimming.resultR2, extractUMIs.fastqR2, p.right]) else fastqR2,
            modules = dragmapModules,
            dragmapHashTable = dragmapHashTable
        }    
        
        call addReadGroups {
            input: 
            outputBam = runDragmap.outputBam
        }
    }

    call bwaMem.bamMerge as bamMerge {
        input:
        bams = addReadGroups.outputReadGroupBam,
        outputFileNamePrefix = outputFileNamePrefix
    }

    call bwaMem.indexBam as indexBam { 
        input: 
        inputBam = bamMerge.outputMergedBam
    }

    if (doTrim) {
        call bwaMem.adapterTrimmingLog as adapterTrimmingLog {
            input:
            inputLogs = select_all(adapterTrimming.log),
            outputFileNamePrefix = outputFileNamePrefix,
            numChunk = numChunk,
            singleEnded = if (defined(fastqR2)) then false else true
        }
    }

    meta {
        author: "Xuemei Luo and Muna Mohamed"
        email: "xuemei.luo@oicr.on.ca, mmohamed@oicr.on.ca"
        description: "This workflow aligns sequence data provided as fastq files using Dragmap (an open source Dragen mapper/aligner). The alignment is completed using a hash table of a reference genome. The workflow imports the bwaMem workflow, to provide options to remove 5' umi sequences and trim off 3' sequencing adapters prior to alignment. Using the bwaMem workflow, this workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single BAM file to decrease workflow runtime. Read-group information must be provided."
        dependencies: [
            {
                name: "dragmap/1.2.1",
                url: "https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.2.1.tar.gz"
            },
            {
                name: "samtools/1.9",
                url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
            },
            {
                name: "picard/3.1.0",
                url: "https://github.com/broadinstitute/picard/archive/refs/tags/3.1.0.tar.gz"
            },
            {
                name: "cutadapt/1.8.3",
                url: "https://cutadapt.readthedocs.io/en/v1.8.3/"
            },
            {
                name: "slicer/0.3.0",
                url: "https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz"
            },
            {
                name: "python/3.7",
                url: "https://www.python.org"
            },      
            {
                name: "barcodex-rs/0.1.2",
                url: "https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz"
            },
            {
                name: "rust/1.2",
                url: "https://www.rust-lang.org/tools/install"
            },
            { 
                name: "gsi software modules: samtools/1.9 dragmap/1.2.1",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            { 
                name: "gsi hg38 modules: dragmap-hash-table-hg38/p12",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            {
                name: "gsi hg19 modules: dragmap-hash-table-hg19/p13",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            },
            {
                name: "gsi mm10 modules: dragmap-hash-table-mm10/p6",
                url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
            }
        ]
        output_meta: {
            dragmapBam: "Merged BAM file aligned to genome with Dragmap",
            dragmapIndex: "Output index file for the merged BAM file",
            log: "A summary log file for adapter trimming",
            cutAdaptAllLogs: "A file with the logs from the adapter trimming of each fastq chunk"
        }
    }

    output {
        File dragmapBam = bamMerge.outputMergedBam
        File dragmapIndex = indexBam.outputBai
        File? log = adapterTrimmingLog.summaryLog
        File? cutAdaptAllLogs = adapterTrimmingLog.allLogs
    }

}


task runDragmap {

    input {
        File read1s
        File? read2s
        String dragmapHashTable
        String? addParam
        String modules
        Int jobMemory = 50
        Int timeout = 96
    }

    parameter_meta {
        read1s: "Fastq file for read 1"
        read2s: "Fastq file for read 2"
        dragmapHashTable: "The hash table prepared by Dragmap using the appropriate reference genome. Used to align sample with Dragmap"
        addParam: "Additional Dragmap parameters"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for the job"
        timeout: "Hours until task timeout"
    }

    String resultBam = "~{basename(read1s)}.bam"
    String tmpDir = "tmp/"

    command <<<
        set -euo pipefail
        mkdir -p ~{tmpDir}
        
        dragen-os \
            -r ~{dragmapHashTable} \
            -1 ~{read1s} \
            ~{if (defined(read2s)) then "-2 ~{read2s}" else ""} \
            ~{addParam} \
        | \
        samtools view -b - \
        | \
        samtools sort -O bam -T ~{tmpDir} -o ~{resultBam} - 
    >>>

    runtime {
        modules: "~{modules}"
        memory:  "~{jobMemory} GB"
        timeout: "~{timeout}"
    }  
    
    output {
        File outputBam = "~{resultBam}"
    }

    meta {
        output_meta: {
            outputBam: "Output BAM file aligned to the appropriate genome. Does not have read-group information."
        }
    }

}


task addReadGroups {

    input {
        File outputBam
        String readGroups
        String modules = "picard/3.1.0"
        Int jobMemory = 12
        Int timeout = 5
    }

    parameter_meta {
        outputBam: "Output BAM file aligned to the appropriate genome. Does not have read-group information."
        readGroups: "The read-group information to be added into the BAM file header"
        modules: "Required environment modules to run the task"
        jobMemory: "Memory allocated for the job"
        timeout: "Hours until task timeout"
    }
    
    String fileNamePrefix = basename(outputBam, ".bam")
    String resultReadGroupBam = "~{fileNamePrefix}_readgroup.bam"

    command <<<
        set -euo pipefail

        # An array containing potential read-group fields
        fieldNames=("ID:" "LB:" "PL:" "PU:" "SM:" "CN:" "DS:" "DT:" "FO:" "KS:" "PG:" "PI:" "PM:")

        # Removes the start of the readGroups string (@RG\\t) and replaces tabs with commas        
        strippedString=$(echo -e ~{readGroups} | sed -E "s/@RG\\t//" | sed -E "s/\\t/,/g")
        
        # Declares an associative array and appends the values of the fields present in readGroups string
        declare -A fieldsArray
        
        for field in "${fieldNames[@]}"; do
            if [[ ${strippedString} == *${field}* ]]; then 
                fieldsArray[${field}]=$(echo "${strippedString}" | awk -F "${field}" '{print $2}' | cut -d ',' -f 1) 
            fi
        done

        java -jar $PICARD_ROOT/picard.jar AddOrReplaceReadGroups \
            I= ~{outputBam} \
            O= ~{resultReadGroupBam} \
            $( [[ -v fieldsArray["ID:"] ]] && echo "RGID=${fieldsArray["ID:"]}" || : ) \
            $( [[ -v fieldsArray["LB:"] ]] && echo "RGLB=${fieldsArray["LB:"]}" || : ) \
            $( [[ -v fieldsArray["PL:"] ]] && echo "RGPL=${fieldsArray["PL:"]}" || : ) \
            $( [[ -v fieldsArray["PU:"] ]] && echo "RGPU=${fieldsArray["PU:"]}" || : ) \
            $( [[ -v fieldsArray["SM:"] ]] && echo "RGSM=${fieldsArray["SM:"]}" || : ) \
            $( [[ -v fieldsArray["CN:"] ]] && echo "RGCN=${fieldsArray["CN:"]}" || : ) \
            $( [[ -v fieldsArray["DS:"] ]] && echo "RGDS=${fieldsArray["DS:"]}" || : ) \
            $( [[ -v fieldsArray["DT:"] ]] && echo "RGDT=${fieldsArray["DT:"]}" || : ) \
            $( [[ -v fieldsArray["FO:"] ]] && echo "RGFO=${fieldsArray["FO:"]}" || : ) \
            $( [[ -v fieldsArray["KS:"] ]] && echo "RGKS=${fieldsArray["KS:"]}" || : ) \
            $( [[ -v fieldsArray["PG:"] ]] && echo "RGPG=${fieldsArray["PG:"]}" || : ) \
            $( [[ -v fieldsArray["PI:"] ]] && echo "RGPI=${fieldsArray["PI:"]}" || : ) \
            $( [[ -v fieldsArray["PM:"] ]] && echo "RGPM=${fieldsArray["PM:"]}" || : )
    >>>

    runtime {
        modules: "~{modules}"
        memory:  "~{jobMemory} GB"
        timeout: "~{timeout}"
    }  
    
    output {
        File outputReadGroupBam = "~{resultReadGroupBam}"
    }

    meta {
        output_meta: {
            outputReadGroupBam: "Output BAM file aligned to the appropriate genome with read-group information in the header"
        }
    }

}
