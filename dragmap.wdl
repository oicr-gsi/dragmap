version 1.0


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
        String readGroups
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
        readGroups: "The read-group information to be added into the BAM file header"
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

    # Validating the read-group information
    call readGroupCheck {  
        input: 
        readGroups = readGroups 
    }

    if (numChunk > 1) {
        call countChunkSize {
            input:
            fastqR1 = fastqR1,
            numChunk = numChunk,
            numReads = numReads
        }
    
        call slicer as slicerR1 { 
            input: 
            fastqR = fastqR1,
            chunkSize = countChunkSize.chunkSize
        }
        if (defined(fastqR2)) {
            # workaround for converting File? to File
            File fastqR2File = select_all([fastqR2])[0]
            call slicer as slicerR2 {
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
            call extractUMIs { 
                input:
                fastq1 = p.left,
                fastq2 = p.right
            }
        }

        if (doTrim) {
            call adapterTrimming { 
                input:
                fastqR1 = select_first([extractUMIs.fastqR1, p.left]),
                fastqR2 = if (defined(fastqR2)) then select_first([extractUMIs.fastqR2, p.right]) else fastqR2
            }
        }

        call runDragmap { 
            input: 
            read1s = select_first([adapterTrimming.resultR1, extractUMIs.fastqR1, p.left]),
            read2s = if (defined(fastqR2)) then select_first([adapterTrimming.resultR2, extractUMIs.fastqR2, p.right]) else fastqR2,
            modules = dragmapModules,
            dragmapHashTable = dragmapHashTable
        }
    }

    call bamMerge {
        input:
        outputBams = runDragmap.outputBam,
        outputFileNamePrefix = outputFileNamePrefix
    }

    call addReadGroups {
        input:
        mergedBam = bamMerge.outputMergedBam,
        readGroups = readGroups
    }

    if (doTrim) {
        call adapterTrimmingLog {
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
        description: "This workflow aligns sequence data provided as fastq files using Dragmap (an open source Dragen mapper/aligner). The alignment is completed using a hash table of a reference genome. The workflow borrows functions from the bwaMem workflow, to provide options to remove 5' umi sequences and trim off 3' sequencing adapters prior to alignment. Using these functions, this workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single BAM file to decrease workflow runtime. Read-group information must be provided."
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
        File dragmapBam = addReadGroups.outputReadGroupBam
        File dragmapIndex = addReadGroups.outputReadGroupBai
        File? log = adapterTrimmingLog.summaryLog
        File? cutAdaptAllLogs = adapterTrimmingLog.allLogs
    }

}

task readGroupCheck { 

    input { 
        String readGroups 
        Int jobMemory = 1
        Int timeout = 1 
    } 

    parameter_meta { 
        readGroups: "The read-group information to be added into the BAM file header" 
        jobMemory: "Memory allocated for this job" 
        timeout: "Hours before task timeout" 
    } 
    
    command <<< 
        set -euo pipefail 

        fieldNames=("ID=" "LB=" "PL=" "PU=" "SM=" "CN=" "DS=" "DT=" "FO=" "KS=" "PG=" "PI=" "PM=") 
        
        # Split the string into an array 
        IFS=, read -ra readFields <<< ~{readGroups}

        for field in "${readFields[@]}"; do 
            tag=${field:0:3}
            validTag=false
            for name in "${fieldNames[@]}"; do
                if [ "$tag" == "$name" ]; then
                    validTag=true
                fi
            done 
            if ! $validTag; then
                # Redirect error message to stderr
                echo "Invalid tag: '$tag'" >&2  
                exit 1
            fi
        done 
    >>> 

    runtime { 
        memory: "~{jobMemory} GB" 
        timeout: "~{timeout}" 
    } 

} 


# Function copied from the bwaMem workflow
task countChunkSize {

    input {
        File fastqR1
        Int numChunk
        Int? numReads
        String modules = "python/3.7"
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        numChunk: "Number of chunks to split fastq file"
        numReads: "Number of reads"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail

        if [ -z "~{numReads}" ]; then
            totalLines=$(zcat ~{fastqR1} | wc -l)
        else totalLines=$((~{numReads}*4))
        fi
        
        python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunk})*4))"
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }
    
    output {
        String chunkSize =read_string(stdout())
    }

    meta {
        output_meta: {
            chunkSize: "output number of lines per chunk"
        }
    }    
   
}


# Function copied from the bwaMem workflow
task slicer {

    input {
        File fastqR         
        String chunkSize
        String modules = "slicer/0.3.0"
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR: "Fastq file"
        chunkSize: "Number of lines per chunk"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail
        slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    } 
    
    output {
        Array[File] chunkFastq = glob("*.fastq.gz")
    }

    meta {
        output_meta: {
            chunkFastq: "Output fastq chunks"
        }
    } 
  
}


# Function copied from the bwaMem workflow
task extractUMIs {

    input {
        String umiList = "umiList"
        String outputPrefix = "extractUMIs_output"
        File fastq1
        File? fastq2
        String pattern1 = "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"
        String pattern2 = "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"
        String modules = "barcodex-rs/0.1.2 rust/1.45.1"
        Int jobMemory = 24
        Int timeout = 12
    }

    parameter_meta {
        umiList: "Reference file with valid UMIs"
        outputPrefix: "Specifies the start of the output files"
        fastq1: "FASTQ file containing read 1"
        fastq2: "FASTQ file containing read 2"
        pattern1: "UMI RegEx pattern 1"
        pattern2: "UMI RegEx pattern 2"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<
        set -euo pipefail

        barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
        --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
        ~{if (defined(fastq2)) then "--pattern2 '~{pattern2}' --r2-in ~{fastq2} " else ""}

        cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt

        tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
        echo "{$(sort -i tmp.txt)}" > new.txt
        tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{jobMemory} GB"
        timeout: "~{timeout}"
    }

    output {
        File fastqR1 = "~{outputPrefix}_R1.fastq.gz"
        File? fastqR2 = "~{outputPrefix}_R2.fastq.gz"
        File discardR1 = "~{outputPrefix}_R1.discarded.fastq.gz"
        File? discardR2 = "~{outputPrefix}_R2.discarded.fastq.gz"
        File extractR1 = "~{outputPrefix}_R1.extracted.fastq.gz"
        File? extractR2 = "~{outputPrefix}_R2.extracted.fastq.gz"
        File umiCounts = "~{outputPrefix}_UMI_counts.json"
        File extractionMetrics = "~{outputPrefix}_extraction_metrics.json"
    }

    meta {
        output_meta: {
            fastqR1: "Read 1 fastq file with UMIs extracted",
            fastqR2: "Read 2 fastq file with UMIs extracted",
            discardR1: "Reads without a matching UMI pattern in read 1",
            discardR2: "Reads without a matching UMI pattern in read 2",
            extractR1: "Extracted reads (UMIs and any spacer sequences) from read 1",
            extractR2: "Extracted reads (UMIs and any spacer sequences) from read 2",
            umiCounts: "Record of UMI counts after extraction",
            extractionMetrics: "Metrics relating to extraction process"
        }
    }

}


# Function copied from the bwaMem workflow
task adapterTrimming {

    input {
        File fastqR1
        File? fastqR2
        String modules = "cutadapt/1.8.3"
        Boolean doUMItrim = false
        Int umiLength = 5
        Int trimMinLength = 1
        Int trimMinQuality = 0
        String adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        String adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
        String? addParam
        Int jobMemory = 16
        Int timeout = 48  
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        fastqR2: "Fastq file for read 2"
        doUMItrim: "If true, do umi trimming"
        umiLength: "The number of bases to trim when doUMItrim is true. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end"
        trimMinLength: "Minimum length of reads to keep"
        trimMinQuality: "Minimum quality of read ends to keep"
        adapter1: "Adapter sequence to trim from read 1"
        adapter2: "Adapter sequence to trim from read 2"
        modules: "Required environment modules"
        addParam: "Additional cutadapt parameters"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
   
    Array[File] inputs = select_all([fastqR1,fastqR2])
    String resultFastqR1 = "~{basename(fastqR1, ".fastq.gz")}.trim.fastq.gz"
    String resultFastqR2 = if (length(inputs) > 1) then "~{basename(inputs[1], ".fastq.gz")}.trim.fastq.gz" else "None"
    String resultLog = "~{basename(fastqR1, ".fastq.gz")}.log"
    
    command <<<
        set -euo pipefail

        cutadapt -q ~{trimMinQuality} \
                -m ~{trimMinLength} \
                -a ~{adapter1} \
                -o ~{resultFastqR1} \
                ~{if (defined(fastqR2)) then "-A ~{adapter2} -p ~{resultFastqR2} " else ""} \
                ~{if (doUMItrim) then "-u ~{umiLength} -U ~{umiLength} " else ""} \
                ~{addParam} \
                ~{fastqR1} \
                ~{fastqR2} > ~{resultLog}
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    } 
    
    output { 
        File resultR1 = "~{resultFastqR1}"
        File? resultR2 = "~{resultFastqR2}"
        File log =  "~{resultLog}"     
    }

    meta {
        output_meta: {
            resultR1: "Output fastq read 1 after trimming",
            resultR2: "Output fastq read 2 after trimming",
            log: "Output adapter trimming log"
        }
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

    command <<<
        set -euo pipefail

        dragen-os \
            -r ~{dragmapHashTable} \
            -1 ~{read1s} \
            ~{if (defined(read2s)) then "-2 ~{read2s}" else ""} \
            ~{addParam} \
        | \
        samtools view -b -o ~{resultBam} -
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


task bamMerge{

    input {
        Array[File] outputBams
        String outputFileNamePrefix
        Int jobMemory = 32
        String modules = "samtools/1.9"
        Int timeout = 72
    }

    parameter_meta {
        outputBams: "Input BAM files"
        outputFileNamePrefix: "Prefix for output file"
        jobMemory: "Memory allocated for the job"
        modules: "Required environment modules"
        timeout: "Hours until task timeout"    
    }

    String resultMergedBam = "~{outputFileNamePrefix}_merged.bam"
    String tmpDir = "tmp/"
    
    command <<<
        set -euo pipefail
        mkdir -p ~{tmpDir}

        samtools merge -O bam - ~{sep=" " outputBams} \
        | \
        samtools sort -O bam -T ~{tmpDir} -o ~{resultMergedBam} - 
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File outputMergedBam = "~{resultMergedBam}"
    }

    meta {
        output_meta: {
            outputMergedBam: "Output merged, sorted BAM aligned to genome"
        }
    } 
     
}


task addReadGroups {

    input {
        File mergedBam
        String readGroups
        String modules = "picard/3.1.0"
        Int jobMemory = 12
        Int timeout = 5
    }

    parameter_meta {
        mergedBam: "Output BAM file aligned to the appropriate genome. Does not have read-group information."
        readGroups: "The read-group information to be added into the BAM file header"
        modules: "Required environment modules to run the task"
        jobMemory: "Memory allocated for the job"
        timeout: "Hours until task timeout"
    }
    
    String fileNamePrefix = basename(mergedBam, "_merged.bam")
    String resultReadGroupBam = "~{fileNamePrefix}.bam"
    String resultReadGroupBai = "~{fileNamePrefix}.bai"

    command <<<
        set -euo pipefail

        # An array containing potential read-group fields
        fieldNames=("ID=" "LB=" "PL=" "PU=" "SM=" "CN=" "DS=" "DT=" "FO=" "KS=" "PG=" "PI=" "PM=")

        # Declares an associative array and appends the values of the fields present in readGroups string
        declare -A fieldsArray
        
        for field in "${fieldNames[@]}"; do
            if [[ ~{readGroups} == *${field}* ]]; then 
                fieldsArray[${field}]=$(echo ~{readGroups} | awk -F "${field}" '{print $2}' | cut -d ',' -f 1) 
            fi
        done

        java -jar $PICARD_ROOT/picard.jar AddOrReplaceReadGroups \
            CREATE_INDEX=true \
            I= ~{mergedBam} \
            O= ~{resultReadGroupBam} \
            $( [[ -v fieldsArray["ID="] ]] && echo "RGID=${fieldsArray["ID="]}" || : ) \
            $( [[ -v fieldsArray["LB="] ]] && echo "RGLB=${fieldsArray["LB="]}" || : ) \
            $( [[ -v fieldsArray["PL="] ]] && echo "RGPL=${fieldsArray["PL="]}" || : ) \
            $( [[ -v fieldsArray["PU="] ]] && echo "RGPU=${fieldsArray["PU="]}" || : ) \
            $( [[ -v fieldsArray["SM="] ]] && echo "RGSM=${fieldsArray["SM="]}" || : ) \
            $( [[ -v fieldsArray["CN="] ]] && echo "RGCN=${fieldsArray["CN="]}" || : ) \
            $( [[ -v fieldsArray["DS="] ]] && echo "RGDS=${fieldsArray["DS="]}" || : ) \
            $( [[ -v fieldsArray["DT="] ]] && echo "RGDT=${fieldsArray["DT="]}" || : ) \
            $( [[ -v fieldsArray["FO="] ]] && echo "RGFO=${fieldsArray["FO="]}" || : ) \
            $( [[ -v fieldsArray["KS="] ]] && echo "RGKS=${fieldsArray["KS="]}" || : ) \
            $( [[ -v fieldsArray["PG="] ]] && echo "RGPG=${fieldsArray["PG="]}" || : ) \
            $( [[ -v fieldsArray["PI="] ]] && echo "RGPI=${fieldsArray["PI="]}" || : ) \
            $( [[ -v fieldsArray["PM="] ]] && echo "RGPM=${fieldsArray["PM="]}" || : )
    >>>

    runtime {
        modules: "~{modules}"
        memory:  "~{jobMemory} GB"
        timeout: "~{timeout}"
    }  
    
    output {
        File outputReadGroupBam = "~{resultReadGroupBam}"
        File outputReadGroupBai = "~{resultReadGroupBai}"
    }

    meta {
        output_meta: {
            outputReadGroupBam: "Output BAM file aligned to the appropriate genome with read-group information in the header",
            outputReadGroupBai: "Output index file for the merged BAM file with read-group information"
        }
    }

}


# Function copied from the bwaMem workflow
task adapterTrimmingLog {

    input {
        Array[File] inputLogs
        String outputFileNamePrefix
        Int numChunk
        Boolean singleEnded = false
        Int jobMemory = 12
        Int timeout = 48
    }
    parameter_meta {
        inputLogs:  "Input log files"
        outputFileNamePrefix: "Prefix for output file"
        numChunk: "Number of chunks to split fastq file"
        singleEnded: "True if reads are single ended"
        jobMemory: "Memory allocated indexing job"
        timeout: "Hours before task timeout"
    }

    String allLog = "~{outputFileNamePrefix}.txt"
    String log = "~{outputFileNamePrefix}.log"

    command <<<
        set -euo pipefail
        awk 'BEGINFILE {print "###################################\n"}{print}' ~{sep=" " inputLogs} > ~{allLog}

        totalBP=$(cat ~{allLog} | grep "Total basepairs processed:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

        bpQualitytrimmed=$(cat ~{allLog} | grep "Quality-trimmed:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g; s/ (.*)//' | awk '{x+=$1}END{print x}')
        percentQualitytrimmed=$(awk -v A="${bpQualitytrimmed}" -v B="${totalBP}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

        bpTotalWritten=$(cat ~{allLog} | grep "Total written (filtered):" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g; s/ (.*)//' | awk '{x+=$1}END{print x}')
        percentBPWritten=$(awk -v A="${bpTotalWritten}" -v B="${totalBP}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

        echo -e "This is a cutadapt summary from ~{numChunk} fastq chunks\n" > ~{log}

        if ! ~{singleEnded} ; then
          totalRead=$(cat ~{allLog} | grep "Total read pairs processed:" | cut -d":" -f2 | sed 's/ //g; s/,//g' | awk '{x+=$1}END{print x}')
          adapterR1=$(cat ~{allLog} | grep " Read 1 with adapter:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentAdapterR1=$(awk -v A="${adapterR1}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
          adapterR2=$(cat ~{allLog} | grep " Read 2 with adapter:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentAdapterR2=$(awk -v A="${adapterR2}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          shortPairs=$(cat ~{allLog} | grep "Pairs that were too short:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentShortPairs=$(awk -v A="${shortPairs}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          pairsWritten=$(cat ~{allLog} | grep "Pairs written (passing filters): " | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentpairsWritten=$(awk -v A="${pairsWritten}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          bpR1=$(cat ~{allLog} | grep -A 2 "Total basepairs processed:" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
          bpR2=$(cat ~{allLog} | grep -A 2 "Total basepairs processed:" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

          bpQualitytrimmedR1=$(cat ~{allLog} | grep -A 2 "Quality-trimmed:" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
          bpQualitytrimmedR2=$(cat ~{allLog} | grep -A 2 "Quality-trimmed:" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

          bpWrittenR1=$(cat ~{allLog} | grep -A 2 "Total written (filtered):" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
          bpWrittenR2=$(cat ~{allLog} | grep -A 2 "Total written (filtered):" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

          echo -e "Total read pairs processed:\t${totalRead}" >> ~{log}
          echo -e "  Read 1 with adapter:\t${adapterR1} (${percentAdapterR1}%)" >> ~{log}
          echo -e "  Read 2 with adapter:\t${adapterR2} (${percentAdapterR2}%)" >> ~{log}
          echo -e "Pairs that were too short:\t${shortPairs} (${percentShortPairs}%)" >> ~{log}
          echo -e "Pairs written (passing filters):\t${pairsWritten} (${percentpairsWritten}%)\n\n" >> ~{log}
          echo -e "Total basepairs processed:\t${totalBP} bp" >> ~{log}
          echo -e "  Read 1:\t${bpR1} bp" >> ~{log}
          echo -e "  Read 2:\t${bpR2} bp" >> ~{log}
          echo -e "Quality-trimmed:\t${bpQualitytrimmed} bp (${percentQualitytrimmed}%)" >> ~{log}
          echo -e "  Read 1:\t${bpQualitytrimmedR1} bp" >> ~{log}
          echo -e "  Read 2:\t${bpQualitytrimmedR2} bp" >> ~{log}
          echo -e "Total written (filtered):\t${bpTotalWritten} bp (${percentBPWritten}%)" >> ~{log}
          echo -e "  Read 1:\t${bpWrittenR1} bp" >> ~{log}
          echo -e "  Read 2:\t${bpWrittenR2} bp" >> ~{log}

        else 
          totalRead=$(cat ~{allLog} | grep "Total reads processed:" | cut -d":" -f2 | sed 's/ //g; s/,//g' | awk '{x+=$1}END{print x}')
          adapterR=$(cat ~{allLog} | grep "Reads with adapters:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentAdapterR=$(awk -v A="${adapterR}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          shortReads=$(cat ~{allLog} | grep "Reads that were too short:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentShortReads=$(awk -v A="${shortReads}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          ReadsWritten=$(cat ~{allLog} | grep "Reads written (passing filters): " | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentreadsWritten=$(awk -v A="${ReadsWritten}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')                 

          echo -e "Total reads processed:\t${totalRead}" >> ~{log}
          echo -e "Reads with adapters:\t${adapterR} (${percentAdapterR}%)" >> ~{log}
          echo -e "Reads that were too short:\t${shortReads} (${percentShortReads}%)" >> ~{log}
          echo -e "Reads written (passing filters):\t${ReadsWritten} (${percentreadsWritten}%)\n\n" >> ~{log}
          echo -e "Total basepairs processed:\t${totalBP} bp" >> ~{log}
          echo -e "Quality-trimmed:\t${bpQualitytrimmed} bp (${percentQualitytrimmed}%)" >> ~{log}
          echo -e "Total written (filtered):\t${bpTotalWritten} bp (${percentBPWritten}%)" >> ~{log}
        fi
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        timeout: "~{timeout}"
    }
  
    output {
        File summaryLog = "~{log}"
        File allLogs = "~{allLog}"
    }

    meta {
        output_meta: {
            summaryLog: "a summary log file for adapter trimming",
            allLogs: "a file containing all logs for adapter trimming for each fastq chunk"
        }
    }

}