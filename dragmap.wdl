version 1.0

import "imports/pull_bwaMem.wdl" bwaMem


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
      ##UPDATE THE PATHS
    "hg19": {"samtools/1.9 dragmap/1.2.1 dragmap-hash-table-hg19/p13", "PATH"},
    "hg38": {"samtools/1.9 dragmap/1.2.1 dragmap-hash-table-hg38/p12", "PATH"},
    "mm10": {"samtools/1.9 dragmap/1.2.1 dragmap-hash-table-mm10/p6", "PATH"}
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
            File fastqR2_ = select_all([fastqR2])[0]
            call bwaMem.slicer as slicerR2 {
                input:
                fastqR = fastqR2,
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
    }

    call bwaMem.bamMerge as bamMerge {
        input:
        bams = runDragmap.outputBam,
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

    ##Do we need to include all of the dependencies that bwaMem has and ours? Would I need to remove the bwa?
    meta {
        author: "Xuemei Luo and Muna Mohamed"
        email: "xuemei.luo@oicr.on.ca, mmohamed@oicr.on.ca"
        description: "This workflow aligns sequence data provided as fastq files against a genomic reference using bwa (burrows-wheeler-aligner).  Prior to alignment, there are options to remove 5' umi sequence and to trim off 3' sequencing adapter. Readgroup information to be injected into the bam header needs to be provided.  The workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single bam file.  This decreases the workflow run time.  Optional bwa mem parameters can be provided to the workflow."
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
        String readGroups
        String modules
        String dragmapHashTable
        ##Is this too many? Based on my runs
        Int jobMemory = 50
        Int timeout = 48
    }

    parameter_meta {
        read1s: "Fastq file for read 1"
        read2s: "Fastq file for read 2"
        ## Unsure if I can keep this, and need to ask Mei if this effects downstream modules
        readGroups: "The readgroup information to be injected into the bam header"
        dragmapHashTable: "The hash table prepared by Dragmap using the appropriate reference genome. Used to align sample with Dragmap"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for the job"
        timeout: "Hours until task timeout"
    }

    # --RGID arg (=1)                                 Read Group ID
    # --RGSM arg (=none)                              Read Group Sample
    #   --num-threads arg (=32)                 Worker threads for mapper/aligner (default = maximum available on system)
    # --output-directory arg                          Output directory
    # --output-file-prefix arg                        Output filename prefix
    String resultBam = "~{basename(read1s)}.bam"
    String tmpDir = "tmp/"

    # Have to do if/else statement for if single read or not
    command <<<
        set -euo pipefail
        mkdir -p ~{tmpDir}
        
        dragen-os \
            -r ~{dragmapHashTable} \
            -1 ~{read1s} \
            -2 ~{read2s} \
            --RGSM arg ~{readGroups} \
        | \
        samtools view -b -
        | \
        samtools sort -O bam -T ~{tmpDir} -o ~{resultBam} - 
    >>>

    runtime {
        modules: "~{modules}"
        memory:  "~{jobMemory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }  
    
    output {
        File outputBam = "~{resultBam}"
    }

    meta {
        output_meta: {
            outputBam: "Output Bam file aligned to the appropriate genome"
        }
    }

}
