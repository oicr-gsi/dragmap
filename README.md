# dragmap

This workflow aligns sequence data provided as fastq files using Dragmap (an open source Dragen mapper/aligner). The alignment is completed using a hash table of a reference genome. The workflow borrows functions from the bwaMem workflow, to provide options to remove 5' umi sequences and trim off 3' sequencing adapters prior to alignment. Using these functions, this workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single BAM file to decrease workflow runtime. Read-group information must be provided.

## Overview

## Dependencies

* [dragmap 1.2.1](https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.2.1.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [picard 3.1.0](https://github.com/broadinstitute/picard/archive/refs/tags/3.1.0.tar.gz)
* [cutadapt 1.8.3](https://cutadapt.readthedocs.io/en/v1.8.3/)
* [slicer 0.3.0](https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz)
* [python 3.7](https://www.python.org)
* [barcodex-rs 0.1.2](https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz)
* [rust 1.2](https://www.rust-lang.org/tools/install)
* [gsi software modules: samtools 1.9 dragmap 1.2.1](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg38 modules: dragmap-hash-table-hg38 p12](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg19 modules: dragmap-hash-table-hg19 p13](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi mm10 modules: dragmap-hash-table-mm10 p6](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run dragmap.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqR1`|File|Fastq file for read 1
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38, mm10
`readGroups`|String|The read-group information to be added into the BAM file header


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR2`|File?|None|Fastq file for read 2
`numChunk`|Int|1|Number of chunks to split fastq file [1, no splitting]
`doUMIextract`|Boolean|false|If true, UMI will be extracted before alignment [false]
`doTrim`|Boolean|false|If true, adapters will be trimmed before alignment [false]
`numReads`|Int?|None|Number of reads


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`readGroupCheck.jobMemory`|Int|1|Memory allocated for this job
`readGroupCheck.timeout`|Int|1|Hours before task timeout
`countChunkSize.modules`|String|"python/3.7"|Required environment modules
`countChunkSize.jobMemory`|Int|16|Memory allocated for this job
`countChunkSize.timeout`|Int|48|Hours before task timeout
`slicerR1.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR1.jobMemory`|Int|16|Memory allocated for this job
`slicerR1.timeout`|Int|48|Hours before task timeout
`slicerR2.modules`|String|"slicer/0.3.0"|Required environment modules
`slicerR2.jobMemory`|Int|16|Memory allocated for this job
`slicerR2.timeout`|Int|48|Hours before task timeout
`extractUMIs.umiList`|String|"umiList"|Reference file with valid UMIs
`extractUMIs.outputPrefix`|String|"extractUMIs_output"|Specifies the start of the output files
`extractUMIs.pattern1`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 1
`extractUMIs.pattern2`|String|"(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"|UMI RegEx pattern 2
`extractUMIs.modules`|String|"barcodex-rs/0.1.2 rust/1.45.1"|Required environment modules
`extractUMIs.jobMemory`|Int|24|Memory allocated for this job
`extractUMIs.timeout`|Int|12|Time in hours before task timeout
`adapterTrimming.modules`|String|"cutadapt/1.8.3"|Required environment modules
`adapterTrimming.doUMItrim`|Boolean|false|If true, do umi trimming
`adapterTrimming.umiLength`|Int|5|The number of bases to trim when doUMItrim is true. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end
`adapterTrimming.trimMinLength`|Int|1|Minimum length of reads to keep
`adapterTrimming.trimMinQuality`|Int|0|Minimum quality of read ends to keep
`adapterTrimming.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|Adapter sequence to trim from read 1
`adapterTrimming.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|Adapter sequence to trim from read 2
`adapterTrimming.addParam`|String?|None|Additional cutadapt parameters
`adapterTrimming.jobMemory`|Int|16|Memory allocated for this job
`adapterTrimming.timeout`|Int|48|Hours before task timeout
`runDragmap.addParam`|String?|None|Additional Dragmap parameters
`runDragmap.jobMemory`|Int|50|Memory allocated for the job
`runDragmap.timeout`|Int|96|Hours until task timeout
`bamMerge.jobMemory`|Int|32|Memory allocated for the job
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.timeout`|Int|72|Hours until task timeout
`addReadGroups.modules`|String|"picard/3.1.0"|Required environment modules to run the task
`addReadGroups.jobMemory`|Int|12|Memory allocated for the job
`addReadGroups.timeout`|Int|5|Hours until task timeout
`adapterTrimmingLog.jobMemory`|Int|12|Memory allocated indexing job
`adapterTrimmingLog.timeout`|Int|48|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`dragmapBam`|File?|Merged BAM file aligned to genome with Dragmap
`dragmapIndex`|File?|Output index file for the merged BAM file
`log`|File?|A summary log file for adapter trimming
`cutAdaptAllLogs`|File?|A file with the logs from the adapter trimming of each fastq chunk


## Commands
 This section lists command(s) run by dragmap workflow
 
 * Running dragmap
 
 === Ensures the read-group information is valid, and in the correct format prior to running the rest of the workflow ===.
 
 ```
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
 ```
 
 === Parallelizes the alignment by splitting the fastq files into chunks. Subsequent steps will be run on the fastq chunks (Optional) ===.
 
 ```
         if [ -z "~{numReads}" ]; then
             totalLines=$(zcat ~{fastqR1} | wc -l)
         else totalLines=$((~{numReads}*4))
         fi
         
         python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunk})*4))"
         slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
 ```
 
 === Trims off the UMI bases (Optional) ===.
 
 ```
         barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
         --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
         ~{if (defined(fastq2)) then "--pattern2 '~{pattern2}' --r2-in ~{fastq2} " else ""}
 
         cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt
 
         tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
         echo "{$(sort -i tmp.txt)}" > new.txt
         tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
 ```
 
 === Trims off adapter sequence (Optional) ===.
 
 ```
         cutadapt -q ~{trimMinQuality} \
                 -m ~{trimMinLength} \
                 -a ~{adapter1} \
                 -o ~{resultFastqR1} \
                 ~{if (defined(fastqR2)) then "-A ~{adapter2} -p ~{resultFastqR2} " else ""} \
                 ~{if (doUMItrim) then "-u ~{umiLength} -U ~{umiLength} " else ""} \
                 ~{addParam} \
                 ~{fastqR1} \
                 ~{fastqR2} > ~{resultLog}
 ```
 
 === Align to reference using Dragmap ===.
 
 ```
         dragen-os \
             -r ~{dragmapHashTable} \
             -1 ~{read1s} \
             ~{if (defined(read2s)) then "-2 ~{read2s}" else ""} \
             ~{addParam} \
         | \
         samtools view -b -o ~{resultBam} -
 ```
 
 === Merge parallelized alignments (if the fastq was split) and sort the BAM file ===.
 
 ```
         mkdir -p ~{tmpDir}
 
         samtools merge -O bam - ~{sep=" " bams} \
         | \
         samtools sort -O bam -T ~{tmpDir} -o ~{resultMergedBam} - 
 ```
 
 === Assigns read-groups and index BAM file ===.
 
 ```
         # An array containing potential read-group fields
         fieldNames=("ID=" "LB=" "PL=" "PU=" "SM=" "CN=" "DS=" "DT=" "FO=" "KS=" "PG=" "PI=" "PM=")
 
         # Declares an associative array and appends the values of the fields present in readGroups string
         declare -A fieldsArray
         
         for field in "${fieldNames[@]}"; do
             if [[ ~{readGroups} == *${field}* ]]; then 
                 fieldsArray[${field}]=$(echo ~{readGroups} | awk -F "${field}" '{print $2}' | cut -d ',' -f 1) 
             fi
         done
 
         java -jar picard.jar AddOrReplaceReadGroups \
             CREATE_INDEX=true \
             I= ~{outputBam} \
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
 ```
 
 === Merges parallelized adapter trimming logs ===.
 
 ```
         COMMANDS NOT SHOWN, see WDL for details
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
