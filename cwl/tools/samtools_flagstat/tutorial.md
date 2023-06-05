


# Samtools Flagstat Tool Translation

## Introduction

This tutorial demonstrates translation of a basic `samtools flagstat` tool from CWL to Nextflow using `janis translate`. <br>

**Source Tool**

The CWL tool used in this tutorial - [samtools_flagstat](https://github.com/genome/analysis-workflows/blob/master/definitions/tools/samtools_flagstat.cwl) -  is taken from the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) [analysis-workflows](https://github.com/genome/analysis-workflows) repository. 

This resource stores publically available analysis pipelines for genomics data. <br>
It is a fantastic piece of research software, and the authors thank MGI for their contribution to open-source research software. 

The underlying software run by this tool - [samtools_flagstat](http://www.htslib.org/doc/samtools-flagstat.html) - displays summary information for an alignment file. 

**Tutorial Outcomes**

In this tutorial we will:
- Install the required software
- Translate the CWL using `janis translate`
- Make manual adjustments to the translation if necessary
- Run the nextflow using sample input data to validate our nextflow code

After completing this short tutorial, you will be familiar with using `janis translate` to migrate workflow tools in CWL to Nextflow.

Other tutorials exist to demonstrate migration from WDL / CWL / Galaxy -> Nextflow in this repository, including full workflow migrations with multiple tasks. 

**Installation**

To begin, make sure you have [nextflow](https://nf-co.re/usage/installation), [docker](https://docs.docker.com/engine/install/), and [janis translate](https://janis.readthedocs.io/en/latest/index.html) installed. <br>
The links above contain installation instructions. 

<br>

## Janis Translate

To translate a workflow,  we use `janis translate`.

```
janis translate --from <src> --to <dest> <filepath>
```

The `--from` specifies the workflow language of the source file(s), and `--to` specifies the destination we want to translate to. 

In our case, we want to translate CWL -> Nextflow, and our source CWL file is located at `source/samtools_flagstat.cwl` relative to this document.

<br>

*using pip*

To translate `samtools_flagstat.cwl` to nextflow, we can write the following in a shell:
```
janis translate --from cwl --to nextflow ./source/samtools_flagstat.cwl
```

*using docker (linux bash)*

If the janis translate docker container is being used, we can write the following:
```
docker run -v $(pwd):/home janis translate --from cwl --to nextflow ./source/samtools_flagstat.cwl
```

<br>

You will see a folder called `translated` appear, and a nextflow process called `samtools_flagstat.nf` will be present inside. 

<br>

## Manual Adjustments

The `translated/samtools_flagstat.nf` file should be similar to the following: 

```
nextflow.enable.dsl=2

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

    input:
    path bam

    output:
    path "${bam[0]}.flagstat", emit: flagstats

    script:
    def bam = bam[0]
    """
    /usr/local/bin/samtools flagstat \
    ${bam} \
    > ${bam}.flagstat \
    """

}
```

We can see that this nextflow process has a single input, a single output, and calls `samtools flagstat` on the input `bam`. 

We can also see that a container image is available for this tool. In the next section we will run this process using some sample data and the specified container. 

This translation is correct for the `samtools_flagstat.cwl` file and needs no adjusting. <br>
Have a look at the source CWL file to see how they match up. 

> Note: <br>
> `def bam = bam[0]` in the script block is used to handle datatypes with secondary files. <br>
> The `bam` input is an indexed bam type, so requires a `.bai` file to also be present in the working directory alongside the `.bam` file. <br><br>
> For this reason the `bam` input is supplied as an Array with 2 files - the `.bam` and the `.bai`. <br>
> Here the `def bam = bam[0]` is used so that `bam` refers to the `.bam` file in that Array. 

<br>

## Running Samtools Flagstat as a Workflow


**Collecting Process Outputs**

Let's add a `publishDir` directive to our translated process so that we can capture the outputs of this process.

```
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    publishDir "./outputs"                                   
```

Nextflow allows us to capture the outputs created by a process using the `publishDir` directive seen above. 

<br>

**Setting up nextflow.config**

To run this process, we will set up a `nextflow.config` file and add some lines to the top of our process definition to turn it into a workflow.

Create a new file called `nextflow.config` in the `translated` folder alongside `samtools_flagstat.nf`. 

Copy and paste the following code into your `nextflow.config` file: 

```
nextflow.enable.dsl=2
docker.enabled = true

params {

    bam = [
        '../../../../sample_data/cwl/2895499223_sorted.bam',
        '../../../../sample_data/cwl/2895499223_sorted.bai',
    ]

}
```
<br>

This tells nextflow how to run, and sets up an input parameter for our indexed bam input.

The `bam` parameter is a list which provides paths to the `.bam` and `.bai` sample data we will use to test the nextflow translation. From here, we can refer to the indexed bam input as `params.bam` in other files.

> NOTE<br>
> `nextflow.enable.dsl=2` ensures that we are using the dsl2 nextflow syntax which is the current standard. <br>
> `docker.enabled = true` tells nextflow to run processes using docker. Our `samtools_flagstat.nf` has a directive with the form `container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"` provided, so it will use the specified image when running this process. 

<br>

**Creating Workflow & Passing Data** 

Now that we have the `nextflow.config` file set up, we will add a few lines to `samtools_flagstat.nf` to turn it into a workflow. 

Copy and paste the following lines at the top of `samtools_flagstat.nf`:

```
ch_bam = Channel.fromPath( params.bam ).toList()

workflow {
    SAMTOOLS_FLAGSTAT(ch_bam)
}
```

The first line creates a nextflow `Channel` for our `bam` input and ensures it is a list. <br>
The `Channel.toList()` part collects our files into a list, as both the `.bam` and `.bai` files must be passed together. <br>
The `params.bam` global variable we set up previously is used to supply the 
paths to our sample data.

The new `workflow {}` section declares the main workflow entry point. <br>
When we run this file, nextflow will look for this section and run the workflow contained within. 

In our case, the workflow only contains a single task, which runs the `SAMTOOLS_FLAGSTAT` process defined below the workflow section. The single `SAMTOOLS_FLAGSTAT` input is being passed data from our `ch_bam` channel we declared at the top of the file. 

<br>

**Running Our Workflow**

Ensure you are in the `translated/` working directory, where `nextflow.config` and `samtools_flagstat.nf` reside. 

```
cd translated/
```

To run the workflow using our sample data, we can now write the following command: 
```
nextflow run samtools_flagstat.nf
```

Nextflow will automatically check if there is a `nextflow.config` file in the working directory, and if so will use that to configure itself. Our inputs are supplied in `nextflow.config` alongside the dsl2 & docker config, so it should run without issue. 

Once completed, we can check the `./outputs` folder to view our results. 

If everything went well, there should be a file called `2895499223_sorted.bam.flagstat` with the following contents:

```
2495 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
1 + 0 supplementary
0 + 0 duplicates
2480 + 0 mapped (99.40% : N/A)
2494 + 0 paired in sequencing
1247 + 0 read1
1247 + 0 read2
2460 + 0 properly paired (98.64% : N/A)
2464 + 0 with itself and mate mapped
15 + 0 singletons (0.60% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

```

If needed, you can check the `./final` folder which contains the files we created in this tutorial.  

<br>

### Conclusion

In this tutorial we explored how to translate a simple CWL tool to a Nextflow process. 

A tutorial for a more complex CWL tool is available in the `cwl/tools/gatk_haplotype_caller` folder. 

A tutorial for a CWL workflow translation to Nextflow is available in the `cwl/workflows/align_sort_markdup` folder. 