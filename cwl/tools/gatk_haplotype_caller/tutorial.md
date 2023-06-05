


# GATK HaplotypeCaller Tool Translation

## Introduction

This tutorial demonstrates translation of a `gatk HaplotypeCaller` tool from CWL to Nextflow using `janis translate`. <br>

**Source Tool**

The CWL tool used in this tutorial is taken from the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) [analysis-workflows](https://github.com/genome/analysis-workflows) repository. 

This resource stores publically available analysis pipelines for genomics data. <br>
It is a fantastic piece of research software, and the authors thank MGI for their contribution to open-source research software. 

The tool using in this tutorial - [samtools_flagstat](https://github.com/genome/analysis-workflows/blob/master/definitions/tools/samtools_flagstat.cwl) - displays summary information for an alignment file. 


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

In our case, we want to translate CWL -> Nextflow, and our source CWL file is located at `source/gatk_haplotype_caller.cwl` relative to this document.

<br>

*using pip*

To translate `gatk_haplotype_caller.cwl` to nextflow, we can write the following in a shell:
```
janis translate --from cwl --to nextflow ./source/gatk_haplotype_caller.cwl
```

*using docker (linux bash)*

If the janis translate docker container is being used, we can write the following:
```
docker run -v $(pwd):/home janis translate --from cwl --to nextflow ./source/gatk_haplotype_caller.cwl
```

<br>

You will see a folder called `translated` appear, and a nextflow process called `gatk_haplotype_caller.nf` will be present inside. 

<br>

## Manual Adjustments

The `translated/gatk_haplotype_caller.nf` file should be similar to the following: 

```
nextflow.enable.dsl=2

process GATK_HAPLOTYPE_CALLER {
    
    container "broadinstitute/gatk:4.1.8.1"

    input:
    path bam
    path reference
    path dbsnp_vcf
    val gvcf_gq_bands
    val intervals
    val emit_reference_confidence
    val contamination_fraction
    val max_alternate_alleles
    val ploidy
    val read_filter

    output:
    tuple path("output.g.vcf.gz"), path("*.tbi"), emit: gvcf

    script:
    def bam = bam[0]
    def dbsnp_vcf = dbsnp_vcf[0] != null ? "--dbsnp ${dbsnp_vcf[0]}" : ""
    def reference = reference[0]
    def gvcf_gq_bands_joined = gvcf_gq_bands.join(' ')
    def intervals_joined = intervals.join(' ')
    def contamination_fraction = contamination_fraction != params.NULL_VALUE ? "-contamination ${contamination_fraction}" : ""
    def max_alternate_alleles = max_alternate_alleles != params.NULL_VALUE ? "--max_alternate_alleles ${max_alternate_alleles}" : ""
    def ploidy = ploidy != params.NULL_VALUE ? "-ploidy ${ploidy}" : ""
    def read_filter = read_filter != params.NULL_VALUE ? "--read_filter ${read_filter}" : ""
    """
    /gatk/gatk --java-options -Xmx16g HaplotypeCaller \
    -R ${reference} \
    -I ${bam} \
    -ERC ${emit_reference_confidence} \
    -GQB ${gvcf_gq_bands_joined} \
    -L ${intervals_joined} \
    ${dbsnp_vcf} \
    ${contamination_fraction} \
    ${max_alternate_alleles} \
    ${ploidy} \
    ${read_filter} \
    -O "output.g.vcf.gz" \
    """

}
```

We can see that this nextflow process has a multiple inputs, single output, and calls `gatk HaplotypeCaller` using the input data we supply to the process.  

We can also see that a container image is available for this tool. In the next section we will run this process using some sample data and the specified container. 

This translation is correct for the `gatk_haplotype_caller.cwl` file and needs no adjusting. <br>
Have a look at the source CWL file to see how they match up. 

> Notes on translation: <br><br>
> **(1)** `def bam = bam[0]`<br><br>
> This pattern is used in the script block to handle datatypes with secondary files. <br>
> The `bam` input is an indexed bam type, so requires a `.bai` file to also be present in the working directory alongside the `.bam` file. <br>
> For this reason the `bam` input is supplied as a list with 2 files - the `.bam` and the `.bai`. <br>
> `def bam = bam[0]` is used so when we reference `${bam}` in the script body, we are refering to the `.bam` file in that list. <br><br>
> **(2)** `def dbsnp_vcf = dbsnp_vcf[0] != null ? "--dbsnp ${dbsnp_vcf[0]}" : ""`<br><br>
> This is performing the same job as item (1), with some extra details. <br>
> The `dbsnp_vcf` input is an indexed vcf datatype. <br>
> Similar to `bam`, it is passed as a list containing a `.vcf.gz` file and a `.vcf.gz.tbi` file. <br>
> It differs because it is an *optional* input, and needs the `--dbsnp` prefix. <br><br>
> Here we use a *ternary operator* to check if the input is **null**. <br>
> The format is `cond_check ? cond_true : cond_false`.<br>
> The `dbsnp_vcf[0] != null` here checks if the `.vcf.gz` file is present. <br>
> If present, it templates a string we can use in the script body to provide the required argument: `--dbsnp ${dbsnp_vcf[0]}`<br>
> If absent, it templates an empty string. <br>
> This ensures that when we use `${dbsnp_vcf}` in the script body, it will be formatted correctly for either case.<br><br>
> **(3)** `def intervals_joined = intervals.join(' ')` <br><br>
> Templates our `intervals` list of strings to a single string.<br>
> Each item is joined by a space: eg ["hello", "there"] -> "hello there". <br><br>
> **(4)** `def ploidy = ploidy != params.NULL_VALUE ? "-ploidy ${ploidy}" : ""` <br><br>
> Same as **(2)** except uses a different check for null value. <br>
> Nextflow doesn't like **null** values to be passed to process inputs, so we use 2 different tricks to make **optional** inputs possible.<br><br>
> For `val` inputs we set up a `NULL_VALUE` param in `nextflow.config` which we use as a placeholder.  <br>
> We will do this in the following section. <br><br>
> For `path` inputs (ie files and directories) we set up a **null** file rather than passing **null** directly. <br>
> This ensures that the file is staged correctly in the working directory when an actual filepath is provided. 

<br>

## Running GATK HaplotypeCaller as a Workflow


**Collecting Process Outputs**

Let's add a `publishDir` directive to our translated process so that we can capture the outputs of this process.

```
    container "broadinstitute/gatk:4.1.8.1"
    publishDir "./outputs"                                   
```

Nextflow allows us to capture the outputs created by a process using the `publishDir` directive seen above. 

<br>

**Setting up nextflow.config**

To run this process, we will set up a `nextflow.config` file and add some lines to the top of our process definition to turn it into a workflow.

Create a new file called `nextflow.config` in the `translated` folder alongside `gatk_haplotype_caller.nf`. 

Copy and paste the following code into your `nextflow.config` file: 

```
docker.enabled = true
nextflow.enable.dsl=2


params {
    NULL_VALUE = 'NULL'

    bam = [
        '../../../../sample_data/cwl/2895499223_sorted.bam',
        '../../../../sample_data/cwl/2895499223_sorted.bam.bai',
    ]
    dbsnp = [
        '../../../../sample_data/cwl/chr17_test_dbsnp.vcf.gz',
        '../../../../sample_data/cwl/chr17_test_dbsnp.vcf.gz.tbi',
    ]
    reference = [
        '../../../../sample_data/cwl/chr17_test.fa',
        '../../../../sample_data/cwl/chr17_test.dict',
        '../../../../sample_data/cwl/chr17_test.fa.fai',
    ]
    gvcf_gq_bands = NULL_VALUE
    intervals = ["chr17"]
    emit_reference_confidence = 'GVCF'
    contamination_fraction = NULL_VALUE
    max_alternate_alleles = NULL_VALUE
    ploidy = NULL_VALUE
    read_filter = NULL_VALUE
}
```
<br>

This tells nextflow how to run, and sets up the sample data as inputs.

*file inputs*

The `bam` parameter is a list which provides paths to the `.bam` and `.bai` sample data we will use to test the nextflow translation. From here, we can refer to the indexed bam input as `params.bam` in other files. The `dbsnp` and `reference` params follow this same pattern. 

*non-file inputs*

We also set up a `NULL_VALUE` param which we use as a *placeholder* for a null value. <br> 
In this case we are providing null values for the `gvcf_gq_bands`, `contamination_fraction`, `max_alternate_alleles`, `ploidy` and `read_filter` inputs as they are all optional.

<br>

> NOTE<br>
> `nextflow.enable.dsl=2` ensures that we are using the dsl2 nextflow syntax which is the current standard. <br>
> `docker.enabled = true` tells nextflow to run processes using docker. Our `gatk_haplotype_caller.nf` has a directive with the form `container "broadinstitute/gatk:4.1.8.1"` provided, so it will use the specified image when running this process. 

<br>

**Creating Workflow & Passing Data** 

Now that we have the `nextflow.config` file set up, we will add a few lines to `gatk_haplotype_caller.nf` to turn it into a workflow. 

Copy and paste the following lines at the top of `gatk_haplotype_caller.nf`:

```
ch_bam = Channel.fromPath( params.bam ).toList()
ch_dbsnp = Channel.fromPath( params.dbsnp ).toList()
ch_reference = Channel.fromPath( params.reference ).toList()

workflow {

    GATK_HAPLOTYPE_CALLER(
        ch_bam,
        ch_reference,
        ch_dbsnp,
        params.intervals,
        params.gvcf_gq_bands,
        params.emit_reference_confidence,
        params.contamination_fraction,
        params.max_alternate_alleles,
        params.ploidy,
        params.read_filter,
    )

}
```

The first 3 lines create nextflow `Channels` for our `bam`, `dbsnp`, and `reference` inputs and ensures they are lists. <br>
The `Channel.toList()` aspect collects our files into a list, as the primary & secondary files for these datatypes must be passed together.

The `params.bam`, `params.dbsnp` and `params.reference` global variables we set up previously are used to supply the paths to our sample data for these channels.

The new `workflow {}` section declares the main workflow entry point. <br>
When we run this file, nextflow will look for this section and run the workflow contained within. 

In our case, the workflow only contains a single task, which runs the `GATK_HAPLOTYPE_CALLER` process defined below the workflow section. We call `GATK_HAPLOTYPE_CALLER` by feeding inputs in the correct order, using the channels we declared at the top of the file, and variables we set up in the global `params` object. 

<br>

**Running Our Workflow**

Ensure you are in the `translated/` working directory, where `nextflow.config` and `gatk_haplotype_caller.nf` reside. 
```
cd translated/
```

To run the workflow using our sample data, we can now write the following command: 
```
nextflow run gatk_haplotype_caller.nf
```

Nextflow will automatically check if there is a `nextflow.config` file in the working directory, and if so will use that to configure itself. Our inputs are supplied in `nextflow.config` alongside the dsl2 & docker config, so it should run without issue. 

Once completed, we can check the `./outputs` folder to view our results. <br>
If everything went well, the `./outputs` folder should contain 2 files: 
- `output.g.vcf.gz`
- `output.g.vcf.gz.tbi`.

If needed, you can check the `./final` folder which contains the files we created in this tutorial as reference.  

<br>

### Conclusion

In this tutorial we explored how to translate the `gatk_haplotype_caller` CWL tool to a Nextflow process. 

A tutorial for a CWL workflow translation to Nextflow is available in the `cwl/workflows/align_sort_markdup` folder. 