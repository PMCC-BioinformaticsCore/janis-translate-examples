
[//]: <> (PREAMBLE)

<img src="../media/melbioinf_logo.png" width="350"> <img src="../media/PRIMARY_A_Vertical_Housed_RGB.png" width="150">

# Janis Translate - Migrating CWL / Galaxy to Nextflow

Anticipated workshop duration when delivered to a group of participants is **3.5 hours**.  

For queries relating to this workshop, contact Melbourne Bioinformatics (bioinformatics-training@unimelb.edu.au).

<br>

## Overview

### Topic

* [ ] Genomics
* [ ] Transcriptomics
* [ ] Proteomics
* [ ] Metabolomics
* [ ] Statistics and visualisation
* [ ] Structural Modelling
* [x] Basic skills

### Skill level

* [ ] Beginner  
* [ ] Intermediate  
* [x] Advanced  

The workshop is conducted in a Unix environment.<br>
Command line experience is required. <br>
Prior experience with developing, running, and troubleshooting Nextflow workflows is strongly recommended.

### Description

Bioinformatics workflows are critical for reproducibly transferring methodologies between research groups and for scaling between computational infrastructures. Research groups currently invest a lot of time and effort in creating and updating workflows; the ability to translate from one workflow language into another can make them easier to share, and maintain with minimal effort. For example, research groups that would like to run an existing Galaxy workflow on HPC, or extend it for their use, might find translating the workflow to Nextflow more suitable for their ongoing use-cases. 

Janis is a framework that provides an abstraction layer for describing workflows, and a tool that can translate workflows between existing languages such as CWL, WDL, Galaxy and Nextflow. Janis aims to translate as much as it can, leaving the user to validate the workflow and make small manual adjustments where direct translations are not possible. Originating from the Portable Pipelines Project between Melbourne Bioinformatics, the Peter MacCallum Cancer Centre, and the Walter and Eliza Hall Institute of Medical Research, this tool is now available for everyone to use.

This workshop provides an introduction to Janis and how it can be used to translate Galaxy and CWL based tools and workflows into Nextflow. Using hands-on examples we’ll step you through the process and demonstrate how to optimise, troubleshoot and test the translated workflows.

*Section 1* covers migration of CWL tools / workflows to Nextflow. <br>
*Section 2* covers migration of Galaxy tool wrappers / workflows to Nextflow.

-------------------------------

### Learning Objectives

By the end of the workshop you should be able to:
- Recognise the main aspects and benefits of workflow translation
- Use Janis to translated Galaxy / CWL tools & workflows to Nextflow
- Configure Nextflow to run translated tools & workflows
- Troubleshoot translated Nextflow tools & workflow 
- Adjust the translated Nextflow tools / workflows & complete missing translations manually

-------------------------------

### Required Software
- Janis
- Nextflow
- Singularity or Docker
- IDE of your choosing (we use VS Code in this workshop)

### Required Data
* The sample data can be downloaded via zenodo at the following link: [link]

-------------------------------

### Author Information

Written by: Grace Hall  
Melbourne Bioinformatics, University of Melbourne

Created/Reviewed: May 2023

-------------------------------

[//]: <> (/PREAMBLE)

-------------------------------------------------------------------
-------------------------------------------------------------------

[//]: <> (SETUP_LINK)

Refer back to the [setup instructions](#Setup) if required.

[//]: <> (/SETUP_LINK)

-------------------------------------------------------------------
-------------------------------------------------------------------

[//]: <> (SOFTWARE_INSTALLATION)

Make sure you have the following installed:
- [Nextflow](https://nf-co.re/usage/installation)
- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) or [Docker](https://docs.docker.com/engine/install/)
- [Janis](https://janis.readthedocs.io/en/latest/index.html)

[//]: <> (/SOFTWARE_INSTALLATION)

-------------------------------------------------------------------
-------------------------------------------------------------------

[//]: <> (VSC_EXTENSIONS)

VS Code can be customised using extensions. 

We will install the Nextflow extension which provides syntax highlighting and quick shortcuts for common code snippets.

[//]: <> (/VSC_EXTENSIONS)

-------------------------------------------------------------------
-------------------------------------------------------------------

[//]: <> (REMOTE_SSH_EXTENSION)

[//]: <> (/REMOTE_SSH_EXTENSION)

-------------------------------------------------------------------
-------------------------------------------------------------------
