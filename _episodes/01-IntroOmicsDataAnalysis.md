---
title: "Introduction to Omics Data Analysis"
teaching: 10
exercises: 0
questions:
- "What is bioinformatics and biostatistics?"
- "What are the similarities between bioinformatics and biostatistics?"
- "How can I combine R and BASH scripts to automate my data analysis workflow?"
- "What will we be covering in this workshop?"
objectives:
- "Learn the key features of the bioinformatics and biostatistics fields."
- "Discover the similarities between bioinformatics and biostatistics."
- "Determine what coding exercises we will be performing in this workshop."
keypoints:
- "Biostaticians apply statistical methods to the analysis of biological data sets."
- "Bioinformaticians use algorithms and computer programming to alayze omics data"
- "There are a huge variety of omics technologies that generate data for bioinformatics analysis."
- "Scrpting pipelines can be used to automate complex bioinformatics analysis workflows."
---

## What is Biostatistics?

Biostatistics is the [application of statistical methods][biostats] to the designing, analyzing, and interpreting of big biological data such as genetic sequences. This requires the development and application of computational algorithms in order to conduct different complex statistical analyses.

![Bioinformatics vs Biostats](../fig/bioinfoStats.jpeg){: width="800" }
*[Image source][bioinfoStats]*

The goal of biostatistical analyses are to answer questions and solve problems in a variety of biological fields, ranging from medicine to agriculture. For example, biostaticians are able to use mathematics study the determining factors that impact the health of people in order to arrive at conclusions about disorders and diseases.

> ## Checklist
>
> Some of the common skills required of a biostatician include:
> - experimental design
> - collecting data
> - visualizing data
> - interpreting results
> - computing
> - problem solving
> - critical thinking
> - mathematical statistics
>   - hypothesis testing
>   - multivariate analysis
>   - analysis of correlated data
>   - parameter estimation
>   - probability theory
>   - modeling
{: .checklist}


## What is Bioinformatics?

Bioinfotmatics is a interdisciplinary field that combines techniques and knowledge from both computer science and biology. It is a computational field that involves the analysis of complex genetic data. This commonly includes DNA, RNA, or protein sequence data. Bioinformatics data is generated through various omics technologies used to analyze different types of biological molecules. Omics technologies include genetic, transcriptomic, proteomic, and metabolomic data. 

![Genomics Overview](../fig/Overview-of-different-omics-sciences-such-as-genomics-transcriptomics-and-proteomics.png){: width="800" }
*[Image source][omicsInfo]*

Bioinformaticians develop algorithms and use computer programming to investigate the patterns in omics data. For example, it is possible to determine the function of newly sequenced genes by comparing their protein sequences with a database of functionally annotated sequences. More powerfully, different omics data can be combined with other biological data to find complex and subtle relationships.

![What is Bioinformatics?](../fig/Bioinformatics-graphic-01-983x640.png){: width="800" }
*[Image source][bioinfoInfo]*

> ## Checklist
>
> Some of the common skills requires of bioinformaticians include:
> - mining (collecting) data
> - visualizing data
> - computing
> - problem solving
> - critical thinking
> - omics data analysis
>   - genomics
>   - transcriptomics
>   - proteomics
>   - metabolomics
{: .checklist}

> ## Discussion
>
> What are some of the similarities between the fields of bioinformatics and biostatistics?
>
>> ## Solution
>>
>> Some of the similarities between the two fields are:
>> - colelcting data
>> - visualizing data
>> - computing
>> - problem solving
>> - critical thinking
>{: .solution}
{: .discussion}


## Data Analysis Pipelines

Together bioinformatics and biostatistics are important interdisciplinary fields for those conducting research in a wide variety of biological sciences. A typical bioinformatics or biostatistics data analysis workflow will have multiple steps, which may include:
1. experimental design
2. data collection
3. quality control
4. data preparation
5. statistical analysis of omics data
6. data visualization
7. interpretation of results

The following graphic depicts a general worklfow for preparing, analyzing, and visualzing results from transcriptomic RNA sequencing data. 

![RNA Sequencing Data Analysis Pipeline](../fig/fbinf-01-693836-g001.jpeg){: width="800" }
*[Image source][omicsWorkflow]*

In this workshop you will gain experience in performing different aspects of bioinformatics and biostatistical data analysis by working through a workflow for analyzing gene expression data. The steps in the bioinformatics workflow for this workshop are:
1. data collection - SRA toolkit
2. quality control - fastqc
3. data preparation
- convert genomic data format - gffread
- align transcriptomic data - hisat2
- quantify transcript alignments - featureCounts
4. basic statistical analysis and visualization (exact tests) - edgeR
5. advanced statistical analysis and visualization (generalized linear models) - edgeR


[bioinfoInfo]: https://www.genomicseducation.hee.nhs.uk/education/core-concepts/what-is-bioinformatics/
[omicsInfo]: https://www.researchgate.net/figure/Overview-of-different-omics-sciences-such-as-genomics-transcriptomics-and-proteomics_fig1_333003279
[biostats]: https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_biostatisticsbasics/bs704_biostatisticsbasics_print.html
[bioinfoStats]: https://cgm.sjtu.edu.cn/summer_school/
[omicsWorkflow]: https://www.frontiersin.org/articles/10.3389/fbinf.2021.693836/full

{% include links.md %}
