---
title: "Bioinformatics Workflow Pipelining"
teaching: 10
exercises: 15
questions:
- "How can I use R and BASH scripts to automate my data analysis process?"
objectives:
- "Develop BASH scripting skills to efficiently integrate R scripts into a coherent data analysis pipeline."
keypoints:
- "TBD"
---

## Pipelining & The Benefits of Scripting

The primary goal of any bioinformatics analysis workflow is to give meaning to large and complex biological data sets. This often requires the development of code that is generalized and can be automated to run on multiple similar data sets, which may have small differences in their structure or content.

![What is Shell Scripting?](../fig/what_is_shell_scripting.jpeg){: width="500" }
*[Image source][scriptingBenefits]*

The analysis of biological data using bioinformatics or biostatistical analysis is often a process that requires multiple steps. Furthermore, these steps often need to be completed in a specific order to achieve a final result. This process involves the output of a previous step being used as input to the subsequent step. This flow of outputs, to inputs, to a final result can be depicted as water flowing from one pipe to the next.

![Pipelineing Depiction](../fig/noun-pvc-pipes-147592.png){: width="500" }

By creating a pipeline of scripts, it is possible to automate much of the bioinformatics workflow. This means making modular scripts to perform different steps of the analysis process. 




## BASH Scripting

So, we can use BASH scripting to create modular pieces of code for use in bioinformatics data analysis pipelines. BASH scripts are text files that have the **.sh** file extension. These are text files that you can use to save the lines of BASH code that you want the interpreter componenet of the computer operating system to execute (run).

![The Interpreter Operating System Component](../fig/interpreter.png){: width="800" }
*[Image source][interpreterComp]*



[scriptingBenefits]: https://www.interviewbit.com/shell-scripting-interview-questions/ 
[experienceCite]: https://www.abstract.tech/blog/let-people-die-in-virtual-reality
[interpreterComp]: https://www.geeksforgeeks.org/difference-between-assembler-and-interpreter/ 

{% include links.md %}
