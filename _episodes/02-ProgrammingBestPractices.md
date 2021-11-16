---
title: "Programming Best Practices"
teaching: 20
exercises: 10
questions:
- "What are the benefits of writing programs?"
- "What are the most helpful programming techniques?"
- "How can I get started with writing a program?"
- "How can I tackle programming mistakes?"
objectives:
- "Be able to write pseudocode to describe the steps of a program in a plain language."
- "Become familiar with methods for writing modular and understandable programs."
- "Be able to break down an overly complex piece of code into smaller, more readily understandable components."
- "Be able to write helpful and simple comments and documents for programs."
keypoints:
- "Use programs to accomplish complex or repetative tasks."
- "Write programs that can be understood by others."
- "Take the time to plan how you will write a program."
- "Make small changes and plan for mistakes."
- "Collaborate with others whenever possible."
- "Always include informative documents for your programs."
---

## Introduction
In this section we will learn about some of the common best practices in programming, which are easy to implememnt into your personal programming process. We will also explore approaches to solving problems and where to begin with designing algorithms.

## How to Be a Good Programmer
The development of custom software programs has become increasingly necessary in biological research. Scientists are often required to create their own programs to analyze data and create publishable results. It is therefore very important that we consider [techniques][goodProgrammer] for improving the reproducibility and reliability of code. 

![Good Programmer Characteristics](../fig/codingmindset-sm.png)
*Image source: https://mitcommlab.mit.edu/broad/commkit/coding-mindset/*

> ## Checklist
>
> These are programming techniques that have been found to be helpful and effective in a variety of research settings.
>
> - Use programs to accomplish complex or repetative tasks
> - Write programs that can be understood by others
> - Take the time to plan how you will write a program
> - Make small changes and plan for mistakes
> - Collaborate with others whenever possible
> - Always include informative documents for your programs and data
> - Carefully structure and track your raw and calculated data
{: .checklist}

## Ways to Approach Programming Tasks
Throughout any programming undertaking we should be thinking about our problem solving thought process. This means that you will need to think critically about how you approach solving problems with programming. Often you will find that there are many routes to the same solution, and which route you take may depend on your intended user or available tools.

> ## Checklist
>
> These are steps you can take to approach solving a problem.
>
> 1. Understand the problem
> 2. Create a plan to solve the problem
> 3. Implement the plan
> 4. Reflect on the results
{: .checklist}

The first step for approaching problem solving requires us to break down the problem before we can begin creating a solution plan. There are a few techniques we can use to help break down a problem before coding:

1. Determine the inputs
2. Determine the outputs
3. Test a simple example
4. Test a complex example

Now, let's put these steps into practice. Keep in mind that the number of steps a task or problem is broken into may depend on the skills of the intended user or available tools.

> ## Challenge
>
> Write an algorithm in pseudocode to complete the task of getting dressed for the day, while considering the:
> - Current weather
> - Clothes available to you
>
> To keep things simple, assume that:
> - You are currently wearing pajamas
> - You will wear only one top and one bottom clothing item
> - You will be outside all day
> - The weather will not change
> 
>> ## Solution
>>
>> In order to determine how to write an algorithm for getting dressed for the day, we should consider the following steps to breaking down a problem.
>>
>> **First,** determine the inputs
>> 1. Current weather
>> 2. Available clothing
>>
>> **Second,** determine the outputs
>> 1. The clothes that you will be wearing for the day
>> 2. The order in which the clothing should be put on
>>
>> **Third,** test a simple example by specifying sample inputs:
>> 1. The weather outside is cold
>> 2. You have access to a pair of pants and a shirt
>>
>> *Our simple algorithm might then be:*
>> 1. Walk to where your clothes are kept
>> 2. Take off pajamas
>> 3. Take out the the pants and shirt
>> 4. Put on pants
>> 5. Put on shirt
>>
>> **Fourth,** test a complex example
>> Let's try out a more complex example by generalizing the inputs:
>> 1. Assume you have a way to check the current weather
>> 2. Assume you have a closet with all types of clothing
>>
>> *Our more complex algorithm might then be:*
>> 1. Check the weather
>> 2. Walk to wear you clothes are kept
>> 3. If the temperature is less than 75 degrees fahrenheit, then
>> - Put on pants
>> - Put on shirt
>> 4. If the temperature is greater than 75 degrees fahrenheit, then
>> - Put on shorts
>> - Put on tank top
>>
>> Note that one way to generalize your algorithm is to use conditional statements, such as the "if" statements in the above example algorithm. Remember that conditonal statements are used in programming to handle descisions, and they have two parts: hypothesis (if) and conclusion (then). So, the outcome of a conditional statement depends on the state of the inputs at that moment.
> {: .solution}
{: .challenge}

After devising a plan for a solution to a problem or task, it is a good idea to stop and think carefully about the plan. This is particularly important for debugging and fixing any errors. 

Some questions you can ask yourself at this point include:
- Is my solution comprehensive?
- Did I make any mistakes?
- How can errors or incorrect outputs arise?
- What can I do next?

Considering our simple solution to the previous challenge of writing an algorithm for getting dressed, there remain other ways that the algorithm may be written to be more comprehensive. For example, what if the intended user or audience is a young child? Then it may be necessary to further break down the steps of the algorithm to meet the needs of the user.

For example, let's re-write step 4 of the simple algorithm from our solution to the above challenge. Instead of leaving the step to "Put on pants", we might break down the step as follows:
1. Hold pants
2. Open waistband
3. Insert right leg into right leg hole of pants
4. Insert left leg into left leg hole of pants
5. Pull left pant leg up so the left foot comes through it
6. Pull pants up from waitsband

[goodProgrammer]: https://mitcommlab.mit.edu/broad/commkit/coding-mindset/

{% include links.md %}
