---
title: "Programming Basics"
teaching: 20
exercises: 10
questions:
- "What is the difference between programming and coding?"
- "What are algorithms and how are they developed?"
- "What is pseudocode and how can it be used?"
- "What are the most commonly used programming commands?"
objectives:
- "Become familiar with the key concepts of programming, coding, algorithms, and pseudocode."
- "Be able to design algorithms to accomplish basic everyday tasks."
- "Become familiar with logical expressions commonly used in programming."
keypoints:
- "Programming is the process of creating instructions or set of related activities to achieve a task or goal."
- "Coding is the process of transforming the set of instructions for a process into a written language that a computer can interpret."
- "Algorithms are the set of step-by-step instructions that explain how to solve a given problem."
- "Pseudocode is the set of instructions for an algorithm written in a plain language."
- "Boolean algebra uses mathematical expressions that are evaluated to one of two values: true or false."
- "Conditional Statements are used in programming to handle descisions, and they have two parts: hypothesis (if) and conclusion (then)."
---

## Introduction
Before we begin learning about how to write helpful programs for data analysis, it is important that we consider fundamental concepts and best practices in programming. 

## Programming vs Coding 
While sometimes used interchangeably, [programming and coding][codingProgramming] have different definitions. 

![This is an image](/fig/What-Are-The-Main-Uses-of-Java-in-The-World-2-1.jpeg)

Based on your personal experiences, let's discuss our current understanding of these important concepts.

> ## Discussion
>
> What is programming?
>
>> ## Solution
>>
>> Programming is the process of creating instructions or set of related activities to achieve a task or goal.
> {: .solution}
>
> What is coding?
>
>> ## Solution
>>
>> Coding is the process of transforming the set of instructions for a process into a written language that a computer can interpret.
> {: .solution}
{: .discussion}

So although programming and coding have different meanings, they are related. The goal of coding is to create the code that acts as a set of computer instructions for a small part of a project. The goal of programming on the other hand, is to produce programs that are complete and ready to use software products.

## Pseudocode, Code, and Algorithms... Oh My!
Although the differences seem small, there are important distinctions that we can make between the concepts of pseudocode, code, and algorithms.

Everyone has some experienece with algorithms in their day-to-day life. For example, if you have ever cooked or done some task that requires you to follow instructions with a sequence of steps. 

> ## Discussion
>
> What are algorithms?
>
>> ## Solution
>>
>> Algorithms are the set of step-by-step instructions that explain how to solve a given problem.
> {: .solution}
{: .discussion}

Algorithms need to be represented by some form of language in order to be understood and shared with others. The process of writing pseudocode can be tremendously helpful for figuring out how to start developing code to solve a problem, or implement an algorithm.

> ## Discussion
>
> What is pseudocode?
>
>> ## Solution
>>
>> Pseudocode is the set of instructions for an algorithm written in a plain language.
> {: .solution}
{: .discussion}

As a first step before you begin developing an algorithm or writing any code, it is a good idea to write out the steps in a plain language. Let's look at an example of pseudocode for a simple algorithm to make tea:

1. Remove a teabag from the package
2. Put the teabag in a cup
3. Boil some water
4. Add the hot water to the cup
5. Allow the tea to steep for 5 minutes
6. Remove the teabag

> ## Challenge
>
> Write your own pseudocode for an algorithm to make buttered toast.
>> ## Solution
>>
>> 1. Take a slice of bread from the package
>> 2. Place the bread in the toaster
>> 3. Allow the bread to toast for 5 minutes
>> 4. Remove the toasted bread from the toaster
>> 5. Put the toasted bread on a plate
>> 6. Open the container of butter
>> 7. Grab a knife by the handle
>> 8. Dip the kife blade into the butter
>> 9. Apply the butter to the toasted side of the bread
> {: .solution}
{: .challenge}

The primary advantage to using pseudocode in your programming process is that it improves the readability of your algorithms. By first writing algorithms for programs in a plain language, it allows you to break down a complex problem into smaller and more manageable pieces for coding. Furthermore, it gives you the chance to easily identify the most complex and potentially troublesome portions for code development.

## Programming with Logic
A fundamental concept of computer programming, Boolean logic is the mathematical logic underlying Boolean algebra. In Boolean algebra mathematical expressions are evaluated to one of two values: true or false. Since an expression may only take on one of two values, Boolean logic is considered "two valued logic".

Note that an expression is a combination of logical operands and operators. In Boolean logic the operands are statements that can be proven true or false, and the operators are the logical AND, OR and NOT.

> ## Discussion
>
> What are some examples of boolean expressions?
>
> Hint: Use comparison operators (<, >, =, >=, <=, !=) to add complexity to your expressions!
>
>> ## Solution
>>
>> 1. It is raining and it is Spring
>> 2. I have green eyes and I am not 4 feet tall
>> 3. My cat is hungry or my cat is cute
>> 4. The temperature is < 32 degrees Fahrenheit and it is snowing 
> {: .solution}
{: .discussion}

We can combine boolean expressions with control statements to set how programs will complete a task. Control statments allow you to make have flexible outcomes by selecting which section of codes are executed. 

The three primary types of control statements are: 
- Sequential statmenetes are executed in the default ordering
- Iterative statements control the number of times a block of code is executed
- Conditional statements control which block of code is executed, and which are not

Iterative statements are defined using either FOR or WHILE loops.

> ## Discussion
>
> What does the following FOR loop output?
>
> ~~~
> for (x in 1:5) {
>   print(x)
> }
> ~~~
>
> Pseudocode:
> 1. For each x in the sequence 1, 2, 3, 4, 5 
> 2. Print the value of x
>
>> ## Solution
>>
>> The FOR loop outputs the current value of x after each iteration:
>> - 1
>> - 2
>> - 3
>> - 4
>> - 5
> {: .solution}
{: .discussion}

> ## Discussion
>
> What does the following WHILE loop output?
>
>> x <- 1
>> while (x < 3) {
>>   print(x)
>>   x <- i + 1
>> }
> {: .language-r}
>
> Pseudocode:
> 1. Set x equal to 1
> 2. While x is less than 3 
> 3. Print the value of x
> 4. Incrememnt the value of x by 1
>
>> ## Solution
>>
>> The WHILE loop outputs the current value of x after each iteration:
>> - 1
>> - 2
>> - 3
> {: .solution}
{: .discussion}

Conditional statements are defined using combinations of the IF... THEN... ELSE format. The most simple form of conditional statement is the IF... THEN form.

> ## Discussion
>
> What does the following IF... THEN conditional statement output?
>
>> x <- 3
>> if (x < 6) {
>>   print(x)
>> }
> {: .language-r}
>
> Pseudocode:
> 1. Set x equal to 3
> 2. If x is less than 6 
> 3. Print the value of x
>
>> ## Solution
>>
>> The IF... THEN statement outputs the value of x if it is less than 6:
>>
>> 3
> {: .solution}
{: .discussion}

The next type of conditional statement adds a level of complexity with the IF... THEN... ELSE format.

> ## Discussion
>
> What does the following IF... THEN... ELSE conditional statement output?
>
>> x <- 7
>> if (x < 6) {
>>   print(x)
>> } else {
>> 	print("x is larger than or equal to 6")
>> }
> {: .language-r}
>
> Pseudocode:
> 1. Set x equal to 7
> 2. If x is less than 6 
> 3. Print the value of x
> 4. Else print "x is larger than or equal to 6"
>
>> ## Solution
>>
>> The IF... THEN... ELSE statement outputs the value of x if it is less than 6, else it prints a message:
>>
>> x is larger than or equal to 6
> {: .solution}
{: .discussion}

The final version of the conditional statement further combines the IF... THEN... ELSE tags to make a statememnt with many alternative outcomes.

> ## Discussion
>
> What does the following IF... THEN... ELSE conditional statement output?
>
>> x <- 7
>> if (x = 6) {
>>   print("x is equal to 6")
>> } else if (x > 6) {
>> 	print("x is larger than 6")
>> } else if (x < 6) {
>> 	print("x is less than 6")
>> }
> {: .language-r}
>
> Pseudocode:
> 1. Set x equal to 7
> 2. If x is equal to 6 
> 3. Print "x is equal to 6"
> 4. Else if x is larger than 6
> 5. print "x is larger than 6"
> 6. Else if x is less than 6
> 7. print "x is less than 6"
>
>> ## Solution
>>
>> The IF... THEN... ELSE statement outputs the value of x if it is less than 6, else it prints a message:
>>
>> x is larger than 6
> {: .solution}
{: .discussion}

[codingProgramming]: https://techbiason.com/coding-vs-programming/

{% include links.md %}

