### Data Provided

File cell-count.csv contains cell count information for various immune cell populations of each patient sample. There are five populations: b_cell, cd8_t_cell, cd4_t_cell, nk_cell, and monocyte. Each row in the file corresponds to a biological sample.

The file also includes sample metadata such as sample_id, indication, treatment, time_from_treatment_start, response, and gender.

### Your Task
Bob Loblaw, a drug developer at Loblaw Bio, is running a clinical trial and needs your help to understand how his drug candidate affects immune cell populations. Your job is to:

Design a Python program that meets Bob’s analytical needs, as outlined in Parts 1-4 below.

Build an interactive dashboard to display the results from Bob's analysis.

### Part 1: Data Management

Using the data provided in cell-count.csv, your first task is to:

Design a relational database schema (using SQLite) that models this data effectively.

Create a Python script named "load_data.py" in the root directory of your repository that:

Initializes the database with your schema.

Loads all rows from cell-count.csv.

Requirements:
The script must be named `load_data.py` and located in the root directory (not in subdirectories like `src/`).
 - When executed with `python load_data.py`, it should create a SQLite database file (`.db` extension) in the repository root.
- The script should be executable directly without command-line arguments or module-style execution (`python -m`).

### Part 2: Initial Analysis - Data Overview

Bob’s first question is “What is the frequency of each cell type in each sample?” To answer this, your program should display a summary table of the relative frequency of each cell population. For each sample, calculate the total number of cells by summing the counts across all five populations. Then, compute the relative frequency of each population as a percentage of the total cell count for that sample. Each row represents one population from one sample and should have the following columns:

sample: the sample id as in column sample in cell-count.csv

total_count: total cell count of sample

population: name of the immune cell population (e.g. b_cell, cd8_t_cell, etc.)

count: cell count

percentage: relative frequency in percentage

### Part 3: Statistical Analysis

As the trial progresses, Bob wants to identify patterns that might predict treatment response and share those findings with his colleague, Yah D’yada. Using the data reported in the summary table, your program should provide functionality to:

Compare the differences in cell population relative frequencies of melanoma patients receiving miraclib who respond (responders) versus those who do not (non-responders), with the overarching aim of predicting response to the treatment miraclib. Response information can be found in column "response", with value "yes" for responding and value "no" for non-responding. Please only include PBMC samples.

Visualize the population relative frequencies comparing responders versus non-responders using a boxplot of for each immune cell population.

Report which cell populations have a significant difference in relative frequencies between responders and non-responders. Statistics are needed to support any conclusion to convince Yah of Bob’s findings. 

### Part 4 Data Subset Analysis: 

Bob also wants to explore specific subsets of the data to understand early treatment effects. Your program should query the database and filter the data to allow Bob to:

Identify all melanoma PBMC samples at baseline (time_from_treatment_start is 0) from patients who have been treated with miraclib. 

Among these samples, extend the query to determine:

How many samples from each project

How many subjects were responders/non-responders 

How many subjects were males/females


### answer this question
Considering Melanoma males of all sample and treatment types, what is the average number of B cells for responders at time=0? Use two decimals (XXX.XX).

### Submission Requirements

Please submit your solution as a GitHub repository link.

Your project should include:

Your Python program with all accompanying files

Any input or output files generated

A README.md with:

Any instructions needed to run your code and reproduce the outputs (We will run your code using GitHub Codespaces).

An explanation of the schema used for the relational database, with rationale for the design and how this would scale if there were hundreds of projects, thousands of samples and various types of analytics you’d want to perform.

A brief overview of your code structure and an explanation of why you designed it the way you did.

A link to the dashboard.

A Makefile in the root directory. We will use this to automatically grade your submission using GitHub Codespaces. Your Makefile must implement the following three targets exactly as named:
make setup: Installs all necessary dependencies for your project (e.g., from a requirements.txt, environment.yml, or pyproject.toml).
make pipeline: Executes your entire data pipeline sequentially from start to finish without any manual intervention. When our grader runs this command, it should initialize the database, load the data (Part 1), and generate all required output tables and plots (Parts 2-4). (Note: You may use pure Python, bash scripts, Snakemake, or any other orchestration tool, as long as make pipeline triggers the complete execution).
make dashboard: Starts the local server for your interactive dashboard.