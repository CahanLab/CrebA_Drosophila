This directory contain all the analysis of integrating CrebA ChIP-seq data with scRNA-seq results. The entire pipeline is in the snakemake file inside `workflow` folder. You can find the requirements for the Python packages in the `python_requirements.txt`. 

The CLI tools used are below. You can also install them by installing the entire [MEME Suite](https://web.mit.edu/meme_v4.11.4/share/doc/install.html). 
| Package     | Version     |
| ----------- | ----------- |
| bedtools    | 2.31.1      |
| dreme       | 5.5.4       |
| tomtom      | 5.5.4.      |


Before running the analysis code, please download the zip file from [here](https://cnobjects.s3.us-east-1.amazonaws.com/drosophila_2023/CrebA_Drosophila/input.zip) and unzip the content into a directory called `input`. 
