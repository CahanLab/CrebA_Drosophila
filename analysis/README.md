### Analysis for scRNA-seq

- Current code: *main.R* calls our current analysis scripts in the correct order. You can find all the R scripts in `R/` folder. Before running the analysis, please download the zip file from [here](https://cnobjects.s3.us-east-1.amazonaws.com/drosophila_2023/CrebA_Drosophila/accessory_data.zip) and unzip the content into a directory called `accessory_data`.

Here are the packages required for the main analysis. 

| Package     | Version     |
| ----------- | ----------- |
| R           | 4.1.2       |
| Seurat      | 4.3.0       |
| harmony     | 0.1.1       |
| presto      | 1.0.0       |
| fgsea       | 1.20.0      |
| ggplot2     | 3.4.4       |
| dplyr       | 1.1.4       |
| cowplot     | 1.1.1       |
| dplyr       | 1.1.4       |
| magrittr    | 2.0.3       |
| Matrix      | 1.5.3       | 
| stringr     | 1.5.0       | 
