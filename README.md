# CDI-FMT-project

Welcome to the repository for the Clostridium Difficile Infection Fecal Microbiota Transplantation (CDI-FMT) study. This repository hosts the Python code utilized in the analysis for the CDI-FMT paper.

## Python Modules Overview

This project is structured into several modules, each designed to perform specific tasks related to the data analysis presented in the CDI-FMT paper.

### Helper
- **Functionality**: Contains basic modules to read in the OTU table and convert it into an operable Pandas DataFrame.

### Preprocessing
- **Dependencies**: Utilizes modules from the Helper.
- **Purpose**: Converts files into Pandas DataFrames for further analysis.

### Boxplot
- **Features**: Includes modules to generate boxplots. 
- **Specifics**: Currently, the function is tailored specifically for plotting the Simpson Diversity Index (SDI).

### Volcano Plot
- **Description**: Implements volcano plotting for all features.
- **Purpose**: This module acts as a feature selection step, identifying significant features for further analysis.

## Getting Started

To replicate the analysis:
1. Ensure you have Python installed along with the necessary libraries (Pandas, Matplotlib, etc.).
2. Clone this repository.
3. Run the scripts in the order specified by their dependencies.

For any issues or further inquiries, please open an issue in this repository.
