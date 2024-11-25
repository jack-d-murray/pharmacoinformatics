# Pharmacoinformatics Respository
This page shares the code and data associated with the submitted manuscript, *"Establishing a Pharmacoinformatics Repository of Approved Medicines: A Database to Support Drug Product Development"*, The manuscript is currently under review. 

[Explore the database without any code using the Streamlit App.](https://pharmacoinformatics.streamlit.app/)

## Background
The need for a machine-readable drug product library which consolidates many aspects of formulation design and performance remains largely unmet. This study presents a scripted, reproducible approach to database curation and explores its potential to streamline oral 
medicine development. Web scraping was used to retrieve a European Public Assessment Report Product (EPAR) Information File for all drug products containing a small molecule drug substance (<1000 g mol<sup>-1</sup>). Regular expressions isolated relevant information. External 
curated databases, including PubChem and ChEMBL, were integrated through programmatic access. Analysis and modelling was performed on the final database.

## Project Structure
- <b>notebooks</b>:
  - <b>database_construction.ipynb</b>: Construction of the database, including web scraping, regular expression matching, and applying curation pipelines.
  - <b>chemical_space.ipynb</b>: Analysing the chemical space of drugs captured by our database.
  - <b>excipient_selection.ipynb</b>: Exploratory hypothesis testing and association rule learning applied to formulation data.
  - <b>fa_model.ipynb</b>: Constructing binary classifiers of human oral fraction absorbed.
  - <b>constructiontools.py</b>: Module containing functions used in database_construction to enhance readability.
  - <b>curation_guidelines.pdf</b>: General rules for curating excipients. Note that decisions were also made on a case-by-case basis. 
- <b>pharmacoinformatics_database</b>: Final database consisting of five `.csv` file linked by drug, excipient, and product identifiers. Includes `data_dictionary.csv`.
- <b>results</b>: Output directory for files describing the results of analysis.    
- <b>csv_files</b>: Intermediate `.csv` files used for manual curation steps.
- <b>figures</b>: Directory for project figures and visualisations.
- <b>requirements.txt</b>: Provided to list the versions of all packages/libraries used.

## Usage
This repository is designed to share code and data associated with a scientific manuscript and is not a contemporaneous source of drug product information. For information regarding drug products, please visit the website of the relevant regulator. 

## Contributing
Contributions are always welcome, please feel free to reach out with any suggestions or if you notice any errors in data curation.

## Authors and Contact
- Jack D. Murray<sup>1</sup> (jackmurray at ucc.ie)
- Harriet Bennett-Lenane<sup>1</sup>
- Patrick J. O'Dwyer<sup>1</sup>
- Brendan T. Griffin<sup>1</sup>

1. School of Pharmacy, University College Cork, Ireland

## Reference
Manuscript submitted and under review in *Molecular Pharmaceutics* Virtual Special Issue: Computational Methods in Drug Delivery

## License
This project is licensed under the [MIT License](https://opensource.org/licenses/MIT). See the [LICENSE](LICENSE) file for more details.
