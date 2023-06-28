# CroMaSt: A workflow for domain family curation through cross-mapping of structural instances between protein domain databases

CroMaSt (<span style="color:red">**Cro**</span>ss <span style="color:red">**Ma**</span>pper of domain <span style="color:red">**St**</span>ructural instances) is an automated iterative workflow to clarify domain definition by cross-mapping of domain structural instances between domain databases. CroMaSt (for Cross-Mapper of domain Structural instances) will classify all structural instances of a given domain type into 3 different categories (core, true and domain-like). 


## Requirements
1. [Conda](https://docs.conda.io/projects/conda/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. [Kpax](http://kpax.loria.fr/download.php)  
Download and install conda (or Miniconda) and Kpax by following the instructions from their official site.


## Get it running 
(Considering the requirements are already met)

1. Clone the repository and change the directory

```
git clone https://gitlab.inria.fr/capsid.public_codes/CroMaSt.git
cd CroMaSt
```

2. Create the conda environment for the workflow
```
conda env create --file yml/environment.yml
conda activate CroMaSt
```

3. Change the path of variables in paramter file
```
sed -i 's/\/home\/hdhondge\/CroMaSt\//\/YOUR\/PATH\/TO_CroMaSt\//g' yml/CroMaSt_input.yml 
```

4. Create the directory to store files from PDB and SIFTS (if not already)
```
mkdir PDB_files SIFTS
```

5. Download the source input data
```
cwl-runner Tools/download_data.cwl yml/download_data.yml
```

## Basic example

### 1. First, we will run the workflow for the KH domain with family identifiers `RRM_1` and `RRM` in Pfam and CATH, respectively.
Run the workflow -

```
cwl-runner --parallel  --outdir=Results/  CroMaSt.cwl yml/CroMaSt_input.yml
```

### 2.  Once the iteration is complete, check the `new_param.yml` file from the `outputdir` (Results), if there is any family identifier in either `pfam` or `cath`; run the next iteration using following command (Until there is no new families explored by workflow) -

```
cwl-runner --parallel  --outdir=Results/  CroMaSt.cwl Results/new_param.yml
```
  
### **Extra:** Start the workflow with multiple families from one or both databases  
If you would like to start the workflow with multiple families from one or both databases, then simply add a comma in between two family identifiers. 
```
pfam: ['PF00076', 'PF08777']
cath: ['3.30.70.330']
```

- **Pro Tip**: Don't forget to give different path to `--outdir` option while running the workflow multiple times or at least move the results to some other location after first run.

## Run the workflow for protein domain of your choice  
### 1. You can run the workflow for the domain of your choice by simply changing the family identifers in `yml/CroMaSt_input.yml` file.

Simply replace the following values of family identifiers (for pfam and cath) with the family identifiers of your choice in `yml/CroMaSt_input.yml` file. 
```
pfam: ['PF00076']
cath: ['3.30.70.330']
```



## Data files used in current version are as follows:
**Files in Data directory can be downloaded as follows**:

1. File used from Pfam database: [pdbmap.gz](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/pdbmap.gz)  

2. File used from CATH database: [cath-domain-description-file.txt](ftp://orengoftp.biochem.ucl.ac.uk:21/cath/releases/latest-release/cath-classification-data/cath-domain-description-file.txt)  

3. Obsolete entries from RCSB PDB
[obsolete_PDB_entry_ids.txt](https://data.rcsb.org/rest/v1/holdings/removed/entry_ids)  


CATH Version - 4.3.0 (Ver_Date - 11-Sep-2019) [FTP site](ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/)  
Pfam Version - 33.0 (Ver_Date - 18-Mar-2020) [FTP site](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/)  

## Reference
```
Poster - 
1. Hrishikesh Dhondge, Isaure Chauvot de Beauchêne, Marie-Dominique Devignes. CroMaSt: A workflow for domain family curation through cross-mapping of structural instances between protein domain databases. 21st European Conference on Computational Biology, Sep 2022, Sitges, Spain. ⟨hal-03789541⟩

```

## Acknowledgements
This  project  has  received  funding  from  the  Marie  Skłodowska-Curie Innovative Training Network (MSCA-ITN) RNAct supported by European Union’s Horizon 2020 research and innovation programme under granta greement No 813239.

