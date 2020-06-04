# AbDesing_for_enzymes
Scripts and data to run AbDesign as described in Tools for protein science 2020. Below I'll describe the different steps, data needed and example command lines to genearte a repertoire of structures for the GH10 xylanase family.   
You are required to install python >3.5 and Rosetta (free for academic users).  

## 1 Structures and Segmentation scheme
The list of xylanase structures was taken from [here](http://www.cazy.org/GH10_structure.html) and downloaded from the PDB.  All structures must be aligned to the template structure, which can be done using the **alignto** command in PyMOL.

## 2 Backbone fragments databases
The segmentation points in this step are:
| Segment | start | end |
| :---: | :---: | :---: |
| 1 | 19 | 47 |
| 2-4 | 44 | 189 |
| 5-6 | 183 | 256 |
| 7-8 | 246 | 324 |

The following command is used to extract a fragment structure corresponding to the fragment starting with *start_res* and ends with *end_res* in the *source* structure (i.e. the template). This command should be used for each protein in the family for each of the segments.  
For example, extracting the first fragment from structure X.pdb would use the following command:
```bash 
ROSETTA_SCRIPTS @backbone_database/flags_cod -s X.pdb -out:path:pdb  -out:path:score scores pdbs -parser:protocol backbone_database/cut_out_domain.xml -parser:script_vars source=data/3w24_template.pdb.gz start_res=19 end_res=47 
```
The fragment will be outputted to a directory named pdbs (needs to be created by you before using this command).
