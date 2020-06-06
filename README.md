# AbDesing_for_enzymes
Scripts and data to run AbDesign as described in Tools for protein science 2020. Below I'll describe the different steps, data needed and example command lines to genearte a repertoire of structures for the GH10 xylanase family.   
To run the following, you are required to install:
  * python >3.5 
  * Rosetta (free for academic users)
  * pymol 

## 1 Structures and Segmentation scheme
The list of xylanase structures was taken from [here](http://www.cazy.org/GH10_structure.html) and downloaded from the PDB.  All structures must be aligned to the template structure, which can be done using the **alignto** command in PyMOL.

## 2 Backbone fragments databases
### 2.1 Initial crude fragment extraction
The segmentation points in this step are (numbering corrspond to template structure pdbID: 3w24):
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
### 2.2 Refinmenet of fragment alignment to template
In this step, fragments from the same segments are refined according to the template sructure to all have the exact same ends of the fragment to better assemble later (details in the paper).
In this step the template's fragments start & end at the following positions: (in different systems, make sure the fragments do not overlap)
| Segment | start | end |
| :---: | :---: | :---: |
| 1 | 19 | 45 |
| 2-4 | 46 | 186 |
| 5-6 | 187 | 249 |
| 7-8 | 250 | 320 |

In the command below, the fragment *x_fragment1.pdb* in being matched to the fragment 1 from the template
```bash 
pymol -c ./backbone_database/find_best_match.py data/blade1_template.pdb x_fragment1.pdb
```
The resulting pdb file will be located at *pdbs/fbm_x_fragment1.pdb*.
