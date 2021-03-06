# AbDesign for enzymes
Scripts and data to run AbDesign as described in Tools for protein science 2021.   
Below I'll describe the different steps, data needed and example command lines to generate a repertoire of structures for the GH10 xylanase family. Steps 2-4 are using PDB accession number 4pud as an example but should be repeated for each of the structures of the chosen protein family.
To run the following, you are required to install:
  * python > 3.5 
  * Rosetta (free for academic users)
  * PyMOL 

## 1 Structures and Segmentation scheme
The list of xylanase structures was taken from [here](http://www.cazy.org/GH10_structure.html). The structures could be downloaded using the following command:   
```bash
./utils/get_pdb.py 3w24 4pud #space separated list of pdb IDs 
```
All structures must be aligned to the template structure, which can be done using the **alignto** command in PyMOL. Then, the segmentation points should be selected.   
A few suggestions for selecting segmentation points:
   1. different segments should correspond to modular fragments- in TIM-barrel for example it would be different beta-alpha blades 
   2. the cut points should be perfectly aligned among **all** fragments in order to generate continuous backbones
   3. If two segments have many interaction among them, consider combining them into a single larger segment
   4. A large segment could be breaken up into two smaller ones in order to generate more possible diversity

## 2 Structures idealization
To generate an idealized version of a protein structure, use (change the *-s* path to your structure):
```bash 
ROSETTA_SCRIPTS @idealization/flags -s idealization/4pud.pdb -out:prefix ideal_ 
```
The idealized version is ideal_4pud.pdb.gz.   
**Time estimation**: 30min on average per structure.

## 3 Backbone fragments databases
### 3.1 Initial crude fragment extraction
*Segmentation point selection*: After we define the segmentation scheme in step 1, here we start by a rough segmentation (to be further refined in following steps). Say we want a segment starting at residue X and ending at residues Y, we will start from X-3 and Y+3.   

The segmentation points in this step are (numbering corresponds to template structure pdbID: 3w24):
| Segment | start | end |
| :---: | :---: | :---: |
| 1 | 19 | 47 |
| 2-4 | 44 | 189 |
| 5-6 | 183 | 256 |
| 7-8 | 246 | 324 |

The following command is used to extract a fragment structure corresponding to the fragment starting with *start_res* and ends with *end_res* in the *source* structure (i.e. the template). This command should be used for each protein in the family for each of the segments.  
For example, extracting the first fragment from structure 4pud.pdb would use the following command:
```bash 
ROSETTA_SCRIPTS @backbone_database/flags_cod -s backbone_database/4pud.pdb -out:prefix blade1_  -parser:protocol backbone_database/cut_out_domain.xml -parser:script_vars source=template_data/3w24_template.pdb.gz start_res=19 end_res=47
```
The output pdb of the fragment will be located at: pdbs/blade1_4pud.pdb.gz (you need to create the pdbs folder before running the command).     
**Time estimation**: 30s or less per structure.

### 3.2 Refinement of fragment alignment to template
In this step, fragments from the same segment are refined according to the template structure, such that they all have the exact same ends of the fragment to better assemble later (details in the paper).  
  
*Segmentation point selection*: Here we are using the exact residues for the start and end of the segment. Make sure the segments are not overlapping (e.g. a segment cannot start at residue 47 if the previous one ended at 48).  
In this step the template's fragments start & end at the following positions: 
| Segment | start | end |
| :---: | :---: | :---: |
| 1 | 19 | 45 |
| 2-4 | 46 | 186 |
| 5-6 | 187 | 249 |
| 7-8 | 250 | 320 |

In the command below, the fragment from step 2.1 is being matched to the fragment 1 from the template
```bash 
pymol -c ./backbone_database/find_best_match.py template_data/blade1_template.pdb backbone_database/blade1_4pud.pdb.gz
```
The resulting pdb file will be located at *pdbs/fbm_blade1_4pud.pdb*.   
**Time estimation**: a few seconds per structure.

## 4 PSSM
First, you need to generate a PSSM for each of the structures of your protein family. We used PSSM generated by [PROSS](https://pross.weizmann.ac.il/step/pross-terms/) (details in the paper).  
The following command extracts a sub-PSSM, matching a fragment of the protein:
```bash
./pssm/cut_pssm_for_fragment.py backbone_database/fbm_blade1_4pud.pdb pssm/4pud.pssm
```
The output pssm is *fbm_blade1_4pud.pssm*.     
**Time estimation**: A full PSSM creation can take up to 1 hour (due to a BLAST search). Extracting a sub-PSSM takes a few seconds.

## 5 Torsion database
This step extracts the torsion angles from each of the fragment. Before running the command, we have to prepare a few files:
  1. **Directories**:  we need to create to the following directories to hold the output: ```mkdir pdbs scores  db```
  2. **Structures names**: rename the structures of the fragments to include only the pdbID, e.g.: 4pud.pdb or 4pud.pdb.gz
  3. **pdb_profile_match**: a file mapping names of pdbs. See example at *torsions_database/pdb_profile_match*
  4. **flags_pssm**: mapping pdb fragments to their PSSM files. See the format at *torsions_database/flags_pssm*
  5. **flags**: When running with a different protein family, change the path to the appropriate template structure, the catalytic residues to keep conformation and all other paths to the correct ones. 
  6. **splice_out.xml**: When running with a different protein family, change the segments section of the Splice mover to match your naming and segmentation scheme. The frm1 & frm2 tags correspond to the start and end of the protein, respectively, which are kept constant in all designs (i.e. residues before and after the first and last segment) 
  7. **template pdb file**: The splice mover in Rosetta writes to the PDB file the source pdb of the different segment. The input template structure should have at the bottom of the file the following lines (should correspond to the segments section in the xml file, see 6 above):   
     ```
     ##Begin comments##   
     segment_frm1 3w24_template   
     segment_blade1 3w24_template   
     segment_blade2_4 3w24_template      
     segment_blade5_6 3w24_template   
     segment_blade7_8 3w24_template   
     segment_frm2 3w24_template   
     ##End comments##  
  
*Segmentation point selection*: The segmentation points here are **one residue into the segment** relative to step 3.2.
| Segment | start | end |
| :---: | :---: | :---: |
| 1 | 20 | 44 |
| 2-4 | 47 | 185 |
| 5-6 | 188 | 248 |
| 7-8 | 251 | 319 |  

To generate the torsion database file for 4pud.pdb (i.e. the fragment's structure from step 3.2) use:
```bash
ROSETTA_SCRIPTS @torsions_database/flags -out:prefix 4pud_ -parser:script_vars source=torsions_database/4pud.pdb db=db/blade1_4pud.db start_res=20 end_res=44 current_segment=blade1
```
The torsion angles of the fragment will be generated at *db/blade1_4pud.db*. The fragment in the context of the template is located at *pdbs/4pud_3w24_template.pdb.gz* (not needed for later steps, but useful for debugging).   

For further explanations please [see](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/SpliceOut).     
**Time estimation**: 6H on average per structure.

### 5.1 Torsion database for N/C-termini regions
The above working scheme works great in case there is no important information encoded in the tail regions of the protein. In case these regions are needed (i.e starting from residue 19 is not an option for you), slight adjustments are needed.   
  1. Skip steps 3.1 and 3.2. Instead, cut manually the N/C-terminus regions needed out of an idealized pdb (step 2). Use segmentation points as described in step 3.2. See an example of a full-size tail at backbone_tails/4pud.pdb.
  2. Follow step 4 to extract a sub PSSM for the tail region. 
  3. **As in step 5**: create **Directories**, rename **Structures names**, prepare **pdb_profile_match** file and a **flags_pssm file**. Make sure all paths in **flags** file point to the appropriate directories.
  4. **splice_out_Ntail.xml**: change the segments section of the SpliceOutTail mover to match your naming and segmentation scheme. As opposed to step 5.0, frm1 and/or frm2 tags are not needed, as there are no “left behind” residues at the tails. “splice_out_Ntail.xml” is adjusted for an N-terminus tail. For a C-terminus, change the tail_segment="n" to tail_segment="c".  
  5. **template pdb file**: as in step 5.0, the input template structure should include a “comment” section at the bottom. Make sure to remove “frm1” and/or “frm2” as the current segmentation requires. 

An example of a command line:   
```
ROSETTA_SCRIPTS @backbone_tails/flags -out:prefix 4pud_ -parser:script_vars source=backbone_tails/4pud.pdb db=db/blade1_4pud.db current_segment=blade1 from_res=44
```

“from_res” is used both in N-and in C-terminus cases to indicate the single defined segmentation point. 

As before, the torsion angles of the fragment will be generated at db/blade1_4pud.db 
The last three digits every line in the database indicate the segmentation points chosen. In case of an N-terminus tail, change the third digit from the end, to the desired segmentation residue number. (Due to a bug, the default output is 1). I.e, “... 188 130.796 175.255 SER 1 0 0 4pud”, should become “...188 130.796 175.255 SER 44 0 0 4pud”. There is no need for such an adjustment in the case of a C-terminus tail.  


## 6 Assembly of backbones
Here we will generate a new backbone by combining different fragments. Input files to prepare:
  1. **Directories**:  we need to create to the following directories to hold the output: ```mkdir pdbs scores```
  2. **pdb_profile_match & flags_pssm**: as previous step
  3. **flags**: as before, change paths if using a different family
  4. **Torsion database**: for each segment, combine all torsions calculated in the previous step to a single file, each line is the torsion of a single fragment of this segment. For example see: *backbone_assembly/blade1.db*
  5. **splice_in.xml**: When running with a different protein family, change the segments section of the Splice mover to match your naming and segmentation scheme. The frm1 & frm2 tags correspond to the start and end of the protein, respectively, which are kept constant in all designs (i.e. residues before and after the first and last segment)
  6. **template pdb file**: see above at *Torsion database*
  
  Example command to generate a new backbone (change the pdbID in the entries to generate a backbone from different fragments):
  ```bash
ROSETTA_SCRIPTS @backbone_assembly/flags -s template_data/3w24_template.pdb.gz -out:prefix 4pud_4qdmB_1xyzA_1e5nB_ -parser:script_vars entry_blade1=4pud entry_blade2_4=4qdmB entry_blade5_6=1xyzA entry_blade7_8=1e5nB
```
### 6.1 Assembly of backbones that include N/C-termini regions
Follow the instructions as written in step 6.0. The only adjustment needed is in the .xml file where a SpliceIn mover needs to be changed to SpliceInTail. An example can be seen at backbone_tails/backbone_assembly/splice_inTail.xml

Command line to run spiceInTail: 
  ```bash
ROSETTA_SCRIPTS @backbone_tails/backbone_assembly/flags -s backbone_tails/3w24_template.pdb.gz -out:prefix 4pud_4qdmB_1xyzA_1e5nB_ -parser:script_vars entry_blade1=4pud entry_blade2_4=4qdmB entry_blade5_6=1xyzA entry_blade7_8=1e5nB 
```   
**Time estimation**: up to 20m per chimera.

## 7. Frequent debugging issues
  1. Misalignment of segments (overlaps, or discontinues backbone formation). 
  2. The *comments* section in the template pdb do not match the segments section in the .xml file.
  3. *flags_pssm* file includes a newline or a whitespace. Make sure to include full paths to pssm files.
  4. There is a difference in capitalization of the “Segment” section between spliceOut and SpliceIn movers. See example torsions_database/splice_out.xml vs. backbone_assembly/splice_in.xml    
  
  For any bugs/questions please write to: rosalie.lipsh@weizmann.ac.il
