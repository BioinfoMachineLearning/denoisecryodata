# A Labeled Dataset for AI-based Cryo-EM Map Enhancement

An open-source dataset for cryo-EM density map denoising comprising 650 high-resolution experimental maps paired with three types of generated label maps: regression maps capturing idealized density distributions, binary classification maps distinguishing structural elements from background, and atom-type classification maps.

## Dataset Download
To keep the data files of denoise dataset permanent, we published all data to the Harvard Dataverse (https://doi.org/10.7910/DVN/CI0J2B), an online data management and sharing platform with a permanent Digital Object Identifier number for each dataset. 

## Description of the dataset
Dataset Link : https://doi.org/10.7910/DVN/CI0J2B .

The dataset can be accessed using the above download link. The dataset follows the format described below, and please read the subsequent section describing what these data mean.

    │── 6bco
        │-- 6bco.mrc
        |-- 6bco.pdb1
        |-- 6bco_classification_situs.mrc
        |-- 6bco_classification_types_situs.mrc
        |-- 6bco_regression_situs.mrc
        |-- 6bco_situs_simulated.mrc
 
    │── 7ki6
        │-- 7ki6.mrc
        |-- 7ki6.pdb1
        |-- 7ki6_classification_situs.mrc
        |-- 7ki6_classification_types_situs.mrc
        |-- 7ki6_regression_situs.mrc
        |-- 7ki6_situs_simulated.mrc

    │── 7o3h
        │-- 7o3h.mrc
        |-- 7o3h.pdb1
        |-- 7o3h_classification_situs.mrc
        |-- 7o3h_classification_types_situs.mrc
        |-- .
        |-- .

    .
    .
    .

As shown in the example data format above, each individual directory for the cryo-EM density map (directory name is the protein PDB Code) provides the following data files:

- ``6bco.mrc`` : The deposited experimental cryo-EM density map. The filename is the PDB code of the cryo-EM density map.
- ``6bco.pdb1`` : The atomic biological assembly of the macromolecule.
- ``6bco_situs_simulated.mrc`` : The simulated cryo-EM density map of ``6bco.pdb1`` generated using Situs package.
- ``6bco_classification_situs.mrc`` : The classification label map, neighboring voxels are assigned a value of 2, central atomic positions (value 1) and the background (value 0).
- ``6bco_classification_types_situs.mrc`` : The atom type classification label map, values to voxels based on the corresponding atom type: Cα (1), Cβ (2), carbonyl carbon (3), oxygen (4), and nitrogen (5).
- ``6bco_regression_situs.mrc`` : This map contains continuous density values derived from the simulated cryo-EM density map. At each voxel position corresponding to an atom in the PDB structure, we assign the density value from the simulated map, providing a noise-free target for regression-based learning.


# Programs to generate the dataset
```
python3 generate_labels_all_conv.py
```

# Evaluating map quality with Phenix mtriage
The below runs phenix.mtriage to compute the FSC scores. The webpage of Phenix mtriage is available here: https://phenix-online.org/documentation/reference/mtriage.html
```
python3 evaluate_phenix_mtriage_exp.py
python3 evaluate_phenix_mtriage_reg.py
python3 evaluate_phenix_mtriage_sim.py
```
# Programs to extract the information from Phenix mtriage output
The FSC scores reported in the paper are generated using Phenix mtriage. Use the below to extract information in ``.csv`` format
```
python3 parse_get_resolution_from_log_file_exp.py
python3 parse_get_resolution_from_log_file_reg.py
python3 parse_get_resolution_from_log_file_sim.py
```


# Contact Information
If you have any question, feel free to open an issue or reach out to us: [ngzvh@missouri.edu](ngzvh@missouri.edu), [chengji@missouri.edu](chengji@missouri.edu).

## Citing this work
If you use the code or data in this package, please cite:

```bibtex
@article{giri2025labeled,
  title={A Labeled Dataset for AI-based Cryo-EM Map Enhancement},
  author={Giri, Nabin and Chen, Xiao and Wang, Liguo and Cheng, Jianlin},
  journal={Computational and Structural Biotechnology Journal},
  year={2025},
  publisher={Elsevier}
}

```
