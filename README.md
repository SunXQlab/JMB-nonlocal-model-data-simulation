# SPATIAL-TUMOR-HETEROGENEITY  
**Mathematical Dissection of Tumor Phenotypic Heterogeneity in a Spatial Data-Informed Nonlocal Reaction-Diffusion Model**  
accompanying code for the *Journal of Mathematical Biology* submission  
Li F., Li H., Lin L., Sun X.*  

---

## 🧪 What the codes do

| Script | Input | Output | Purpose |
|--------|-------|--------|---------|
| `interpolation.m` | deconvolution results: `spot.csv`, `st_ct_proportion.csv` | `phenotype.csv`, `malignant.csv`, `endothelial.csv`, `vascular_pos.csv` | Natural-neighbour interpolation of ST-deconvolution data → smooth 2-D fields for **initial data** (Fig. 5a) |
| `ADEI.cpp` | the above `*.csv` files| `tumor.csv`, `oxygen.csv` | Alternating-Direction Explicit-Implicit (ADEI) solver for the coupled non-local reaction-diffusion system (eq. 1.1) (Fig. 5b)|
| `MoranI.py` | any 2-D csv (e.g. `malignant.csv`) | `MoranIndex.csv` | Moran’s Index of tumour-cell number, vascular distribution & local mean phenotypic state (Fig. 6) |
