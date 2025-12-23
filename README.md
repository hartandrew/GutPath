# GutPath: A Single-Cell Atlas of the Ileum during Diverse Infections

## 📌 Introduction
**GutPath** is a comprehensive single-cell atlas of the mouse ileum and draining mesenteric lymph node (MLN) during six immunologically distinct acute infections and provides a much-needed resource for the immunology community — enabling studies of ileal biology across infectious and inflammatory conditions. This atlas captures **over 500,000 cells from the ileum and MLN** during infection or colonization by fungal (*Candida albicans*), bacterial (*Yersinia pseudotuberculosis* and Segmented filamentous bacteria *Savagella spp*), viral (Murine norovirus CR6 strain), protozoan (*Cryptosporidium parvum*), and helminth (*Nippostrongylus brasiliensis*) organisms. While single-cell technologies have advanced mucosal immunology, there is an unmet need for large-scale profiling of the ileum across a wide range of infections to understand how diverse microbes reshape epithelial and immune responses.

The **ileum** plays a central role in the absorption of vitamin B12, bile acids, fat-soluble vitamins (A, D, E, K), water, sodium, and chloride. The ileum is also a complex immunological site where tolerance to microbial commensals and food antigens must be balanced by appropriate immunological responses to harmful pathogens. The **intestinal epithelium** — composed of enterocytes, goblet cells, Paneth cells, tuft cells, and others — is forefront as both a critical regulator of nutrient absorption and as a barrier to pathogen entry. Using GutPath, we demonstrate **pathogen-specific enterocyte transcriptional programs**, including a novel **pro-inflammatory enterocyte state** associated with *Yersinia* infection and localizing to pyogranuloma structures rich in immune infiltrate.

### 🌐 Public Resources

GutPath is **publicly accessible** in several formats and can be easily explored as a Cell x Gene object hosted through [GutPath.org](https://gutpath.org). GutPath can also be utilized as a reference for cell annotation through cell state marker profiles (DataS1 and DataS2 of our manuscript) or by accessing our the data directly as a downloadable seurat object on [GutPath.org](https://gutpath.org).

---

## Citation
If you find our work useful, please cite us. 
This work is currently under peer review and is available as a preprint 

## Methods

### Infections
We infected mice with six distinct pathogens, all orally administered (except subcutaneous *Nippostrongylus* L3 larvae), and collected tissues during early infection:

| Pathogen | Route | Timepoint |
|--------|-------|---------|
| Naive | NA | NA |
| _Cryptosporidium parvum_ | Oral gavage | 4 dpi |
| _Yersinia pseudotuberculosis_ | Oral gavage | 5 dpi |
| _Nippostrongylus brasiliensis_ | Subcutaneous L3 larvae | 5 dpi |
| _Candida albicans_ | Oral gavage | 7 dpi |
| Segmented filamentous bacteria (SFB) | Oral gavage | 7 dpi |
| Murine norovirus CR6 (MNV) | Oral gavage | 7 dpi |

Tissues were collected from:
- Ileum epithelial layer (**IEC**)
- Ileum lamina propria (**LP**)
- Draining mesenteric lymph node (**MLN**)

### Sample Preparation & CITE-seq
The ileal draining mesenteric lymph node and the distal 3rd of the small intestine isolated. Peyer's patches were excised and discarded. Tissue dissociation using enzymatic and mechanical methods were performed to generate single cell suspensions. CITE-seq libraries containing single cell 10X 3' RNA and data for ~100 Abs from TotalSeq B Mouse kit (Biolegend) were constructed and sequenced (NEXTSeq 2000).


## 📁 Repository Files

### Raw Data
- `TotalSeq_B_Universal_Cocktail_v1_140_Antibodies.xlsx`  
  → Antibody clones used in CITE-seq.

### Analysis Scripts
- `Scripts/Processing_00_CellRangerMetrics.Rmd`  
  → Visualize the raw Cell Ranger output quality control statistics (Figure S1A & Figure S1B).
- `Scripts/Processing_1.r`
  → Load samples into R, filter, and merge.
- `Scripts/Analysis_1_Clustering.r`
  → Normalizes and scales data, performs PCA, clustering, integration, cell cycle scoring, doublet removal, and WNN dual modality analyses.
- `Scripts/Analysis_2_Annotation.r`
  → Abbreviated script of automated and manual annotation of the GutPath atlas.
- `Scripts/Analysis_3_Exploration_FigurePlotting.r`  
  → Explores and graphs the annotated MLN and Ileum data sets.
  → Figure 2A, Figure 2B, Figure S2A, Figure S2B, Figure 2D, Figure S2G, Figure 2C, Figure S2F, Figure S2D, Figure S2E, Figure S2C, Figure S1D, Figure S3A, Figure 2E, Figure 3A, Figure 3B
- `Scripts/Analysis_4a_MarkerGenes.r`
  → Differential Expression based marker gene creation for GutPath data.
  → Data S1.
- `Scripts/Analysis_4b_python_nsForest.ipynb`
  → Random forest-based minimal and sufficient marker gene analysis for GutPath data.
  → Data S2.
- `Scripts/Analysis_5_DEG_Enrichment.rmd`
  → MAST based DEG analysis was performed and the results analyzed. Next, UCell enrichment analysis of gene sets was performed
  → Figure S5A, Figure 4A, Figure 4B, Figure 4C, Figure 3C , Figure 3D, Figure 3F, Figure 3H, Figure S4E, Figure 3E, Figure 3G , Figure S4A, Figure S4B, Figure S4C, Figure S4D
- `Scripts/DEG_MAST.r`  
  → Accompanies Analysis_5_DEG_Enrichment.rmd for MAST DEG analysis
- `Scripts/Hallmark_UCell_enrichment.r`
  → Accompanies Analysis_5_DEG_Enrichment.rmd for gene set enrichment analysis
- `Scripts/Analysis_6_Enterocyte_lineage.rmd`
  → Focuses on enterocyte lineage within GutPath. Pseudotime is performed as well as pseudotime modeled DEG analysis via Tradeseq. Investigations of the distinct *Yersinia* subset of enterocytes are performed
  → Figure S6C, Figure 6A, Figure 4D, #Figure S5C, Figure 4F, #Figure 4G, Figure 4H, Figure 6B, Figure S6A, Figure S6B , Figure 6C, Figure 6E, Figure 6D, Figure S5B, Figure 6F, Figure S6D
- `Scripts/Analysis_7_scCellfie.ipynb`
  → Metabolic activity for single cells was inferred using scCellfie
  → Figure 5A, Figure 5B, Figure 5C
- `Scripts/Analysis_8_CellChat_Nippostrongylus.R`  
  → Cell Cell Communication analysis is performed using Cell Chat
  → Figure 5E Figure 5F Figure 5D
- `Scripts/Proseg_Process`
  → Description fo the Xenium analysis and data processing (cell segmentation by Proseg) following raw data generation
- `Scripts/Analysis_9_Xenium.R`
  → Xenium data is analyzed, clustered, annotated, and expression evaluated. Cell Chat analysis is performed.
  → Figure 7A, Figure S7C, Figure S7D, Figure 7B, Figure 7C, Figure 7D, figure 7E, #Figure 7F, Figure 7G, Figure 7H, Figure 7I Figure 7J
- `Scripts/Analysis10_Xenium_UCell_Yp.R`
  → Secondary script needed to assist with Analysis_9_Xenium.R

---

## 🧩 CELLxGENE Fields

Downsampled data are available for direct exploration utilizing CellxGene VIP and can be accessed via [GutPath.org](https://gutpath.org)

| Field Name | Description | 
| :--- | :--- | 
| **Cell Cycle Phase** | The cell cycle phase (G1, S, or G2/M) assigned to the cell by Seurat's cell cycle scoring function. | 
| **Annotation 1 Lineages** | The annotation of cells at low resolution defining overarching cell lineages. | 
| **Annotation 2 Cell Types** | The annotation of cells at intermediate resolution defining common cell types. | 
| **Annotation 3 Cell States** | The annotation of cells at fine resolution defining individual cell states. | 
| **Infection** | The organism the animal was infected with at the time of sample collection. | 
| **Mouse** | The individual mouse ID. Useful for comparing biological replicates within an infection model. | 
| **Tissue** | The anatomical niche represented by the sample (duodenum or ileum & intestinal epithelial cells or lamina propria) (mesenteric lymph node). | 
| **Sample Type** | The anatomical niche isolation of the data (intestinal epithelial cells or lamina propria) (mesenteric lymph node). | 
| **NCount_RNA** | A measure of the transcript counts per cell. | 
| **nFeature_RNA** | The number of unique genes identified per cell. | 
| **Percent Mitochondrial** | The proportion of reads mapping to mitochondrial genes. Primarily for quality control. | 
| **Percent Ribosomal** | The proportion of reads mapping to ribosomal genes. Primarily for quality control. | 
| **Sample Name** | The unique identifier for the mouse/tissue sample combination. | 


## 📡 Data Access

- Raw data: [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) (Accession: PRJNA1378118)
- Processed Seurat objects: [GutPath.org](https://gutpath.org)
- Interactive visualization: [GutPath.org](https://gutpath.org)

---


## 📧 Contact

For questions or collaboration:  

📧 hartand@upenn.edu  
🌐 [GutPath.org](https://gutpath.org)  

> *This atlas enables research into mucosal immunology and epithelial plasticity during infection.*
