 Single-Cell RNA-seq Analysis of NSCLC

From Raw Data to Biological Insights (GSE131907)

This repository provides a complete, beginner-friendly workflow for performing single-cell RNA sequencing (scRNA-seq) analysis using publicly available non-small cell lung cancer (NSCLC) data.

The tutorial demonstrates how raw UMI count matrices can be transformed into meaningful biological insights using reproducible methods implemented in R.

⸻

🎯 Project Goals

The aim of this project is to:

• reconstruct a cellular atlas of NSCLC
• identify tumor vs normal transcriptional differences
• explore immune and epithelial heterogeneity
• perform differential expression analysis
• interpret functional pathways associated with tumor progression

⸻

📂 Dataset

GSE131907 – NSCLC Lung Cancer Atlas

Includes:
	•	Tumor lung tissue
	•	Normal lung tissue
	•	Lymph node samples
	•	Effusion samples

Total cells analyzed: ~180,000+
Genes profiled: ~29,000

⸻

🛠 Tools Used

• R
• Seurat
• ggplot2
• gprofiler2
• patchwork

The workflow avoids heavy dependencies and focuses on stable, reproducible methods.

⸻

🔬 Analysis Workflow

1️⃣ Data Preprocessing

✔ Load UMI matrix
✔ Match cell barcodes with annotation
✔ Create Seurat object
✔ Compute QC metrics:
	•	nFeature_RNA
	•	nCount_RNA
	•	percent mitochondrial reads

⸻

2️⃣ Normalization & Dimensionality Reduction

✔ SCTransform normalization
✔ Identify highly variable genes
✔ PCA for feature reduction
✔ UMAP for visualization
✔ Graph-based clustering

⸻

3️⃣ Cell Type Annotation

Using provided metadata:
	•	Epithelial cells
	•	Myeloid cells
	•	T/NK cells
	•	Fibroblasts
	•	Endothelial cells

⸻

4️⃣ Sub-Atlas Construction

Focused analyses were conducted on:

🧫 Epithelial Cells
• Tumor vs normal comparison
• Identification of heterogeneous epithelial states
• Transcriptional programs associated with tumor progression

🧬 Myeloid Cells
• Monocyte/macrophage heterogeneity
• Tumor-associated immune remodeling

⸻

5️⃣ Differential Expression Analysis

Comparisons performed:
• Tumor vs Normal (Epithelial)
• Tumor vs Normal (Myeloid)

Visualizations generated:
✔ Volcano plots
✔ Violin plots
✔ DEG tables

⸻

6️⃣ Functional Enrichment Analysis

Pathway enrichment revealed:

• Immune activation pathways
• Inflammatory signaling
• Epithelial plasticity
• Tumor microenvironment remodeling

⸻

📊 Figures Generated

• UMAP atlas by tissue origin
• Cell composition plots
• Epithelial sub-atlas
• Myeloid clustering
• Volcano plots
• Functional enrichment plots

All figures are produced in publication-ready format.

⸻

🧠 Key Biological Insights

The analysis demonstrates:

• Tumor epithelial cells exhibit distinct transcriptional programs
• Myeloid populations show tumor-associated activation signatures
• Immune and inflammatory pathways are enriched in tumor states
• Cellular heterogeneity reflects tumor microenvironment dynamics

⸻

🚀 Who Is This Tutorial For?

✔ Beginners in single-cell analysis
✔ Researchers transitioning from bulk RNA-seq
✔ Students learning reproducible genomics workflows
✔ Anyone interested in tumor microenvironment analysis
