# EPICv2-methylation-analysis
This is a pipeline for the analysis of EPICv2 microarray data using R 4.4.1. 


It has two phases: preprocess and analysis:

## Preprocessing
We start from the IDATs files plus the SampleSheet. 
1. The data is imported and annotated using EPICv2 annotation files.
2. Quality control of the samples
3. Normalization
4. Probes filtering
   
## Analysis
1. MDS analysis: Evaluation of batch effects and assesment of the impact of age
2. DMPs analysis
3. DMRs analysis
4. Functional enrichment analysis

