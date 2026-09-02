

# all marker --------------------------------------------------------------
Epithelial_cell <- c("WFDC2","EPCAM", "SLPI", "EHF", "TMC5")#"CXADR", developing neuron system

Hepatocytes_feature <- c("APOC3", "APOA2", "ALB", "SAA1", "ORM1", "HP")
Hepatocytes_feature <- c(Hepatocytes_feature, "TF", "CYP3A4", "CYP2E1")#, "ASS1", "APOE")

Neuron_feature <- c( "GPC5", "NRXN1", "PCDH9", #astrocyte #G6PAM endothelial marker
                     "RBFOX1", "CADM2", "CSMD1", "KCNIP4")#, #exitatory neu
#"NRXN3", "CNTNAP2", "ROBO2") # inhibitory neu

Fibroblast_feature <- c("PDGFRA", "VIM", "COL4A1", "FN1", "COL1A1", "COL1A2", "COL5A1",  "LUM", "FBLN1") # VCAN
smooth_muscle <- c("TAGLN", "MYH11","ACTA2","CNN1", "DES", "PRUNE2")

pericyte_feature <- c("RGS5", "GJC1", "ADCY3",toupper(c("Higd1b", "Cox4i2", "Notch3", "pdgfrb", "mcam", "Cspg4", "Kcnj8","Abcc9", "Vtn ")))# ADCY3, neuronal primary cilia and obesity

Chondrocyte_feature <- c("IBSP", "SP7", "IFITM5", "SATB2", "INSC")


Endothelial_feature <- c("RAMP2","VWF", "ST6GALNAC3", "NOTCH1", "FLT1")


Erythroid_marker <- c("HBA2", "HBA1","HBB", "ANK1")
Macrophage_feature <- c("C1QA","C1QC","LYZ", "CD86", "CD68","CD14", "FCGR3A","FCGR2B")
focused_B_cell_feature <- c("CD19","MS4A1","BANK1", "TNFRSF13C")
plasma_feature <- c("CD38","MZB1","PDK1","IGHG1")
T_cell_feature <- c("CD247", "CD3E", "CD3D", "CD2", "CD7") 


all_marker <- list(Hepatocytes_feature, Neuron_feature,Epithelial_cell, Fibroblast_feature,smooth_muscle, pericyte_feature,Chondrocyte_feature,  Endothelial_feature, 
                   Macrophage_feature, plasma_feature, focused_B_cell_feature, 
                   T_cell_feature, Erythroid_marker)
names(all_marker) <- c("Hepatocyte", "Neuron", "Epithelial", "Fibroblast","Smooth muscle", "Pericyte","Chondrocyte",   "Endothelial", "Macrophage", "Plasma cell","B cell", "T cell", "Erythroblst")



# Yangyi's marker ---------------------------------------------------------
Macrophage_feature <- c("C1QA","C1QC","LYZ","FCER1A","EREG","CD68","CD14","FCGR3A", "DHRS9", "CD80",
                        "CD86", "MS4A1", "CD19")## pro-inflammatory
M1_feature <- c("IL23A", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86", "IL1A", "IL1B", "IL6", "CCL5",
                "IRF5", "IRF1", "CD40", "IDO1", "KYNU", "CCR7", "CXCL3", "PTGS2", "CD14")
## anti-inflammatory
M2_feature <- c("IL4R", "CCL4", "CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB",
                "VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14",
                "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASLG", "TNFSF12", "TNFSF8", "CD276", "VTCN1", "MSR1", "FN1", "IRF4",
                "ARG1", "FLT1")
phagocytosis_feature <- c("MRC1", "CD163","MERTK","C1QB")
angiogenesis_feature <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN",
                          "HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1",
                          "TNFAIP6","TYMP","VAV2","VCAN","VEGFA")
TAM_feature <- c("FCN1", "C1QC", "SPP1")
Monocyte_feature <- c("S100A9","S100A8","LYZ", "KCNMA1", "PLXDC2", "SAT1", "LRMDA", "RBM47", "CD14")  ## CD14

focused_Macrophage_feature <- unique(c("CD86","C1QA","C1QB","C1QC",
                                       "IRF1","KYNU","CD14",
                                       "CTSB","CTSC","CTSD","MSR1","MRC1","CD274","CD163",
                                       "SPP1","EZH2","CD68","MERTK",
                                       "CD44","CXCR4","E2F3","EZH2","ITGAV","SPP1","TYMP","VCAN","VEGFA"))
focused_M1_M2_feature <- c("IRF1","KYNU","CD14",
                           "CTSB","CTSC","CTSD","MSR1","MRC1","CD274","CD163")
focused_phago_angio_feature <- c("MRC1", "CD163","MERTK","C1QB",
                                 "CD44","CXCR4","E2F3","EZH2","ITGAV","SPP1","TYMP","VCAN","VEGFA")

NK_feature <- c("KLRB1","KLRG1","XCL2","XCL1","NCR1","KLRD1","NKG7","FGFBP2","GZMB", "NCAM1")
Fibroblast_feature <- c("DCN","LUM")
Endothelial_feature <- c("RAMP2","VWF")
Epithelial_feature <- c("XBP1","EPCAM","PCA3","KLK3","KRT18","KRT8")
Mast_feature <- c("TPSAB1","CPA3")
pDC_feature <- c("IRF7","IRF4", "CD80", "CD86")
mDC_feature <- c("LYZ","FCER1A","CD1C", "CD80", "CD86")

Treg_feature <- c("CTLA4","IL2RA","FOXP3","CD3D","CD4","CCR7","SELL")
Th17_feature <- c("KLRB1","IL7R","RORC","IL17A","IL17F","CD3D","CD4")
Th1_feature <- c("IL7R","IL2","TNF","CD4")
Naive_Th_feature <- c("IL7R","CCR7","SELL")
CTL_feature <- c("IL7R","IL2","TNF","IFNG","XCL2","XCL1","GZMB","CD8B","CD8A", "CD3G", "CD247","CD3D","CD3E", "FOXP3", "PTPRC")
CD8_effector_feature <- c("GZMB","CD8B","CD8A", "CD3E", "CD8", "FOXP3","CD3G", "CD247", "PTPRC")
# CD4 T: 'IL17A','IL17F', are removed since no cell expresses it
CD4_T_cell_feature <- c('KLRB1', 'IL7R', 'SELL', 'IL2', 'CTLA4', "CD3G", "CD247",'CD3D', 'CD3E', 'CD4', 'RORC', 'TNF', 'FOXP3', 'IL2RA', 'CCR7', 'PTPRC', 'CD8A')
CD8_T_cell_feature <- unique(CTL_feature, CD8_effector_feature)

functional_T_cell_feature <- c("CD247","IL7R", "SELL","FOXP3","IL2RA", "TIGIT", "PDCD1","CTLA4","PRF1","GZMA")
focused_naive_T_cell_feature <- c("CD247","IL7R", "SELL")
focused_eff_reg_exh_T_cell_feature <- c("PRF1","GZMA","FOXP3","IL2RA", "TIGIT", "PDCD1","CTLA4")
main_T_cell_feature <- c('CD3D', 'CD3E', "CD3G", 'CD4', "CD8A","CD8B")
focused_T_cell_feature <- unique(c( main_T_cell_feature, focused_eff_reg_exh_T_cell_feature,functional_T_cell_feature))

Plasma_cell_feature <- c("MZB1","XBP1","CD79B", "FCGR3A", "FCGR2A", "FCGRT", "FCGR1A", "FCGR3B", "FCGR2B", "CD40LG")
# "MS4A1"=CD20, "CD40LG" = IgM, "MME"=CD10, "FCER2"=CD23, "CR2"=CD21
B_cell_feature <- c("MZB1","XBP1","CD79B", "FCGR3A", "FCGR2A", "FCGRT", "FCGR1A", "FCGR3B", "FCGR2B", "CD40LG","FCGR3A", "FCGR2A", "FCGRT", "FCGR1A", "FCGR3B", "FCGR2B", "CD40LG",
                    "MS4A1","CD79B","CD79A","BANK1","GPR183","TCL1A", "CD19", "CD27", "CD80", "CD86","CR2","MME","DNTT","ALKBH3",
                    "DNTT","MME","CD19","MS4A1","CD22","HLA-DRA","CD24","CD40LG","CR2","FCRL2","FCER2","SDC1","CD80","CD84","CD86",
                    "FAS","IL10","CXCL13","IL1A","IL1B","CD1D","IL2RA","BSG","HAVCR1","CD5","CD27","CTLA4","PDCD1","TNFRSF13C","TNFRSF17","TAC1")
B_development_feature <- c("DNTT","MME","CD19","MS4A1","CD22","HLA-DRA","CD24","CD40LG","CR2","FCRL2","FCER2","SDC1","CD80","CD84","CD86",
                           "FAS","IL10","CXCL13","IL1A","IL1B","CD1D","IL2RA","BSG","HAVCR1","CD5","CD27","CTLA4","PDCD1","TNFRSF13C")
focused_B_cell_feature <- c("MZB1","XBP1","CD38","SDC1","BSG",
                            "CD44","MS4A1","CD19","FCGR2B","BANK1","HLA-DRA","CD79B","CD79A",
                            "CD86","TNFRSF13C",
                            "CD22","CD24",
                            "FCRL2","IL2RA","CD27")
focused_plasma_cell_feature <- c("MZB1","XBP1","CD38","SDC1","BSG")
focused_activated_B_feature <- c("CD44","MS4A1","CD19","FCGR2B","BANK1","HLA-DRA","CD79B","CD79A","CD86","TNFRSF13C")
focused_naive_memory_B_feature <- c("CD22","CD24","FCRL2","IL2RA","CD27")

all_feature <- unique(c(Macrophage_feature, Monocyte_feature, CD4_T_cell_feature, CD8_T_cell_feature,
                        NK_feature, Fibroblast_feature, Plasma_cell_feature, B_cell_feature, Endothelial_feature,
                        Epithelial_feature, Mast_feature, pDC_feature, mDC_feature))


# B cell ------------------------------------------------------------------
bcell_marker <- list(plasma = c("FKBP11", "JCHAIN", "CD38", "XBP1", "MZB1"), 
                     immunoglobulin = c("IGLC3", "IGLC2", "IGKC", "IGHG2", "IGHG3", "IGHG1", "IGHA1"), 
                     naiveB = c("SIPA1L1", "AC120193.1", "IGHD", "ADAM28", "IL4R", "COL19A1"),
                     memory.naive = c("LTB", "CD83", "HLA−DQB1", "HLA−DQA1", "BANK1", "MARCH1", "LINC00926", "COTL1", "ARHGAP24", "MS4A1"), 
                     preB = c("NSMCE1", "NEIL1", "BCL7A", "RUBCNL", "HCK", "FAM129C", "TCL1A", "PCDH9", "ACSM3", "CCDC191"),
                     proB = c("RAG1", "CD9", "PSD3", "ERG", "VPREB1", "DNTT", "AKAP12", "LINC01013", "ARPP21", "SLC8A1−AS1"), 
                     Erythroid =  c("HBA2", "HBA1","HBB", "ANK1"))

# macrophage --------------------------------------------------------------
macro_marker <- list(monocyte = c("CD14", "FCGR3A" , "CD68", "C1QC", "C1QA", "HLA-DPA1"),
                     angio_TAM = c("VCAN", "FCN1","MARCO",  "OLR1", "S100A9", "S100A8"), 
                     TRM = c("F13A1", "LYVE1", "CD163L1", "CD163", "MRC1"), 
                     TAM = c("SPP1", "CD83", "CD109", "TM4SF19"))

# endothelial -------------------------------------------------------------
 #tip cell #https://www.nature.com/articles/s41422-022-00615-z

endo_marker <- list(tip_like = c("SLC45A4", "CXCR4", "IGFBP3", "ESM1"),
                    arterial = c("GJA5", "FBLN5", "LTBP4", "SOX17"), 
                    vein = c("ACKR1", "SELP", "VCAM1", "NR2F2"), 
                    lymphatic = c("PKHD1L1", "STON2", "PDGFC", "PROX1"), 
                    liver_sinusoidal = c("ATRNL1", "OIT3", "FGF23", "DNASE1L3")) 

# Epithelial --------------------------------------------------------------
epi_marker <- list(basal = c("KRT5", "KRT15", "VAV3", "TP63", "TG"),
                   luminal = c("LMO7", "CLDN3", "ALDH1A3", "ALDH1A2"),
                   ciliated = c("ADCY2", "NR4A1", "PTPRN2"),
                   aveoliar = c("ITK", "DOCK2", "ABCA3", "SFTPB"),
                   seminal_vesicle = c("SEMG1", "SEMG2", "PATE4", "MUC6", "PGC",
                                       "PAEP", "SERPINA5", "CYP4F8"))

# T cell ------------------------------------------------------------------
tcell_marker <- list(Treg = c("CTLA4" ,"IL2RA", "FOXP3", "RTKN2","IKZF2"), 
                     CD8_effector = c("CD8A", "CD8B", "CCL5", "GZMK"),
                     CD4_helper = c( "IL2", "TNF", "IL4R", "CD4"),
                     Naive = c( "BACH2", "RGS10", "SERINC5", "LEF1", "NELL2", "SELL", "TCF7", "IL7R"),
                     NK_cell = c( "GZMB","NCAM1", "GNLY","NKG7",  "KLRF1", "FCGR3A", "FCER1G"))



# fibroblast --------------------------------------------------------------
fibro_marker <- list(smooth_muscle = c("ACTA2", "MYH11", "KCNJ8", "CNN1"), 
                     chondrocyte = c("IBSP", "SP7"), 
                     skeletal_muscle = c("IGF1", "C3", "ABCA9"), 
                     resting_fibroblast = c("S100A4", "CFD", "DPT"), 
                     inflam = c("FAP", "COL1A1", "POSTN","CTHRC1"), 
                     adipose =c("ABCC9", "RGS5", "GJC1", "ADCY3","COX4I2", "LPL"))

# Cell-type marker gene signatures
# Source: Nat Commun 2023 (41467-2023-36325), Supplementary Data 3
# Parsed from 41467_2023_36325_MOESM6_ESM.xlsx
# NOTE: gene symbols are as-published (GRCh37/pre-HGNC-update in places,
#       e.g. DARC->ACKR1, C19orf10->MYDGF, KIAA0226->RUBCN). Intersect with
#       rownames(seu) before scoring.
 
# ---- Option A: Tumor signature only ----
prostate_tumor_signature <- c("AMACR", "CACNA1D", "PCA3", "ERG", "FABP5", "COL9A2", "GCNT1", "PHGR1")
 
# ---- Option B: full named list of all cell types ----
marker_signatures <- list(
  `Mast cells` = c("TPSAB1", "CPA3", "HPGDS", "RGS13", "IL1RL1", "KIT", "GATA2", "FOSB", "RP11-354E11.2", "VWA5A", "SLC18A2", "HDC", "CAPG"),
  `Endothelial cells` = c("RAMP2", "TM4SF1", "RNASE1", "EGFL7", "RAMP3", "PLVAP", "AQP1", "ECSCR", "FKBP1A", "AC011526.1", "EMP1", "DARC", "VWF", "EMCN"),
  `Pericytes` = c("RGS5", "ACTA2", "MYH11", "MT1M", "FRZB", "MT1A", "NDUFA4L2", "PPP1R14A", "MYLK", "PHLDA1"),
  `Fibroblasts` = c("DCN", "LUM", "PTN", "IGF1", "APOD", "COL1A2", "FBLN1", "MEG3", "CXCL12"),
  `PDC` = c("IRF7", "IRF4", "LILRA4", "PPP1R14B", "SOX4", "TSPAN13", "KIAA0226", "PTCRA", "RAB11FIP1", "CXCR3", "IL3RA"),
  `B cells` = c("MS4A1", "CD79B", "BTG1", "VPREB3", "BANK1", "CD79A", "IGLL5"),
  `Plasma cells` = c("SEC11C", "XBP1", "PRDX4", "SPCS2", "SSR3", "SDF2L1", "C19orf10", "MANF", "TMEM258", "DNAJB9"),
  `Macrophage` = c("C1QA", "C1QC", "C1QB", "GPR34", "CSF1R", "MS4A4A"),
  `Monocytes` = c("S100A9", "FCN1", "S100A8", "CSTA", "EREG", "NEAT1", "TKT", "THBS1", "VCAN", "TSPO", "LYZ"),
  `mDC` = c("CD1C", "PKIB", "INSIG1", "CLEC10A", "C15orf48", "PPA1"),
  `NK` = c("NKG7", "GNLY", "KLRD1", "KLRB1", "FGFBP2", "PRF1"),
  `CTL` = c("CD8A", "CD8B", "GZMH", "GZMA", "PTPRC"),
  `Naive Th` = c("CCR7", "SELL"),
  `Th1` = c("CD4", "IL2", "TNF"),
  `Treg` = c("TIGIT", "CTLA4", "SOD1", "TNFRSF4", "TNFRSF18", "RTKN2", "FOXP3"),
  `Epithelial Club` = c("SCGB1A3", "WFDC2", "LCN2", "MMP7", "KRT4", "TACSTD2", "SCGB3A1"),
  `Epithelial Hillock` = c("KRT13", "S100A16", "S100A14", "KRT19"),
  `Epithelial Basal` = c("TP63", "KRT14", "KRT5"),
  `Epithelial Luminal` = c("KLK4", "KLK3", "KLK2", "ACPP", "AR"),
  `Tumor` = c("AMACR", "CACNA1D", "PCA3", "ERG", "FABP5", "COL9A2", "GCNT1", "PHGR1"),
  `Th17` = c("IL17A", "IL17F", "RORC", "CD4", "CCL20", "CCR6", "RORA"),
  `CD8+ effector` = c("CD8A", "CD8B", "IFNG", "GZMK"),
  `CD56bright NK` = c("KLRC1", "CD44", "COTL1", "XCL1", "XCL2", "TBX21", "EOMES"),
  `CD56dim NK` = c("GZMB", "FGFBP2", "PRF1", "FCGR3A", "TBX21"),
  `NKT` = c("CD44", "CD8A", "CD3D", "NKG7")
)

Epithelial_subtype <- list(
  `Club`    = c("SCGB1A3", "WFDC2", "LCN2", "MMP7", "KRT4", "TACSTD2", "SCGB3A1"),
  `Hillock` = c("KRT13", "S100A16", "S100A14", "KRT19"),
  `Basal`   = c("TP63", "KRT14", "KRT5"),
  `Luminal` = c("KLK4", "KLK3", "KLK2", "ACPP", "AR")
)

Stromal_subtype <- list(
  `Endothelial` = c("RAMP2", "TM4SF1", "RNASE1", "EGFL7", "RAMP3", "PLVAP", "AQP1",
                    "ECSCR", "FKBP1A", "AC011526.1", "EMP1", "DARC", "VWF", "EMCN"),
  `Pericyte`    = c("RGS5", "ACTA2", "MYH11", "MT1M", "FRZB", "MT1A", "NDUFA4L2",
                    "PPP1R14A", "MYLK", "PHLDA1"),
  `Fibroblast`  = c("DCN", "LUM", "PTN", "IGF1", "APOD", "COL1A2", "FBLN1", "MEG3", "CXCL12")
)

Myeloid_subtype <- list(
  `Macrophage` = c("C1QA", "C1QC", "C1QB", "GPR34", "CSF1R", "MS4A4A"),
  `Monocyte`   = c("S100A9", "FCN1", "S100A8", "CSTA", "EREG", "NEAT1", "TKT",
                   "THBS1", "VCAN", "TSPO", "LYZ"),
  `mDC`        = c("CD1C", "PKIB", "INSIG1", "CLEC10A", "C15orf48", "PPA1"),
  `pDC`        = c("IRF7", "IRF4", "LILRA4", "PPP1R14B", "SOX4", "TSPAN13",
                   "KIAA0226", "PTCRA", "RAB11FIP1", "CXCR3", "IL3RA"),
  `Mast`       = c("TPSAB1", "CPA3", "HPGDS", "RGS13", "IL1RL1", "KIT", "GATA2",
                   "FOSB", "RP11-354E11.2", "VWA5A", "SLC18A2", "HDC", "CAPG")
)

Tcell_subtype <- list(
  `CTL`           = c("CD8A", "CD8B", "GZMH", "GZMA", "PTPRC"),
  `CD8_effector`  = c("CD8A", "CD8B", "IFNG", "GZMK"),
  `Naive_Th`      = c("CCR7", "SELL"),
  `Th1`           = c("CD4", "IL2", "TNF"),
  `Th17`          = c("IL17A", "IL17F", "RORC", "CD4", "CCL20", "CCR6", "RORA"),
  `Treg`          = c("TIGIT", "CTLA4", "SOD1", "TNFRSF4", "TNFRSF18", "RTKN2", "FOXP3"),
  `NK`            = c("NKG7", "GNLY", "KLRD1", "KLRB1", "FGFBP2", "PRF1"),
  `CD56bright_NK` = c("KLRC1", "CD44", "COTL1", "XCL1", "XCL2", "TBX21", "EOMES"),
  `CD56dim_NK`    = c("GZMB", "FGFBP2", "PRF1", "FCGR3A", "TBX21"),
  `NKT`           = c("CD44", "CD8A", "CD3D", "NKG7")
)

Bcell_subtype <- list(
  `B`      = c("MS4A1", "CD79B", "BTG1", "VPREB3", "BANK1", "CD79A", "IGLL5"),
  `Plasma` = c("SEC11C", "XBP1", "PRDX4", "SPCS2", "SSR3", "SDF2L1",
               "C19orf10", "MANF", "TMEM258", "DNAJB9")
)
