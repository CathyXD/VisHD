#!/bin/bash/

#USAGE: this script runs spaceranger count on 3.1.2 for HD assays
# sh /home/jcomben/Batch_scripts/run_VisHD.sh [output directory] [fastq directory] [CytAssist alignment scan filename] [reference]
# eg ${1} = /home/sweng/VisHD/raw
# eg ${2} = /pipeline/Runs/NextSeq/251120_VH01624_389_222CWN5NX/ProjectFolders/Project_Sirui-Weng
# eg ${3} = CAVG10520_2025-11-13_12-54-16_2025-11-13_12-25-15_H1-CWVHCK2_A1_LUT-245-20.tif
# eg ${4} = refdata-gex-mm10-2020-A OR refdata-gex-GRCh38-2024-A
#
# **REQUIREMENTS:
# Cytassist alignment scans must be in [output directory]/scans/...
# FASTQs must be named identically to the CytAssist identifier
#
# typical CytAssist libraries (~100M reads) take 2.5 hours and 16GB to run on 4 cores
# eg. sinteractive -c 4 -p prod_med --time 0-04:00 --mem 24G

#specify inputs
outdir=${1}
fastqpath=${2}
cytScan=${3}
ref=${4}

export PATH=/home/sweng/SpaceRanger/spaceranger-4.1.0:$PATH

#references in static path
refsource="/home/sweng/SpaceRanger"
probesource="/home/sweng/SpaceRanger"

#automatically assign probeset based on reference
if [[ ${ref} = "refdata-gex-mm10-2020-A" ]]; then
	probeset="Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv"
elif [[ ${ref} = "refdata-gex-GRCh38-2024-A" ]]; then
	probeset="Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv"
fi

refpath=${refsource}/${ref}
probepath=${probesource}/${probeset}

#extract sample info from CytAssist alignment scan
#eg. CAVG10520_2024-12-13_12-21-46_2024-12-13_11-51-56_H1-84N3HPG_D1_MM03E1
if [[ ${cytScan} =~ (H1-[-[:alnum:]_]+) ]]; then slideref=${BASH_REMATCH[1]}
fi

if [[ ${slideref} =~ (H1-[-[:alnum:]]+) ]]; then slideID=${BASH_REMATCH[1]}
fi

if [[ ${slideref} =~ H1-[-[:alnum:]]+_([AD][12]) ]]; then wellID=${BASH_REMATCH[1]}
fi

if [[ ${slideref} =~ H1-[-[:alnum:]]+_[AD][12]_(.*) ]]; then sampleID=${BASH_REMATCH[1]}
fi

echo "sample name: ${sampleID}"
echo "reference: ${ref}"
echo "probeset: ${probeset}"
echo "fastqs: ${fastqpath}/Sample_${sampleID}"

#set up working directories
cd ${outdir}

#run space ranger count
spaceranger count --id=${sampleID} \
--create-bam=false \
--transcriptome=${refpath} \
--probe-set=${probepath} \
--fastqs=${fastqpath}/Sample_${sampleID} \
--sample=${sampleID} \
--cytaimage=${outdir}/scans/${cytScan} \
--image=${outdir}/images/${sampleID}_40x.ome.tif \
--slide=${slideID} \
--area=${wellID}
