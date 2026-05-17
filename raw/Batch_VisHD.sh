#!/bin/bash/
#USAGE: this script creates a new job for each CytAssist slide scan in specified path and runs spaceranger-3.1.2 as per /home/jcomben/Batch_scripts/Run_VisHD_Space_Ranger.sh
# sh /home/jcomben/Batch_scripts/Batch_VisHD.sh [output directory] [fastq directory] [reference] [optional: sample IDs]
# eg ${1} = /home/sweng/VisHD
# eg ${2} = /pipeline/Runs/NextSeq/251120_VH01624_389_222CWN5NX/ProjectFolders/Project_Sirui-Weng
# eg ${3} = refdata-gex-GRCh38-2024-A
# eg ${4} = "sample1,sample2" (optional - comma-separated list of sample IDs to process)
#
# **REQUIREMENTS
# SpaceRanger 3.1.2 must be manually installed as it's not set up on the cluster
# This MUST be run from the login node to kick off downstream jobs
# FASTQs must be named identically to the CytAssist identifier
# CytAssist alignment scans must be in ${outdir}/scans/...
echo "--------------usage---------------"
echo "this script batch submits SpaceRanger jobs to the cluster for each CytAssist scan in specified path"
echo 
echo "syntax:"
echo "sh /home/jcomben/Batch_scripts/Batch_VisHD.sh [output directory] [fastq directory] [reference] [optional: sample IDs]"
echo
echo "-----------requirements-----------"
echo "fastqs must be named identically to the CytAssist identifier"
echo "CytAssist scans must be in {outdir}/scans/..."
echo "40x images must be in {outdir}/images/..."
echo "The images must end with *_40x.ome.tif"
echo
echo
echo "available references:"
echo "refdata-gex-GRCh38-2024-A"
echo
outdir=${1}
fastqpath=${2}
ref=${3}
filter_samples=${4}   # optional: comma-separated sample IDs e.g. "sample1,sample2"
part=rhel_short
cores=32
mem=50GB
time=0-10:00

# Parse the optional sample filter into an array
declare -a sample_filter=()
if [[ -n "${filter_samples}" ]]; then
    IFS=',' read -ra sample_filter <<< "${filter_samples}"
    # Trim whitespace from each element
    sample_filter=("${sample_filter[@]// /}")
    echo "----------sample filter-----------"
    echo "Restricting to samples: ${sample_filter[*]}"
fi

echo "-----job resources to request-----"
echo "partition: ${part}"
echo "wall-time: ${time}"
echo "threads: ${cores}"
echo "memory: ${mem}"
echo "---------checking inputs----------"
mkdir -p ${outdir}
if [ -d ${outdir} ];
	then
		echo "READY outdir: ${outdir}"
	else
		echo "NOT READY outdir unavailable"
fi
#specify the location that contains the Sample_${sampleID} directories
if [[ -z $(find ${fastqpath} -maxdepth 1 -type d -name "Sample_*") ]];
	then
		echo "NOT READY directory does not contain "Sample_[ID]" subdirectories"
	else
		echo "READY fastqdir: ${fastqpath}"
fi
#whitelist of available references with both STAR ref and probeset defined
case "${ref}" in
	refdata-gex-GRCh38-2024-A)
    echo "READY reference: ${ref}"
    ;;
  *)
    echo "NOT READY reference unavailable"
    ;;
esac
echo

# Helper function: check if a sampleID is in the filter list (returns 0 if match or no filter set)
in_sample_filter() {
    local id="${1}"
    if [[ ${#sample_filter[@]} -eq 0 ]]; then
        return 0   # no filter — accept all
    fi
    for s in "${sample_filter[@]}"; do
        if [[ "${s}" == "${id}" ]]; then
            return 0   # match found
        fi
    done
    return 1   # not in filter
}

#Cytassist scans must be in ${outdir}/scans/
ready_samples=()
echo "Detected samples:"
for cytScan in $(find "${outdir}/scans" -type f -name "*CAVG*.tif" -exec basename -- {} \;)
	do
		if [[ ${cytScan} =~ H1-[-[:alnum:]]+_[AD][12]_(.*).tif ]];
			then sampleID=${BASH_REMATCH[1]}
		fi

		# Skip samples not in the filter list (if a filter was provided)
		if ! in_sample_filter "${sampleID}"; then
			echo "SKIPPED ${sampleID} (not in sample filter)"
			continue
		fi

		if test -n "$(find ${fastqpath} -mindepth 2 -maxdepth 2 -type f -name "${sampleID}*.fastq.gz" -print -quit)";
			then
				echo "READY ${sampleID}"
				ready_samples+=("${sampleID}|${cytScan}")
			else
				echo "MISMATCH ${sampleID}"
	fi
done

# Warn if any requested filter samples were never detected
if [[ ${#sample_filter[@]} -gt 0 ]]; then
    for s in "${sample_filter[@]}"; do
        found=0
        for entry in "${ready_samples[@]}"; do
            if [[ "${entry%%|*}" == "${s}" ]]; then
                found=1; break
            fi
        done
        if [[ ${found} -eq 0 ]]; then
            echo "WARNING: requested sample '${s}' was not found / not ready"
        fi
    done
fi

echo "----------------------------------"
echo 
read -p "Proceed with above inputs? (y/n)" -n 1 -r
#make job log files go into the output directory rather than /home
cd ${outdir}
# -p ${part}
#submits a new cluster job for each CytAssist scan found in above list.
#typical CytAssist libraries (~100M reads) take 2.5 hours and 16GB to run on 4 cores
if [[ $REPLY =~ ^[Yy]$ ]];
	then
		echo
		for entry in "${ready_samples[@]}"
			do
				sampleID="${entry%%|*}"
				cytScan="${entry##*|}"
				echo "starting job for ${sampleID}"
				sbatch -c ${cores} --time ${time} --partition ${part} --output ${outdir}/slurmout/${sampleID}.%j.out --mem ${mem} --mail-type=ALL --mail-user=sirui.weng@petermac.org --wrap="sh /home/sweng/VisHD/raw/Run_VisHD_Space_Ranger.sh ${outdir} ${fastqpath} ${cytScan} ${ref}"
				sleep 1
		done
		echo
	else
		echo
		echo "aborted"
		echo
fi
