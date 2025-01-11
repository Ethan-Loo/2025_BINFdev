# Author: S Sevilla
# Date: 1/8/25
# Description: Run script which deploys the nextflow workflow

# Input args
profile=$1 
outDir=$2
projectID=$3
samplesheet=$4

# Print usage if no arguments provided
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <profile> <outDir> <projectID>"
    echo "Example: $0 standard /path/to/output PROJECT001"
    exit 1
fi

# housekeeping
if [[ -z $outDir || -z $profile || -z $projectID ]]; then echo "All variables are required: profile outDir projectID"; exit; fi
if [[ ! -d $outDir ]]; then mkdir -p $outDir; fi

# run the workflow
nextflow run main.nf \
	-profile $profile \
	--outdir $outDir \
 	--projectID $projectID \
	--input $samplesheet \
	-resume

# Check if the workflow ran successfully
if [[ $? -eq 0 ]]; then
	echo "Workflow completed successfully"
else
	echo "Workflow failed"
fi