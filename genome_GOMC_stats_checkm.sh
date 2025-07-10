#!/bin/bash

# Default values
input_dir="../Bins"
output_file="Genome_Statistics_$(date +%Y%m%d).csv"
file_type="auto"
checkm_file=""

# Usage information
usage() {
    echo "Usage: $0 [-i input_directory] [-o output_file] [-t file_type] [-c checkm_file] [-h]"
    echo "Options:"
    echo "  -i  input directory (default: $input_dir)"
    echo "  -o  output file path (default: $output_file)"
    echo "  -t  file organization type: 'flat' or 'nested' (default: auto-detect)"
    echo "  -c  CheckM results file (required for quality assessment)"
    echo "  -h  display this help message"
    exit 0
}

# Parse command-line arguments
while getopts ":i:o:t:c:h" opt; do
    case $opt in
        i) input_dir="$OPTARG" ;;
        o) output_file="$OPTARG" ;;
        t) file_type="$OPTARG" ;;
        c) checkm_file="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

# Validate file type
if [[ "$file_type" != "auto" && "$file_type" != "flat" && "$file_type" != "nested" ]]; then
    echo "Invalid file type: $file_type. Must be 'auto', 'flat' or 'nested'" >&2
    exit 1
fi

# Validate CheckM file
if [[ -z "$checkm_file" ]]; then
    echo "Error: CheckM results file is required for quality assessment. Use -c option." >&2
    exit 1
elif [[ ! -f "$checkm_file" ]]; then
    echo "CheckM file not found: $checkm_file" >&2
    exit 1
fi

# Create output directory if it doesn't exist
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"

# ==== Load CheckM data ====
declare -A checkm_data
echo "Loading CheckM quality data from: $checkm_file"
while IFS=$'\t' read -r bin_id _ _ _ completeness contamination _; do
    # Skip header lines
    [[ "$bin_id" == "Bin Id" || "$completeness" == "Completeness" ]] && continue
    
    # Clean numerical values
    completeness="${completeness//[^0-9.]/}"
    contamination="${contamination//[^0-9.]/}"
    
    # Store in associative array
    checkm_data["$bin_id"]="$completeness,$contamination"
done < "$checkm_file"
echo "Loaded quality data for ${#checkm_data[@]} genomes"

# ==== Function to calculate quality classification ====
classify_quality() {
    local completeness="$1"
    local contamination="$2"
    
    # Calculate Quality Score (QS)
    local qs=$(echo "$completeness - 5 * $contamination" | bc -l)
    
    # Classify according to MIMAG standards
    if (( $(echo "$completeness >= 90 && $contamination <= 5" | bc -l) )); then
        echo "near-complete,$qs"
    elif (( $(echo "$completeness >= 70 && $contamination <= 10" | bc -l) )); then
        echo "high-quality,$qs"
    elif (( $(echo "$completeness >= 50 && $contamination <= 10" | bc -l) )); then
        echo "medium-quality,$qs"
    else
        echo "low-quality,$qs"
    fi
}

# ==== Statistics calculation function ====
calc_stats() {
    awk '
    BEGIN { total=0; max=0; min=1e9; n=0; seqlen=0 }
    /^>/ { 
        if (seqlen) { 
            total += seqlen
            len[n++] = seqlen 
            if(seqlen>max) max=seqlen
            if(seqlen<min) min=seqlen
        }
        seqlen=0
        next
    }
    { seqlen += length($0) }
    END {
        if (seqlen) {
            total += seqlen
            len[n++] = seqlen
            if(seqlen>max) max=seqlen
            if(seqlen<min) min=seqlen
        }
        if(n==0) { print "0,0,0,0,0,0"; exit }
        asort(len, sorted)
        for(i=n; i>=1; i--) rev[n-i+1] = sorted[i]
        sum=0; half=total/2
        for(i=1; i<=n; i++) {
            sum += rev[i]
            if(sum >= half) {
                N50=rev[i]
                L50=i
                break
            }
        }
        print total","n","max","min","N50","L50
    }'
}

# ==== Generate bin identifier ====
generate_bin_id() {
    local file_path="$1"
    local use_type="$2"
    
    if [[ "$use_type" == "nested" ]]; then
        local dir_name=$(basename "$(dirname "$file_path")")
        local base_name=$(basename "$file_path")
        base_name="${base_name%.*}"  # Remove one extension
        base_name="${base_name%.*}"  # Remove second extension for .fa.gz
        echo "${dir_name}_${base_name}"
    else
        local base_name=$(basename "$file_path")
        base_name="${base_name%.*}"  # Remove one extension
        base_name="${base_name%.*}"  # Remove second extension for .fa.gz
        echo "$base_name"
    fi
}

# ==== Directory processing function ====
process_dir() {
    local target_dir="$1"
    [ ! -d "$target_dir" ] && echo "Directory not found: $target_dir" >&2 && return

    # Auto-detect file organization
    local detected_type="flat"
    if find "$target_dir" -maxdepth 1 -type d | grep -q .; then
        detected_type="nested"
    fi
    
    # Use user-specified type or auto-detected
    local use_type="$file_type"
    if [[ "$use_type" == "auto" ]]; then
        use_type="$detected_type"
    fi

    # Enable recursive globbing
    shopt -s globstar
    for file in "$target_dir"/**/*.{fa,fa.gz,fna,fna.gz}; do
        [ -f "$file" ] || continue
        
        # Generate filename based on organization type
        if [[ "$use_type" == "nested" ]]; then
            dir_name=$(basename "$(dirname "$file")")
            base_name=$(basename "$file")
            new_name="${dir_name}_${base_name}"
        else
            new_name=$(basename "$file")
        fi
        
        # Generate bin identifier for CheckM matching
        bin_id=$(generate_bin_id "$file" "$use_type")
        
        # Process file
        if [[ "$file" == *.gz ]]; then
            md5_val=$(md5sum "$file" | awk '{print $1}')
            stats=$(gzip -dc "$file" | calc_stats)
        else
            md5_val=$(md5sum "$file" | awk '{print $1}')
            stats=$(calc_stats < "$file")
        fi
        
        # Initialize quality variables
        local completeness="NA"
        local contamination="NA"
        local qs="NA"
        local quality_class="NA"
        
        # If CheckM data is available for this bin
        if [[ -n "${checkm_data[$bin_id]}" ]]; then
            IFS=',' read -r completeness contamination <<< "${checkm_data[$bin_id]}"
            IFS=',' read -r quality_class qs <<< "$(classify_quality "$completeness" "$contamination")"
        fi
        
        # Output with quality info
        echo "$new_name,$md5_val,$stats,$completeness,$contamination,$qs,$quality_class" >> "$output_file"
    done
    shopt -u globstar
}

# Main execution
echo "fasta_file_name,fasta_file_md5,total_size(bp),sequences,largest_seq(bp),smallest_seq(bp),N50(bp),L50,completeness(%),contamination(%),QS,quality_class" > "$output_file"
process_dir "$input_dir"

echo "Processing complete! Results saved to: $output_file"
