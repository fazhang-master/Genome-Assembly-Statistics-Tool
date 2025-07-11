#!/bin/bash

# Default values
input_dir="../Bins"
output_file="Genome_Statistics_$(date +%Y%m%d).csv"
file_type="auto"

# Usage information
usage() {
    echo "Usage: $0 [-i input_directory] [-o output_file] [-t file_type] [-h]"
    echo "Options:"
    echo "  -i  input directory (default: $input_dir)"
    echo "  -o  output file path (default: $output_file)"
    echo "  -t  file organization type: 'flat' or 'nested' (default: auto-detect)"
    echo "  -h  display this help message"
    exit 0
}

# Parse command-line arguments
while getopts ":i:o:t:h" opt; do
    case $opt in
        i) input_dir="$OPTARG" ;;
        o) output_file="$OPTARG" ;;
        t) file_type="$OPTARG" ;;
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

# Create output directory if it doesn't exist
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"

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
        # Protection against empty files
        if(n==0) { print "0,0,0,0,0,0"; exit }
        # N50/L50 calculation
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
    for file in "$target_dir"/**/*.{fa,fa.gz}; do
        [ -f "$file" ] || continue  # Skip non-files
        
        # Generate filename based on organization type
        if [[ "$use_type" == "nested" ]]; then
            dir_name=$(basename "$(dirname "$file")")
            base_name=$(basename "$file")
            new_name="${dir_name}_${base_name}"  # Format: dirname_original.fa
        else
            new_name=$(basename "$file")  # Format: original.fa
        fi
        
        # Process file
        if [[ "$file" == *.gz ]]; then
            md5_val=$(md5sum "$file" | awk '{print $1}')
            stats=$(gzip -dc "$file" | calc_stats)
        else
            md5_val=$(md5sum "$file" | awk '{print $1}')
            stats=$(calc_stats < "$file")
        fi
        
        # Output with new filename
        echo "$new_name,$md5_val,$stats" >> "$output_file"
    done
    shopt -u globstar  # Disable recursive globbing
}

# Main execution
echo "fasta_file_name,fasta_file_md5,total_size(bp),sequences,largest_seq(bp),smallest_seq(bp),N50(bp),L50" > "$output_file"
process_dir "$input_dir"

echo "Processing complete! Results saved to: $output_file"