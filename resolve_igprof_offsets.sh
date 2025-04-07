#!/bin/bash

# Usage: ./resolve_igprof_offsets.sh <igprof_report> > resolved_report.txt

# Cache for binary base addresses
declare -A base_addresses

# Resolve a binary+offset into symbol name and location
resolve_address() {
    local bin="$1"
    local offset_dec="$2"

    # Get base address if not cached
    if [[ -z "${base_addresses[$bin]}" ]]; then
        if [[ ! -f "$bin" ]]; then
#            echo "!! Missing binary: $bin" >&2
            base_addresses[$bin]=0
            return 1
        fi
        base_addr_hex=$(readelf -Wl "$bin" 2>/dev/null | awk '/LOAD/ {print $3; exit}')
        base_addr=$(($base_addr_hex))
        base_addresses[$bin]=$base_addr
    else
        base_addr=${base_addresses[$bin]}
    fi

    # Compute full virtual address
    full_addr=$((base_addr + offset_dec))
    full_addr_hex=$(printf "0x%x" "$full_addr")

    # Resolve via addr2line
    result=$(addr2line -e "$bin" -f -C "$full_addr_hex" | head -n 1 2>/dev/null)
    echo "$result"
}

# Process the report line-by-line
while IFS= read -r line; do
    if [[ "$line" =~ \@\{([^+]+)\+([0-9]+)\} ]]; then
        bin="${BASH_REMATCH[1]}"
        offset="${BASH_REMATCH[2]}"

        resolved=$(resolve_address "$bin" "$offset")
        if [ "$?" -eq 0 ]; then
            # Replace the matched @{...} with the resolved info
            line="${line//@\{${bin}+${offset}\}/$resolved}"
        fi
    fi
    echo "$line"
done < "$1"
