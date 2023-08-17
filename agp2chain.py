#!/usr/bin/env python3
import argparse


def agp_to_chain(agp_file, chain_file):
    with open(agp_file, 'r') as agp, open(chain_file, 'w') as chain:
        current_chain_id = 1
        for line in agp:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            object_id, object_beg, object_end, part_number, component_type = fields[:5]

            # If it's a component line
            if component_type not in ['N', 'U']:
                component_id, component_beg, component_end, orientation = fields[5:9]
                block_size = int(component_end) - int(component_beg) + 1
                chain.write(
                    f"chain {current_chain_id} {object_id} {object_end} + {object_beg} {object_end} {component_id} {component_end} {orientation} {component_beg} {component_end}\n"
                    )
                chain.write(f"{block_size}\n\n")
                current_chain_id += 1
            # If it's a gap line
            else:
                gap_length = int(fields[5])
                chain.write(
                    f"chain {current_chain_id} {object_id} {object_end} + {object_beg} {object_end} gap {gap_length} + 0 {gap_length}\n"
                    )
                chain.write(f"{gap_length}\n\n")
                current_chain_id += 1


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="Convert AGP to CHAIN")
    parser.add_argument("agp_file", help="AGP file")
    parser.add_argument("chain_file", help="CHAIN file")

    args = parser.parse_args()
    agp_to_chain(args.agp_file, args.chain_file)

if __name__ == "__main__":
    main()
