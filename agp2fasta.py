#!/usr/bin/env python
"""
create assembly from agp file and input fasta file
it hande only W and N component type

"""

import argparse


def read_fasta(fasta_file):
    """
    read fasta file and return dictionary with contig name as key and sequence as value
    """
    # first arch sequence is list of lines
    fasta_dict = {}
    with open(fasta_file, "r") as f:
        for line in f:
            if line[0] == ">":
                contig_name = line.strip()[1:]
                fasta_dict[contig_name] = []
                continue
            else:
                fasta_dict[contig_name].append(line.strip())
    # concatenate list to str
    fasta_dict = {k: "".join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def create_scaffold(agp_list, input_fasta, scaffold_size):
    """
    create dictionary of bytearrays with scaffolds
    """
    scaffold_dict = {}
    for scaffold in scaffold_size:
        scaffold_dict[scaffold] = []
    for agp in agp_list:
        if agp.component_type == "W":
            # check in replace bytearrays is same length:
            if agp.end - agp.start == agp.component_end - agp.component_start:
                component_seq = input_fasta[agp.component_id][
                                agp.component_start:agp.component_end]
                if agp.orientation == "-":
                    component_seq = reverse_complement(component_seq)

                scaffold_dict[agp.seqid].append(component_seq)
            else:
                raise ValueError(
                    "component {} and scaffold part {} are not same "
                    "length".format(agp.component_id, agp.seqid)
                    )
        elif agp.component_type == "N":
            # fill with N character
            scaffold_dict[agp.seqid].append("N" * agp.gap_length)
    # concatenate to strings
    scaffold_dict = {k: "".join(v) for k, v in scaffold_dict.items()}
    # check if scaffolds are correct size:
    for scaffold in scaffold_dict:
        if len(scaffold_dict[scaffold]) != scaffold_size[scaffold]:
            print("incorrect size of scaffold {}".format(scaffold))
            print(scaffold, len(scaffold_dict[scaffold]), scaffold_size[scaffold])
            print(
                "------------------------------"
                )  # raise ValueError("scaffold  {} size is not correct".format(scaffold))
        else:
            print("correct size of scaffold {}".format(scaffold))
    return scaffold_dict


class AGP:
    """
    class to store agp file
    """

    def __init__(
            self, seqid, start, end, partno, component_type, component_id,
            component_start, component_end, orientation
            ):
        self.seqid = seqid
        self.start = int(start) - 1
        self.end = int(end)
        self.partno = int(partno)
        self.component_type = component_type
        # agp is one based, convert to zero based
        if component_type == "W":
            self.component_id = component_id
            self.component_start = int(component_start) - 1
            self.component_end = int(component_end)
            self.orientation = orientation
            self.gap_length = None
            self.gap_type = None
            self.linkage = None
            self.linkage_evidence = None
        elif component_type == "N":
            self.component_id = None
            self.component_start = None
            self.component_end = None
            self.orientation = None
            self.gap_length = int(component_id)
            self.gap_type = component_start
            self.linkage = component_end
            self.linkage_evidence = orientation


def read_agp(agp_file):
    """
    read agp file using csv module
    column of agp file are:
        seqid
        start
        end
        partno
        component_type
        component_id
        component_start
        component_end
        orientation
    """
    ags_list = []
    with open(agp_file, "r") as f:
        for line in f:
            items = line.strip().split()
            print(items)
            ags_list.append(AGP(*items))
    return ags_list


def get_scaffold_size(ags_list):
    """

    :param ags_list:
    :return:
    dictionary with scaffold name as key and size as value
    """
    scaffold_size = {}
    for ags in ags_list:
        if ags.seqid not in scaffold_size:
            scaffold_size[ags.seqid] = ags.end
        else:
            if ags.end > scaffold_size[ags.seqid]:
                scaffold_size[ags.seqid] = ags.end
    return scaffold_size


def reverse_complement(seq):
    """
    reverse complement sequence
    """
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join([comp[base] for base in seq[::-1]])


def validate_scaffold(scaffolds_dict, agp_list, input_fasta):
    """
    validate scaffold size
    """
    for agp in agp_list:
        if agp.component_type == "W":
            scaffold_part = scaffolds_dict[agp.seqid][agp.start:agp.end]
            component_seq = input_fasta[agp.component_id][
                            agp.component_start:agp.component_end]
            if agp.orientation == "-":
                component_seq = reverse_complement(component_seq)
            if scaffold_part != component_seq:
                raise ValueError(
                    "scaffold part {} is not equal to component {}".format(
                        agp.seqid, agp.component_id
                        )
                    )
        if agp.component_type == "N":
            scaffold_part = scaffolds_dict[agp.seqid][agp.start:agp.end]
            if scaffold_part != "N" * agp.gap_length:
                raise ValueError("scaffold part {} is not equal to N".format(agp.seqid))
    return True


def main():
    """
    main function

    """
    # get arguments from command line
    parser = argparse.ArgumentParser(
        description="""This script is used to create assembly from agp file and input 
        fasta file.""", formatter_class=argparse.RawTextHelpFormatter, )
    parser.add_argument(
        "-a", "--agp", default=None, required=True, help="agp file", type=str,
        action='store'
        )
    parser.add_argument(
        "-f", "--fasta", default=None, required=True, help="fasta file", type=str,
        action='store'
        )
    parser.add_argument(
        "-o", "--output", default=None, required=True, help="output file name", type=str,
        action='store'
        )
    args = parser.parse_args()

    # read fasta file
    input_fasta = read_fasta(args.fasta)
    agp_list = read_agp(args.agp)
    # create assembly
    scaffold_size = get_scaffold_size(agp_list)
    print(scaffold_size)
    scaffold_dict = create_scaffold(agp_list, input_fasta, scaffold_size)
    # write assembly to file

    validate_scaffold(scaffold_dict, agp_list, input_fasta)

    with open(args.output, "wt") as f:
        for scaffold in scaffold_dict:
            # break scaffold into 80 character lines
            f.write(">{}\n".format(scaffold))
            for i in range(0, len(scaffold_dict[scaffold]), 80):
                f.write("{}\n".format(scaffold_dict[scaffold][i:i + 80]))


if __name__ == '__main__':
    main()
