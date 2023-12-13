"""
Mask initial bases from alignment FASTA
"""
import argparse
from augur.io import open_file, read_sequences, write_sequences
import Bio
import Bio.SeqIO
from Bio.Seq import Seq


def mask_terminal_gaps(seq):
    L = len(seq)
    seq_trimmed = seq.lstrip("-")
    left_gaps = L - len(seq_trimmed)
    seq_trimmed = seq_trimmed.rstrip("-")
    right_gaps = L - len(seq_trimmed) - left_gaps
    return "N" * left_gaps + seq_trimmed + "N" * right_gaps


def mask_all_gaps(seq):
    return seq.replace("-", "N")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of alignment")
    parser.add_argument(
        "--mask-terminal-gaps",
        action="store_true",
        help="fill all terminal gaps with N as they likely represent missing data",
    )
    parser.add_argument(
        "--mask-all-gaps", action="store_true", help="fill all gaps with N"
    )
    parser.add_argument(
        "--mask-from-beginning", type=int, help="number of bases to mask from start"
    )
    parser.add_argument(
        "--mask-from-end", type=int, help="number of bases to mask from end"
    )
    parser.add_argument(
        "--mask-sites", nargs="+", type=int, help="list of sites to mask"
    )
    parser.add_argument("--mask-site-file", help="list of sites to mask")
    parser.add_argument(
        "--output", required=True, help="FASTA file of output alignment"
    )
    args = parser.parse_args()

    begin_length = 0
    if args.mask_from_beginning:
        begin_length = args.mask_from_beginning
    end_length = 0
    if args.mask_from_end:
        end_length = args.mask_from_end

    mask_sites = []
    if args.mask_site_file:
        with open_file(args.mask_site_file, "r") as f:
            mask_sites = [int(x) for x in f.read().split()]

    with open_file(args.output, "w") as outfile:
        for record in read_sequences(args.alignment):
            seq = str(record.seq).upper()
            if args.mask_terminal_gaps:
                seq = mask_terminal_gaps(seq)
            if args.mask_all_gaps:
                seq = mask_all_gaps(seq)

            start = "N" * begin_length
            middle = (
                seq[begin_length:-end_length] if end_length > 0 else seq[begin_length:]
            )
            end = "N" * end_length
            seq_list = list(start + middle + end)
            if mask_sites:
                for site in mask_sites:
                    if seq_list[site - 1] != "-":
                        seq_list[site - 1] = "N"
            record.seq = Seq("".join(seq_list))
            write_sequences(record, outfile)
