import argparse
import os
import pysam
yaml = __import__('yaml')

def make_letter_bam(letter, read_seqs, out_dir, region, fasta_path, fg='C', bg='A'):
    chrom, span = region.split(':')
    start_pos, end_pos = [int(_) for _ in span.split('-')]
    bam = f'{out_dir}/{letter}.fg{fg}.bg{bg}.bam'
    ref_fasta = pysam.FastaFile(fasta_path)
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": r, "LN": ref_fasta.get_reference_length(r)} for r in ref_fasta.references]
    }
    ref_id = ref_fasta.references.index(chrom)
    with pysam.AlignmentFile(bam, "wb", header=header) as outf:
        for i, seq in enumerate(read_seqs[::-1]):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = ref_id
            a.reference_start = start_pos
            a.mapping_quality = 60
            a.cigar = [(0, len(seq))]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            outf.write(a)
    os.system(f"samtools index {bam}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--letters_dir', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--region_yaml', required=True)
    parser.add_argument('--fg', default='C')
    parser.add_argument('--bg', default='A')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    stretch_starts = yaml.safe_load(open(args.region_yaml))
    region = stretch_starts[args.bg][70]

    for file in os.listdir(args.letters_dir):
        if file.endswith('.txt'):
            dst_letter = file.split('.')[0]
            with open(os.path.join(args.letters_dir, file)) as f:
                read_seqs = [line.strip() for line in f if line.strip()]
            make_letter_bam(dst_letter, read_seqs, args.out_dir, region, args.fasta, fg=args.fg, bg=args.bg)

if __name__ == '__main__':
    main()

