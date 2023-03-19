import argparse
import vcf
from spliceai.utils import Annotator, get_delta_scores
from collections import namedtuple
Record = namedtuple('Record', ['chrom', 'pos', 'ref', 'alts'])

def process_vcf(in_vcf:str, ann: Annotator, distance: int, mask: int, out_vcf:str):
    reader = vcf.Reader(filename=in_vcf)
    writer = vcf.Writer(open(out_vcf, 'w'), reader)
    for row in reader:
        if row.is_indel and row.INFO.get('SpliceAI'):
            record = Record(chrom=row.CHROM, pos=row.POS, ref=row.REF, alts=row.ALT)
            scores =  get_delta_scores(record=record, ann=ann, dist_var=distance, mask=mask)
            if scores:
                row.INFO['SpliceAI'] = ','.join(scores)
        writer.write_record(row)

def run_main_cmd(args):
    ann = Annotator(args.reference, args.builder.replace('hg', 'grch'))
    process_vcf(in_vcf=args.input, ann=ann, distance=args.distance, mask=args.mask, out_vcf=args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SpliceAI Wrapper')
    parser.add_argument('--input', '-i', required=True, help='path to the input VCF file')
    parser.add_argument('--output', '-o', required=True, help='path to the output VCF file')
    parser.add_argument('--reference', '-r', required=True, help='path to the reference genome fasta file')
    parser.add_argument('--builder', '-b', choices=['grch37', 'hg19', 'grch38', 'hg38'], default='grch38', help='genomer version, defaults to grch38')
    parser.add_argument('--distance', '-d', default=50, type=int, help='maximum distance between the variant and gained/lost splice site, defaults to 50')
    parser.add_argument('--mask', '-m', default=0, type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and unannotated acceptor/donor loss, defaults to 0')
    parser.set_defaults(func=run_main_cmd)
    args = parser.parse_args()
    args.func(args)