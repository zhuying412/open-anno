#! /usr/bin/env python3
import argparse
import re
from vcf import Reader as VcfReader



def format_variant(chrom: str, pos: int, ref: str, alt: str) -> str:
    start = pos
    if len(ref) > 1 or len(alt) > 1 and ref != alt:
        if ref.startswith(alt) or ref.endswith(alt):
            if ref.startswith(alt):
                start = start + len(alt)
            ref = ref.replace(alt, '', 1)
            alt = ''
        elif alt.startswith(ref) or alt.endswith(ref):
            start = start + len(ref) - 1 if alt.startswith(ref) else start - len(alt) + len(ref)
            alt = alt.replace(ref, '', 1)
            ref = ''
        else:
            ref_rev, alt_rev, substr, stop, index = ref[::-1], alt[::-1], '', False, 0
            while index < len(ref) and index < len(alt):
                if ref_rev[index] != alt_rev[index]:
                    stop = True
                if ref_rev[index] == alt_rev[index] and not stop:
                    substr = ref_rev[index] + substr
                index += 1
            ref = re.sub(r'%s$' % substr, '', ref)
            alt = re.sub(r'%s$' % substr, '', alt)
            substr, stop, index = '', False, 0
            while index < len(ref) and index < len(alt):
                if ref[index] != alt[index]:
                    stop = True
                if ref[index] == alt[index] and not stop:
                    substr += ref[index]
                index += 1
            ref = re.sub(r'^%s' % substr, '', ref)
            alt = re.sub(r'^%s' % substr, '', alt)
            start += len(substr) - 1 if len(substr) and not ref else len(substr)
    end = start + len(ref) - 1 if ref else start
    return f'{chrom}\t{start}\t{end}\t{ref if ref else "-"}\t{alt if alt else "-"}'


def main(args):
    reader = VcfReader(filename=args.vcf)
    with open(args.anno_input, 'w') as writer:
        for row in reader:
            chrom = str(row.CHROM)
            if (chrom.startswith('chr') and len(chrom) > 5) or (not chrom.startswith('chr') and len(chrom) > 2):
                continue
            depth = row.INFO['VD']
            filter = ','.join(row.FILTER) or 'PASS'
            qual = row.QUAL
            gt = '/'.join(row.genotype(reader.samples[0]).gt_alleles)
            for i, alt in enumerate(row.ALT):
                alt = str(alt)
                if re.findall(r'[^ATGC]', alt):
                    if alt != '<DEL>':
                        continue
                    alt = '-'
                raw = f'{row.CHROM}:{row.POS}:{row.REF}:{alt}'
                vaf = row.INFO['AF'][i] if isinstance(row.INFO['AF'], list) else row.INFO['AF']
                variant = format_variant(chrom, int(row.POS), str(row.REF), alt)
                writer.write(f'{variant}\tDEPTH={depth};VAF={vaf};FILTER={filter};QUAL={qual};GT={gt};RAW={raw}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert VCF to AnnoInput')
    parser.add_argument('--vcf', '-i', required=True, help='input file, VCF file')
    parser.add_argument('--anno_input', '-o', required=True, help='output file, anno input file')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)

