#! /usr/bin/env python3
import argparse
import vcf

def calc_baf(args):
    output_file = open(args.outputfile, 'w')
    output_file.write("chromosome\tstart\tend\tlocus\tbaf\n")
    vcf_file = args.inputfile
    with open(vcf_file, 'r') as vcf_input_file: 
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        sampleid = vcf_reader.samples[0]
        for record in vcf_reader:
            chrom = record.CHROM
            end = record.POS
            start = end - 1
            locus = "{}:{}-{}".format(chrom, start, end)
            dp = record.genotype(sampleid)['DP']
            variant_call = record.genotype(sampleid).is_variant

            if dp >= args.mindepth: 
                if variant_call == True:
                    ad = record.genotype(sampleid)['AD']
                    baf = (ad[1] / dp) * 100
                elif variant_call == False:
                    try:
                        ad = [record.genotype(sampleid)['AD']]
                        baf = ((dp - ad[0]) / (dp)) * 100
                    except:  # Bypass to remove positions without AD
                        continue
                else: # is_variant status can be None
                    continue

            else:
                continue

            if len(ad) <= 2:  # skip multiallelic sites
                output_file.write("{chrom}\t{start}\t{end}\t{locus}\t{baf:.2f}\n".format(
                    chrom=chrom,
                    start=start,
                    end=end,
                    locus=locus,
                    baf=baf
                    ))

    output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help='input VCF file (genotyped VCF)')
    parser.add_argument('outputfile', help='output filename')
    parser.add_argument('--mindepth', default= 15, type=int, help='Threshold for minimum depth (DP) of SNV (default = 15)')
    args = parser.parse_args()
    calc_baf(args)