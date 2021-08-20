#! /usr/bin/env python3
import sys
import argparse
import vcf

def calc_baf(args):
    if args.outputfile:
        output_file = open(args.outputfile, 'w')

    with open(args.inputfile, 'r') as vcf_input_file:
        vcf_reader = vcf.Reader(vcf_input_file)
        if 'fileformat' not in vcf_reader.metadata: # Check if true VCF file
            sys.exit("Input file {} is not a correct VCF file. "\
                "Field \"fileformat\" was not detected in VCF".format(args.inputfile)
                )
        if len(vcf_reader.samples) > 1: # Check is VCF is not a multisample VCF
            sys.exit("Single sample VCF support only. Input file {} is a multisample VCF".format(args.inputfile))

        sampleid = vcf_reader.samples[0]

        header = ("#track type=igv name=BAF_{} color=0,100,0 altColor==0,100,0 "\
            "graphType=points windowingFunction=none maxHeightPixels=100 viewLimits=0,100"\
            "\nchromosome\tstart\tend\tlocus\tbaf".format(sampleid)
            )
        if args.outputfile:
            output_file.write("{}\n".format(header))
        else:
            print(header)
        
        for record in vcf_reader:  # Parse all variant positions, and determine BAF if possible.
            chrom = record.CHROM
            end = record.POS
            start = end - 1
            locus = "{}:{}-{}".format(chrom, start, end)

            if 'DP' in record.genotype(sampleid).data._asdict():
                dp = record.genotype(sampleid)['DP']
            else: # Skip variant position as BAF could not be determined because DP is missing
                continue

            variant_call = record.genotype(sampleid).is_variant  # Returns True if a variant-call, and thus not a reference-call

            if dp >= args.mindepth and variant_call == True: 
                """ Calculate ad and baf for variant-calls if QC is correct """
                ad = record.genotype(sampleid)['AD']
                baf = (ad[1] / dp) * 100

            elif dp >= args.mindepth and variant_call == False and 'AD' in record.genotype(sampleid).data._asdict():
                """ Calculate ad and baf for reference-calls if QC is correct """
                ad = [record.genotype(sampleid)['AD']]
                baf = ((dp - ad[0]) / (dp)) * 100

            else:
                continue  

            if len(ad) <= 2:  # skip multiallelic sites
                variant_line = "{chrom}\t{start}\t{end}\t{locus}\t{baf:.2f}".format(
                    chrom=chrom,
                    start=start,
                    end=end,
                    locus=locus,
                    baf=baf
                    )

                if args.outputfile:
                    output_file.write("{}\n".format(variant_line))
                else:
                    print(variant_line)

    if args.outputfile:
        output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help='input VCF file (genotyped VCF). VCF must be uncompressed')
    parser.add_argument('-o', '--outputfile', help='output filename. Without argument, output will be printed in stdout')
    parser.add_argument('--mindepth', default= 15, type=int, help='Threshold for minimum depth (DP) of SNV (default = 15)')
    args = parser.parse_args()
    calc_baf(args)
