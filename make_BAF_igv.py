#! /usr/bin/env python3
import sys
import argparse
import vcf

def calc_baf(args):
    output_file = open(args.outputfile, 'w')
    vcf_file = args.inputfile
    with open(vcf_file, 'r') as vcf_input_file: 
        vcf_reader = vcf.Reader(vcf_input_file)
        if len(vcf_reader.samples) > 1:
            sys.exit("Single sample VCF support only. Input file {} is a multisample VCF".format(args.inputfile))
        sampleid = vcf_reader.samples[0]
        output_file.write("#track type=igv name=BAF_{} color=0,100,0 altColor==0,100,0 graphType=points windowingFunction=none maxHeightPixels=100 viewLimits=0,100\n".format(sampleid))
        output_file.write("chromosome\tstart\tend\tlocus\tbaf\n")
        for record in vcf_reader:  # Parse all variant positions, and determine BAF if possible.
            chrom = record.CHROM
            end = record.POS
            start = end - 1
            locus = "{}:{}-{}".format(chrom, start, end)
            try:
                dp = record.genotype(sampleid)['DP']
            except:  # Skip ref/ref cals without AD field (this does exist)
                print("Locus {} has no DP field, skipping variant position as BAF could not be determined".format(locus))
                continue
            variant_call = record.genotype(sampleid).is_variant  # Returns True if a variant-call, and thus not a reference-call
            if dp >= args.mindepth:
                if variant_call == True: # Calculate ad and baf for variant-calls
                    ad = record.genotype(sampleid)['AD']
                    baf = (ad[1] / dp) * 100
                elif variant_call == False: # Calculate ad and baf for reference-calls
                    try:
                        ad = [record.genotype(sampleid)['AD']]
                        baf = ((dp - ad[0]) / (dp)) * 100
                    except:  # Skip ref/ref cals without AD field (this does exist)
                        print("Locus {} has no AD field, skipping variant position as BAF could not be determined".format(locus))
                        continue
                else: # In case variant status is not known (no call or ./.)
                    continue
            else:  # If threshold if not met, skip position
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
