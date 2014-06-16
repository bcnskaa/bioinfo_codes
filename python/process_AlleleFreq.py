#########################
#
# Copyright (c) 2013, SK Woolf <bcnskaa@gmail.com>.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#########################
import sys, getopt

def main(argv):
    vcf_fn = ''

    try:
        opts, args = getopt.getopt(argv, "i:o")
    except getopt.GetoptError:
        #print "process_AlleleFreq.py -i VCF-FRQ-INFILE -o PROCESSED-VCF-FRQ-OUTFILE"
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-i':
            vcf_fn = arg
            out_fn = vcf_fn + ".processed_freq"

    #print "Reading from ", vcf_fn

    vcfs = open(vcf_fn).read()
    vcfs = vcfs.splitlines()

    #print 'Number of vcf read: ', len(vcfs)


    with open(out_fn, "wt") as f:
        f.write("CHROM\tPOS\tN_ALLELES\tN_CHR\tALLELE_A\tALLELE_A_FREQ\tALLELE_B\tALLELE_B_FREQ\n");
        for vcf in vcfs:
            if vcf[0] != 'C':
                items = vcf.split("\t")
                
                allele_A = items[4].split(":")
                allele_B = items[5].split(":")
                f.write(items[0] + "\t" + items[1] + "\t" + items[2] + "\t" + items[3] + "\t" + allele_A[0] + "\t" + allele_A[1] + "\t" +allele_B[0] + "\t" +allele_A[1] + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
