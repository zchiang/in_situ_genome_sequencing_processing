# Author: Zachary Chiang, Buenrostro Lab, Harvard University
# Converts ex situ sequencing barcodes to colorspace 

import sys

# Cy5, FITC, Cy3, Texas red (PGP1_001)
#dibase2color = {"AA":2, "AC":3, "AG":4, "AT":1, "CA":3, "CC":2, "CG":1, "CT":4, "GA":4, "GC":1, "GG":2, "GT":3, "TA":1, "TC":4, "TG":3, "TT":2}

# Cy5, Cy3, FITC, Texas red (EMBRYO_002)
dibase2color = {"AA":3, "AC":2, "AG":4, "AT":1, "CA":2, "CC":3, "CG":1, "CT":4, "GA":4, "GC":1, "GG":3, "GT":2, "TA":1, "TC":4, "TG":2, "TT":3}

for line in sys.stdin:

        # find barcode and reverse

	column = line.rstrip().split()
	seq = column[2]
	seq = seq[::-1]

        if True: # can add a QC check here

		cs_read = ""

                # N
                
                cs_read = cs_read + str(dibase2color[seq[3:5]]) # cycle 1, L
		cs_read = cs_read + str(dibase2color[seq[8:10]]) # cycle 2, LL
		cs_read = cs_read + str(dibase2color[seq[13:15]]) # cycle 3, LLL
		cs_read = cs_read + str(dibase2color[seq[18:20]]) # cycle 4, LLLL
                
                # N-1

		cs_read = cs_read + str(dibase2color[seq[2:4]]) # cycle 5,L
		cs_read = cs_read + str(dibase2color[seq[7:9]]) # cycle 6, LL
		cs_read = cs_read + str(dibase2color[seq[12:14]]) # cycle 7, LLL
                cs_read = cs_read + str(dibase2color[seq[17:19]]) # cycle 8, LLLL

                # N-2

		cs_read = cs_read + str(dibase2color[seq[1:3]]) # cycle 9, L
		cs_read = cs_read + str(dibase2color[seq[6:8]]) # cycle 10, LL
		cs_read = cs_read + str(dibase2color[seq[11:13]]) # cycle 11, LLL
		cs_read = cs_read + str(dibase2color[seq[16:18]]) # cycle 12, LLLL

                # N-3
		
                cs_read = cs_read + str(dibase2color[seq[0:2]]) # cycle 13, L
		cs_read = cs_read + str(dibase2color[seq[5:7]]) # cycle 14, LL
		cs_read = cs_read + str(dibase2color[seq[10:12]]) # cycle 15, LLL
		cs_read = cs_read + str(dibase2color[seq[15:17]]) # cycle 16, LLLL

                # N-4

		cs_read = cs_read + str(dibase2color[seq[4:6]]) # cycle 17, LL
		cs_read = cs_read + str(dibase2color[seq[9:11]]) # cycle 18, LLL
		cs_read = cs_read + str(dibase2color[seq[14:16]]) # cycle 19, LLLL
	
                # add to correct column

		if len(column) == 5:
			column.append(cs_read)
		else:
			column[5] = cs_read

		print "\t".join(column)

