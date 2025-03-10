#!/usr/bin/env python3
import sys
from PyPDF2 import PdfFileReader
def pdf_list_anchors(fh,ofh):
    reader = PdfFileReader(fh)
    destinations = reader.getNamedDestinations()
    for name in destinations:
        ofh.write(name+"\n")
f = open(sys.argv[2],'w')
pdf_list_anchors(open(sys.argv[1],'rb'),f)

