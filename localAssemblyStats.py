#!/usr/bin/env python
import argparse
import os
import sys
import string 
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO
from LocalAssembly import LocalAssembly
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("localAssemblyDirectories", help="list of all of the directoreis with local assemblies" )
parser.add_argument("--psvGraphURL", help="url where the psv graph can be found", default="https://eichlerlab.gs.washington.edu/help/mvollger/psvGraphs/" )
parser.add_argument("--psvGraphLoc", help="place where all the psv graphs will be linked to", default="/net/eichler/vol2/home/mvollger/public_html/psvGraphs/" )
parser.add_argument("--out", help="where to place the table of generated results", default="/net/eichler/vol2/home/mvollger/public_html/localAssemblyStats.tsv" )
args = parser.parse_args()
localAssemblyDirectories = args.localAssemblyDirectories
psvGraphLoc = os.path.abspath( args.psvGraphLoc.strip() )
psvURL = args.psvGraphURL

def main():
    out = ""
    header = ""
    LAS = open(localAssemblyDirectories)
    
    for idx, directory in enumerate(LAS):
        temp = LocalAssembly(directory, psvGraphLoc, psvURL)
        out += str(temp)
        if(idx == 0 ):
            header = temp.getHeader()
        #if(idx == 5): break
    
    #    
    # create df
    #
    out = header + out
    dataStream = StringIO(out)
    df = pd.read_csv(dataStream, sep="\t", header = 0)
    # sort the way I want it sorted
    df = df.sort_values(['Status', 'region_in_falcon', "CC_ID"], ascending=[0, 1, 1])

    # 
    # write to tsv
    #
    df.to_csv(args.out, sep='\t', index = False)
    
    #
    # write to excel 
    #
    # add hyperlinks
    df['psvURL'] = '=HYPERLINK("' + df['psvURL'].astype(str) + '")'
    #print(df) 
    # setup the writer
    writer = pd.ExcelWriter(args.out + ".xlsx")
    df.to_excel(writer,'Sheet1', index = False)
    worksheet = writer.sheets['Sheet1']
    
    # set width based on column name width
    for idx, letter in enumerate(string.uppercase[:len(list(df))]):
        length = len(list(df)[idx]) + 5
        worksheet.column_dimensions[letter].width = length
    
    # save excel
    writer.save()

main()



