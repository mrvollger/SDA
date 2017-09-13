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
    dfs = [] 
    for idx, directory in enumerate(LAS):
        LA = LocalAssembly(directory, psvGraphLoc, psvURL).asPD()
        dfs.append(LA)
        print(type(LA), idx)
        #if(idx == 20):
        #    break
        #out += str(temp)
        #if(idx == 0 ):
        #    header = temp.getHeader()
        #if(idx == 5): break
    
    #    
    # create df
    #
    df = pd.concat(dfs, axis=0, ignore_index=True)
    # get a list of columns
    cols = list(df)
    length = len(cols)
    # move the column to head of list using index, pop and insert
    cols.insert(0, cols.pop(cols.index('numOfPSVs')))
    cols.insert(0, cols.pop(cols.index('copies_in_reference')))
    cols.insert(0, cols.pop(cols.index('number_of_CC_groups')))
    cols.insert(0, cols.pop(cols.index('CC_ID')))
    cols.insert(0, cols.pop(cols.index('Status')))
    cols.insert(0, cols.pop(cols.index('region_in_falcon')))
    # move the cc matrix to the back:
    counter = 0
    while(True):
        key = str(counter) + "."
        if(key == "0."):
            key = "copy=>"
        if(key not in cols):
            print("done changing cols", key)
            break
        cols.insert(length-1, cols.pop(cols.index(key)))
        counter += 1

    print(cols)
    # reorder the columns 
    df = df.loc[:, cols]


    #out = header + out
    #print(out)
    #dataStream = StringIO(out)
    #df = pd.read_csv(dataStream, sep="\t", header = 0)
    # sort the way I want it sorted
    #df = df.sort_values(['Status', 'region_in_falcon', "CC_ID"], ascending=[0, 1, 1])

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



