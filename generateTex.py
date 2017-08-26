#!/usr/bin/env python
import re
import glob
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--dir", help="dir that has all the files to be included", default="")
parser.add_argument("--ext", help="extenshion of files to be included", default="pdf")
parser.add_argument("--regex",help="regular expression to replace the use of dir and ext", default=None)
parser.add_argument("--caption", help="name of caption file in same dir as png", default=None)
parser.add_argument("--out", help="name of tex file, will have .tex added to it", default="booklet")
args = parser.parse_args()

mydir = args.dir
ext = args.ext
captionFile = args.caption 
myglob = mydir + "/*." + ext
if(args.regex is not None):
    myglob = args.regex
outfile = open(args.out + ".tex", "w+")


header='''\\documentclass{article}
\\usepackage{graphicx}
\\usepackage{fancyvrb}
\\usepackage[margin=0.15in]{geometry}
\\usepackage{hyperref}

\\pdfminorversion=5

\\begin{document}
\\tableofcontents
\\newpage
'''

footer='''
\\end{document}
'''

figureFormat='''
\\begin{{figure}}[!ht]
\\centering
\\includegraphics[width=0.8\\textwidth,height=0.8\\textwidth]{{{{{}}}{}}}
\\end{{figure}}
'''
#\\caption[LoF entry]{{ {} }}

mainText=""


def get_caption(filenmae):
    rtn = ""
    if(captionFile is None):
        return("")
    path = os.path.dirname(filenmae)    
    myfile = os.path.join(path, captionFile)
    if(not os.path.exists(myfile)):
        return("")

    f = open(myfile)
    rtn += "\n\\begin{Verbatim}[fontsize=\\scriptsize]\n" 
    rtn += f.read()
    rtn += "\\end{Verbatim}\n"
    #rtn += "\\newpage \n"
    #rtn  = re.sub("_", "-", rtn)
    #rtn = re.sub("\n", r"\\\\\hspace{\\textwidth}", rtn) 
    return(rtn)


for idx, filename in enumerate(glob.glob(myglob)):
    filename=os.path.abspath(filename)
    picname, ext = os.path.splitext(filename)
    caption = get_caption(filename)
    
    if(len(caption) > 0):
        mainText += "\\section{{ {} }}\n".format(caption.split("\n")[2].split(":")[1] )
    else:
        mainText += "\\section{{ {} }}\n".format( idx )
    mainText +=  figureFormat.format( picname, ext)#, caption)
    mainText +=  caption 
    mainText +=  "\\newpage \n" 


outfile.write(header+mainText+footer)
outfile.close()



splitMain = mainText.split("\\newpage")
counter = 0
running = ""
for idx, text in enumerate(splitMain):
    running += text + "\\newpage \n" 
    if(  (idx+1) % 50 == 0  ):
        f = open("booklets/booklet{}.tex".format(counter), "w+")
        f.write(header+running+footer)
        f.close()
        running = ""
        counter += 1









