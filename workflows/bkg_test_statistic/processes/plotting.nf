
process CrossSection{

 input:
  path gibuu_files;

 output:
  path "*.png", emit: plt;
  path "*.txt", emit: txt;

 script:
 """
#!/usr/bin/python3
import sys
import re
import glob
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def extract_xs(files, process):
    dict={'process': process,'e' : [], 'xs': []}

    for file in files:
        f = open(file,'r')
        lines = f.readlines()
        dict['e'] .append(float(lines[1].split()[0]))
        dict['xs'].append(float(lines[1].split()[1]))
    return dict

def xsDicts2Ascii(dictionaries, filename):
    for i in range(1,len(dictionaries)):
        if (dictionaries[i]['e']!=dictionaries[i-1]['e']):
            sys.exit("Energies are not the same")

    header='#'
    table=[]
    table.append(dictionaries[0]['e'])
    for d in dictionaries:
        table.append(d['xs'])
        header += '\\t' + d['process']
    with open(filename, "a") as f:
        f.write(header+"\\n")
        np.savetxt(f,np.column_stack(table))

xs_files = sorted("${gibuu_files}".split())

regexp = re.compile('xs_(\\w+)_(\\w+)_(.+)_Z(\\d+)A(\\d+)_(\\w+).dat')

processes = defaultdict(list)

flavor=''
interaction=''
z=''
a=''
for file in xs_files:
    basename=os.path.basename(file)
    if (regexp.match(basename)!=None):
        fields = regexp.split(basename)
        flavor = fields[1]
        interaction = fields[2]
        z = fields[4]
        a = fields[5]
        if (fields[6]!='numbers' and fields[6]!='highRES'):
            processes[fields[6]].append(file)

output_ascii = "xs_"+flavor+"_"+interaction+"_Z"+z+"A"+a+".txt"
output_plt   = "xs_"+flavor+"_"+interaction+"_Z"+z+"A"+a+".png"

xs_dicts = []
for key in processes:
    d = extract_xs(processes[key], key)
    xs_dicts.append(d)
    plt.plot(np.array(d['e']),np.array(d['xs'])/np.array(d['e']),
             marker = '.', label=key, linestyle='--')

xsDicts2Ascii(xs_dicts, output_ascii)

plt.ylabel('Cross section/whatever_units_gibuu_produces')
plt.xlabel('Energy/GeV')
plt.grid(linestyle='--', which='both')
plt.yscale('log')
plt.legend()
plt.savefig(output_plt)

 """
}
