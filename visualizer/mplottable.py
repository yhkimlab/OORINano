#!/home/joonho/anaconda3/bin/python
''' read file and plot '''

import argparse
import re
import sys
import numpy as np
from .mplot2d import mplot_table_line
from .pltcommon import *
import csv

def fwhite2table(inf, icx=1):
    '''
    convert white space to table
    if column labels, stored in ytitle
    '''
    x=[]
    y2d=[]
    ytitle = []
    with open(inf,"r") as f:
        lines=f.readlines()
        table_lines = []     # y as 2d [ [y1], [y2], ...] this needs to be transposed
        ### if 1 file, get line-labels here
        for i, line in enumerate(lines):
            items = line.strip().split()
            if i == 0 and re.search('[a-zA-Z]', line):
                ytitle.extend(items)
                continue
            table_lines.append(items)   # 2D list

    for i, row in enumerate(table_lines):
        #print(f"use icx {icx}")
        x.append(row.pop(0))
        y2d.append(row)
    return x, y2d, ytitle

def fcsv2table(inf):
    '''
    convert csv to table
    '''
    fields=[]
    rows=[]
    with open(inf, 'r') as f:
        csvreader = csv.reader(f)
        fields = next(csvreader)
        for row in csvreader:
            rows.append(row)
        print("Total no. of rows: %d"%(csvreader.line_num))
    print('Field names are:' + ', '.join(field for field in fields))
    #print(f'{rows}')
    for row in rows:
        for col in row:
            print("%10s"%col, end=' ')
        print("")

    return fields, rows


def plot_dostable(inf,ix,icy,title,xlabel,ylabel,ylegend_in,colors,lvertical=None,orb=None, xlim=None):
    '''
    inf : table format
    fmt : white space or csv
    if file has column title, it will be ylabel
    dict_vert: 'center', 'max'
    '''
    ### modify get title
    if not title:
        title = inf.split('.')[0]
    ### define file format
    fnamelist = inf.split('.')
    if len(fnamelist) == 2:
        if fnamelist[1] == 'csv':
            fmt = 'csv'
        else:
            fmt = 'white'
    else:
        print(f"can't find extension using single . in {inf}")
        sys.exit(11)
    ### scan file list

    if fmt == 'white':
        x, y2d, ylegend = fwhite2table(inf)
        y2 = np.array(y2d).T
    elif fmt == 'csv':
        fields, rows = fcsv2table(inf)
        ylegend=fields[1:]
        tab = np.array(rows).T
        x   = tab[0,:]
        y2  = tab[1:,:]
    
    ### change string to value
    #print(f"{y2}")

    ### find x is int/float/string
    for i in range(10):
        if re.search('[a-zA-Z]', x[i]):
            print(f"xvalue {x[i]}")
            xl = 'string'
            break
        elif re.search('\.', x[i]):
            xl = 'float'
            break
        xl = 'int'
    #print(f"xvalue type {xl}")
    if xl == 'float':
        x = [ float(i) for i in x ]
    elif xl == 'int':
        x = [ int(i) for i in x ]

    ### change y value to float
    y2value = [ [ float(y) for y in ys ] for ys in y2 ]
    #print(f"{y2value}")
    #print(f"size: x {len(x)}, y: {len(ylegend)}, shape of data {np.array(y2value).shape}")

    ### var: x, y2, ylegend
    ys = []
    if not ylegend:
        ylegend=ylegend_in
    if not colors and ylegend:
        colors = []
        for yl in ylegend:
            colors.append(getcolor_orbital(yl))

    if not icy:
        icy =  list(range(len(y2value)))
    for i in icy:
        ys.append(y2value[i-2][:])
        
    #print(f"title {title} xlabel {xlabel} ylabel {ylabel}")
    #print(f"ylegends {ylegend}")
    ### x will be used for just len(x)
    mplot_table_line(x, ys, title=title, xlabel=xlabel, ylabel=ylabel, legend=ylegend, colors=colors,lvert=lvertical, orb=orb, xlim=xlim)
    return 0

def main():
    parser = argparse.ArgumentParser(description='Drawing files of table')

    parser.add_argument('inf', help='read table from file')
    parser.add_argument('-l', '--level', action='store_true', help='turn on to draw energy levels')
    parser.add_argument('-icx', '--icolumn_x', type=int, default=0, help='column index of X')
    g_file=parser.add_argument_group('Files', description="get input files")
    g_file.add_argument('-icy', '--icolumn_y', nargs='*', type=int, help='column indices of Y')
    g_file.add_argument('-ys', '--y_scale', default=[1], nargs="+", help='scale factor for Y [value|str|str-], use for str- for "-"')
    #g_file.add_argument('-ys', '--y_scale', nargs="*", help='scale factor for Y [value|str|str-], use for str- for "-"')
    g_twin = parser.add_argument_group('Twin-X', description='to plot using two y-axes')
    g_twin.add_argument('-tx', '--twinx', action="store_true", help='using two y-axes with twin x ticks')
    g_twin.add_argument('-icy2', '--second_iy', default=[2], nargs="+", type=int, help='designate the index of y for 2nd y-axis')
    g_twin.add_argument('-yl2', '--second_yl', default='G (eV)', help='input left y-axis title')
    #parser.add_argument('-j', '--job', help='job of qcmo|ai|gromacs')
    parser.add_argument('-t', '--title', default='PDOS', help='title of figure would be filename')
    parser.add_argument('-xl', '--xlabel', default=r'E - E\$_F\$ [eV]', help='X title, label in mpl')
    parser.add_argument('-yl', '--ylabel', default='DOS', help='Y title, label in mpl')
    parser.add_argument('-yls', '--ylabels', nargs='*', help='Y labels for legend')
    parser.add_argument('-c', '--colors', nargs='*', help='Y label for legend')
    parser.add_argument('-s', '--save', action='store_true', help='Save figure')
    args = parser.parse_args()

    ### n columns in 1 file, twinx
    draw_table(args.inf,args.level,args.icolumn_x,args.icolumn_y,args.title,args.xlabel,args.ylabel,args.ylabels,args.save,args.y_scale, args.colors, args.twinx, args.second_iy, args.second_yl)
    return 0

if __name__=='__main__':
	main()	


