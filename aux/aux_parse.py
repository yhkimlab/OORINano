'''
auxiliary functions for parsing
'''
import re


def parse_line(line):
    '''
    read line and return list
    '''
    li = line.strip().split()
    return li

def parse_lines(lines):
    '''
    read line block return 2d list
    '''
    li2d=[]
    for line in lines:
        li = parse_line(line)
        li2d.append(li)

    return li2d