'''
auxiliary functions
    parsing
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

def convert_time2human(sec):
    day = sec // (24 * 3600)
    sec %= 24 * 3600
    hour = sec // 3600
    sec %= 3600
    mini = sec // 60
    sec %= 60
    st = f"{day}d - {hour:d}:{mini:02d}:{sec:02d}"
    return st

def check_file(fname, st):
    with open(fname, 'r') as f:
        for line in f.readlines():
            if st in line:
                return True
    return False


