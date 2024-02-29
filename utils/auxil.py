'''
auxiliary functions
    string to number
    parsing
    time2human
    check_file
    
'''
import re, inspect, os, glob, sys
import numpy as np

### change string (read input file) to value
re_num = re.compile("^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$")   # string for all kinds of number

def isnumber(str):
    '''
    str     string
    return  number or not
    '''
    if re_num.match(str):
        return True
    else:
        return False

def isint(value):
    '''
    value   input string is value
    return  try to convert string to integer
    '''
    flag = True
    try:
        int(value)
    except ValueError:
        flag = False
    return flag


### 
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

def get_digits_4str(word, string):
    ''' obtain digits after key-word '''
    if word in string:
        substr = re.split(word,string)[1]
        mat = re.match('\d+', substr)
        digits = mat.group()
    else:
        digits = None
    return digits


def whereami(rank=1):
    '''
    rank can be [0,1]
    '''
    return inspect.stack()[rank][3]

def list2str_old(li, delimit=None):
    '''
    list to string
        delimit delimiter between list elements, default=''
    '''

    if delimit == None:
        st = "".join(str(x) for x in li)
    else:
        st = delimit.join(str(x) for x in li)
    return st

def list2str(lis, decimal = 0, head=None, end=None):
    '''
    list to string
    delimit     delimiter between list elements, default=''
    decimal     change float to decimal place in case float

    '''
    ### change float to the formmated decimal placement
    if type(decimal) == int:
        lis = list(np.around(np.array(lis), decimal))

    #if delimit == None:
    #    st = "".join(str(x) for x in lis)
    #else:
    #    st = delimit.join(str(x) for x in lis)

    st = "   ".join(['{:10.2f}'.format(i) for i in lis])
    if end:
        st += '\n'

    if head:
        headst = '{:^10}'.format(head)
        st = headst + st
    return st

def flist2str(lis, delimit='  ', head=None, decimal=2):
    '''
    float list to string
        delimit delimiter between list elements, default=''
    formatted   for decimal number
    '''
    if decimal == 2:
        st = delimit.join('%10.2f' % x for x in lis)
    else:
        print("Exit: more format is required")
        sys.exit(111)

    if head:
        headst = '%10s' % head
        st = headst + st
    return st


def list2dict(li):
    it = iter(li)
    dic = dict(zip(it,it))
    return dic

def print_list(li):
    st = list2str(li)
    print(st)
    return st

def list2_format(li, decimal=2):
    return list(np.around(np.array(li), decimal))

def headlist2str(li, head, Lwrite=True):

    ### formated list string
    #st = str([ f'{x:7.2f}' for x in li])

    st = f"{head:^10s}" + str(li)
    if Lwrite:
        st += '\n'
    return st

def f_ext(fname):
    return fname.split('.')[-1]
def f_root(fname):
    return fname.split('/')[-1].split('.')[0]
fname_root = f_root
fname_ext  = f_ext
def fname_decom(fname):
    fname_parts = fname.split('.')
    if len(fname_parts) < 2:
        print("Error: %s has not dot in file name" % fname)
        exit(1)
    suff = fname_parts.pop()
    fn = '_'.join(fname_parts)
    return fn, suff        

fname_parsing = fname_decom

def find_file(dname, f_root):
    files = os.listdir(dname)
    for f in files:
        if f.endswith(f_root):
            return f
        elif f.startswith(f_root):
            return f
    return None

def search_dirs(dir_pre, fname):
    dirs =  glob.glob(f"{dir_pre}*")
    list_dir=[]
    #print(fname, dirs)
    for d in dirs:
        if os.path.isfile(f"{d}/{fname}"):
            list_dir.append(d)
    return list_dir            
        

def expand_dim_str(lstring):
    """ expand 1D string to 2D string """
    a = np.array(lstring)
    dim = len(a.shape)
    new_2d=[]
    if dim == 1:
        ### ele is string
        for ele in a:
            new_list=ele.split()
            new_2d.append(new_list)
    return new_2d

def f_parsing(fname, sep='num'):
    with open(fname, 'r') as f:
        lines = f.readlines()
        i=0
        y=[]
        for line in lines:
            if re.search('[a-zA-Z]', line):
                pass
            else:
                i += 1
                line = line.strip()
                ele = line.split()
                if i == 1:
                    n = len(ele)
                    for j in range(n):
                        y.append([])         # declare list as many as the number of column
                for j in range(n):
                    y[j].append(float(ele[j]))
    for j in range(n):
        if j > 0:
            if len(y[j-1]) != len(y[j]):
                print("in zip file to column, Error: length of colums are different")
                sys.exit(1)
    return y

def fname_index(pre, suff, cwd):
    files = os.listdir(cwd)
    for f in files:
        if re.match(pre, f) and re.search(suff, f):
            pass
    return 0
