'''
auxiliary functions
    parsing
    time2human
    check_file
    
'''
import re, inspect, os, glob
import numpy as np

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

def list2str(li, delimit=None):
    '''
    list to string
        delimit delimiter between list elements, default=''
    '''

    if delimit == None:
        st = "".join(str(x) for x in li)
    else:
        st = delimit.join(str(x) for x in li)
    return st

def list2dict(li):
    it = iter(li)
    dic = dict(zip(it,it))
    return dic

def print_list(li):
    st = list2str(li)
    print(st)
    return st

def list2_2f(li, decimal=2):
    return list(np.around(np.array(li), decimal))


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
