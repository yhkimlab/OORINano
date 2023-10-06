'''
    ln -s env_yours.py env.py
        link softly your envfile to env.py
'''

from __future__ import print_function

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
