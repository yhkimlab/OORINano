import re, sys

def getcolor_orbital(yl):
    '''
    plot color depending on legend
    '''
    ### as for dos plot
    if yl == 's':
        return 'r'
    elif yl == 'p':
        return 'b'
    elif yl == 'd':
        return 'm'
    elif re.search('t', yl, re.IGNORECASE) :
        return 'lightgray'
    ### add more color for legend
    else:
        print("legend cannot define color")
        sys.exit(11)


