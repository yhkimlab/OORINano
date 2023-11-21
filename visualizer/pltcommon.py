import re, sys

def getcolor_orbital(yl):
    '''
    plot color depending on legend
    '''
    ### as for dos plot
    if re.match('s', yl):
        return 'r'
    elif re.match('p', yl):
        return 'b'
    elif re.match('d', yl):
        return 'm'
    elif re.search('t', yl, re.IGNORECASE) :
        return 'lightgray'
    ### add more color for legend
    else:
        print(f"legend cannot define color for {yl}")
        sys.exit(11)


