import re, sys
from ..utils.auxil import whereami

def getcolor_orbital(yl):
    '''
    plot color depending on legend
    '''
    ### as for dos plot
    if re.match('s', yl) and not re.match('sum', yl):  # end with 's'
        return 'r'
    elif re.match('p', yl):
        return 'b'
    elif re.match('d', yl):
        return 'm'
    elif re.match('sum', yl, re.IGNORECASE) :
        return 'lightgray'
    ### add more color for legend
    else:
        print(f"legend cannot define color for {yl} in {__name__}")
        sys.exit(11)


