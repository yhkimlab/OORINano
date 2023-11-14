"""
    2d plot library
    draw_2d
    fdraw   where f == "file"
    barplot
"""    
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import sys, re
import scipy
import numpy as np
from .pltcommon import *
from ..utils.auxil import *

#plt.switch_backend('agg')
size_title = 20
size_label = 18
#size_tick = 15
text_x = 0.75
text_y = 0.8
text_twinx_x = 0.8
text_twinx_y = 0.9

### functions
#   my_font
#   common_figure: will be moved to myplot_default and deleted
#   common_figure_after: needed after plt.plot()
#   draw_dots_two: option: twinx
#   xtitle_font
#   mplot_twinx
#   mplot_nvector
#   draw_histogram
#

#from myplot_default import *

def my_font(pack='amp'):
    if pack == 'amp':
        font={  'family': 'normal',
                'weight': 'bold',
                'size'  : 22,
                }
    else:
        print("package: %s is not included" % pack)
        sys.exit(0)
    mpl.rc('font', **font)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    return
### choice between "import myplot_default"|call common_figurea

def lighten_color(color, ncolors):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    print(f"h l s {c[0]} {c[1]} {c[2]}")
    lmin = 0.2
    lmax = c[1]*1.9
    dl  = (lmax-lmin)/(ncolors-1)
    print(f"min, max, dl {lmin} {lmax} {dl}")
    color_rgb=[]
    for i in range(ncolors):
        print(f"lightness {lmin + dl*i}")
        color_rgb.append(colorsys.hls_to_rgb(c[0], lmin + dl*i, c[2]))
    amount=0.5 
    lightness = c[1]*0.3 + amount * (1 - c[1])
    print(f" lightness: {lightness }")
    #return colorsys.hls_to_rgb(c[0], lightness, c[2])
    return color_rgb
    #return colorsys.hls_to_rgb(c[0], amount * (1 - c[1]), c[2])

def lighten_2color(color, ncolors):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys

    color1=color[0]
    color2=color[1]
    nlight = int(ncolors/2)
    ndark  = ncolors -nlight

    c = colorsys.rgb_to_hls(*mc.to_rgb(color1))
    print(f"h l s {c[0]} {c[1]} {c[2]}")
    lmin = 0.25
    lmax = 0.75
    dl  = (lmax-lmin)/(ncolors-1)*2
    print(f"min, max, dl {lmin} {lmax} {dl}")
    color_rgb=[]
    color_light=[0.5, 0.8]
    for i in range(nlight):
        print(f"lightness {lmin + dl*i}")
        color_rgb.append(colorsys.hls_to_rgb(c[0], color_light[i], c[2]))
    c = colorsys.rgb_to_hls(*mc.to_rgb(color2))
    print(f"h l s {c[0]} {c[1]} {c[2]}")
    color_light=[0.8, 0.6, 0.3]
    for i in range(ndark):
        print(f"darkness {lmax - dl*i}")
        color_rgb.append(colorsys.hls_to_rgb(c[0], color_light[i], c[2]))
    #amount=0.5 
    #lightness = c[1]*0.3 + amount * (1 - c[1])
    #print(f" lightness: {lightness }")
    #return colorsys.hls_to_rgb(c[0], lightness, c[2])
    return color_rgb
    #return colorsys.hls_to_rgb(c[0], amount * (1 - c[1]), c[2])

def common_figure(ctype='dark', ncolor=4, Ltwinx=False):
    '''
    ctype   darken to change intensity
            cycle to use designated color turns
    '''
    if ctype == 'darken':
        import matplotlib.colors as ms
        import colorsys
    else:
        from cycler import cycler

    fig = plt.figure(figsize=(10,6))
    ax = plt.axes()
    mpl.rcParams.update({'font.size':12})
    #ax.tick_params(axis='both', which='major', labelsize=25)
    #ax.tick_params(axis='x', labelsize=30)
    if ctype == 'darken':
        pass   
    else:
        if Ltwinx:
            print(f"make twinx axis")
            ax2 = ax.twinx()
            custom_cycler = (cycler(color=['r','m', 'orange']))
            custom_cycler2 = (cycler(color=['b', 'g']))
            ax.set_prop_cycle(custom_cycler)
            ax2.set_prop_cycle(custom_cycler2)
            return fig, ax, ax2
        else:
            if ncolor == 2:
                custom_cycler = (cycler(color=['r','g']))                                          # Figure 8(b)
                custom_cycler = (cycler(color=['darkcyan','b']))                                   # Figure S16(a)
                custom_cycler = (cycler(color=['r','darkcyan']))                                    # Figure S16(b)
            elif ncolor == 3:
                custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))       # Figure 8(a)
                custom_cycler = (cycler(color=['r','m','g','b']))
            elif ncolor == 4:
                custom_cycler = (cycler(color=['orange','m','g','b'])+ cycler(lw=[1,1,1,2]))       # Figure 8(a)
                custom_cycler = (cycler(color=['r','m','g','b']))
            else:
                custom_cycler = (cycler(color=['r','m','g','b']))
            ax.set_prop_cycle(custom_cycler)
            return fig, ax
    return fig, ax
    
### draw_dots_two was upgraded to twinx
def draw_2subdots(y, h, title, suptitle, Ltwinx=None, escale=1.0,Colors=['r','b','o'], Ldiff=True):
    '''
    this makes error in serial plotting
    '''
    fig, ax = common_figure()
    escale = 1.0
    nlen = len(y)
    h_conv = np.array(h) * escale        # escale = my_chem.ev2kj
    y_conv = np.array(y) * escale
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())
    max_res = abs(max(diff, key=abs))
    #max_res = max(diff, key=abs)
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text_pos_x = nlen*0.85                  # 0.85, 0.2
    text_pos_y = max(y_conv)*0.2
    text="E_rms(test) = {:7.3f}\nE_maxres = {:7.3f}".format(rmse, max_res)

    if Colors:  color = Colors.pop(0)       #'tab:' + Colors.pop(0)
    else:       color='tab:green'
    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    #mpl.rcParams.update({'font.size':22})
    plt.title(title, fontsize=20)
#    plt.suptitle(suptitle, x=0.5, y=0.92, va='top', fontsize=18)
    plt.suptitle(suptitle, fontsize=10)
    if escale == 1.0:
        ax.set_ylabel('PE(eV)', color='b', fontsize=15)
    elif escale == my_chem.ev2kj:
        ax.set_ylabel('PE(kJ/mol)', fontsize=15)
    plt.xlabel('data', fontsize=15)
    #ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='y', labelcolor='b', labelsize=10)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    if Ltwinx:
        ax2=ax.twinx()
        ax2.set_ylabel("Difference(eV)", color='g')
        #ax2.plot(x, ys[1,:], 'o-', color=color)
        ax2.tick_params(axis='y', labelcolor='g', labelsize=10)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    p1  = ax.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = ax.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    print(y_conv)
    
    if Ltwinx:
        if Ldiff:
            p3, = ax2.plot(range(nlen), diff, c='g', label='difference')
            plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
        else:
            plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_twinx_x, text_twinx_y, text, fontsize=10, transform=ax.transAxes)
    else:
        #p3, = plt.plot(range(nlen), diff, c='g', label='difference')
        #plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
        plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_x, text_y, text, fontsize=10, transform=ax.transAxes)
    plt.plot(range(nlen), ones)

    plt.show()
    return rmse, max_res
### draw_dots_two was upgraded to twinx
### still used by amp_run.py
def draw_amp_twinx(y, h, title, suptitle, natom=1, Ltwinx=None, escale=1.0,Colors=['r','b','o'], Ldiff=True):
    '''
    this makes error in serial plotting
    '''
    fig, ax = common_figure()
    escale = 1.0
    nlen = len(y)
    h_conv = np.array(h) * escale        # escale = my_chem.ev2kj
    y_conv = np.array(y) * escale
    ymin = min(y_conv)
    ymax = max(y_conv)
    y_width = ymax - ymin
    diff =  np.subtract(h_conv,y_conv)
    rmse = np.sqrt((diff**2).mean())/natom                      # eV/atom
    max_res = abs(max(diff, key=abs))/natom                     # eV/atom
    #print("{:10.3f} {:10.3f}".format(rmse,max_res))
    ### input text inside figure
    text_pos_x = nlen*0.85                  # 0.85, 0.2
    text_pos_y = max(y_conv)*0.2
    text="E_rms(test) = {:7.3f} eV/atom\nE_maxres   = {:7.3f} eV/atom".format(rmse, max_res)

    if Colors:  color = Colors.pop(0)       #'tab:' + Colors.pop(0)
    else:       color='tab:green'
    ones = np.zeros((len(y_conv)))
    #my_font('amp')
    #mpl.rcParams.update({'font.size':22})
    plt.title(title+'\n', fontsize=size_title)
    #suptitle=suptitle+'\n'
    plt.suptitle(suptitle, x=0.5, y=0.96, va='top', fontsize=18)
    #plt.suptitle(suptitle, fontsize=size_tick)
    if escale == 1.0:
        ax.set_ylabel('PE(eV)', color='b', fontsize=size_label)
        ax.set_ylim(ymin-1, ymax+0.25)
    elif escale == my_chem.ev2kj:
        ax.set_ylabel('PE(kJ/mol)', fontsize=size_tick)
    plt.xlabel('data', fontsize=size_label)
    #ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='y', labelcolor='b', labelsize=size_tick)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    if Ltwinx:
        ax2=ax.twinx()
        ax2.set_ylabel("Difference(eV)", color='g')
        ax2.set_ylim(-0.001, 0.01)
        #ax2.plot(x, ys[1,:], 'o-', color=color)
        ax2.tick_params(axis='y', labelcolor='g', labelsize=size_tick)
    #plt.scatter(x, y, 'r', x, y_bar, 'b')
    #print(y_conv, h_conv)
    p1  = ax.scatter(range(nlen), y_conv, c='r', marker='o', label='true value')
    p2  = ax.scatter(range(nlen), h_conv, c='b', marker='^', label='hypothesis')
    if Ltwinx:
        if Ldiff:
            p3, = ax2.plot(range(nlen), diff, c='g', label='difference')
            plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.2))
        else:
            plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_twinx_x, text_twinx_y, text, fontsize=size_tick, transform=ax.transAxes)
    else:
        #p3, = plt.plot(range(nlen), diff, c='g', label='difference')
        #plt.legend([p1,p2,p3],['true value', 'hypothesis', 'difference'],loc=(0.0, 0.1))
        plt.legend([p1,p2],['true value', 'hypothesis'],loc=(0.0, 0.1))
        plt.text(text_x, text_y, text, fontsize=10, transform=ax.transAxes)
    plt.plot(range(nlen), ones)

    plt.show()
    return rmse, max_res





def xtitle_font(tit):
    st = "\'{}\', fontsize=20".format(tit)
    print(st)
    return st


### twinx1: used for md.ene, normal data file
def mplot_twinx(x, y, iy_right, title=None, xlabel=None, ylabel='E [eV]', legend=None, Lsave=False, Colors=None):
    '''
    called from "amp_plot_stat.py"
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    len(ylabel) == 2
    '''
    if iy_right:
        fig, ax, ax2 = common_figure(ctype='cycle', ncolor = len(y), Ltwinx=True)
    else:
        fig, ax = common_figure(ctype='cycle', ncolor = len(y))
    #fig, ax = plt.subplots(figsize=(12,8))
    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    ### function: mplot_twinx
    plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=25)
    if isinstance(ylabel, str): ylabel1 = ylabel2 = ylabel
    elif isinstance(ylabel, list):
        ylabel1 = ylabel[0]
        ylabel2 = ylabel[1]
    ### try autocolor
    #plt.ylabel(ylabel1, fontsize=25, color='r')
    plt.ylabel(ylabel1, fontsize=25)
    #ax.tick_params(axis='y', colors='r')
    ax.tick_params(axis='y')
    #ax.set_ylim(-2,2)
    #ax.xaxis.set_major_locator(plt.NullLocator())
    #print(f"x, y shape:: {np.array(x).shape} {np.array(y).shape} and ylabel {ylabel} in {whereami()}")
    #ax2.set_ylabel(ylabel2, fontsize=25, color='g')
    ax2.set_ylabel(ylabel2, fontsize=25)
    ax2.set_ylim(-2,2)
    pls=[]
    print(f"{iy_right} {len(ys)} {whereami()}")
    for i in range(len(ys)):
        if i in iy_right: 
            #if Colors:  color = Colors.pop(i)       #'tab:' + Colors.pop(0)
            #else:       color='tab:green'
            #plt.yticks(color='g')
            #p2, = ax2.plot(x, ys[i,:], '-', color='g', label=legend[i])
            p2, = ax2.plot(x, ys[i,:], '-', label=legend[i])
            pls.append(p2)
        else:
            #ax2.tick_params(axis='y')
            #if Colors:  color = Colors.pop(i)       #color = 'tab:' + Colors.pop(0)
            #else:       color = 'tab:red'
            #print(f"shape of x, ys[i] = {np.array(x).shape} {ys[i,:].shape}")
            #plt.yticks(color='r')
            #p1, = ax.plot(x, ys[i,:], '-', color='r',  label=legend[i])
            p1, = ax.plot(x, ys[i,:], '-', label=legend[i])
            pls.append(p1)
    #plt.legend(pls, Ylabels, loc=2)
    ax.legend(loc=2)            # 2
    ax2.legend(loc=4)           # 1
    plt.legend()
    #common_figure_after()
    #x_ticks = ['PP', 'PPP', 'PNP', 'PNP-bridged']
    #plt.xticks(x_ticks)
    #ax.set_xticklabels(x_ticks)
    #plt.locator_params(axis='x', nbins=10)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0

def make_double(x, y):
    xnew = []
    ynew = []
    line_pairs=[]
    for i, (x1, y1) in enumerate(zip(x, y)):
        xnew.append(x1)
        xnew.append(x1)
        ynew.append(y1)
        ynew.append(y1)
        if y1 is not None:
            line_pairs.append([[2*i, 2*i+1],[y1,y1]])
    return xnew, np.array(ynew), line_pairs


def mplot_levels(x, ys, title=None, xlabel=None, ylabel=None, legend=None,Lsave=False, colors=None):
    '''
    plot table 
    x   from 1st column
    y   other columns [0,1,2...] in  [[column],[size],...[size]]
    '''
    fig, ax = common_figure()
    ys = np.array(ys)
    if len(x) != ys.shape[1]:
        print(f"error in shape: x {len(x)}, ys.shape {ys.shape}")
    
    plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=35)
    plt.ylabel(ylabel, fontsize=35)
    
    if not colors:
        colors = ["blue", "red", "green", "orange", "purple"]
    else:
        colors = colors
    ### make levels: x=list, ys=array
    for i, y in enumerate(ys):
        xd, yd, linepair = make_double(x, y)
        print(f"Double: xd {xd} yd {yd}")
        xs = np.arange(len(xd))
        print(f"{np.max(xs)}")
        series = yd.astype(np.double)
        mask = np.isfinite(series)
        plt.plot(xs[mask], series[mask], '--', color=colors[i], label=legend[i])
        ### overdraw thick level
        for xl, yl in linepair: 
            plt.plot(xl,yl, '-', color=colors[i], lw=5 )
    #ax.yaxis.set_major_locator(ticker.MultipleLocator())
    #ax.set_xlim(-0.5, np.max(xs)+0.5)
    plt.legend(loc=1)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0


def auto_nvector(x,y):
    fig, ax = plt.subplots()
    for i in range(len(y)):
        plt.plot(x,y[i])
    plt.show()
    return 0
### most1
def mplot_table_line(x, y, dx=1.0, title=None, xlabel=None, ylabel=None, legend=None,Lsave=False, colors=None, lvert=None, v_legend=None, orb=None, xlim=None):
    '''
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[size],[size],...[size]]
    '''
    if not colors:
        fig, ax = common_figure(ncolor = len(legend))
    else:
        fig, ax = common_figure()

    #plt.locator_params(axis='x', nbins=10)
    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    #print(f"x={len(x)} y={ys.shape} in {whereami()}()")
    plt.title(title)
    xlabel = 'E - E$_F$ [eV]'
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    print(f"xlabel: {xlabel} ")
    if ys.ndim == 1:
        plt.plot(x, y, '-')
        #plt.scatter(x, y)
    elif ys.ndim == 2:
        for i in range(len(ys)):
            if re.search('t', legend[i]):
                d = scipy.zeros(len(ys[i,:]))
                print(f"shape d {np.array(d).shape}")
                ax.fill_between(x, ys[i,:], where=ys[i,:]>=d, color=colors[i])
                plt.plot(x,ys[i,:],  label='TDOS' , color=colors[i])
            else:
                plt.plot(x,ys[i,:],  label=legend[i] , color=colors[i])

    else:
        print(f"Error:: obscure in y-dim {ys.ndim}")
    
    ylim = ax.get_ylim()    # after draw dos, get_ylim
    ### plot vertical lines
    if lvert:
        lstyle=['-', '--']
        for i, key in enumerate(lvert.keys()):
            label = f"{orb}-{key}"
            plt.axvline(x=lvert[key], color=getcolor_orbital(orb), linestyle=lstyle[i], label=label)
    
    ### ADD LEGEND
    plt.legend(loc=1)               # locate after plot
    #plt.xlim([-20.0, 15.0])
    #ax.set_xlim([-25, 15])
    ax.set_ylim([0.0, ylim[1]])     # set ymin=0 after get ymax
    if xlim:
        ax.set_xlim(xlim)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    plt.show()
    if Lsave:
        fig.savefig("pdos.png", dpi=150)
    return 0


def mplot_nvector_v1(x, y, dx=1.0, Title=None, Xtitle=None, Ytitle=None, Ylabels=None, Lsave=False, Colors=None, Ltwinx=None):
    '''
    call with x=[] and y=[ [...
    x:: [] or [size]
    y:: [size] or [[multi],[multi],...size]
    '''
    if Ltwinx:
        ax2 = ax.twinx()
        ax2.set_ylabel("Kinetic energy(kJ/mol)")
        ax2.tick_params(axis='y', labelcolor='g', labelsize=10)

    ys = np.array(y)
    if len(x) != 0:
        xsize = len(x)
    else:
        xsize = ys.shape[0]
        x=range(xsize)
    #print("hmm: {}".format(xsize))
    print(f"{x} :: {ys}")
    #plt.xticks(np.arange(min(x), max(x)+1, int(max(x)/dx)))
    #plt.xticks(np.arange(min(x), max(x)+1))
    #if tag=='x-sub':
    #    #xlabels = [item.get_text() for item in ax.get_xticklabels()]
    #    xlabels = tag_value
    #    ax.set_xticklabels(xlabels)

    plt.title(Title)
    if Xtitle:
        plt.xlabel(Xtitle, fontsize=15)
    plt.ylabel(Ytitle, fontsize=15)
    #ax.xaxis.set_major_locator(plt.NullLocator())
    #print(f"x, y shape:: {np.array(x).shape} {y.shape}")
    if ys.ndim == 1:
        plt.plot(x, y, 'bo-')
        #plt.scatter(x, y)
    elif ys.ndim == 2:
        for i in range(len(Ylabels)):
            #plt.plot(x,ys[i,:], 'o-', label=Ylabels[i] )
            plt.plot(x,ys[i,:], label=Ylabels[i] )
    else:
        print(f"Error:: obscure in y-dim {ys.ndim}")
    ### ADD LEGEND
    plt.legend(loc=2)                   # locate after plot
    common_figure_after()              # comment to remove legend box 
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
    return 0

def draw_histogram(y, nbin, Lsave, fname):
    n, bins, patches = plt.hist(y, nbin)
    if Lsave:
        plt.savefig(fname, dpi=150)
    plt.show()

    return 0

def _mplot_2f(f1, f2, Lsave, figname, Title, Xtitle, Ytitle):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
    
def _f_draw(fname, dp, Lsave, figname, Title, Xtitle, Ytitle):
    x1=[]
    y1=[]
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    with open(fname, 'r') as f:
        for line in f:
            #if re.search("[a-zA-Z]", line):
            #    continue
            xy=line.split()
            if not xy[0]:
                del xy[0]
            xvalue=round(float(xy[0]), dp)
            yvalue=round(float(xy[1]), dp)

            x1.append(xvalue)
            y1.append(yvalue)
    print(x1, y1)
    #draw_2d(x1, y1, Lsave, fname )
    plt.plot(x1, y1)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0

def mplot_vector_one(x, y, Title=None, Xtitle=None, Ytitle=None ,Lsave=False):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
plot_line = mplot_vector_one

### this does not make an error in serial plotting
def mplot_vector_two(x, y, Title=None, Xtitle=None, Ytitle=None ,Lsave=False):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    for i in range(len(x)):
        print(x[i], y[i])
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0


def _mplot_2c(x, y, Lsave, figname, Title, Xtitle, Ytitle):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0

def _mplot_3c(x, y, y2, figname, Title, Xtitle, Ytitle, Lsave):
    plt.xlabel(Xtitle)
    plt.ylabel(Ytitle)
    plt.title(Title)
    plt.plot(x, y )
    plt.plot(x, y2)
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0
def _mplot_2f3c(x, y1, y2,f1, f2, figname, Title, Xtitle, Ytitle, Lsave):
    #handles, labels = ax.get_legend_handels_labels()
    #ax.legend(handles, labels)
    plt.xlabel(Xtitle, fontsize=font_size)
    plt.ylabel(Ytitle, fontsize=font_size)
    plt.title(Title, fontsize=font_size)
    plt.plot(x, y1, label=f1)
    plt.plot(x, y2, label=f2)
    plt.legend()
    plt.show()
    if Lsave:
        plt.savefig(figname, dpi=150)
            
    return 0


def barplot2(hash, name, title):
    """Makes a barplot of the fingerprint about the O atom."""
    fp = descriptor.fingerprints[hash][0]
    fig, ax = pyplot.subplots()
    ax.bar(range(len(fp[1])), fp[1])
    ax.set_title(title)
    ax.set_ylim(0., 2.)
    ax.set_xlabel('fingerprint')
    ax.set_ylabel('value')
    fig.savefig(name)

def barplot_y(ys, name=None, xlabel=None, ylabel=None, title=None):
    """Makes a barplot of the fingerprint about the O atom."""
    fig, ax = plt.subplots()
    ax.bar(range(1,len(ys)+1), ys)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    if name:
        fig.savefig(name)
