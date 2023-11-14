### Inside doscar directory, run
../pdoscar.py -j pdanpl
    -j  pdanpl  (default)
        pd      pdos calculation
        an      analysis such as dosmax, doscenter, fermi_abundance
        pl      plot option
    necessary option:
        pd  [pdos calculation]
            -a      atom list
        an  [orbital analysis]
            -l      angular momentum for analysis
        pl  [plot option]
            -po     plot option: designate orbital to plot [default: 'a',all]
    E.G.::
        pdos calculation
            python vas_pdos.py -j pd    -a 43
        orbital analysis
            python vas_pdos.py -j pdan  -a 43   -l p
            python vas_pdos.py -j an            -l p
        plot orbital
            python vas_pdos.py -j pdanpl -a 43  -l p -po a
            python vas_pdos.py -j anpl          -l p -po a
            python vas_pdos.py -j pl                 -po a

### Explanation for each job ['pd'(pdos), 'an'(orbital analysis), 'pl'(plot)]
(1) Abtain pdos:
	    E.g.:
		python pdoscar.py -a 43
		python pdoscar.py -j plods -i DOSCAR -a atomlist -o split
	    Options for pdos
		-j  pd        calculate pdos from DOSCAR
		    an        analysis for the previous cal_pdos by reading 'SUM_ATOM.dat'
		    pl        plot with options
		-i  'DOSCAR' (default)
		-a  atom list index starts from 1
		-os spin option for spin = 2
		    None      up and down spin summed: format(ene, s, p, d, s+p+d)
		    split     separate up/down spin: format(ene, s_up, s_down, ... 9 colm)
		    polar     down spin has - value
(2)	Analysis of PDOS
	    E.g.:
		pdoscar.py -j anal -l 1
	    Options for analysis
		-j  an        to get l-center, DOSmax, Fermi abundance
		-p  PDOS file name 'SUM_ATOM.dat' (default)
		-l  angular quantum number for analysis
		-oa options for analysis
		    max       energy for max DOS of l-orbital
		    center    average E of l-orbital
		    fa        E_fermi abundance
(3)	Plot orbital
	    E.g.:
		pdoscar.py -j anal -l p -po t s p d
	    Options for plot
		-j  pl
		-po plot orbital: t s p d, t in first for clean line

