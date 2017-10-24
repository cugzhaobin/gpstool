# Written by He Kefeng, Institute of Seismology, CEA. 2017-06-09
# Modified by Zhao Bin, Institute of Seismology, CEA. 2017-09-26

# import the libs
import re
import sys
import numpy as np
from numpy import sin, cos, log, exp, pi
import matplotlib.pyplot as plt
import coord_time as cts
from scipy.optimize import curve_fit


def main(baselinefile = '', brksfile='', fitflag=(1,1,1)):
    '''
    Input:
        baselinefile:  baseline time series file
        brksfile    :  breaks in time series
        fitflag     :  flag for plotting figure
                       (1,1,1) means long-term velocity, period and offsets
    '''
    # read baseline file
    baseline = ReadBaseline(baselinefile)
    # get file name
    pattern = re.compile(r'(?:.*\\|.*/|^)(\w{9})(?=\.)')
    # get station-station name
    sitename = pattern.match(baselinefile).group(1) 
    # get dirname
    dirname = baselinefile[:baselinefile.find(sitename)]
    #
    dirname = './' if dirname == '' else dirname 
#   print(sitename, dirname)
    brkinfo = Readbrks(brksfile, sitename)
    onebl = OnePosfitinfo(baseline, brkinfo)
    fig = onebl.posfitting(name=sitename, fitflag=fitflag)
    fig.savefig(dirname+sitename+'.png')
    plt.show(fig)
    return 

class OnePosfitinfo:
    '''
    # Class of fitting time series
    '''
    def __init__(self, baseline, brkinfo):
        self.baseline = baseline
        self.brkinfo = brkinfo
        
    def Getpmenu(self):
        '''
        obtain parameter list
        '''
        # make sure the breaks time is unique
        t = self.t
        pmenu = ['contant', 'velocity'] + ['cycle']*4
        if len(self.brkinfo) != 0:
            pmenu.extend(['break']*len(self.brkinfo))
        self.pmenu = pmenu 
        return pmenu

    def Setpoptinit(self, pmenu):
        '''
        init bounds for parameters
        '''
        length = len(pmenu)
        pinit = [1]*length
        return pinit 
    
    def Setpoptbound(self, pmenu):
        '''
        setting bounds for parameters
        '''
        length = len(pmenu)
        plower = [-np.inf]*length
        pupper = [np.inf]*length
        pbound = [tuple(plower), tuple(pupper)]
        return pbound 
    
    @property
    def t(self):
        #convert to decimal year
#       YMD = [[int(tt[0:4]), int(tt[4:6]), int(tt[6:8])]for tt in self.baseline['YMD']]
#       Time = cts.YMD2DotY(YMD)
        Time = self.baseline[:,0]
        return Time
       
    def Getstepindex(self,pmenu):
        '''
        set index of breaks
        '''
        t = self.t
        stepindex = []
        for i,name in enumerate(pmenu):
            if name == 'break':
                stepindex.append(i)
        return stepindex
        
    def  full_filter(self):
        t = self.t
        brks = self.brkinfo
        t0 = np.mean(t)
        def pos_filter(t, *p):
            '''
            fitting time series
            '''
            y, index = 0, 0
            y += Line(t, t0, *p[0:2])
            index += 2
            y += Cycle(t, *p[index:4+index])
            index += 4
            # if we have breaks
            if len(brks) != 0:
                y += Breaks(t, brks, p[index:])
            return y
        return pos_filter   
    
    def Getfitpopt(self):
        '''
        get poptE, perrE, poptN, perrN
        '''
        data = self.baseline        
        # get time
        t = self.t
        # 
        pmenu = self.Getpmenu()
        pinit = self.Setpoptinit(pmenu)
        pbounds = self.Setpoptbound(pmenu)
        popt, pcov = curve_fit(self.full_filter(), data[:,0], data[:,1], p0=pinit, bounds=pbounds, maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        self.popt, self.perr = popt, perr
        return popt, perr
    
    def Showfig(self, fig):
        '''
        Show figure
        '''
        fig.show()
        return 
    
    def GetImg(self, ax0, name):
        '''
        ax0:t, posdata, simt, simE
        ax1:t, posdata, simt, simN
        EQname: earthquke name
        '''
        t, pos, simt, sim = ax0 
        delta = pos[0]
        pos = (pos-delta)*1e3
        sim = (sim-delta)*1e3
        fig, axes = plt.subplots(nrows=1, ncols=1)
        axes.plot(simt, sim, color='r', label='fitted curve')
        delta = np.max(sim)-np.min(sim)
        axes.set_ylim([np.min(sim)-delta, np.max(sim)+delta])
        axes.set_xlim([np.min(simt), np.max(simt)])
        axes.scatter(t, pos, label='data')
        axes.set_xlabel('Time (year)')
        axes.set_ylabel('Length (mm)')
        axes.legend(framealpha=1, shadow=True, loc='best')
        axes.get_xaxis().get_major_formatter().set_useOffset(False)
        fig.suptitle(name, fontsize='large', y=0.98)
        return fig
        
    def posfitting(self, fitflag=[1, 1, 1], name = ''):
        '''
        [liner, cycle, break]
        eqname: for coseismic and postseismic 
       
        '''
        # obtain time
        t = self.t
        t_mean = np.mean(t)        
        pmenu = self.Getpmenu()
        fitting = self.full_filter()
        popt, perr = self.Getfitpopt()
        stepindex = self.Getstepindex(pmenu)
        # fitting time list
        tend = t[-1]+0.05
        simt = np.arange(t[0]-0.05, tend, 0.05)
        sim = fitting(simt, *popt)        

        allbrk = self.brkinfo
        allbrk_seq = sorted(zip(allbrk, range(len(allbrk))), key = lambda i: i[0], reverse=True)
        allbrk = [i[0] for i in allbrk_seq]
        brk_seq = [i[1] for i in allbrk_seq]
        insertseq = []
        vel = popt[1]        
        #obtain break time
        #calcualte model results
        for i in allbrk:
            insertseq.append(len(simt[simt<=i]))
        #
        pos = self.baseline[:,1]
        pmenu = np.asarray(pmenu)
        pflagsim, pflagpos = np.zeros(pmenu.shape), np.ones(pmenu.shape)
        if fitflag[1] == 1:
            pflagsim = pflagsim + (pmenu == 'cycle')
            pflagpos = pflagpos - (pmenu == 'cycle')
            
        temp = 0
        if fitflag[2] == 1:
            pflagsim = pflagsim + (pmenu == 'break')
            pflagpos = pflagpos - (pmenu == 'break')
            temp = 1
        if fitflag[0] == 1:
            pflagsim = pflagsim + (pmenu == 'contant') + (pmenu == 'velocity')
            pflagpos = pflagpos - (pmenu == 'contant') - (pmenu == 'velocity') 
            pos = pos - fitting(t, *(popt*pflagpos))
            sim = fitting(simt, *(popt*pflagsim))                
        else:
            sim = fitting(simt, *(popt*pflagsim))
            pos = pos - fitting(t, *(popt*pflagpos)) 
            vel = 0
        # save data with no insert values
        sim, simt =sim.tolist(), simt.tolist()         
        for i in range(len(insertseq)):
            if allbrk[i] > simt[0]:
                sim.insert(insertseq[i], sim[insertseq[i]-1] + vel*(allbrk[i] - simt[insertseq[i]-1]))
                simt.insert(insertseq[i], allbrk[i])
                sim.insert(insertseq[i]+1, sim[insertseq[i]] + temp*popt[stepindex[brk_seq[i]]])
                simt.insert(insertseq[i]+1, allbrk[i])
        # plot figure
        self.simt, self.sim = np.asarray(simt), np.asarray(sim)
        fig = self.GetImg([t, pos, self.simt, self.sim], name)
        return fig

    
        
#########################################################################
###               define functions for linear fitting terms           ###
#########################################################################

def Break(t, t0, param):
    '''
    define single break
    '''
    return param * np.heaviside(t - t0, 0)

def Breaks(t, Brks, params):
    '''
    define break terms
    '''
    temp = 0
    for i in range(len(Brks)):
        temp += Break(t, Brks[i], params[i])
    return temp

def Line(t, t0, a, b):
    '''
    define linear term
    '''
    return a + b*(t - t0)


def Cycle(t, c, d, e, f):
    '''
    define periodic term
    '''
    return (c*sin(2*pi*t) + d*cos(2*pi*t) + e*sin(4*pi*t) +
            f*cos(4*pi*t))


#########################################################################
###                   define reading break files                      ###
#########################################################################
# 

def Readbrks(brksfile, sitename):
    '''
    read eq_rename file, 
    '''
    brks = []
    name = ''
    try:
        with open(brksfile, 'rt') as fin:
            for Line in fin:
                if re.match(r'^\s*$', Line): continue
                if Line[0] == '#': continue
                if Line[0] == '*': continue
                if Line[0] == ' ':
                    if Line.find('break') != -1:
                        Linelist = Line.split()
                        if Linelist[1] == sitename:
                            if len(Linelist) == 7:
                                brks.append(cts.YMD2DotY([int(i) for i in Linelist[2:5]], [int(i) for i in Linelist[5:]]+[0])[0])
                            elif len(Linelist) == 3:
                                brks.append(float(Linelist[2]))
#       print brks
        return brks
    except IOError as ioerr:
        print("FIle error: " + str(ioerr))

def ReadBaseline(basefile):
    '''
    read baseline files
    '''
#   data = np.genfromtxt(basefile, names=('YMD', 'DiffY', 'Val'), dtype='S8,f8,f8')
    data = np.genfromtxt(basefile)
    return data
    
    
if __name__ == "__main__": 
    if len(sys.argv) >= 2:
        name = sys.argv[1]
        main(name, brksfile='brkfile', fitflag=(1,1,0))
    else:
#       main('AHAQ_JSLS.txt', brksfile='brkfile', fitflag=(1,1,0))
        print sys.argv[0] + ' <baseline> [brkfile]'
