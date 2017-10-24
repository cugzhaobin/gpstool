# plot time series of date processing quality assurance file
# Written by Zhao Bin, Institute of Seismology, CEA. Sep 23, 2017.
# If you have any questions please drop lines to cugzhaobin@163.com
import sys, datetime
import numpy as np
import matplotlib.pylab as plt

def read_qa(qafile):
     data=np.genfromtxt(qafile, names=('YMD', 'MP1', 'MP2', 'RMS', 'A', 'B','DECYR'), dtype='S8,f8,f8,f8,f8,f8,f8', usecols=[12,5,6,9,10,11,14])
     return data

def YMD2datetime(strYMD):
    '''
    '''
    fmtdt = []
    for i in range(len(strYMD)):
        year  = int(strYMD[i][0:4])
        month = int(strYMD[i][4:6])
        day   = int(strYMD[i][6:8])
        dt    = datetime.date(year, month, day)
        fmtdt.append(dt)
    return fmtdt

def zero2nan(dat):
    '''
    '''
    for i in range(len(dat)):
        if dat[i] == 0.0:
            dat[i] = 'NaN'
    return dat


if __name__ == "__main__": 
    if len(sys.argv) >= 2:
        qafile = sys.argv[1]
        data   = read_qa(qafile)
        strYMD = data['YMD']
        decyr  = data['DECYR']
        site   = qafile.split('.')[0]
        pdfile = site+'_qa.pdf'
        
        
        mp1    = zero2nan(data['MP1'])
        mp2    = zero2nan(data['MP2'])
        rms    = zero2nan(data['RMS'])
        a      = zero2nan(data['A'])
        b      = zero2nan(data['B'])

        plt.figure(figsize=(10,13))
#        plt.subplot(3,2,1)
#        plt.title(site+': Processing Quality Assurance')
#        plt.scatter(decyr, mp1, s=2, c='r', linewidths=0)
#        plt.xlabel('Date (year)')
#        plt.ylabel('MP1 (m)')
#        plt.subplot(3,2,2)
#        plt.scatter(decyr, mp2, s=2, c='g', linewidths=0)
#        plt.xlabel('Date (year)')
#        plt.ylabel('MP2 (m)')
#        plt.subplot(3,2,3)
#        plt.scatter(decyr, a, s=2, linewidths=0)
#        plt.xlabel('Date (year)')
#        plt.ylabel('A (mm)')
#        plt.subplot(3,2,4)
#        plt.scatter(decyr, b, s=2, linewidths=0)
#        plt.xlabel('Date (year)')
#        plt.ylabel('B (mm)')
#        ax = plt.subplot(3,1,3)
#        ax.get_xaxis().get_major_formatter().set_useOffset(False)
#        plt.scatter(decyr, rms, s=2, linewidths=0)
#        plt.xlabel('Date (year)')
#        plt.ylabel('RMS (mm)')        
        ax = plt.subplot(5,1,1)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.scatter(decyr, mp1, s=2, linewidths=0)
        plt.xlabel('Date (year)')
        plt.ylabel('MP1 (m)')
        plt.ylim([0,1])
        plt.title('Data Processing Quality for '+site)

        ax = plt.subplot(5,1,2)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.scatter(decyr, mp2, s=2, linewidths=0)
        plt.xlabel('Date (year)')
        plt.ylabel('MP2 (m)')
        plt.ylim([0,1])

        ax = plt.subplot(5,1,3)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.scatter(decyr, a, s=2, linewidths=0)
        plt.xlabel('Date (year)')
        plt.ylabel('A (mm)')

        ax = plt.subplot(5,1,4)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.scatter(decyr, b, s=2, linewidths=0)
        plt.xlabel('Date (year)')
        plt.ylabel('B (mm)')

        ax = plt.subplot(5,1,5)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.scatter(decyr, rms, s=2, linewidths=0)
        plt.xlabel('Date (year)')
        plt.ylabel('RMS (mm)')
        
        plt.subplots_adjust(hspace=0.6, wspace=0.6)
        plt.savefig(pdfile, format='pdf')
        
        
    else:
        print(sys.argv[0] + ' <qa file>')
        print(' Please download from ftp://ftp.cgps.ac.cn/by_site/')
        sys.exit()
        
