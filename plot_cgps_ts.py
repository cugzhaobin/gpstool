# plot cgps time series.
# please download time series data from ftp://ftp.cgps.ac.cn/products/position/gamit/
# now we only support ????.ios_gamit_detrend.neu file
# By Zhao Bin, Institute of Seismology, CEA. Sep 23, 2017
# If you have any questions, please drop lines to cugzhaobin@163.com

# import libs
import sys
import numpy as np
import gpstime as gpst
import matplotlib.pyplot as plt
import os
from numpy import sin, cos, pi

# functions
def obs_min_linear(obs_decyr, obs, vel, offset, ref_decyr):
    '''
    raw time series - linear term
    '''
    for i in range(0, len(obs)):
        obs[i] = obs[i] - (vel*(obs_decyr[i]-ref_decyr) + offset)
    return obs

def mod_add_linear(obs_decyr, obs, vel, offset, ref_decyr):
    '''
    mod time series + linear term
    '''
    for i in range(0, len(obs)):
        obs[i] = obs[i] + (vel*(obs_decyr[i]-ref_decyr) + offset)
    return obs

def obs_min_period(obs_decyr, obs, a_sin, a_cos, sa_sin, sa_cos):
    '''
    raw time series - period term
    '''
    for i in range(0, len(obs)):
        obs[i] = obs[i] - a_sin*sin(2*pi*obs_decyr[i]) - a_cos*cos(2*pi*obs_decyr[i]) - sa_sin*sin(4*pi*obs_decyr[i]) - sa_cos*cos(4*pi*obs_decyr[i])
    return obs

def mod_add_period(obs_decyr, obs, a_sin, a_cos, sa_sin, sa_cos):
    '''
    mod time series + period term
    '''
    for i in range(0, len(obs)):
        obs[i] = obs[i] + a_sin*sin(2*pi*obs_decyr[i]) + a_cos*cos(2*pi*obs_decyr[i]) + sa_sin*sin(4*pi*obs_decyr[i]) + sa_cos*cos(4*pi*obs_decyr[i])
    return obs

def obs_min_break(obs_decyr, obs, brk_date, brk_offset):
    '''
    raw time series - break/eqoffset
    '''
    for i in range(0, len(obs)):
        if obs_decyr[i] >= brk_date:
            obs[i] = obs[i] - brk_offset
    return obs_decyr, obs

def mod_add_break(obs_decyr, obs, brk_date, brk_offset):
    '''
    mod time series + break/eqoffset
    '''
    for i in range(0, len(obs)):
        if obs_decyr[i] >= brk_date:
            obs[i] = obs[i] + brk_offset
    return obs_decyr, obs

def obs_min_eqlog(obs_decyr, obs, log_date, log_amp, log_tau):
    '''
    raw time series - eqlog term
    '''
    for i in range(0, len(obs)):
        if obs_decyr[i] >= log_date:
            obs[i] = obs[i] - log_amp*np.log(1+((obs_decyr[i]-log_date)*365.0/log_tau))
    return obs_decyr, obs

def mod_add_eqlog(obs_decyr, obs, log_date, log_amp, log_tau):
    '''
    mod time series + eqlog term
    '''
    for i in range(0, len(obs)):
        if obs_decyr[i] >= log_date:
            obs[i] = obs[i] + log_amp*np.log(1+((obs_decyr[i]-log_date)*365.0/log_tau))
    return obs_decyr, obs

def obs_min_eqexp(obs_decyr, obs, exp_date, exp_amp, exp_tau):
    '''
    raw time series - eqexp term
    '''
    for i in range(0, len(obs)):
        if obs_decyr[i] >= exp_date:
            obs[i] = obs[i] + exp_amp*np.exp(1-(-(obs_decyr[i]-exp_date)*365.0/exp_tau))
    return obs_decyr, obs

def mod_add_eqexp(obs_decyr, obs, exp_date, exp_amp, exp_tau):
    '''
    mod time series + eqexp term
    '''
    for i in range(0, len(obs)):
        if obs_decyr[i] >= exp_date:
            obs[i] = obs[i] + exp_amp*np.exp(1+(-(obs_decyr[i]-exp_date)*365.0/exp_tau))
    return obs_decyr, obs


def parse_resfile(resfile, plotopt):

    plotopt = str(plotopt)
    # check file exit
    if os.path.exists(resfile) == False:
        print ' '+resfile+" does not exist, please check!"
        sys.exit()  

    # read in the block file
    fid = open(resfile)

    # read all lines into content
    content = fid.readlines()

    # close the file
    fid.close()

    # init parameters
    # break term
    brk_date    = []
    brk_noffset = []
    brk_eoffset = []
    brk_uoffset = []

    # exp term
    exp_date    = []
    exp_n       = []
    exp_e       = []
    exp_u       = []
    exp_tau     = []

    # log term
    log_date    = []
    log_n       = []
    log_e       = []
    log_u       = []
    log_tau     = []

    # eqoffset term
    eq_date = []
    eq_noffset = []
    eq_eoffset = []
    eq_uoffset = []

    # period term
    n_sa_sin = 0.0
    e_sa_sin = 0.0
    u_sa_sin = 0.0
    n_sa_cos = 0.0
    e_sa_cos = 0.0
    u_sa_cos = 0.0
    n_a_sin = 0.0
    e_a_sin = 0.0
    u_a_sin = 0.0
    n_a_cos = 0.0
    e_a_cos = 0.0
    u_a_cos = 0.0

    # offsets
    noffset = 0.0
    eoffset = 0.0
    uoffset = 0.0

    # velocity
    nvel    = 0.0
    evel    = 0.0
    uvel    = 0.0

    # decyear
    obs_n_decyr = []
    obs_e_decyr = []
    obs_u_decyr = []

    mod_n_decyr = []
    mod_e_decyr = []
    mod_u_decyr = []

    # obs list
    obs_n = []
    obs_e = []
    obs_u = []

    # mod list
    mod_n = []
    mod_e = []
    mod_u = []

    # error list
    err_n = []
    err_e = []
    err_u = []

    # begin the job
    print ' Begin parse the file '+resfile
    for i in range(1, len(content)):
        # siteid
        if content[i][0:14] == "4-character ID":
            siteid = content[i].split(":")[1].strip()
        # first epoch
        if content[i][0:11] == "First Epoch":
            fepoch = content[i].split()[3]
        # last epoch
        if content[i][0:14] == "Last Epoch    ":
            eepoch = content[i].split()[3]
        # XYZ
        if content[i][0:3] == 'XYZ':
            xyz  = content[i].split()[4:7]
            xyz  = map(eval, xyz)
        # NEU
        if content[i][0:3] == 'NEU':
            neu  = content[i].split()[4:7]
            neu  = map(eval, neu)
        # intercept
        if content[i][0:8] == "Offsets ":
            offset  = content[i].split(":")[1]
            noffset =  float(offset.split()[0])
            eoffset =  float(offset.split()[3])
            uoffset =  float(offset.split()[6])
        # long term rate
        if content[i][0:6] == "Rates ":
            vel   =  content[i].split(":")[1]
            nvel  =  float(vel.split()[0])
            evel  =  float(vel.split()[3])
            uvel  =  float(vel.split()[6])
        # annual sesaonal term
        if content[i][0:14] == "Periodic : Cos" and content[i].split()[3] == "365.24" :
            n_a_cos = float(content[i].split()[5])
            e_a_cos = float(content[i].split()[8])
            u_a_cos = float(content[i].split()[11])
        if content[i][0:14] == "Periodic : Sin" and content[i].split()[3] == "365.24" :
            n_a_sin = float(content[i].split()[5])
            e_a_sin = float(content[i].split()[8])
            u_a_sin = float(content[i].split()[11])
        # semi-annual sesaonal term
        if content[i][0:14] == "Periodic : Cos" and content[i].split()[3] == "182.62" :
            n_sa_cos = float(content[i].split()[5])
            e_sa_cos = float(content[i].split()[8])
            u_sa_cos = float(content[i].split()[11])
        if content[i][0:14] == "Periodic : Sin" and content[i].split()[3] == "182.62" :
            n_sa_sin = float(content[i].split()[5])
            e_sa_sin = float(content[i].split()[8])
            u_sa_sin = float(content[i].split()[11])
        # earthquake offset
        if content[i][0:6] == "OffEq ":
            date = content[i].split()[2]
            year = int(date[0:4])
            month= int(date[4:6])
            day  = int(date[6:8])
            decyr = gpst.ymd_to_decyrs(year, month, day) 
            eq_date.append(decyr)
            eq_noffset.append(float(content[i].split()[4]))
            eq_eoffset.append(float(content[i].split()[7]))
            eq_uoffset.append(float(content[i].split()[10]))
        # postseismic 
        if content[i][0:6] == "EqLog ":
            date = content[i].split()[2]
            year = int(date[0:4])
            month= int(date[4:6])
            day  = int(date[6:8])
            decyr = gpst.ymd_to_decyrs(year, month, day) 
            log_date.append(decyr)
            log_n.append(float(content[i].split()[4]))
            log_e.append(float(content[i].split()[7]))
            log_u.append(float(content[i].split()[10]))
            log_tau.append(float(content[i].split()[16]))
        # break term
        if content[i][0:6] == "Break ":
            date = content[i].split()[2]
            year = int(date[0:4])
            month= int(date[4:6])
            day  = int(date[6:8])
            decyr = gpst.ymd_to_decyrs(year, month, day) 
            print decyr
            brk_date.append(decyr)
            brk_noffset.append(float(content[i].split()[4]))
            brk_eoffset.append(float(content[i].split()[7]))
            brk_uoffset.append(float(content[i].split()[10]))
        # postseismic term
        if content[i][0:3] == "Exp":
            date = content[i].split()[2]
            year = int(date[0:4])
            month= int(date[4:6])
            day  = int(date[6:8])
            decyr = gpst.ymd_to_decyrs(year, month, day) 
            exp_date.append(decyr)
            exp_n.append(float(content[i].split()[4]))
            exp_e.append(float(content[i].split()[7]))
            exp_u.append(float(content[i].split()[10]))
            exp_tau.append(float(content[i].split()[16]))
        # read the data
        if content[i][0] == " ":
            line = content[i].split()
            # north
            if float(line[9]) < 5:
                obs_n_decyr.append(float(line[2]))
                obs_n.append(float(line[4]))
                err_n.append(float(line[8]))
            # east
            if float(line[10]) < 5:
                obs_e_decyr.append(float(line[2]))
                obs_e.append(float(line[5]))
                err_e.append(float(line[10]))
            # veritcal
            if float(line[12]) < 15:
                obs_u_decyr.append(float(line[2]))
                obs_u.append(float(line[6]))
                err_u.append(float(line[12]))
    #
    print ' Finished parse file ...'


    # get start dec year
    fyr = int(fepoch[0:4])
    fmo = int(fepoch[4:6])
    fdy = int(fepoch[6:8])
    start_decyr = gpst.ymd_to_decyrs(fyr, fmo, fdy)

    # get stop decimal year
    eyr = int(eepoch[0:4])
    emo = int(eepoch[4:6])
    edy = int(eepoch[6:8])
    stop_decyr = gpst.ymd_to_decyrs(eyr, emo, edy)

    # compute reference epoch
    ref_decyr = 0.5*(start_decyr+stop_decyr)

    # get the model data
    print ' create obs-data and mod-data using plot option: '+plotopt

    # init the model data
    for decyr in np.arange(start_decyr, stop_decyr+0.002739726, 0.002739726):
        mod_n_decyr.append(decyr)
        mod_n.append(0.0)
        mod_e_decyr.append(decyr)
        mod_e.append(0.0)
        mod_u_decyr.append(decyr)
        mod_u.append(0.0)


    # linear term
    print ' work for long term ...'
    if plotopt[4] == "1":
        obs_n = obs_min_linear(obs_n_decyr, obs_n, nvel, noffset, ref_decyr)
        obs_e = obs_min_linear(obs_e_decyr, obs_e, evel, eoffset, ref_decyr)
        obs_u = obs_min_linear(obs_u_decyr, obs_u, uvel, uoffset, ref_decyr)
    elif plotopt[4] == "0":
        mod_n = mod_add_linear(mod_n_decyr, mod_n, nvel, noffset, ref_decyr)
        mod_e = mod_add_linear(mod_e_decyr, mod_e, evel, eoffset, ref_decyr)
        mod_u = mod_add_linear(mod_u_decyr, mod_u, uvel, uoffset, ref_decyr)

    # period term
    print ' work for period term ...'
    if plotopt[0] == "1":
        obs_n = obs_min_period(obs_n_decyr, obs_n, n_a_sin, n_a_cos, n_sa_sin, n_sa_cos)
        obs_e = obs_min_period(obs_e_decyr, obs_e, e_a_sin, e_a_cos, e_sa_sin, e_sa_cos)
        obs_u = obs_min_period(obs_u_decyr, obs_u, u_a_sin, u_a_cos, u_sa_sin, u_sa_cos)
    if plotopt[0] == "0":
        mod_n = mod_add_period(mod_n_decyr, mod_n, n_a_sin, n_a_cos, n_sa_sin, n_sa_cos)
        mod_e = mod_add_period(mod_e_decyr, mod_e, e_a_sin, e_a_cos, e_sa_sin, e_sa_cos)
        mod_u = mod_add_period(mod_u_decyr, mod_u, u_a_sin, u_a_cos, u_sa_sin, u_sa_cos)

    # break term
    print ' work for break term ...'
    if plotopt[3] == "1":
        for i in range(0, len(brk_date)):
            obs_n_decyr, obs_n = obs_min_break(obs_n_decyr, obs_n, brk_date[i], brk_noffset[i])
            obs_e_decyr, obs_e = obs_min_break(obs_e_decyr, obs_e, brk_date[i], brk_eoffset[i])
            obs_u_decyr, obs_u = obs_min_break(obs_u_decyr, obs_u, brk_date[i], brk_uoffset[i])
    if plotopt[3] == "0":
        for i in range(0, len(brk_date)):
            mod_n_decyr, mod_n = mod_add_break(mod_n_decyr, mod_n, brk_date[i], brk_noffset[i])
            mod_e_decyr, mod_e = mod_add_break(mod_e_decyr, mod_e, brk_date[i], brk_eoffset[i])
            mod_u_decyr, mod_u = mod_add_break(mod_u_decyr, mod_u, brk_date[i], brk_uoffset[i])


    # eqoffset term
    print ' work for earthquake offset term'
    if plotopt[2] == "1":
        for i in range(0, len(eq_date)):
            obs_n_decyr, obs_n = obs_min_break(obs_n_decyr, obs_n, eq_date[i], eq_noffset[i])
            obs_e_decyr, obs_e = obs_min_break(obs_e_decyr, obs_e, eq_date[i], eq_eoffset[i])
            obs_u_decyr, obs_u = obs_min_break(obs_u_decyr, obs_u, eq_date[i], eq_uoffset[i])
    if plotopt[2] == "0":
        for i in range(0, len(eq_date)):
            mod_n_decyr, mod_n = mod_add_break(mod_n_decyr, mod_n, eq_date[i], eq_noffset[i])
            mod_e_decyr, mod_e = mod_add_break(mod_e_decyr, mod_e, eq_date[i], eq_eoffset[i])
            mod_u_decyr, mod_u = mod_add_break(mod_u_decyr, mod_u, eq_date[i], eq_uoffset[i])

    # eq log term
    print ' work for earthquake postseismic term using log function'
    if  plotopt[1] == "1":
        for i in range(0, len(log_date)):
            obs_n_decyr, obs_n = obs_min_eqlog(obs_n_decyr, obs_n, log_date[i], log_n[i], log_tau[i])
            obs_e_decyr, obs_e = obs_min_eqlog(obs_e_decyr, obs_e, log_date[i], log_e[i], log_tau[i])
            obs_u_decyr, obs_u = obs_min_eqlog(obs_u_decyr, obs_u, log_date[i], log_u[i], log_tau[i])
    if  plotopt[1] == "0":
        for i in range(0, len(log_date)):
            mod_n_decyr, mod_n = mod_add_eqlog(mod_n_decyr, mod_n, log_date[i], log_n[i], log_tau[i])
            mod_e_decyr, mod_e = mod_add_eqlog(mod_e_decyr, mod_e, log_date[i], log_e[i], log_tau[i])
            mod_u_decyr, mod_u = mod_add_eqlog(mod_u_decyr, mod_u, log_date[i], log_u[i], log_tau[i])


    # eq exp term
    print ' work for earthquake postseismic term using exp function'
    if  plotopt[1] == "1":
        for i in range(0, len(exp_date)):
            obs_n_decyr, obs_n = obs_min_eqexp(obs_n_decyr, obs_n, exp_date[i], exp_n, exp_tau)
            obs_e_decyr, obs_e = obs_min_eqexp(obs_e_decyr, obs_e, exp_date[i], exp_e, exp_tau)
            obs_u_decyr, obs_u = obs_min_eqexp(obs_u_decyr, obs_u, exp_date[i], exp_u, exp_tau)
    if  plotopt[1] == "0":
        for i in range(0, len(exp_date)):
            mod_n_decyr, mod_n = mod_add_eqexp(mod_n_decyr, mod_n, exp_date[i], exp_n, exp_tau)
            mod_e_decyr, mod_e = mod_add_eqexp(mod_e_decyr, mod_e, exp_date[i], exp_e, exp_tau)
            mod_u_decyr, mod_u = mod_add_eqexp(mod_u_decyr, mod_u, exp_date[i], exp_u, exp_tau)

    # output data
    # observation
    outmatrix   = np.vstack((obs_n_decyr, obs_n, err_n)).T
    np.savetxt('t.no', outmatrix)

    outmatrix   = np.vstack((obs_e_decyr, obs_e, err_e)).T
    np.savetxt('t.eo', outmatrix)

    outmatrix   = np.vstack((obs_u_decyr, obs_u, err_u)).T
    np.savetxt('t.uo', outmatrix)

    # model
    outmatrix   = np.vstack((mod_n_decyr, mod_n)).T
    np.savetxt('t.nm', outmatrix)

    outmatrix   = np.vstack((mod_e_decyr, mod_e)).T
    np.savetxt('t.em', outmatrix)

    outmatrix   = np.vstack((mod_u_decyr, mod_u)).T
    np.savetxt('t.um', outmatrix)

    # breaks
    if len(brk_date) > 0:
        fid = open('brk.n', 'w')
        for i in range(0,len(brk_date)):
            if np.abs(brk_noffset[i]) > 1:
                line = "%10.4f%10.2f\n" %(brk_date[i], -99999)
                fid.write(line)
                line = "%10.4f%10.2f\n" %(brk_date[i],  99999)
                fid.write(line)
                fid.write(">\n")
        fid = open('brk.e', 'w')
        for i in range(0,len(brk_date)):
            if np.abs(brk_eoffset[i]) > 1:
                line = "%10.4f%10.2f\n" %(brk_date[i], -99999)
                fid.write(line)
                line = "%10.4f%10.2f\n" %(brk_date[i],  99999)
                fid.write(line)
                fid.write(">\n")
        fid = open('brk.u', 'w')
        for i in range(0,len(brk_date)):
            if np.abs(brk_uoffset[i]) > 1:
                line = "%10.4f%10.2f\n" %(brk_date[i], -99999)
                fid.write(line)
                line = "%10.4f%10.2f\n" %(brk_date[i],  99999)
                fid.write(line)
                fid.write(">\n")

    # eqoffset
    if len(eq_date) > 0:
        fid = open('eqs.n', 'w')
        for i in range(0,len(eq_date)):
            if np.abs(eq_noffset[i]) > 1:
                line = "%10.4f%10.2f\n" %(eq_date[i], -99999)
                fid.write(line)
                line = "%10.4f%10.2f\n" %(eq_date[i],  99999)
                fid.write(line)
                fid.write(">\n")
        fid = open('eqs.e', 'w')
        for i in range(0,len(eq_date)):
            if np.abs(eq_eoffset[i]) > 1:
                line = "%10.4f%10.2f\n" %(eq_date[i], -99999)
                fid.write(line)
                line = "%10.4f%10.2f\n" %(eq_date[i],  99999)
                fid.write(line)
                fid.write(">\n")
        fid = open('eqs.u', 'w')
        for i in range(0,len(eq_date)):
            if np.abs(eq_uoffset[i]) > 1:
                line = "%10.4f%10.2f\n" %(eq_date[i], -99999)
                fid.write(line)
                line = "%10.4f%10.2f\n" %(eq_date[i],  99999)
                fid.write(line)
                fid.write(">\n")


    res_dict = {'site': siteid, 'eqdt': eq_date, 'brkdt': brk_date}
    return res_dict

if __name__ == '__main__':
##############################################################################
# main program
##############################################################################
    if len(sys.argv) != 3:
        print ' Usage: plot_pbo_res.py <residual file> <10111>'
        print ' option:'
        print '  1       1           1         1      1'
        print '  period  postseismic eqoffset  break  linear'

        sys.exit()

    # get the input file name
    resfile = sys.argv[1]
    plotopt = str(sys.argv[2])
    resdict = parse_resfile(resfile, plotopt)
    site    = resdict['site']

    # read time series and plot figure

    plt.figure(figsize=(10,13))
    # north component
    nobs = np.genfromtxt('t.no')
    nmod = np.genfromtxt('t.nm')
    plt.subplot(3,1,1)
    plt.scatter(nobs[:,0], nobs[:,1], s=2, color='r')
    plt.plot(nmod[:,0], nmod[:,1])
    plt.title('Displacement Time series at '+site)
    for dt in resdict['eqdt']:
        if dt > np.min(nobs[:,0]) and dt < np.max(nobs[:,0]):
            plt.plot([dt,dt], [-1000,1000])
    for dt in resdict['brkdt']:
        if dt > np.min(nobs[:,0]) and dt < np.max(nobs[:,0]):
            plt.plot([dt,dt], [-1000,1000])
    plt.xlim([np.min(nmod[:,0]), np.max(nmod[:,0])])
    if plotopt[4] == "0":
        plt.ylim([np.min(nmod[:,1])*1.5, np.max(nmod[:,1])*1.5])
    else:
        delta = np.max(nmod[:,1])-np.min(nmod[:,1])
        if delta == 0:
            delta = 20
        plt.ylim([np.min(nmod[:,1])-delta, np.max(nmod[:,1])+delta])
    plt.xlabel('Time (year)')
    plt.ylabel('North displacement (mm)')



    # east  component
    eobs = np.genfromtxt('t.eo')
    emod = np.genfromtxt('t.em')
    plt.subplot(3,1,2)
    plt.scatter(eobs[:,0], eobs[:,1], s=2, color='g')
    plt.plot(emod[:,0], emod[:,1])
    for dt in resdict['eqdt']:
        if dt > np.min(eobs[:,0]) and dt < np.max(eobs[:,0]):
            plt.plot([dt,dt], [-1000,1000])
    for dt in resdict['brkdt']:
        if dt > np.min(eobs[:,0]) and dt < np.max(eobs[:,0]):
            plt.plot([dt,dt], [-1000,1000])
    plt.xlim([np.min(nmod[:,0]), np.max(nmod[:,0])])
    if plotopt[4] == "0":
        plt.ylim([np.min(emod[:,1])*1.5, np.max(emod[:,1])*1.5])
    else:
        delta = np.max(emod[:,1])-np.min(emod[:,1])
        if delta == 0:
            delta = 20
        plt.ylim([np.min(emod[:,1])-delta, np.max(emod[:,1])+delta])

    plt.xlabel('Time (year)')
    plt.ylabel('East displacement (mm)')

    # up  component
    uobs = np.genfromtxt('t.uo')
    umod = np.genfromtxt('t.um')
    plt.subplot(3,1,3)
    plt.scatter(uobs[:,0], uobs[:,1], s=2, color='b')
    plt.plot(umod[:,0], umod[:,1])
    for dt in resdict['eqdt']:
        if dt > np.min(uobs[:,0]) and dt < np.max(uobs[:,0]):
            plt.plot([dt,dt], [-1000,1000])
    for dt in resdict['brkdt']:
        if dt > np.min(uobs[:,0]) and dt < np.max(uobs[:,0]):
            plt.plot([dt,dt], [-1000,1000])
    plt.xlim([np.min(nmod[:,0]), np.max(nmod[:,0])])
    if plotopt[4] == "0":
        plt.ylim([np.min(umod[:,1])*3.0, np.max(umod[:,1])*3.0])
    else:
        delta = np.max(umod[:,1])-np.min(umod[:,1])
        if delta == 0:
            delta = 10
        plt.ylim([np.min(umod[:,1])-delta*3, np.max(umod[:,1])+delta*3])
    plt.xlabel('Time (year)')
    plt.ylabel('Up displacement (mm)')

    outfile = resfile.split('.')[0]+'_ts.pdf'
    plt.savefig(outfile)
