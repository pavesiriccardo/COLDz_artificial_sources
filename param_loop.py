#[noisefunc((256,256),30.968+(chn-1)*chanwidth,np.sqrt(2)*noisedic[chn/2]) for chn in range(0,4034)]


#9e-5*3/5.65

#1.e-4
#######
#if need to wait:
import subprocess,os,time,numpy as np


def check_pid(pid):        
    """ Check For the existence of a unix pid. """
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True



#while check_pid(12366):
#	time.sleep(60)


#######

import sys,subprocess,os

processes=[]
#!ulimit -n 
current=0
flux_list_mJy=np.array([.3,1.58571429,0.62142857,1.26428571,0.94285714,0.51428571,1.8,1.05,1.47857143,0.72857143,1.37142857,0.40714286,1.15714286,1.69285714,0.83571429,0.2,2.0,0.1,2.2,3.2,2.7,3.7])


for idx1,flux in enumerate(1e-3*flux_list_mJy):
	for idx2,std_spat in enumerate([0.85 ,  2.12,  3.82]):    #corresponds to 1,2.5,4.5 arcsec real source size
		for idx3,std_freq in enumerate([4.95/2.,4.95,4.95*3./2.]):
			if  (idx1*9+idx2*3+idx3)>-1 and ('art'+str(idx1)+str(idx2)+str(idx3) not in os.listdir('.')):# and (idx1*9+idx2*3+idx3)%3==2:
				#P=subprocess.call(['nohup','python','analyze.py','art'+str(idx1)+str(idx2)+str(idx3),str(flux),str(std_spat),str(std_freq),'&'])
				#sys.argv=['','art'+str(idx1)+str(idx2)+str(idx3),flux,std_spat,std_freq]
				subprocess.call(['nohup','python', 'completeArtificial.py','art'+str(idx1)+str(idx2)+str(idx3),str(flux),str(std_spat),str(std_freq),'&'])
				#execfile('completeArtificial.py')



print 'finished all'


#first run inject from python
#then use:
#nohup /usr/local/casapy-stable-42.0.26945-001-64b/bin/casapy -c param_loop.py &
#to execfile completeArtificial
#last run analyze from python
