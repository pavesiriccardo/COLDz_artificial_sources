from inj_source_class import inj_source

objs=[]
for flux in range(22):
	for spat in range(3):
		for freq in range(3):
			os.chdir('art'+str(flux)+str(spat)+str(freq))
			inp=open('analysis_injected_objects.dat')
			objs+=cPickle.load(inp)
			inp.close()
			os.chdir('..')




#injected point of view

comple=dict()
comple2=dict()

Iflux_bins=np.arange(0,2.5,.01)

for spat_bin in [0,1,2]:
        for freq_bin in [0,1,2]:
                comple[(spat_bin,freq_bin)]=np.zeros(Iflux_bins.shape[0])

tot=dict()
for spat_bin in [0,1,2]:
        for freq_bin in [0,1,2]:
                tot[(spat_bin,freq_bin)]=np.zeros(Iflux_bins.shape[0])



inj_spat_list=[2.00175,4.9926,  8.9961]
inj_freq_list=[5.828625, 11.65725,17.485875]
SNR_thresh=5.
for obj in objs:
	idx2=inj_spat_list.index(obj.inj_spat_fwhm)
	idx3=inj_freq_list.index(obj.inj_freq_fwhm)
	flux=obj.inj_flux
	if obj.SNR>SNR_thresh:
		comple[(idx2,idx3)][np.digitize(flux*1.064*obj.inj_freq_fwhm*.004/34.*3e5,Iflux_bins)-1]+=1  #(std_spat/2.6072)**2*
	tot[(idx2,idx3)][np.digitize(flux*1.064*obj.inj_freq_fwhm*.004/34.*3e5,Iflux_bins)-1]+=1 




#FWHM of the line in km/s std_freq*2.355*.004/34.*3e5
#integrated line is is flux*2.5066*std_freq*.004/34.*3e5
from scipy.optimize import minimize
from scipy.special import erf
myfit=lambda f,d,f0: 1-(   1./(f+d)*np.exp(-f/f0)   )
#myfit=lambda f,k,f0p:erf((f-f0p)*k)

chi2=lambda (d,f0),f,dat: np.nansum((myfit(f,d,f0)-dat)**2)
#chi2=lambda (k,f0p),f,dat: np.nansum((myfit(f,k,f0p)-dat)**2+f0p**2/5)

all_params=[]
labels=['spat=0,freq=0','freq=1','freq=2','spat=1','','','spat=2','','']

for spat_bin in [0,1,2]:
        for freq_bin in [0,1,2]:
                comple2[(spat_bin,freq_bin)]=np.array(comple[(spat_bin,freq_bin)])/np.array(tot[(spat_bin,freq_bin)])
                plt.scatter(Iflux_bins,comple2[(spat_bin,freq_bin)],color=['r','g','b'][freq_bin],marker=['o','s','^'][spat_bin])
		large_enough=comple2[(spat_bin,freq_bin)]>.1
                params=minimize(chi2,[1,.05],args=(Iflux_bins[large_enough],comple2[(spat_bin,freq_bin)][large_enough]),method='Nelder-Mead').x  #,method='Nelder-Mead'
                plt.plot(Iflux_bins,myfit(Iflux_bins,*params),color=['r','g','b'][freq_bin],linestyle=['-','--',':'][spat_bin],linewidth=2,label=labels[3*spat_bin+freq_bin])
                print params
                all_params.append(params)
                plt.ylim(0,1.1)
		plt.xlim(0,2.5)
		if spat_bin==0 and freq_bin==0:
	                plt.axhline(1)
                #plt.xlim((min(Iflux_bins),max(Iflux_bins)))
                tr=raw_input()
		#plt.clf()
                #break
        #break



import matplotlib.pyplot as mpl
mpl.rcParams['font.size'] = 20

plt.ylabel('Completeness',fontsize=20)
plt.xlabel('Integrated Flux (Jy km s$^{-1}$)',fontsize=20)
plt.legend(loc='lower right',fontsize=20,frameon=False)

from matplotlib.ticker import MultipleLocator
ax=plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(.05))
ax.yaxis.set_minor_locator(MultipleLocator(.05))
ax.tick_params(which='major',length=8,width=1.5)
ax.tick_params(which='minor',length=4,width=1.5)


np.reshape([list(x) for x in all_params],(3,3,2))


