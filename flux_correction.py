#Probability of the aperture flux ratio, given the measured size, and SNR. Assumes the prob. for injected size, given the measured, which itself relies on a prior on sizes.
from inj_source_class import inj_source

objs=[]
for flux in range(12):
	for spat in range(3):
		for freq in range(3):
			os.chdir('art'+str(flux)+str(spat)+str(freq))
			inp=open('analysis_injected_objects.dat')
			objs+=cPickle.load(inp)
			inp.close()
			os.chdir('..')



measured_to_injected=[]
for obj in objs:	
	if not obj.isnan and obj.SNR>4:
		if obj.peak_templ[0]>-1:
			flux_meas=obj.aper_flux
		else:
			flux_meas=obj.pix_nosmooth_flux
		if np.isnan(flux_meas):
			flux_meas=0
		measured_to_injected.append((0,obj.SNR,flux_meas*1e-3/obj.inj_flux,obj.inj_flux,obj.inj_spat_fwhm,obj.inj_freq_fwhm,obj.peak_templ[0],obj.peak_templ[1]))





measured_to_injected=np.array(measured_to_injected)


bins=np.arange(4,7,.1)


meas_spat=-1
plt.clf()

flux_ratio_range=np.arange(.5,10,.1)
inj_spat_list=[2.00175,4.9926,  8.9961]
r_given_meas=np.zeros((bins.shape[0],flux_ratio_range.shape[0]-1))
for inj in range(3):
        objects=[[obj[2] for obj in measured_to_injected if bin-.5<obj[1]<bin+.5 and obj[6]==meas_spat and obj[4]==inj_spat_list[inj]] for bin in bins]
        r_given_meas+=np.array([posterior[meas_spat][inj][idbin]*np.histogram(objects[idbin],flux_ratio_range,normed=True)[0] for idbin,bin in enumerate(bins)])



for idbin,bin in enumerate(bins):
        throwaway=plt.plot(flux_ratio_range[:-1],r_given_meas[idbin])
        for ind_ratio_bin,ratio_bin in enumerate(flux_ratio_range):
                if np.nansum(r_given_meas[idbin][:ind_ratio_bin])>1/.1/2.:
                        break
        print 'SNR=',bin,'median is between',flux_ratio_range[ind_ratio_bin-1],'-',ratio_bin

for ind_ratio_bin,ratio_bin in enumerate(flux_ratio_range):
                if np.sum(np.nanmean(r_given_meas,0)[:ind_ratio_bin])>1/.1/2.:
                        break

print 'mean over SNR, median is between',flux_ratio_range[ind_ratio_bin-1],'-',ratio_bin
throwaway=plt.plot(flux_ratio_range[:-1],np.nanmean(r_given_meas,0),linewidth=2,c='k')

################
#all plots
cx=1
for meas_spat in [-1,0,2,4,6,8]:
                plt.subplot(2,3,cx)
                #plt.title('measured:'+str(meas_spat/2)+'\'\'')
                print cx
                cx+=1
                flux_ratio_range=np.arange(0.,10,.1)
                r_given_meas=np.zeros((bins.shape[0],flux_ratio_range.shape[0]-1))
                for inj in range(3):
                        objects=[[obj[2] for obj in measured_to_injected if bin-.5<obj[1]<bin+.5 and obj[6]==meas_spat and obj[4]==inj_spat_list[inj]] for bin in bins]
                        r_given_meas+=np.array([posterior[meas_spat][inj][idbin]*np.histogram(objects[idbin],flux_ratio_range,normed=True)[0] for idbin,bin in enumerate(bins)])
                for idbin,bin in enumerate(bins):
                        if 10<=idbin<20:    #this and the 10:20 two lines down is to look at 5<SNR<6
                                throwaway=plt.plot(flux_ratio_range[:-1],r_given_meas[idbin])
                throwaway=plt.plot(flux_ratio_range[:-1],np.nanmean(r_given_meas[10:20],0),linewidth=2,c='k')
                thro2=plt.plot(flux_ratio_range[:-1],my_lognorm(flux_ratio_range[:-1],*fit_lognorm(flux_ratio_range[:-1],np.nanmean(r_given_meas[10:20],0))),linewidth=2,c='r')
                print 'bestfit params:',fit_lognorm(flux_ratio_range[:-1],np.nanmean(r_given_meas[10:20],0))
                #plt.xlim((0,5))
		plt.xlabel('Flux-factor')
		#plt.ylabel('Probability density')
		plt.annotate(str(meas_spat/2)+'\'\'',xy=(.1,.1),xytext=(.05,.9),textcoords='axes fraction',fontsize=16)
                plt.gca().get_yaxis().set_visible(False)
                if cx==2 or cx==5:
                        plt.ylabel('Probability density')
                if cx==4 or cx==7:
                        plt.xlim((0,5))
                else:
                        plt.xlim((0,4.9))
		if cx<5:
			plt.gca().set_xticklabels([])
			plt.xlabel('')

plt.subplots_adjust(wspace=0,hspace=0)

####################
#Fit a lognormal to the black lines above...
def my_lognorm(x,mu,sigma):
        return 1/(x*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))

def fit_lognorm(x,y):
        chi2=lambda (mu,sigma):np.nansum((my_lognorm(x,mu,sigma)-y)**2)
        import scipy.optimize as op
        result=op.minimize(chi2, [0.,1.],method='Nelder-Mead')
        return result['x']
