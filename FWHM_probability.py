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

############################

#size probability

meas_and_inject=[]
for obj in objs:	
	if not obj.isnan and obj.SNR>4:
		meas_and_inject.append((0,obj.SNR,obj.inj_flux,obj.inj_spat_fwhm,obj.inj_freq_fwhm,obj.peak_templ[0],obj.peak_templ[1]))




meas_and_inject=np.array(meas_and_inject)
bins=np.arange(4,7,.1)
likelihood=[]
inj_spat_list=[2.00175,4.9926,  8.9961]
inj_freq_list=[5.828625, 11.65725,17.485875]

for inj_freq in range(3):
        counts={4:np.zeros(bins.shape[0]),8:np.zeros(bins.shape[0]),12:np.zeros(bins.shape[0]),16:np.zeros(bins.shape[0]),20:np.zeros(bins.shape[0])}
        for idxx,meas_freq in enumerate([4,8,12,16,20]):
                objects=[[obj[1] for obj in meas_and_inject if bin-.05<obj[1]<bin+.05 and obj[6]==meas_freq and obj[4]==inj_freq_list[inj_freq]] for bin in bins]  #and obj[4]==0
                counts[meas_freq]=np.array([1.*len(objects[idx]) for idx,bin in enumerate(bins)])
        for idbin,bin in enumerate(bins):
                norm=0
                for idxx,meas_freq in enumerate([4,8,12,16,20]):
                        norm+=counts[meas_freq][idbin]
                for idxx,meas_freq in enumerate([4,8,12,16,20]):
                        counts[meas_freq][idbin]=counts[meas_freq][idbin]/norm
        likelihood.append(counts)



Prior_injected=np.array([.34,.33,.33])

posterior=dict()
for idxx,meas_freq in enumerate([4,8,12,16,20]):
        post=np.einsum('ab,a->ab',np.array([likelihood[inj][meas_freq] for inj in [0,1,2]]),Prior_injected)
        posterior[meas_freq]=post/np.sum(post,0)



import matplotlib.colors as colors
import matplotlib.cm as cmx

jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=7)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
cNorm  = colors.Normalize(vmin=0, vmax=2)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

for idxx,meas_freq in enumerate([4,8,12,16,20]):
        plt.subplot(2,3,idxx+1)
        for idinj in range(3):
                plt.plot(bins[1:],posterior[meas_freq][idinj][1:],c=scalarMap.to_rgba(idinj))
                plt.ylim(0,1)
                plt.xlabel('SNR')
                plt.text(4.2,.85,str(meas_freq)+' chans',fontsize=16)
                if idxx==0 or idxx==3:
                        plt.ylabel('Probability')
                else:
                        #plt.gca().get_yaxis().set_visible(False)
			plt.gca().set_yticklabels([])
                        plt.xlim((bins[1],7))
                #plt.title('measured:'+str(meas_freq)+' chans')
		if idxx<3:
			plt.gca().set_xticklabels([])
			plt.xlabel('')
		else:
			plt.ylim((0,.99))

plt.subplots_adjust(wspace=0,hspace=0)
###################################


