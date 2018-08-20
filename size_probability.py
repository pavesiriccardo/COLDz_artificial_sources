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
inp=open('all_obj.dat')
objs=cPickle.load(inp)
inp.close()

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

for inj_spat in range(3):
        counts={-1:np.zeros(bins.shape[0]),0:np.zeros(bins.shape[0]),2:np.zeros(bins.shape[0]),4:np.zeros(bins.shape[0]),6:np.zeros(bins.shape[0]),8:np.zeros(bins.shape[0]),10:np.zeros(bins.shape[0])}
        for idxx,meas_spat in enumerate([-1,0,2,4,6,8,10]):
                objects=[[obj[1] for obj in meas_and_inject if bin-.05<obj[1]<bin+.05 and obj[5]==meas_spat and obj[3]==inj_spat_list[inj_spat]] for bin in bins]  #and obj[4]==0
                counts[meas_spat]=np.array([1.*len(objects[idx]) for idx,bin in enumerate(bins)])
        for idbin,bin in enumerate(bins):
                norm=0
                for idxx,meas_spat in enumerate([-1,0,2,4,6,8,10]):
                        norm+=counts[meas_spat][idbin]
                for idxx,meas_spat in enumerate([-1,0,2,4,6,8,10]):
                        counts[meas_spat][idbin]=counts[meas_spat][idbin]/norm
        likelihood.append(counts)



Prior_injected=np.array([.88,.1,.02])

posterior=dict()
for idxx,meas_spat in enumerate([-1,0,2,4,6,8,10]):
        post=np.einsum('ab,a->ab',np.array([likelihood[inj][meas_spat] for inj in [0,1,2]]),Prior_injected)
        posterior[meas_spat]=post/np.sum(post,0)



import matplotlib.colors as colors
import matplotlib.cm as cmx

jet = cm = plt.get_cmap('jet')
cNorm  = colors.Normalize(vmin=0, vmax=7)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
cNorm  = colors.Normalize(vmin=0, vmax=2)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

for idxx,meas_spat in enumerate([-1,0,2,4,6,8]):
        plt.subplot(2,3,idxx+1)
        for idinj in range(3):
                plt.plot(bins[1:],posterior[meas_spat][idinj][1:],c=scalarMap.to_rgba(idinj))
                plt.ylim(0,1)
		plt.xlabel('SNR')
		plt.text(4.2,.85,str(meas_spat/2)+'\'\'',fontsize=16)
                #plt.title('measured:'+str(meas_spat/2)+'\'\'')
		if idxx==0 or idxx==3:
                        plt.ylabel('Probability')
                else:
                        #plt.gca().get_yaxis().set_visible(False)
			plt.gca().set_yticklabels([])
                        plt.xlim((bins[1],7))
		if idxx<3:
			plt.gca().set_xticklabels([])
			plt.xlabel('')
		else:
			plt.ylim((0,.99))


plt.subplots_adjust(wspace=0,hspace=0)
###################################


