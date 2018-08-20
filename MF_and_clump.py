import pyfits,numpy as np,os,sys,gc,cPickle

def calc_templ(N1,N2,N3,templ_freq,spat_fwhm):
	fwhm=2.355
	def gauss(x,mean,fw_hm):
		return np.exp(-(x-mean)**2/fw_hm**2*2.7725887222397811)
	if spat_fwhm==0:
		s1=1e-4/fwhm
		s2=1e-4/fwhm
	else:
		s1=spat_fwhm/fwhm
		s2=spat_fwhm/fwhm
	s3=templ_freq/fwhm
	r1=N1/(s1*2*np.pi)
	r2=N2/(s2*2*np.pi)
	r3=N3/(s3*2*np.pi)
	normaliz=np.sqrt(N1*N2*N3/np.sqrt(np.pi)**3/(r1*r2*r3))
	y,x=np.indices((N2,N1))
	center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
	r = np.hypot(x - center[0], y - center[1])
	spat=np.exp(-r**2/2/(r1)**2)
	my_version=np.reshape(np.tile(spat,(int(N3),1)),(int(N3),N1,N2))
	argh=np.array([gauss(x,N3/2,r3*fwhm) for x in range(int(N3))])
	my_version=np.einsum('i,ikl->ikl',argh,my_version)
	my_version*=normaliz
	return np.fft.ifftshift(my_version)



if (sys.argv[1]).endswith('/'):
	newdir=sys.argv[1]
else:
	newdir=sys.argv[1]+'/'


os.chdir(newdir)
#########
#In case I dont want to inject directly into SNR cubes:
#os.system('cp signal_smooth.fits SNR_smooth.fits')
#SNR_smooth=pyfits.open('SNR_smooth.fits',mode='update')
#noise_smooth=pyfits.open('../noise_edge_correct_smooth.fits')
#snr_data=SNR_smooth[0].data/noise_smooth[0].data
#noise_smooth.close()
#SNR_smooth[0].data=np.where(np.isnan(snr_data),0,snr_data)
#del snr_data
#SNR_smooth.close()


#os.system('cp signal_nosmooth.fits SNR_nosmooth.fits')
#SNR_nosmooth=pyfits.open('SNR_nosmooth.fits',mode='update')
#noise_nosmooth=pyfits.open('../noise_edge_correct_nosmooth.fits')
#snr_data=SNR_nosmooth[0].data/noise_nosmooth[0].data
#noise_nosmooth.close()
#SNR_nosmooth[0].data=np.where(np.isnan(snr_data),0,snr_data)
#del snr_data
#SNR_nosmooth.close()
################

noise_set14=np.loadtxt('../noise_set14.dat')
noise_notset14=np.loadtxt('../noise_notset14.dat')
f1_set14=pyfits.open('../f1_set14.fits')
f1_set14_cube=f1_set14[0].data
f2_notset14=pyfits.open('../f2_notset14.fits')
f2_notset14_cube=f2_notset14[0].data



def find_top(peaks,tag,shiftchan,shifty,shiftx):
	#tops.put((np.max(peaks),np.unravel_index(np.argmax(peaks),peaks.shape),tag))
	small_posn=np.unravel_index(np.argmax(peaks),peaks.shape)
	return (np.max(peaks),(small_posn[2]+shiftx,small_posn[1]+shifty,small_posn[0]+shiftchan),tag)
	


def clump_different_templates(combined):
	import operator
	#tag_combined=[(x[0],(x[1][2]+shiftx,x[1][1]+shifty,x[1][0]+shiftchan),x[2]) for x in combined]
	combined.sort(key=operator.itemgetter(0),reverse=True)
	return combined[0]
	


#Load the injected positions
f=pyfits.open('SNR_smooth.fits')
SNR_cube=f[0].data
f=pyfits.open('SNR_nosmooth.fits')
SNR_nosmooth=f[0].data
f.close()

centers3=np.loadtxt('injected_list.txt')
centers3=[[int(obj[0]),int(obj[1]),int(obj[2])] for obj in centers3]


blockL=70   #70
blockdf=70
templFFT=[]
for spat_fwhm in [0,2,4,6,8]:
		for freq_fwhm in [4,8,12,16,20]:
			templFFT.append(calc_templ(blockL+1,blockL+1,blockdf+1,freq_fwhm,spat_fwhm))



def make_SNR_block(cube,obj):
	startchan=np.max([0,obj[2]-blockdf/2])
	endchan=np.min([2015,obj[2]+blockdf/2+1])
	startx=np.max([0,obj[0]-blockL/2])
	endx=np.min([1500,obj[0]+blockL/2+1])
	starty=np.max([0,obj[1]-blockL/2])
	endy=np.min([1500,obj[1]+blockL/2+1])
	initial_block=cube[startchan:endchan,starty:endy,startx:endx]
	return np.pad(initial_block,((startchan-(obj[2]-blockdf/2),(obj[2]+blockdf/2+1)-endchan),(starty-(obj[1]-blockL/2),(obj[1]+blockL/2+1)-endy),(startx-(obj[0]-blockL/2),(obj[0]+blockL/2+1)-endx)),'constant')
	


def M_filter(idx_temp):
			blockMF=np.real(np.fft.ifftn(np.fft.fftn(SNR_block)*templFFT[idx_temp]))
			noise=np.sqrt(noise_set14[idx_temp]**2*f1_set14_cube[obj[2],obj[1],obj[0]]**2+noise_notset14[idx_temp]**2*f2_notset14_cube[obj[2],obj[1],obj[0]]**2)
			blockMF/=noise
			return find_top(blockMF[(blockdf/2-small_size):(blockdf/2+small_size),(blockL/2-small_size):(blockL/2+small_size),(blockL/2-small_size):(blockL/2+small_size)],template_list[idx_temp],shiftchan,shifty,shiftx)


def M_filter_nosmooth(idx_temp): ##DO THIS
			blockMF=np.real(np.fft.ifftn(np.fft.fftn(SNR_block_nosmooth)*templFFT[idx_temp]))
			noise=[0.000150622570289,0.000150903613705,0.000151172867448,0.00015140215036,0.000151593834529][idx_temp]
			#noise=np.sqrt(noise_set14[idx_temp]**2*f1_set14_cube[obj[2],obj[1],obj[0]]**2+noise_notset14[idx_temp]**2*f2_notset14_cube[obj[2],obj[1],obj[0]]**2)
			blockMF/=noise
			return find_top(blockMF[(blockdf/2-small_size):(blockdf/2+small_size),(blockL/2-small_size):(blockL/2+small_size),(blockL/2-small_size):(blockL/2+small_size)],template_list_nosmooth[idx_temp],shiftchan,shifty,shiftx)


template_list=[(0, 4), (0, 8), (0, 12), (0, 16), (0, 20), (2, 4), (2, 8), (2, 12), (2, 16), (2, 20), (4, 4), (4, 8), (4, 12), (4, 16), (4, 20), (6, 4), (6, 8), (6, 12), (6, 16), (6, 20), (8, 4), (8, 8), (8, 12), (8, 16), (8, 20)]
template_list_nosmooth=[(-1, 4), (-1, 8), (-1, 12), (-1, 16), (-1, 20)]

import multiprocessing

reduint=[]
combined_all=[]
for obj in centers3:
	SNR_block=make_SNR_block(SNR_cube,obj)
	small_size=5
	shiftchan,shifty,shiftx=obj[2]-small_size,obj[1]-small_size,obj[0]-small_size
	pool=multiprocessing.Pool(processes=25)
	combined=pool.map(M_filter,range(25))
	pool.close()
	pool.terminate()
	SNR_block_nosmooth=make_SNR_block(SNR_nosmooth,obj)
	pool=multiprocessing.Pool(processes=5)
	combined_nosmooth=pool.map(M_filter_nosmooth,range(5))
	pool.close()
	pool.terminate()
	combined=combined+combined_nosmooth
	combined_all.append(combined)
	reduint_obj=clump_different_templates(combined)
	reduint.append(reduint_obj)


#watch out because some of them are going to have nan if they overlap with the RFI channels, where the SNR is nan



outp=open('reduint.dat','w')
cPickle.dump(reduint,outp)
outp.close()

outp=open('combined.dat','w')
cPickle.dump(combined_all,outp)
outp.close()


'''
SNR=pyfits.open('SNR.fits')
snr_data=SNR[0].data
FTt=np.fft.fftn(snr_data)
SNR.close()
del SNR
del snr_data
tosave=np.memmap('FTSNR_conj',dtype='complex128',mode='w+',shape=(2018,512,512))
tosave[:]=np.conjugate(FTt[:])
del FTt
del tosave
gc.collect()
logoutp=open('log.txt','w')

def M_filter((spat_templ,freq_templ)):
		FTt=np.memmap('FTSNR_conj',dtype='complex128',mode='r',shape=(2018,512,512))
		FTt2=np.memmap('../templates/templ'+str(spat_templ)+'.'+str(freq_templ),dtype='float64',mode='r',shape=(2018,512,512))
		#FTt2=np.load('../templates/templ'+str(spat_templ)+'.'+str(freq_templ)+'.npy')
		product=FTt*FTt2
		del FTt2
		del FTt
		peaks=np.fft.fftn(product)
		del product
		peaks=np.real(peaks)/2018/512/512
		tosave=np.memmap('MFiltered.'+str(spat_templ)+'.'+str(freq_templ),dtype='float64',mode='w+',shape=peaks.shape)
		tosave[:]=peaks[:]
		del tosave
		del peaks
		return 'done '+str(spat_templ)+'.'+str(freq_templ)





import multiprocessing,logging,time,random
argumen=[]
for spat_templ in [0,2,4,6,8,10,12]:#[3,6,7,8.5,10]:
	for freq_templ in [4,8,12,16,20]:
		argumen.append((spat_templ,freq_templ))
		#print >>logoutp, M_filter((spat_templ,freq_templ))
		#logoutp.flush()




#SECOND method

processes=[]
done_ones=[]
next_to_do=argumen
Numproc=len(next_to_do)

if True:#random.random()<.8:
	pool_size=3#4
else:
	pool_size=3



for i in range(pool_size):
	P=(next_to_do[len(next_to_do)-1], multiprocessing.Process(target=M_filter, args=(next_to_do.pop(),)))
	processes.append(P)
	P[1].start()



done=False
while not done:
	time.sleep(10.)
	for process in processes:
		if 'MFiltered.'+str(process[0][0])+'.'+str(process[0][1]) in os.listdir('.'):
			if len(processes)==1:
				process[1].join()
			done_ones.append(process[0])
			processes.remove(process)
			if len(done_ones)==Numproc:
				done=True
				print >>logoutp,done_ones
				logoutp.flush()
				break
			if len(next_to_do)>0:
				next_up=next_to_do.pop()
				P=(next_up, multiprocessing.Process(target=M_filter,args=(next_up,)))
				processes.append(P)
				P[1].start()
		elif not process[1].is_alive():
			print >>logoutp, 'redoing',process[0]
			logoutp.flush()
			processes.remove(process)
			P=(process[0], multiprocessing.Process(target=M_filter, args=(process[0],)))
			processes.append(P)
			P[1].start()


#######




def find_tops_and_clump(peaks,tops,thresh,tag,conv_factor):
	hipoints_chan,hipoints_y,hipoints_x=np.where(peaks>thresh*conv_factor)
	temp_tops=[]
	for obj_id,chan in enumerate(hipoints_chan):
		value=peaks[chan,hipoints_y[obj_id],hipoints_x[obj_id]]
		freq_rang=6				#was 8 before
		spat_rang=4
		startchan=np.max([0,chan-freq_rang])
		endchan=np.min([2017,chan+freq_rang])
		if np.max(peaks[startchan:endchan,hipoints_y[obj_id]-spat_rang:hipoints_y[obj_id]+spat_rang,hipoints_x[obj_id]-spat_rang:hipoints_x[obj_id]+spat_rang])<=value:
			temp_tops.append((value/conv_factor,(hipoints_x[obj_id],hipoints_y[obj_id],chan)))
	tops[tag]=temp_tops



SNR=pyfits.open('SNR.fits')
snr_data=SNR[0].data
combined=dict()
conv_factors=open('conv_factors.txt','w')
for spat_templ in [0,2,4,6,8,10,12]:#[3,6,7,8.5,10]:
	for freq_templ in [4,8,12,16,20]:
		peaks=np.memmap('MFiltered.'+str(spat_templ)+'.'+str(freq_templ),dtype='float64',mode='r',shape=(2018,512,512))
		conv_factor=np.std(peaks[snr_data!=0])
		print >>conv_factors,conv_factor
		find_tops_and_clump(peaks,combined,4,(spat_templ,freq_templ),conv_factor)
		del peaks






outp=open('combined.dat','w')
cPickle.dump(combined,outp)
outp.close()

import operator
def distsq((a,b,c),(d,e,f)):
	return (1.*a-1.*d)**2+(1.*b-1.*e)**2+(1.*c-1.*f)**2



tag_combined=[(x[0],x[1],key) for key in combined.keys() for x in combined[key]]
tag_combined.sort(key=operator.itemgetter(0),reverse=True)
redu=[]
reduced=[]
for obj in tag_combined:
		found=False
		for idx,red in enumerate(reduced):
			if distsq(redu[idx][1],obj[1])<28. and redu[idx][2]!=obj[2]:  #was 28
				red.append(obj)
				if redu[idx][0]<obj[0]:
					redu[idx]=(obj[0],tuple(np.average([x[1] for x in red],axis=0,weights=[x[0]**2 for x in red])),obj[2])
				else:
					redu[idx]=(redu[idx][0],tuple(np.average([x[1] for x in red],axis=0,weights=[x[0]**2 for x in red])),redu[idx][2])
				found=True
		if not found:
				reduced.append([obj])
				redu.append(obj)

reduint=[(x[0],(int(round(x[1][0])),int(round(x[1][1])),int(round(x[1][2]))   ),len(reduced[idx]),x[2]) for idx,x in enumerate(redu)]
reduced_sorted=[x for (y,x) in sorted(zip(reduint,reduced), key=lambda (k,v): operator.itemgetter(0)(v),reverse=True)]
reduint.sort(key=operator.itemgetter(0),reverse=True)



outp=open('reduint.dat','w')
cPickle.dump(reduint,outp)
outp.close()

os.system('rm FTSNR_conj')
import subprocess
print >>logoutp, subprocess.check_output('ls MFiltered* | wc', shell=True)
'''


