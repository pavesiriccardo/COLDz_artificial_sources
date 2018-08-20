
import sys
import pyfits, numpy as np,matplotlib.pyplot as plt,os
fwhm=2.355
chanwidth=0.002


if (sys.argv[1]).endswith('/'):
	newdir=sys.argv[1]
else:
	newdir=sys.argv[1]+'/'



from inj_source_class import inj_source

os.chdir(newdir)


std_spat=float(sys.argv[3])
std_freq=float(sys.argv[4])
flux=float(sys.argv[2])
#NN= (np.pi*3.25*2.9)/np.log(2)/4*4

##############ANALYZE

#outp=open('analysis.log','w')
#print >>outp,flux,std_spat,std_freq


import cPickle
inp=open('reduint.dat')
reduint=cPickle.load(inp)
inp.close()

inp=open('combined.dat')
combined=cPickle.load(inp)
inp.close()


centerint=np.loadtxt('injected_list.txt')
centerint=centerint[:,:3].astype(int)

injected_objects=[inj_source(pos,flux,std_spat,std_freq) for pos in centerint]

for idx,red in enumerate(reduint):
	injected_objects[idx].set_reduint(red,combined[idx])


inp=open('aper_fluxes_smooth.txt')
for obj in injected_objects:
	if not obj.isnan:
		line=inp.readline()
		if line.split()[6]=='no_fit':
			obj.set_aper_flux('no_fit')
		else:
			obj.set_aper_flux(line)


inp.close()

inp=open('1pix_fluxes_smooth.txt')
for obj in injected_objects:
	if not obj.isnan:
		line=inp.readline()
		if line.split()[6]=='no_fit':
			obj.set_1pix_flux_smooth('no_fit')
		else:
			obj.set_1pix_flux_smooth(line)



inp.close()

inp=open('1pix_fluxes_nosmooth.txt')
for obj in injected_objects:
	if not obj.isnan:
		line=inp.readline()
		if line.split()[6]=='no_fit':
			obj.set_1pix_flux_nosmooth('no_fit')
		else:
			obj.set_1pix_flux_nosmooth(line)



inp.close()

outp=open('analysis_injected_objects.dat','w')
cPickle.dump(injected_objects,outp)
outp.close()


'''
for idx,obj in enumerate(centerint):
	coord=(obj[0]+1,1+obj[1],obj[2])
	bin=int(round(fwhm*std_freq))
	cube_narrow=signal[0].section[(coord[2]-bin/2):(coord[2]+bin/2),coord[1]-30:coord[1]+30,coord[0]-30:coord[0]+30]
	startchan=max([0,coord[2]-25])
	endchan=min([2017,coord[2]+25])
	cube=signal[0].section[startchan:endchan,coord[1]-30:coord[1]+30,coord[0]-30:coord[0]+30] 
	y,x=np.mgrid[:60,:60]
	img=np.sum(cube_narrow,0)
	ys, xs = np.mgrid[:10, :10]
	try:
		resfit,covm=curve_fit(gauss2D,(xs,ys),np.array(img[30-5:30+5,30-5:30+5]).flatten(),(5,5,1,1,1,0))
		maxI= np.absolute(fwhm*resfit[3])/2.
		chanwidth=0.002
		minI=np.absolute(fwhm*resfit[4])/2.
		angI=resfit[5]
		c=np.cos(angI*np.pi/180.)
		s=np.sin(angI*np.pi/180.)
		rot=np.array([[c,-s],[s,c]])
		A=np.dot(np.dot(rot,np.diag([1./maxI**2,1./minI**2])),np.transpose(rot))
		fu=A[0,0]*(x-resfit[0]-25)**2+2*A[0,1]*(x-resfit[0]-25)*(y-resfit[1]-25)+A[1,1]*(y-resfit[1]-25)**2
		spectrum=[np.sum(np.where(fu<1,cube[i],0))/NN*2.  for i in range(cube.shape[0])]
		fitspec,covm_spec=curve_fit(gauss,range(cube.shape[0]),spectrum,(len(spectrum)/2,3,1))
		cent_line=(fitspec[0]+startchan)*2*chanwidth+30.965
		cent_line_err=(np.sqrt(covm_spec[0,0]))*2*chanwidth
		Sdeltv=fitspec[2]*np.sqrt(2*np.pi)*fitspec[1]*2*chanwidth/cent_line*3e5
		if not np.all(covm<np.inf):
			covm=np.zeros((6,6))
		if not np.all(covm_spec<np.inf):
			covm_spec=np.zeros((6,6))
		print >>outp,obj,'noise:',centers[idx][3],'fit:',maxI,np.sqrt(covm[3,3])*fwhm/2.,minI,np.sqrt(covm[4,4])*fwhm/2.,cent_line,cent_line_err,fitspec[2]*1e3,1e3*np.sqrt(covm_spec[2,2]),fwhm*fitspec[1]*2*chanwidth/cent_line*3e5,fwhm*np.sqrt(covm_spec[1,1])*2*chanwidth/cent_line*3e5 ,Sdeltv,np.sqrt(2*np.pi)*2*chanwidth/cent_line*3e5*np.sqrt(covm_spec[1,1]/fitspec[1]**2+covm_spec[2,2]/fitspec[2]**2)*fitspec[1]*fitspec[2]
		#single pix spectrum
		peak=[np.argmax(img[30-5:30+5,30-5:30+5])/10,np.argmax(img[30-5:30+5,30-5:30+5])%10]
		spectrum_1pix=[cube[ix][peak[0]+25,25+peak[1]] for ix in range(cube.shape[0])]
		fitspec,covm_spec=curve_fit(gauss,range(cube.shape[0]),spectrum_1pix,(coord[2]-startchan,bin/fwhm,2e-3))
		if not np.all(covm_spec<np.inf):
			covm_spec=np.zeros((6,6))
		cent_line_err=(np.sqrt(covm_spec[0,0]))*2*chanwidth
		Sdeltv=fitspec[2]*np.sqrt(2*np.pi)*fitspec[1]*2*chanwidth/cent_line*3e5
		print >>outp,obj,cent_line,cent_line_err,fitspec[2]*1e3,1e3*np.sqrt(covm_spec[2,2]),fwhm*fitspec[1]*2*chanwidth/cent_line*3e5,fwhm*np.sqrt(covm_spec[1,1])*2*chanwidth/cent_line*3e5 ,Sdeltv,np.sqrt(2*np.pi)*2*chanwidth/cent_line*3e5*np.sqrt(covm_spec[1,1]/fitspec[1]**2+covm_spec[2,2]/fitspec[2]**2)*fitspec[1]*fitspec[2]
	except (RuntimeError,TypeError) as e:
		print >>outp,obj,'no_fit',e
	def equiv((a,b,e),(c,d,f)):
		return (a-c)**2+(b-d)**2+(e-f)**2<30+15
	found=False
	for index_reduint,obj2 in enumerate(reduint):
		if equiv(obj,obj2[1]):
					flux_file=open('fluxes_positive.txt')
					for id_line,line in enumerate(flux_file):
						if id_line==index_reduint:
							print >>outp,line[:-1]
					flux_file.close()
					del flux_file
					flux_file=open('fluxes_positive_1pix.txt')
					for id_line,line in enumerate(flux_file):
						if id_line==index_reduint:
							print >>outp,line[:-1]
					flux_file.close()
					del flux_file
					found=True
	if not found:
		print >>outp,'not found'
	listSNR=[]
	for spat_templ in [0,2,4,6,8,10,12]:
		for freq_templ in [4,8,12,16,20]:
			listSNR.append(toload[(spat_templ,freq_templ)][obj[2],obj[1],obj[0]]/conv_factors[(spat_templ,freq_templ)])
	print >>outp,listSNR
	print >>outp,'\n'
	outp.flush()


outp.close()
'''

print 'finished',newdir,flux,std_spat,std_freq
#os.system('rm -r bin*')
#os.system('rm -r *.image')
#os.system('rm -r *.fits')



