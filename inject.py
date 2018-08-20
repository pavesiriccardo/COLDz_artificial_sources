

import pyfits, numpy as np,matplotlib.pyplot as plt,os,sys
import cPickle
fwhm=2.355
chanwidth=0.002
def gauss(x,mean,std):
		return np.exp(-(x-mean)**2/(2.*std**2))


def mosaicsens(xpix,ypix,freqGHz=34.,sorte=False):
        nSamples=10000
        factor=freqGHz/60*0.5
        pctrs=[(750, 750),(917.42261089116289, 452.48369657494271), (806.96861768041856, 452.42720948105159), (917.27641014690346, 643.54361577204429), (806.91944976081186, 643.48717830145449), (972.55227473819843, 548.05291404475804), (862.1467847912545, 547.96840887824237), (751.74125821542486, 547.93998944819577), (1027.8765912332021, 452.59629406062078), (1027.6333576756713, 643.65611434189077), (696.56242900649295, 643.48680190600885), (1082.9577056577434, 548.19350494100092), (972.35753377858964, 739.11284658633997), (862.04908851816867, 739.02841569859129), (751.7406066611278, 739.00002124865318), (696.5145640484626, 452.42683275461894), (1027.3900902000746, 834.6959976649191), (917.13018901551754, 834.58359802467032), (806.87027498486566, 834.52721018446596), (696.61030063911267, 834.52683411992791), (641.33571740936554, 547.96765577055623), (1082.6659200634335, 739.25331390538986), (586.06054236450245, 452.48256641714437), (586.20544017218651, 643.54248660716166), (862.24445720273411, 356.90848560839879), (751.74190961058332, 356.88004120398051), (641.23934777571037, 356.90773183893964), (530.93025470155146, 548.05140780128795), (641.432110586313, 739.02766325283676), (807.01777873967694, 261.34730471967583), (696.46670576900578, 261.34692766232587), (1137.990200059244, 643.824673938302), (475.84850564656352, 643.6542323830455), (916.98391698884677, 1025.6435023421591), (972.16274537836364, 930.17274167789776), (861.9513684454613, 930.08838508692384), (586.35035818532685, 834.58246985276378), (806.82108309263072, 1025.5871641495185), (531.12369254100315, 739.11134166685088), (420.81530503608366, 739.25105652666957), (751.73995494811072, 930.06001562330835), (641.5285272455709, 930.08763330328316), (861.8536245837787, 1121.148296060634), (476.09046999248693, 834.69411736097834), (696.65818893444293, 1025.5867884161421), (586.49532664011895, 1025.6423751635773), (806.77187406519113, 1216.6670992020192), (531.31717750393295, 930.17123808266251), (751.73930307643855, 1121.1199515893752), (751.74256084660067, 165.82017749913916), (1137.6498863317202, 834.86440903243613), (1027.1467380502311, 1025.7558030183136), (916.83759401061093, 1216.7233877297674), (1192.9741551254456, 739.44981751698367), (1082.3740633877021, 930.31308538970177), (971.96790955885888, 1121.2325783366778), (861.75585693364269, 1312.2081476370672),(640.14301571,   164.84791243)]	
        inverseIncrementRadius_p=float(nSamples-1)/43.
        A=[0 for i in range(58)]
        def poly(x):
                coef=[ 1.0]
                coef.append( -1.300633e-03)
                coef.append( 6.480550e-07)
                coef.append(-1.267928e-10)
                if x<=int(43.*inverseIncrementRadius_p)/inverseIncrementRadius_p:
                        return coef[0]+coef[1]*x**2+coef[2]*x**4+coef[3]*x**6
                else:
                        return 0
        for pnt in range(1,58):
                        A[pnt]=poly(int(np.sqrt((xpix-pctrs[pnt][0])**2+(ypix-pctrs[pnt][1])**2)*factor*inverseIncrementRadius_p)/inverseIncrementRadius_p)
	if sorte:
	        return (np.flipud(np.argsort(A)),np.flipud(np.sort(A)))
	else:
		return np.array(A)


fwhm=2.355
chanwidth=0.002
def gauss(x,mean,std):
		return np.exp(-(x-mean)**2/(2.*std**2))



def convolve(maxI,minI,angI,maxB,minB,angB):
        import numpy as np
        c=np.cos(angI*np.pi/180.)
        s=np.sin(angI*np.pi/180.)
        rot=np.array([[c,-s],[s,c]])
        A=np.dot(np.dot(rot,np.diag([1./maxI**2,1./minI**2])),np.transpose(rot))
        cB=np.cos(angB*np.pi/180.)
        sB=np.sin(angB*np.pi/180.)
        rot=np.array([[cB,-sB],[sB,cB]])
        B=np.dot(np.dot(rot,np.diag([1./maxB**2,1./minB**2])),np.transpose(rot))
        result=np.linalg.inv(np.linalg.inv(A)+np.linalg.inv(B))
        minR,maxR=1/np.sqrt(np.linalg.eig(result)[0])
        angR=-np.arctan(np.linalg.eig(result)[1][1,1]/np.linalg.eig(result)[1][1,0])/np.pi*180.
        if angR<0:
                angR+=180
        return [maxR,minR,angR]



def gauss2D((x,y),xmean,ymean,norm,a,b,theta):
	c=np.cos(theta*np.pi/180.)
	s=np.sin(theta*np.pi/180.)
	rot=np.array([[c,-s],[s,c]])
	A=np.dot(np.dot(rot,np.diag([1./a**2,1./b**2])),np.transpose(rot))
	fu=A[0,0]*(x-xmean)**2+2*A[0,1]*(x-xmean)*(y-ymean)+A[1,1]*(y-ymean)**2
	return norm*np.exp(-fu/2).ravel()


#these are the beam sizes as they came out of clean
inp=open('GNbeams.dat')
GNbeams=cPickle.load(inp)
inp.close()


if (sys.argv[1]).endswith('/'):
	newdir=sys.argv[1]
else:
	newdir=sys.argv[1]+'/'




os.mkdir(newdir)   #REMOVE THIS COMMENT


os.system('cp GNmosaic_nosmooth.fits '+newdir+'signal_nosmooth.fits')  #REMOVE THIS COMMENT
os.system('cp GNmosaic_smooth.fits '+newdir+'signal_smooth.fits')  #REMOVE THIS COMMENT
os.system('cp SNR_cube_smooth.fits '+newdir+'SNR_smooth.fits')  #REMOVE THIS COMMENT
os.system('cp SNR_cube_nosmooth.fits '+newdir+'SNR_nosmooth.fits')  #REMOVE THIS COMMENT
os.chdir(newdir)

###################
#Load pointing noise before and after smoothing
inp=open('../pre_smoothing_noise.dat')
rr_pre=np.array(cPickle.load(inp))
inp.close()
rr_pre=np.transpose(rr_pre)
rtemp=np.zeros((2016,58))
rtemp[:-1,1:]=rr_pre
rr_pre=rtemp
rr_pre[:,0]=1
rr_pre[rr_pre==0]=np.inf
####
rr_post=[]
f=open('../rlist_57.txt','r')
for line in f:
	rr_post.append(map(float,('1 '+line).split()))

f.close()
rr_post=np.array(rr_post)
rr_post[rr_post==0]=np.inf
###################


std_spat=float(sys.argv[3])
std_freq=float(sys.argv[4])
flux=float(sys.argv[2])
#NN= (np.pi*3.38*2.91)/np.log(2)/4*4


f=pyfits.open('signal_nosmooth.fits')
Num=12500
centers=(np.random.rand(Num,2)*1500).astype(int)
centers2=np.zeros((Num,3))
centers2[:,:2]=centers
centers2[:,2]=np.random.rand(Num)*2015
centers2=centers2.astype(int)
centers=np.array([(int(x),int(y),chn) for (x,y,chn) in centers2 if f[0].section[chn,y,x]!=0 ])#[:1000]
#check!!!!!:
#[x for x in centers[:,2] if x<1.5*fwhm*std_freq or x>(4034-1.5*fwhm*std_freq)]
len(centers)
conflict=[]
for obj in centers:
	for obj2 in centers:
			if np.sum(obj==obj2)<3:
				if (obj[0]-obj2[0])**2+(obj[1]-obj2[1])**2<(4*7)**2 and np.absolute(obj[2]-obj2[2])<2*32:
					conflict.append((tuple(obj),tuple(obj2)))

len(conflict)/2 #PAIRS!
for obj in conflict:
	conflict.remove((obj[1],obj[0]))




toremove=set([a for (a,b) in conflict])
toremove=toremove|set([tuple(x) for x in centers if x[2]<1.5*fwhm*std_freq or x[2]>(2015-1.5*fwhm*std_freq)])
len(toremove)
Newcenters=[x for x in centers if tuple(x) not in toremove][:2500]
len(Newcenters)
centers3=[(x,y,chn) for (x,y,chn) in Newcenters]


#If have list of coords already:
#centers3=np.loadtxt('injected_list.txt')
#centers3=[[int(obj[0]),int(obj[1]),int(obj[2])] for obj in centers3]

#######################################################################################################
#Inject in the unsmoothed mosaic, unsing the beam sizes from the original pointings, from CASA

f=pyfits.open('signal_nosmooth.fits',mode='update')
oldimage=np.copy(f[0].data==0)
np.savetxt('injected_list.txt',centers3)


physical_gau=[std_spat*fwhm/2,std_spat*fwhm/2,0.]

df=2*int(1.5*fwhm*std_freq)+1
zs,ys,xs=np.mgrid[:df,:30,:30]

data_cube=f[0].data
for obj in centers3:
	pos=obj
	freq=29.9605+0.004*pos[2]
	temp=np.zeros((df,30,30))
	Alis=mosaicsens(pos[0],pos[1],freq)
	inv_r=1/rr_pre[pos[2],:]
	for pnt in range(1,58):
		if Alis[pnt]>0:
				bmaj,bmin,bpa=list(GNbeams[pnt][pos[2]][:3])
				smaj,smin,spa=convolve(*(physical_gau+[bmaj,bmin,bpa]))
				#print obj,pnt,bmaj,bmin,bpa,smaj,smin,spa,flux*(bmaj*bmin)/(smaj*smin)*1e3
				temp+=Alis[pnt]**2*inv_r[pnt]**2*flux*(bmaj*bmin)/(smaj*smin)*gauss2D((xs+obj[0]-15,ys+obj[1]-15),pos[0],pos[1],1,smaj/fwhm*2,smin/fwhm*2,spa).reshape((df,30,30))*gauss(zs+pos[2]-int(1.5*fwhm*std_freq),pos[2],std_freq)
	temp/=np.sum(Alis**2*inv_r**2)
	data_cube[(pos[2]-int(1.5*fwhm*std_freq)):(pos[2]+1+int(1.5*fwhm*std_freq)),(obj[1]-15):(obj[1]+15),(obj[0]-15):(obj[0]+15)]+=temp




data_cube[oldimage]=0
f[0].data=data_cube
f.close()


##################################################################
#Inject in the smoothed mosaic. Using the beam sizes measured from the smoothed PSF of the individual pointings.
fits=np.loadtxt('../fit_sizes.txt')   #make sure to use the version with the angles as well

f=pyfits.open('signal_smooth.fits',mode='update')
oldimage=np.copy(f[0].data==0)

physical_gau=[std_spat*fwhm/2,std_spat*fwhm/2,0.]

df=2*int(1.5*fwhm*std_freq)+1
zs,ys,xs=np.mgrid[:df,:30,:30]

data_cube=f[0].data
for obj in centers3:
	pos=obj
	freq=29.9605+0.004*pos[2]
	temp=np.zeros((df,30,30))
	Alis=mosaicsens(pos[0],pos[1],freq)
	inv_r=1/rr_post[pos[2],:]
	for pnt in range(1,58):
		if Alis[pnt]>0:
				bmaj,bmin,bpa=list(fits[2015*(pnt-1)+pos[2],2:5])
				smaj,smin,spa=convolve(*(physical_gau+[bmaj,bmin,bpa]))
				#print obj,pnt,bmaj,bmin,bpa,smaj,smin,spa,flux*(bmaj*bmin)/(smaj*smin)*1e3
				temp+=Alis[pnt]**2*inv_r[pnt]**2*flux*(bmaj*bmin)/(smaj*smin)*gauss2D((xs+obj[0]-15,ys+obj[1]-15),pos[0],pos[1],1,smaj/fwhm*2,smin/fwhm*2,spa).reshape((df,30,30))*gauss(zs+pos[2]-int(1.5*fwhm*std_freq),pos[2],std_freq)
	temp/=np.sum(Alis**2*inv_r**2)
	data_cube[(pos[2]-int(1.5*fwhm*std_freq)):(pos[2]+1+int(1.5*fwhm*std_freq)),(obj[1]-15):(obj[1]+15),(obj[0]-15):(obj[0]+15)]+=temp




data_cube[oldimage]=0
f[0].data=data_cube
f.close()

#######################
#Inject in SNR_cube, nosmooth
#I've previously made sure these SNR cubes dont have any nan's, or 1s or infinities
f=pyfits.open('SNR_nosmooth.fits',mode='update')
noise_nosmooth=pyfits.open('../noise_edge_correct_nosmooth.fits')
noise_nosmooth_cube=noise_nosmooth[0].data
noise_nosmooth.close()
oldimage=np.copy(f[0].data==0)

physical_gau=[std_spat*fwhm/2,std_spat*fwhm/2,0.]

df=2*int(1.5*fwhm*std_freq)+1
zs,ys,xs=np.mgrid[:df,:30,:30]

data_cube=f[0].data
for obj in centers3:
	pos=obj
	freq=29.9605+0.004*pos[2]
	temp=np.zeros((df,30,30))
	Alis=mosaicsens(pos[0],pos[1],freq)
	inv_r=1/rr_pre[pos[2],:]
	for pnt in range(1,58):
		if Alis[pnt]>0:
				bmaj,bmin,bpa=list(GNbeams[pnt][pos[2]][:3])
				smaj,smin,spa=convolve(*(physical_gau+[bmaj,bmin,bpa]))
				#print obj,pnt,bmaj,bmin,bpa,smaj,smin,spa,flux*(bmaj*bmin)/(smaj*smin)*1e3
				temp+=Alis[pnt]**2*inv_r[pnt]**2*flux*(bmaj*bmin)/(smaj*smin)*gauss2D((xs+obj[0]-15,ys+obj[1]-15),pos[0],pos[1],1,smaj/fwhm*2,smin/fwhm*2,spa).reshape((df,30,30))*gauss(zs+pos[2]-int(1.5*fwhm*std_freq),pos[2],std_freq)
	temp/=np.sum(Alis**2*inv_r**2)
	noise_cube=noise_nosmooth_cube[(pos[2]-int(1.5*fwhm*std_freq)):(pos[2]+1+int(1.5*fwhm*std_freq)),(obj[1]-15):(obj[1]+15),(obj[0]-15):(obj[0]+15)]
	data_cube[(pos[2]-int(1.5*fwhm*std_freq)):(pos[2]+1+int(1.5*fwhm*std_freq)),(obj[1]-15):(obj[1]+15),(obj[0]-15):(obj[0]+15)]+=(temp/noise_cube)




data_cube[oldimage]=0
f[0].data=data_cube
f.close()
del noise_nosmooth_cube
#################################
#Inject in SNR_cube smooth
f=pyfits.open('SNR_smooth.fits',mode='update')
noise_smooth=pyfits.open('../noise_edge_correct_smooth.fits')
noise_smooth_cube=noise_smooth[0].data
noise_smooth.close()
oldimage=np.copy(f[0].data==0)

physical_gau=[std_spat*fwhm/2,std_spat*fwhm/2,0.]

df=2*int(1.5*fwhm*std_freq)+1
zs,ys,xs=np.mgrid[:df,:30,:30]

data_cube=f[0].data
for obj in centers3:
	pos=obj
	freq=29.9605+0.004*pos[2]
	temp=np.zeros((df,30,30))
	Alis=mosaicsens(pos[0],pos[1],freq)
	inv_r=1/rr_post[pos[2],:]
	for pnt in range(1,58):
		if Alis[pnt]>0:
				bmaj,bmin,bpa=list(fits[2015*(pnt-1)+pos[2],2:5])
				smaj,smin,spa=convolve(*(physical_gau+[bmaj,bmin,bpa]))
				#print obj,pnt,bmaj,bmin,bpa,smaj,smin,spa,flux*(bmaj*bmin)/(smaj*smin)*1e3
				temp+=Alis[pnt]**2*inv_r[pnt]**2*flux*(bmaj*bmin)/(smaj*smin)*gauss2D((xs+obj[0]-15,ys+obj[1]-15),pos[0],pos[1],1,smaj/fwhm*2,smin/fwhm*2,spa).reshape((df,30,30))*gauss(zs+pos[2]-int(1.5*fwhm*std_freq),pos[2],std_freq)
	temp/=np.sum(Alis**2*inv_r**2)
	noise_cube=noise_smooth_cube[(pos[2]-int(1.5*fwhm*std_freq)):(pos[2]+1+int(1.5*fwhm*std_freq)),(obj[1]-15):(obj[1]+15),(obj[0]-15):(obj[0]+15)]
	data_cube[(pos[2]-int(1.5*fwhm*std_freq)):(pos[2]+1+int(1.5*fwhm*std_freq)),(obj[1]-15):(obj[1]+15),(obj[0]-15):(obj[0]+15)]+=(temp/noise_cube)




data_cube[oldimage]=0
f[0].data=data_cube
f.close()
del noise_smooth_cube
###################




