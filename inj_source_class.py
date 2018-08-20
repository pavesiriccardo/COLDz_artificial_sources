import numpy as np
class inj_source(object):
        def __init__(self,inj_posn,flux,std_spat,std_freq):
                self.inj_posn=inj_posn
                self.inj_flux=flux   #this is the peak flux, in the sense integrated spatially, but at the peak channel
                self.inj_spat_fwhm=std_spat*2.355
                self.inj_freq_fwhm=std_freq*2.355
                self.isnan=False
        def set_reduint(self,red,comb):
                self.SNR=red[0]
                if np.isnan(self.SNR):
                        self.isnan=True
			self.aper_maj=np.nan
                	self.aper_min=np.nan
                	self.aper_freq=np.nan
                	self.aper_flux=np.nan
                	self.aper_FWHM=np.nan
                	self.aper_int_flux=np.nan
			self.pix_smooth_freq=np.nan
                	self.pix_smooth_flux=np.nan
               		self.pix_smooth_FWHM=np.nan
                	self.pix_smooth_int_flux=np.nan
			self.pix_nosmooth_freq=np.nan
                	self.pix_nosmooth_flux=np.nan
                	self.pix_nosmooth_FWHM=np.nan
                	self.pix_nosmooth_int_flux=np.nan
                self.reduint_posn=red[1]
                self.peak_templ=red[2]
		self.combined_templates=comb
        def set_aper_flux(self,aper_string):
		if aper_string=='no_fit':
			self.aper_maj=np.nan
                	self.aper_min=np.nan
                	self.aper_freq=np.nan
                	self.aper_flux=np.nan
                	self.aper_FWHM=np.nan
                	self.aper_int_flux=np.nan
		else:
               		splitt=aper_string.split()
                	if (int(splitt[1][1:-1]),int(splitt[2][:-1]),int(splitt[3][:-2]))!= self.reduint_posn:
                	        print 'wrong aper_flux line'
                	self.aper_maj=float(splitt[6])
                	self.aper_min=float(splitt[8])
                	self.aper_freq=float(splitt[10])
                	self.aper_flux=float(splitt[12])
                	self.aper_FWHM=float(splitt[14])
                	self.aper_int_flux=float(splitt[16])
        def set_1pix_flux_smooth(self,pix_smooth_string):
		if pix_smooth_string=='no_fit':
			self.pix_smooth_freq=np.nan
                	self.pix_smooth_flux=np.nan
               		self.pix_smooth_FWHM=np.nan
                	self.pix_smooth_int_flux=np.nan
		else:
                	splitt=pix_smooth_string.split()
                	if (int(splitt[1][1:-1]),int(splitt[2][:-1]),int(splitt[3][:-2]))!= self.reduint_posn:
                        	print 'wrong pix_smooth line'
                	self.pix_smooth_freq=float(splitt[6])
                	self.pix_smooth_flux=float(splitt[8])
                	self.pix_smooth_FWHM=float(splitt[10])
                	self.pix_smooth_int_flux=float(splitt[12])
        def set_1pix_flux_nosmooth(self,pix_nosmooth_string):
		if pix_nosmooth_string=='no_fit':
			self.pix_nosmooth_freq=np.nan
                	self.pix_nosmooth_flux=np.nan
                	self.pix_nosmooth_FWHM=np.nan
                	self.pix_nosmooth_int_flux=np.nan
		else:
	                splitt=pix_nosmooth_string.split()
        	        if (int(splitt[1][1:-1]),int(splitt[2][:-1]),int(splitt[3][:-2]))!= self.reduint_posn:
        	                print 'wrong pix_nosmooth line'
        	        self.pix_nosmooth_freq=float(splitt[6])
        	        self.pix_nosmooth_flux=float(splitt[8])
        	        self.pix_nosmooth_FWHM=float(splitt[10])
        	        self.pix_nosmooth_int_flux=float(splitt[12])
