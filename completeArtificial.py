#for each injected source do 2D gaussian fit and extract size, then spectral gaussian params
#calculate expected SNR
#check if in aper12 selection
#check if in aper1 selection and 
#correct SNR_1pix

#should call with folder_name flux std_spat std_freq



#MUST calculate all the theoretical SNR possible
import pyfits, numpy as np,matplotlib.pyplot as plt,os,sys,subprocess

fwhm=2.355
chanwidth=0.002

std_spat=float(sys.argv[3])
std_freq=float(sys.argv[4])
flux=float(sys.argv[2])

if (sys.argv[1]).endswith('/'):
	newdir=sys.argv[1]
else:
	newdir=sys.argv[1]+'/'


def emailme(subject_my):
	import smtplib
	import email.utils
	from email.mime.text import MIMEText
	smtp = smtplib.SMTP('smtp.gmail.com:587')
	smtp.starttls()
	smtp.login('rp462@cornell.edu', 'HalphaTycho1289')
	msg = MIMEText('This is the body of the message.')
	#msg['To'] = email.utils.formataddr(( 'rp462@cornell.edu'))
	#msg['From'] = email.utils.formataddr(('rp462@cornell.edu'))
	msg['Subject'] = subject_my
	smtp.set_debuglevel(True) # show communication with the server
	try:
		smtp.sendmail('pavesiriccardo3@gmail.com', ['pavesiriccardo3@gmail.com'], msg.as_string())
	finally:
		smtp.quit()




if True or (std_spat==2.66 and std_freq==4.95/2):
	emailme('starting '+newdir)




subprocess.call(['nohup','python', 'inject.py',newdir,str(flux),str(std_spat),str(std_freq),'&'])
subprocess.call(['nohup','python', 'MF_and_clump.py',newdir,'&'])
subprocess.call(['nohup','python', 'flux_extract.py',newdir,'&'])
subprocess.call(['nohup','python', 'analyze.py',newdir,str(flux),str(std_spat),str(std_freq),'&'])
os.system('rm '+newdir+'signal*fits')
os.system('rm '+newdir+'SNR*.fits')





