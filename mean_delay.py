import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy import signal
import random
import pandas as pd
import os

#************************
#	parameters 	#
#************************

#Parameters to be change:
lamb = 0.5500	# wavelength in micro-m
Mbh  = 1.3e8 	# mass of the black hole in solar masses
LLE  = 0.1	# L/LE Eddintong ratio
n    = 0.1	# radiative efficiency
zs   = 0.658	# redshift of the source
alf  = 6	# 6 = Schwarszschild black hole or 1 = Kerr black hole
pix  = 2001	# resolution of the kernel in pixels (always odd numbers)
incd = 30.0	# inclination in degrees
norm = 6.12e13	# pixel scale in cm 
size = np.asarray([0.5,1.0,2.0]) # size of the disk in times R0 (0.5*R0)
posan= np.asarray([0.0,45.0,90.0]) # position angles to be tested


#********************************
#	Magnification map	#
#********************************


in1 = raw_input('Do you want to create a magnification maps using Gerlumph? y/n: ')
if in1 == 'y':
	print('james include the code =p')
else:
	in2 = raw_input('Are you going to analyzed 2 or 4 magnification maps? ')
	in3 = raw_input('do you want to calculate the delay for the different disk and magnification maps y/n ')
	if in3 == 'y':
		if in2 == str(2):
			mapA = pyfits.open(str(raw_input('name for image A map: ')))[0].data
			mapB = pyfits.open(str(raw_input('name for image B map: ')))[0].data
			mapl = np.array(['mapA','mapB'])
			print('coffee? this is going to take "some" minutes =)')
		else:
			mapA = pyfits.open(str(raw_input('name for image A map: ')))[0].data
			mapB = pyfits.open(str(raw_input('name for image B map: ')))[0].data
			mapC = pyfits.open(str(raw_input('name for image C map: ')))[0].data
			mapD = pyfits.open(str(raw_input('name for image D map: ')))[0].data
			mapl = np.array(['mapA','mapB','mapC','mapD'])
			print('coffee? this is going to take "some" minutes =)')

#********************************
#	Accretion disk model	#
#********************************

def sb2d(lamb,Mbh,LLE,n,zs,alf,lmin,lmax,step,pa,inc,norm,size):
	R0 = ( 9.7e15 * ((1.0/(1.0+zs))*(lamb))**(4.0/3.0) * (Mbh/10**9)**(2.0/3.0) * (LLE/n)**(1.0/3.0) )
	Rcod = size*R0
	R0p = ( Rcod/norm )
#	G = 6.67e-11 #m**3 kg**-1 s-2
#	c = 299792458 #m/s
#	Rin = ((alf * G * Mbh*1.989e30)/c**2) * 100 #1msun = 1.989e30kg
#	Rinp = R0p/100.0
	Rin = 0.0
	Rinp = Rin/norm
	rx = np.arange(lmin,lmax,step)
	ry = rx[:,np.newaxis]
	rxn = rx*np.cos(pa) - ry*np.sin(pa) #rotation matrix
	ryn = rx*np.sin(pa) + ry*np.cos(pa)
	R = np.sqrt( (rxn**2)/np.cos(inc)**2 + (ryn**2))
	Rinc = rxn*np.tan(inc)
	xi = (R/R0p)**(3./4.) #* ( 1.0 - np.sqrt(Rinp/R)**(-1.0/4.0) )
	G = ( ( xi * np.exp(xi) )/( np.exp(xi) - 1.0 )**2 ) #G(xi) eq9 T&K18
	return G, Rcod, Rinc, R

#conversions
lmin = -(pix/2.0) 		# centerx = totalpixel/2.0
lmax = (pix/2.0)  		# centery = totalpixel/2.0
inc  = (incd*np.pi)/180.0	# inclination of the disk in radians
c = 299792458*100 		#c in meters
days = (60*60*24) 		#to transform to days
step = 1.0			# steps for the array (usually 1.0)

#************************
#	time delay	#
#************************
dtsel0 = []
dtincsel0 = []
dtlampsel0 = []
if in3 == 'y':
	index = random.sample( range(0,len(mapA.flatten())), 300000 )	#300000 random points to produced the cumulative curve
	for k in range(len(mapl)):
		slist = []	# these 3 empty list will helpful to create the final dataframe list
		palist = []
		inlist = []
		mmap = eval(str(mapl[k]))	#magnification map
#		os.remove(str(mapl[k])+'_val.dat')	
		res = open(str(mapl[k])+'_val.dat',"a")
		res.write('size'+','+'inc'+','+'pa'+','+'dt_m'+','+'dt_dev'+','+'dt_lamp'+','+'dtinc_m'+','+'dtinc_dev'+','+'dtlamp_m'+','+'dtlamp_dev'+'\n')
		for i in range(len(size)):
			for j in range(len(posan)):
				slist.append(size[i])
				palist.append(posan[j])
				inlist.append(incd)
				print('calculating delay for a disk with size',size[i],'*R0, PA:',posan[j],'and convolving with: ', mapl[k])
				G , R0, Rinc, Rlamp = sb2d(lamb,Mbh,LLE,n,zs,alf,lmin,lmax,step,(posan[j]*np.pi)/180.0,inc,norm,size[i])
				R = Rlamp-Rinc
				GRG = (np.mean(R*G)/np.mean(G))/(R0/norm) #<RG>/<G>
				mapG = signal.fftconvolve(mmap,G,mode='same')
				mapRin = signal.fftconvolve(mmap,Rinc*G,mode='same')
				mapRlamp = signal.fftconvolve(mmap,Rlamp*G,mode='same')
				lamp_dt = ( (GRG * (1.0+zs) * R0)/(c) ) / days
				divinc = mapRin/mapG
				divlamp = mapRlamp/mapG
				dtinc = ( ( ( divinc * (1.0 + zs) )/(c/norm) )/days)
				dtlamp = ( ( ( divlamp * (1.0 + zs) )/(c/norm) )/days) - lamp_dt
				dt = dtlamp - dtinc	
				print("saving time delay maps as fits")
				pyfits.writeto('dt'+str(mapl[k])+str(size[i])+'_'+str(posan[j])+'_'+str(incd)+'.fits',dt,clobber=True)
				pyfits.writeto('dtinc'+str(mapl[k])+str(size[i])+'_'+str(posan[j])+'_'+str(incd)+'.fits',dtinc,clobber=True)
				pyfits.writeto('dtlamp'+str(mapl[k])+str(size[i])+'_'+str(posan[j])+'_'+str(incd)+'.fits',dtlamp,clobber=True)
				print("saving 300000 random points")
				dtsel = dt.flatten()[index]
				dtincsel = dtinc.flatten()[index]
				dtlampsel = dtlamp.flatten()[index]
				dtsel0.append(dtsel)
				dtincsel0.append(dtincsel)
				dtlampsel0.append(dtlampsel)
				meddt = np.mean(dt.flatten())
				devdt = np.std(dt.flatten())
				meddtinc = np.mean(dtinc.flatten())
				devdtinc = np.std(dtinc.flatten())
				meddtlamp  = np.mean(dtlamp.flatten())
				devdtlamp  = np.std(dtlamp.flatten())
				print('mean value all of dt:',meddt)
				print('std:', devdt)
				res.write(str(size[i])+','+str(incd)+','+str(posan[j])+','+str(meddt)+','+str(devdt)+','+str(lamp_dt)+','+str(meddtinc)+','+str(devdtinc)+','+str(meddtlamp)+','+str(devdtlamp)+'\n')
		res.close()		
		np.savetxt('dtsel',dtsel0)
		np.savetxt('dtincsel',dtincsel0)
		np.savetxt('dtlampsel',dtlampsel0)

#****************************************
#	cumulative distribution		#
#****************************************

in4 = raw_input('The cumulative plot only work if you saved the time delay maps, choose 3 different sizes for the accretion disk, 3 position angle, and one inclination. Do you want skip this part? y/n')
if in4 == 'n':
	if in3 == 'n':
		dtsel0 = np.loadtxt('dtsel')
		dtincsel0 = np.loadtxt('dtincsel')
		dtlampsel0 = np.loadtxt('dtlampsel')
	n = 15
	bins = np.arange(-10.0,10.0,0.001)
	ig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(ax10,ax11,ax12)) = plt.subplots(4,3,sharex='col',sharey=True,figsize=(20, 10))
	color = ['black','blue','red']
	ax1.hist([dtincsel0[0],dtincsel0[1],dtincsel0[2]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax1.hist([dtincsel0[3],dtincsel0[4],dtincsel0[5]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax1.hist([dtincsel0[6],dtincsel0[7],dtincsel0[8]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax2.hist([dtlampsel0[0],dtlampsel0[1],dtlampsel0[2]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax2.hist([dtlampsel0[3],dtlampsel0[4],dtlampsel0[5]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax2.hist([dtlampsel0[6],dtlampsel0[7],dtlampsel0[8]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax3.hist([dtsel0[0],dtsel0[1],dtsel0[2]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax3.hist([dtsel0[3],dtsel0[4],dtsel0[5]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax3.hist([dtsel0[6],dtsel0[7],dtsel0[8]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax4.hist([dtincsel0[9],dtincsel0[10],dtincsel0[11]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax4.hist([dtincsel0[12],dtincsel0[13],dtincsel0[14]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax4.hist([dtincsel0[15],dtincsel0[16],dtincsel0[17]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax5.hist([dtlampsel0[9],dtlampsel0[10],dtlampsel0[11]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax5.hist([dtlampsel0[12],dtlampsel0[13],dtlampsel0[14]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax5.hist([dtlampsel0[15],dtlampsel0[16],dtlampsel0[17]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax6.hist([dtsel0[9],dtsel0[10],dtsel0[11]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax6.hist([dtsel0[12],dtsel0[13],dtsel0[14]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax6.hist([dtsel0[15],dtsel0[16],dtsel0[17]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax7.hist([dtincsel0[18],dtincsel0[19],dtincsel0[20]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax7.hist([dtincsel0[21],dtincsel0[22],dtincsel0[23]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax7.hist([dtincsel0[24],dtincsel0[25],dtincsel0[26]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax8.hist([dtlampsel0[18],dtlampsel0[19],dtlampsel0[20]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax8.hist([dtlampsel0[21],dtlampsel0[22],dtlampsel0[23]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax8.hist([dtlampsel0[24],dtlampsel0[25],dtlampsel0[26]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax9.hist([dtsel0[18],dtsel0[19],dtsel0[20]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax9.hist([dtsel0[21],dtsel0[22],dtsel0[23]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax9.hist([dtsel0[24],dtsel0[25],dtsel0[26]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax10.hist([dtincsel0[27],dtincsel0[28],dtincsel0[29]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax10.hist([dtincsel0[30],dtincsel0[31],dtincsel0[32]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax10.hist([dtincsel0[33],dtincsel0[34],dtincsel0[35]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax11.hist([dtlampsel0[27],dtlampsel0[28],dtlampsel0[29]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax11.hist([dtlampsel0[30],dtlampsel0[31],dtlampsel0[32]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax11.hist([dtlampsel0[33],dtlampsel0[34],dtlampsel0[35]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax12.hist([dtsel0[27],dtsel0[28],dtsel0[29]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle=':')
	ax12.hist([dtsel0[30],dtsel0[31],dtsel0[32]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='-')
	ax12.hist([dtsel0[33],dtsel0[34],dtsel0[35]],cumulative='True',bins=bins,histtype='step',normed=1,color=color,linestyle='--')
	ax10.set_xlabel("dt( x*sin(i) ) [days]")
	ax1.set_xlim([-1.5,1.5])
	ax1.set_ylim([0.0,1.1])
	ax1.set_xticks([-1.0,-1.0,0.0,1.0])
	ax11.set_xlabel("dt(R) [days]")
	ax2.set_xlim([-1.5,3.0])
	ax2.set_xticks([-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
	ax12.set_xlabel("dt(x*sin(i) - R) [days]")
	ax3.set_xlim([-1.5,3.0])
	ax3.set_xticks([-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])
	plt.savefig('cumulative.pdf')
	plt.show()

#************************
#	table		#
#************************

in5 = raw_input('Do you want to compute a final table with the differences of the mean delay between the images? y/n ')
if in3 == 'n':
	slist = []	# these 3 empty list will helpful to create the final dataframe list
	palist = []
	inlist = []
	for i in range(len(size)):
		for j in range(len(posan)):
			slist.append(size[i])
			palist.append(posan[j])
			inlist.append(incd)
final = pd.DataFrame({'size':slist,'inc':inlist,'pa':palist})
if in5 == 'y':
	if in2 == str(2):
		mapA = pd.read_csv('mapA_val.dat')
		mapB = pd.read_csv('mapB_val.dat')
		final['B-A (m)'] = mapB['dt_m'] - mapA['dt_m']
		final['B-A (dev)'] = np.sqrt((mapB['dt_dev'])**2 + (mapA['dt_dev'])**2)
	else:
		mapA = pd.read_csv('mapA_val.dat')
		mapB = pd.read_csv('mapB_val.dat')
		mapC = pd.read_csv('mapC_val.dat')
		mapD = pd.read_csv('mapD_val.dat')
		final['B-A (m)'] = mapB['dt_m'] - mapA['dt_m']
		final['B-A (dev)'] = np.sqrt((mapB['dt_dev'])**2 + (mapA['dt_dev'])**2)
		final['C-A (m)'] = mapC['dt_m'] - mapA['dt_m']
		final['C-A (dev)'] = np.sqrt((mapC['dt_dev'])**2 + (mapA['dt_dev'])**2)
		final['D-A (m)'] = mapD['dt_m'] - mapA['dt_m']
		final['D-A (dev)'] = np.sqrt((mapD['dt_dev'])**2 + (mapA['dt_dev'])**2)
		final['C-B (m)'] = mapC['dt_m'] - mapB['dt_m']
		final['C-B (dev)'] = np.sqrt((mapC['dt_dev'])**2 + (mapB['dt_dev'])**2)
		final['D-B (m)'] = mapD['dt_m'] - mapB['dt_m']
		final['D-B (dev)'] = np.sqrt((mapD['dt_dev'])**2 + (mapB['dt_dev'])**2)
		final['D-C (m)'] = mapD['dt_m'] - mapC['dt_m']
		final['D-C (dev)'] = np.sqrt((mapD['dt_dev'])**2 + (mapC['dt_dev'])**2)
	final.to_csv('finaltable')
	print('Final Table:')
	print(final)
