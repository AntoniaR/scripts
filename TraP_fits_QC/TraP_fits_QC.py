import sys
import os
import pyfits
import numpy as np
import glob
import getRMS
import tools

hdu=0
if len(sys.argv) != 6:
    print 'python TraP_fits_QC.py <sigmaClip> <sigmaRej> <f> <frequencies> <theory>'
    exit()
sigmaClip = float(sys.argv[1])
sigmaRej = float(sys.argv[2])
f = float(sys.argv[3])
allFreqs=sys.argv[4]
RejTheory=sys.argv[5]

freqs=[]
imageData=[]
for filename in glob.glob('*.fits'):
    data = getRMS.read_data(pyfits.open(filename)[hdu], plane=None)
    hdulist=pyfits.open(filename)
    prihdr=hdulist[0].header
    obsrms=getRMS.rms_with_clipped_subregion(data, sigmaClip, f)
    if 'THRY_RMS' in prihdr.keys():
        theory=prihdr['THRY_RMS']
    else:
        theory=obsrms
    freq=int((float(prihdr['CRVAL3'])/1e6))
    imageData.append([filename, freq, obsrms*1000., obsrms/theory])
    if freq not in freqs:
        freqs.append(freq)

thresholds={}
thresholds2={}

if allFreqs=='T':
    for frequency in freqs:
        noise_avg_log, noise_scatter_log, noise_threshold_log = tools.fit_hist([np.log10(n[2]) for n in imageData if n[1]==frequency], sigmaRej, r'Observed RMS (mJy)', 'rms', frequency)
        noise_avg=10.**(noise_avg_log)
        noise_max=10.**(noise_avg_log+noise_scatter_log)-10.**(noise_avg_log)
        noise_min=10.**(noise_avg_log)-10.**(noise_avg_log-noise_scatter_log)
        print 'Average RMS Noise in images (1 sigma range, frequency='+str(frequency)+' MHz): '+str(frequency)+' MHz): '+str(noise_avg)+' (+'+str(noise_max)+',-'+str(noise_min)+') mJy'
        thresholds[frequency]=([noise_avg,noise_max*sigmaRej,noise_min*sigmaRej*-1.])

        noiserat_avg_log, noiserat_scatter_log, noiserat_threshold_log = tools.fit_hist([np.log10(n[3]) for n in imageData if n[1]==frequency], sigmaRej, r'RMS / Theory', 'ratio', frequency)
        noiserat_avg=10.**(noiserat_avg_log)
        noiserat_max=10.**(noiserat_avg_log+noiserat_scatter_log)-10.**(noiserat_avg_log)
        noiserat_min=10.**(noiserat_avg_log)-10.**(noiserat_avg_log-noiserat_scatter_log)
        print 'Average RMS Noise ratio in images (1 sigma range, frequency='+str(frequency)+' MHz): '+str(frequency)+' MHz): '+str(noiserat_avg)+' (+'+str(noiserat_max)+',-'+str(noiserat_min)+') mJy'
        thresholds2[frequency]=([noiserat_avg,noiserat_max*sigmaRej,noiserat_min*sigmaRej*-1.])
        

frequency='all'
TMPdata = np.array([np.log10(n[2]) for n in imageData])
TMPdata = TMPdata[np.isfinite(TMPdata)]
noise_avg_log, noise_scatter_log, noise_threshold_log = tools.fit_hist(TMPdata, sigmaRej, r'Observed RMS (mJy)', 'rms', frequency)
noise_avg=10.**(noise_avg_log)
noise_max=10.**(noise_avg_log+noise_scatter_log)-10.**(noise_avg_log)
noise_min=10.**(noise_avg_log)-10.**(noise_avg_log-noise_scatter_log)
print 'Average RMS Noise in images (1 sigma range, frequency='+str(frequency)+' MHz): '+str(noise_avg)+' (+'+str(noise_max)+',-'+str(noise_min)+') mJy'
thresholds[frequency]=[noise_avg,noise_max*sigmaRej,noise_min*sigmaRej*-1.]

TMPdata2 = np.array([np.log10(n[3]) for n in imageData])
TMPdata2 = TMPdata2[np.isfinite(TMPdata2)]
noiserat_avg_log, noiserat_scatter_log, noiserat_threshold_log = tools.fit_hist(TMPdata2, sigmaRej, r'RMS / Theory', 'ratio', frequency)
noiserat_avg=10.**(noiserat_avg_log)
noiserat_max=10.**(noiserat_avg_log+noiserat_scatter_log)-10.**(noiserat_avg_log)
noiserat_min=10.**(noiserat_avg_log)-10.**(noiserat_avg_log-noiserat_scatter_log)
print 'Average RMS Noise ratio in images (1 sigma range, frequency='+str(frequency)+' MHz): '+str(frequency)+' MHz): '+str(noiserat_avg)+' (+'+str(noiserat_max)+',-'+str(noiserat_min)+') mJy'
thresholds2[frequency]=([noiserat_avg,noiserat_max*sigmaRej,noiserat_min*sigmaRej*-1.])


goodImg=[]
for image in imageData:
    if RejTheory=='T':
        if allFreqs=='F':
            if image[3]<thresholds2['all'][0]+thresholds2['all'][1] and image[3]>thresholds2['all'][0]+thresholds2['all'][2]:
                goodImg.append(image[0])
            else:
                print 'Bad image:',image[0],image[1],image[2],image[3]
        else:
            if image[3]<thresholds2[image[1]][0]+thresholds2[image[1]][1] and image[3]>thresholds2[image[1]][0]+thresholds2[image[1]][2]:
                goodImg.append(image[0])
            else:
                print 'Bad image:',image[0],image[1],image[2],image[3]
    else:
        if allFreqs=='F':
            if image[2]<thresholds['all'][0]+thresholds['all'][1] and image[2]>thresholds['all'][0]+thresholds['all'][2]:
                goodImg.append(image[0])
            else:
                print 'Bad image:',image[0],image[1],image[2],image[3]
        else:
            if image[2]<thresholds[image[1]][0]+thresholds[image[1]][1] and image[2]>thresholds[image[1]][0]+thresholds[image[1]][2]:
                goodImg.append(image[0])
            else:
                print 'Bad image:',image[0],image[1],image[2],image[3]
        

goodImg = [os.getcwd()+'/'+x for x in goodImg]
f = open('images_to_process.py', 'w')
f.write('images = ')
f.write('[')
f.writelines(["'%s',\n" % good for good in goodImg])
f.write(']\n')
f.write('''#Just for show:
print "******** IMAGES: ********"
for f in images:
    print f
print "*************************"
''')
f.close()
