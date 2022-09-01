import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML
from astropy.io import fits
from astropy import modeling
from astropy.wcs import WCS
from astropy.time import Time
import astropy.units as u
import scipy.ndimage as nd
from scipy.io import readsav
from scipy.interpolate import griddata
import numpy as np
import glob
import math
import warnings

warnings.filterwarnings("ignore")


#routine to match files based on time
#assumes files are organized by date in path
#reads in files on day by day bases
#even if only passing in a single date, date must be a string array i.e. ['20220401']
#returns list of heaaders, imgs as well as an Nx4 array of the indices of headers/images to form mosaics
#if prep is set to true, solohi_prep routine will be applied to images (ONLY USE ON L1)
def match_files(dates, path='', prep=False):
 flist=[]
 #search for files
 for i in dates:
     if len(glob.glob(path+i)) > 0:
        flist=flist+glob.glob(path+i+'/*.fits')

 #sort fits files based on time
 x=flist[0].find('.fits')
 flist=sorted(flist, key=lambda s : s[x-19:] )
 n=len(flist)
 #create ouptut arrays
 imgs=np.zeros([n,1024,1024], dtype=np.float16)
 inx=0   
 hlist=[]         
 for i in flist:
     #read files
     img, hdr = fits.getdata(i, header=True)
     pxs=hdr['pxend2']-hdr['pxbeg2']+(hdr['pxend1']-hdr['pxbeg1'])+2
     #remove any non-fullframe files
     if pxs==3968:
       hlist.append(hdr)
       s=img.shape
       #prep file if needed
       if prep == True:
           img, hdr = solohi_prep(hdr, img=img)
       else:
           img=img/1e-12
       imgs[inx, 0:s[0],0:s[1]]=img
       inx+=1
 dets=[]
 times=[]

 n=len(hlist)

 #get time and detector for each file
 for i in range(n):
    dets.append(hlist[i]['detector'])
    times.append(hlist[i]['date-avg'])
  
 w1 = [i for i, j in enumerate(dets) if j == '1']  
 w2 = [i for i, j in enumerate(dets) if j == '2']
 w3 = [i for i, j in enumerate(dets) if j == '3']
 w4 = [i for i, j in enumerate(dets) if j == '4']  
 
 
 dm, loc = min((dm, loc) for (loc, dm) in enumerate([len(w1),len(w2),len(w3),len(w4)]))
 
 mos_inx=np.zeros([dm, 4])
 t=Time(times)
 #round times to nearest 20 minutes (gives better performance)
 mj=t.mjd
 mj=np.around((mj-mj[0])*1200)
   
 if loc == 0:
     wbase=w1
 elif loc == 1:
     wbase=w2
 elif loc == 2:
     wbase=w3
 else:
     wbase=w4
 #find closest file to "base" file
 for i in range(dm):
     bt=mj[wbase[i]]
     d1, l1 = min((d1, l1) for (l1, d1) in enumerate(abs(mj[w1]-bt)))
     d2, l2 = min((d2, l2) for (l2, d2) in enumerate(abs(mj[w2]-bt)))    
     d3, l3 = min((d3, l3) for (l3, d3) in enumerate(abs(mj[w3]-bt)))     
     d4, l4 = min((d4, l4) for (l4, d4) in enumerate(abs(mj[w4]-bt)))
     mos_inx[i,:]=[w1[l1],w2[l2],w3[l3],w4[l4]]    
     
 return hlist, imgs[0:n,:,:], mos_inx.astype(int)
         
#combine files from each detecctor into mosaic
#filenames can either be string array of fits files or list of headers
#if filenames is list of headers, pass image array with imgs
#if filenames is list of fits files, specify path with path
#returns an output mosaic and heaader
#use prep if imgs/fileames are unprocessed Level 1 data
#time is used to specify output time of output
def shi_mosaic(filenames, prep=False, path='', imgs='', time=''):
 
 #determine how to read in data
 if imgs == '' and prep==False:
     img1, hdr1 = fits.getdata(path+filenames[0], header=True)
     img2, hdr2 = fits.getdata(path+filenames[1], header=True)
     img3, hdr3 = fits.getdata(path+filenames[2], header=True)
     img4, hdr4 = fits.getdata(path+filenames[3], header=True)

     img1=img1.astype(float)
     img2=img2.astype(float)
     img3=img3.astype(float)
     img4=img4.astype(float)
 elif prep == False:
     img1=imgs[0,0:960,0:1024].astype(float)
     img2=imgs[1,0:1024,0:960].astype(float)
     img3=imgs[2,0:960,0:1024].astype(float)
     img4=imgs[3,0:1024,0:960].astype(float)
     hdr1=filenames[0]
     hdr2=filenames[1]
     hdr3=filenames[2]
     hdr4=filenames[3]
     
 if prep==True:
     img1, hdr1 = solohi_prep(path+filenames[0])
     img2, hdr2 = solohi_prep(path+filenames[1])
     img3, hdr3 = solohi_prep(path+filenames[2])
     img4, hdr4 = solohi_prep(path+filenames[3])     
 #create output image
 outimg=np.zeros([2072, 2008], dtype=np.float16, order='C')
 
 #cut extra space from input images 
 outimg[0:1024, 0:960]=img4
 outimg[1048:2008, 0:1024]=img3
 outimg[64:1024, 984:2008]=img1
 outimg[1048:2072, 1048:2008]=img2

 #These rows have no data but can have values introduced by processing, better to zero them out
 outimg[0:64,:]=0
 outimg[2008:2072,:]=0

 #set header up correctly
 outhdr=hdr4.copy()
 #print(time)
 outhdr['naxis1']=2008
 outhdr['naxus2']=2072
 outhdr['detector']='M'
 outhdr['date-obs']=time
 outhdr['date-avg']=time

 #return output image and header
 return outimg, outhdr

#preps solohi data to 'level 2'
#hdr can be a file name or an L1 header
#if hdr is L1 header, img should be correspoonding L1 image
#DO NOT USE IF NOT NECESSARY
def solohi_prep(hdr, img=''):
    
    if img=='':
        img, hdr = fits.getdata(hdr, header=True)
        
    img=img.astype(float)
    
    #each tile has two sets of unique models, for the bias removal and linearity correction
    #check for the specific tile and read in appropriate models
    #for tiles 2 and 4, bias images must also be rotated to match image dimensions
    if int(hdr['detector']) == 1:  
        bcs=readsav('shi_cal/newbmd1.sav')
        ics=readsav('shi_cal/d1lincor.sav')
        lfs=bcs['lfs']
        corfunc=ics['gf2']
        corfunc=np.append(corfunc, np.ones(400)*corfunc[15999])
        brow=np.flip(lfs[0]+lfs[1]*hdr['detect_t'])
        bimage=np.tile(brow, (960, 1))
    elif int(hdr['detector']) == 2:
        bcs=readsav('shi_cal/newbmd2.sav')
        ics=readsav('shi_cal/d2lincor.sav')
        lfs=bcs['lfs']
        corfunc=ics['gf2']
        corfunc=np.append(corfunc, np.ones(400)*corfunc[15999])
        brow=np.flip(lfs[0]+lfs[1]*hdr['detect_t'])
        bimage=np.tile(brow, (960, 1)).transpose()
    elif int(hdr['detector']) == 3:
        bcs=readsav('shi_cal/newbmd3.sav')
        ics=readsav('shi_cal/d3lincor.sav')
        lfs=bcs['lfs']
        corfunc=ics['gf2']
        corfunc=np.append(corfunc, np.ones(400)*corfunc[15999])        
        brow=lfs[0]+lfs[1]*hdr['detect_t']
        bimage=np.tile(brow, (960, 1))
    elif int(hdr['detector']) == 4:
        bcs=readsav('shi_cal/newbmd4.sav')
        ics=readsav('shi_cal/d4lincor.sav')
        lfs=bcs['lfs']
        corfunc=ics['gf2']
        corfunc=np.append(corfunc, np.ones(400)*corfunc[15999])
        brow=lfs[0]+lfs[1]*hdr['detect_t']
        bimage=np.tile(brow, (960, 1)).transpose()
    
    #account for truncation of images
    img=img*2**hdr['iptrunc']

    #subtract bias image and linearity correction. these corrections were determined for unbinned, single exposure
    #the data must therefore account for binning and sums
    img=img-bimage*(hdr['nbin']*hdr['nsumexp'])    
    imcor=corfunc[(img/(hdr['nbin']*hdr['nsumexp'])).astype(int)]
    img=img-imcor*(hdr['nbin']*hdr['nsumexp'])
    
    #Divid be exposure time and binning to get DN/s
    img=img/(hdr['xposure']*hdr['nbin'])
        
    return img, hdr

#use output from match_files and shi_mosaic to create a combined, mosaic data array
#hlist, imgs and mos_inx are outputs directly from match_files
#data can be prepped
def shi_mov_cube(hlist, imgs, mos_inx, prep=False):
    #create outputs
    n=mos_inx.shape
    
    proc_img=np.zeros([n[0], 2072, 2008], dtype=np.float16)
    mos_times=[]
    hdrs=[]
    for i in range(n[0]):
        #determine 'average time' of the full mosaic
        hs=[]
        ts=[]
        ims=np.zeros([4,1024,1024], dtype=np.float16)
        for j in range (4):
            hs.append(hlist[mos_inx[i, j]])
            #print(hs[j]['date-avg'])
            ts.append(hs[j]['date-avg'])
           # print(mos_inx[i,j], hs[j]['date-avg'])
            ims[j,:,:]=imgs[mos_inx[i,j],:,:]
        tm=Time(ts)
        tavg=(tm.max().mjd-tm.min().mjd)/2.+tm.min().mjd
        tat=Time(tavg, format='mjd')
        mos_times.append(tat.isot)
      #  print(tm, tavg, tat.isot, hs[0]['filename'], hs[1]['filename'], hs[2]['filename'], hs[3]['filename'])
        #combine the individual images into the output array
        proc_img[i,:,:], h=shi_mosaic(hs, imgs=ims, prep=prep, time=tat.isot)
        hdrs.append(h)
        
    return proc_img, hdrs

#perform some basic image processing on the data cube
#NOTE: NOT ALL OF THESE CAN BE COMBINED 
#rdiff: set to true for running difference
#ratio: set to true for running ratio
#mednorm: set to true to divide each image by it's median value
#median: takes the median of each pixel in the data cube, creating a 2d 'background' to remove
#smooth: if set to number other than 0, performs a smoothing over that number of pixels.
def process_cube(img_cube, rdiff=False, ratio=False, median=False, smooth=0, mednorm=False):
    s=img_cube.shape
    if mednorm == True:
        for i in range(s[0]):
            img_cube[i,:,:]=img_cube[i,:,:]/np.nanmedian(img_cube[i,:,:])
    if rdiff == True:
        out_img = np.zeros([s[0]-1, s[1], s[2]], dtype=np.float16)

        for i in range(s[0]-1):
            out_img[i,:,:]=img_cube[i+1,:,:]-img_cube[i,:,:]
    if ratio == True:
        out_img = np.zeros([s[0]-1, s[1], s[2]], dtype=np.float16)
        
        for i in range(s[0]-1):
            out_img[i,:,:]=img_cube[i+1,:,:]/img_cube[i,:,:]
    if median == True:
        med=np.median(img_cube, 0)
        out_img=img_cube/med[np.newaxis,:,:]
    #Notice, smoothing is done on each tile to keep the boundary regions from causing problems    
    if smooth != 0:
        for i in range(out_img.shape[0]):
            out_img[i,64:1023, 5:959]=nd.uniform_filter(out_img[i,64:1023, 5:959].astype(float), size=smooth, mode='nearest')
            out_img[i,1050:2002, 5:1023]=nd.uniform_filter(out_img[i,1050:2002, 5:1023].astype(float), size=smooth, mode='nearest')
            out_img[i,69:1023, 983:2002]=nd.uniform_filter(out_img[i,69:1023, 983:2002].astype(float), size=smooth, mode='nearest')
            out_img[i,1047:2002, 1047:2002]=nd.uniform_filter(out_img[i,1047:2002, 1047:2002].astype(float), size=smooth, mode='nearest')
        out_img[~np.isfinite(out_img)]=np.nan
    out_img[np.where(out_img==0)]=np.nan    
    return out_img
                                                              
#display and potentially save output movies. 
#data: can be cube from shi_mov_cube or process_cube
#hdrs: needed if you want to display timestamps on movie or calculate wcs information 
#vrange: specify min/max values for byte scaling the display, if not selected this will be automatically determined
#writefits: if set to true, individual output fits files that form the movie                                                     
#scale: if set to true, the output fits files will be bytescaled
#savepath: place to write fits files
#qcheck: if set to true program will look for files with at least 30% of pixels outside bytescaling range and romethem
#grid: Add grid onto the movie
#wcscor: If combined with grid=true, WCS HPC grid is put on output
#redrawwcs: If false (default) the WCS grid from the first image is used for all. If true the WCS grid from each header will be redrawn
#proj: If combined with grid=true and input data is projected onto HPC grid, will put the grid on the output images
#elongation & latitude: if proj=true and grid=true, these need to be the extent of the coordinates used to calculate images
#moviepath: path to save the output movie file
#moviename: animated movie output
def run_movie(data, vrange=[0,0], hdrs='', writefits=False, savepath='', scale=False, grid=False, qcheck=False, wcscor=False, redrawwcs=False,proj=False, elongation=[-50,0], latitude=[-25,25], moviepath='', moviename='test.mp4'):
    s=data.shape
    #remove nan values
    data[~np.isfinite(data)]=1
    inx=0

    #if needed, establish bytescale
    if vrange[0]==vrange[1]:
        vrange=scale_cube(data)
    
    #if desired by user, check for excessive pixels in each individual frame for bytescaling
    if qcheck == True:
      hdrs2=[]
      print('Checking individual frames for quality...')
      for i in range(s[0]):    
       ba=np.logical_and(data[i,:,:] <= vrange[1], data[i,:,:] >= vrange[0])
       cnt=np.count_nonzero(np.where(ba)[0])
       if cnt/(s[1]*s[2]) > .7:
                if inx==0:
                    data2=data[i,:,:]
                    hdrs2.append(hdrs[i])
                    inx+=1
                else:
                    data2=np.dstack((data2, data[i,:,:]))
                    hdrs2.append(hdrs[i])
                    inx+=1
    
      data=np.transpose(data2, (2,0,1))
      hdrs=hdrs2  
      print(inx, ' of ', s[0], ' frames used')
        
    #if desired, write processed output fits files    
    if writefits == True:
        print('Writing FITS files...')
        for i in range(s[0]):
               times=hdrs[i]['date-avg']
               x=hdrs[i]['filename'].find('solohi-')
               pfilename=hdrs[i]['filename'][0:x+7]+'mft_'+times[0:4]+times[5:7]+times[8:13]+times[14:16]+times[17:19]+'_V00.fits'
               hdrs[i]['filename']=pfilename
               outdat=data[i,:,:].copy()
               if scale == True:
                       outdat=(outdat-vrange[0])/(vrange[1]-vrange[0])*255
                       outdat[outdat < 0]=0
                       outdat[outdat > 255]=255
               htmp=fits.PrimaryHDU(outdat.astype(float), header=hdrs[i])
               htmp.writeto(savepath+pfilename, overwrite=True)
    #set up the animation window
    print('Launching animation...')
    fig=plt.figure()
    ax=plt.subplot()
    if grid==False:
       im=ax.imshow(data[0,:,:], origin='lower', vmin=vrange[0], vmax=vrange[1], cmap='gray')
       plt.xticks([])
       plt.yticks([])
       plt.axis('off')
    elif proj==True:
       im=ax.imshow(data[0,:,:], origin='lower', vmin=vrange[0], vmax=vrange[1], cmap='gray', extent=[elongation[0], elongation[1], latitude[0], latitude[1]])
       ax.set_xlabel('HPC Elongation (Deg)')
       ax.set_ylabel('HPC Latitude (Deg)') 
       ax.grid(color='white')
    elif wcscor==True:    
       wcs=WCS(hdrs[0])
       ax=plt.subplot(projection=wcs)
       im=ax.imshow(data[0,:,:], vmin=vrange[0], vmax=vrange[1], cmap='gray', origin='lower')
       ax.grid(color='white')
       ax.coords[0].set_format_unit(u.deg)
       ax.coords[1].set_format_unit(u.deg)
    else:
       im=ax.imshow(data[0,:,:], origin='lower', vmin=vrange[0], vmax=vrange[1], cmap='gray')
       ax.grid(color='white')

    def init():
        if hdrs != '':
           times=hdrs[0]['date-avg']
           ax.set_title(times)
        else:
           im.set_title='Frame: 0'
        im.set_data(data[0,:,:])
        return im, ax

    def animate(i, im=im, ax=ax):
        if hdrs != '':
           print('Rendering frame ', i+1, 'of', s[0])
           times=hdrs[i]['date-avg']
           
           if grid==True and wcscor==True and redrawwcs==True:
              wcs=WCS(hdrs[i])
              ax=plt.subplot(projection=wcs)
              im=ax.imshow(data[i,:,:], vmin=vrange[0], vmax=vrange[1], cmap='gray', origin='lower')
              ax.grid(color='white')
              ax.coords[0].set_format_unit(u.deg)
              ax.coords[1].set_format_unit(u.deg)
           else:
              im.set_data(data[i,:,:])
           ax.set_title(times)
        else:
           ax.set_title('Frame: '+str(i))
    #animate cube and save the output. ffmpeg arugments can be changed as needed
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=data.shape[0], interval=10, blit=True)
    anim.save(moviepath+moviename, writer=animation.FFMpegWriter(fps=20, extra_args=["-crf", "25", "-s", "864x576", "-vcodec", "libx264"]), dpi=250)

#rudimentary auto-scale function based on gaussian modeling of histograms and standard deviations
def scale_cube(img_cube):
    print('Determining scale...')
    avg=np.nanmean(img_cube[0,:,:])
    var=np.nanstd(img_cube[0,:,:])
    minmax=[avg-var,avg+var]
    cnt, bins = np.histogram(img_cube, range=minmax, bins=20)

    while max(cnt/cnt.sum()) > .4:
        cnt, bins = np.histogram(img_cube, range=[bins[5], bins[-6]], bins=20)
    
    inx=(img_cube > bins[0]) & (img_cube < bins[-1])
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D(amplitude=.25, mean=np.nanmean(img_cube[inx]), stddev=np.nanstd(img_cube[inx]))
    fm=fitter(model, bins[1:], cnt/cnt.sum())    
    minmax=[fm.mean.value-fm.stddev.value*2, fm.mean.value+fm.stddev.value*2]
    print('Byte-scale range: ', minmax)
    return minmax

#projects a SoloHI image onto an HPC coordiante based grid
#display=True will plot the result
#the elongation and latitude ranges as well as the outsize can be speccified by the user
def proj_img(hdr, img, display=False, elongation=[-50,0], latitude=[-25,25], outsize=[1000,1000], vrange=[0,0]):

    #get WCS info from header
    wcs=WCS(hdr)
    
    #establish regular output grid
    x=np.linspace(elongation[0], elongation[1], outsize[0])
    y=np.linspace(latitude[0], latitude[1], outsize[1])
    x2, y2 = np.meshgrid(x, y)

    #convert get pixel coordinates of output grid
    xp, yp = wcs.all_world2pix(x2, y2, 0)

    #interpolate image onto the outut grid
    test=nd.map_coordinates(img.astype(float), [yp, xp], mode='constant', order=1)
    #if desired, plot results
    if display == True:
      if vrange[0]==vrange[1]:
           vrange=scale_cube(data)   
      ax2=plt.subplot()
      ax2.imshow(test, origin='lower', vmin=vrange[0], vmax=vrange[1], cmap='gray', extent=[x[0], x[999], y[0], y[999]])
      ax2.grid(color='white')
    return test

#creates a movie of project SoloHI images
#vrange can be establised or determined automatically
#the elongation and latitude ranges as well as the outsize can be speccified by the
def proj_movie(data, hdrs, vrange=[0,0], elongation=[-50,0], latitude=[-25,25], outsize=[1000,1000], moviepath='', moviename='test.mp4'):
    s=data.shape
    if vrange[0]==vrange[1]:
        vrange=scale_cube(data)
    #set up data cube
    outdat=np.zeros([s[0], outsize[0], outsize[1]], dtype=np.float16)
    #project images
    for i in range(s[0]):
        outdat[i,:,:]=proj_img(hdrs[i], data[i,:,:], elongation=elongation, latitude=latitude, outsize=outsize, vrange=vrange)
    #pass data to run_movie for display/saving
    run_movie(outdat, hdrs=hdrs, vrange=vrange, grid=True, proj=True, elongation=elongation, latitude=latitude, moviepath=moviepath, moviename=moviename)
    
#does the projection for an image into the polar coordinte system used to construct a jmap.
#if display is set to true the output will be displayed
#if cubic is set to true, use cubic interpolation instead of nearest (will be SLOW)
def jframe(hdr, img, display=False, vrange=[0,0], cubic=False):
    
    #use(proj_img to get the image into the HPC coordinates
    pimg=proj_img(hdr, img)
    pimg=np.nan_to_num(pimg)

    #establish the grid of HPC coordinattes
    x=np.linspace(-50,0,1000)
    y=np.linspace(-25,25,1000)
    x2, y2 = np.meshgrid(x, y)
    
    #convert the grid to polar coordinates
    rho=-np.sqrt(x2**2+y2**2)
    theta=(np.arctan2(x2,y2))*180/math.pi
    #get the output image coordinates
    gridx, gridy=np.mgrid[-5:-45:1000j, -45:-135:1000j]
    #interpolate from the original grid to the new
    if cubic != True:
        test=griddata((rho.flatten(), theta.flatten()), pimg.flatten(), (gridx.flatten(), gridy.flatten()), method='nearest')
    else: 
        test=griddata((rho.flatten(), theta.flatten()), pimg.flatten(), (gridx.flatten(), gridy.flatten()), method='cubic')

    gout=np.reshape(test, (1000,1000))
    if display == True:
      if vrange[0]==vrange[1]:
           vrange=scale_cube(data) 
      ax2=plt.subplot()
    #ax2=plt.plot(test)
      ax2.imshow(gout, origin='lower', vmin=vrange[0], vmax=vrange[1], cmap='gray')
    return gout


#generates a data cube using the jframe routine to turn into jmaps
#if cubic is set to true, use cubic interpolation instead of nearest (will be SLOW)
def jmovie(data, hdrs, cubic=False):
   s=data.shape
   jmov=np.zeros([s[0],1000,1000])
   for i in range(s[0]):
       print('Projecting frame: ', i+1, ' of ', s[0])
       jmov[i,:,:]=jframe(hdrs[i], data[i,:,:], cubic=cubic)
   print('Done')
   return jmov

#wrapper to construct the jamps
#jmov must be the output from jmovie
#pa can be any radial position angle from the Sun. Nearly all SoloHI data is off the west limb of the Sun, 
#so Postiion angles between 60-120 will return the best results
#the jmaps are construcuted by taking the median of the strip around the position angle pixel in each row
#the width of that strip is defined by width
#set the scaling with vrange, should probably be a little tighter than for a movie
#outx- the width of the jmap along the time axis. this can be played with to give the desired aspect ratio
#unit- default is degrees, but can also be set for Rsun, AU and KM. Will use the mean spacecraft distance and the
#sin projection to convert elongation into height
#the jmap will determine the starting location for each strip. the width of the output strip will cover the time
#from the current file until the next file
def make_jmap(jmov, hdrs, pa=90, width=10, vrange=[0,0], outx=2000, unit='Deg'):
    if vrange[0]==vrange[1]:
       vrange=scale_cube(data) 
    s=jmov.shape
    #determine the time span of the whole data cube
    t1=Time(hdrs[0]['date-avg']).mjd
    t2=Time(hdrs[s[0]-1]['date-avg']).mjd
    if outx/s[0] < 10:
        outx=s[0]*10
    dsuns=np.zeros(s[0])
    #convert time span to mintes
    tspan=(t2-t1)*24*60
    #get the pixels per minute in the x direction
    ppm=(outx/tspan)
    tcur=t1
    inx=0
    jimg=np.zeros([1000, outx])
    for i in range(s[0]):
        #determine pixel space between files
        if i < s[0]-1:
            tnext=Time(hdrs[i+1]['date-avg']).mjd
            pixw=int(((tnext-tcur)*24*60)*ppm)
            tcur=tnext
        else:
            pixw=outx-1-inx
        #get the spacecraft distance for each file    
        dsuns[i]=hdrs[i]['dsun_obs']
        #fill the output image with the signal from the current file at the right columns
        jimg[:,inx:inx+pixw]=np.transpose(np.tile(np.median(jmov[i,:,int((pa-45)/90*1000-width/2):int((pa-45)/90*1000+width/2)], axis=1), (pixw,1)))
        #increase the index of the current file
        inx+=pixw
    dsun_obs=np.mean(dsuns)
    titleu='m'
    #set the y-axis data of the current file
    extent=[0, tspan/60, 5, 45]
    #if unit is Rs, adjust dsun_obs anad convert elongations
    if unit=='Rs' or unit=='Rsun' or unit=='RS' or unit=='RSUN': 
        dsun_obs=dsun_obs/695700000
        extent=[0, tspan/60, dsun_obs*math.sin(math.radians(5)), dsun_obs*math.sin(math.radians(45))]
        titleu=unit
    #if unit is au, adjust dsun_obs anad convert elongations
    elif unit=='AU' or unit=='au':
        dsun_obs=dsun_obs/149597870700        
        extent=[0, tspan/60, dsun_obs*math.sin(math.radians(5)), dsun_obs*math.sin(math.radians(45))] 
        titleu=unit
    #if unit is km, adjust dsun_obs anad convert elongations    
    elif unit=='KM' or unit=='km' or unit=='Km':
        dsun_obs=dsun_obs/1000         
        extent=[0, tspan/60, dsun_obs*math.sin(math.radians(5)), dsun_obs*math.sin(math.radians(45))]
        titleu='km'
    #if unit isn't recognized keep it in degrees
    elif unit!= 'Deg' or unit != 'Degrees' or unit != 'DEG' or unit !='DEGREES':
        print('Invalid unit, using elongation (degrees)')
    #display the map   
    ax2=plt.subplot()
    ax2.imshow(jimg, origin='lower', vmin=vrange[0], vmax=vrange[1], cmap='gray', extent=extent, aspect='auto')
    #Set up the axes
    ax2.set_xlabel('Hrs after Start Time: '+hdrs[0]['date-avg'])
    ax2.set_ylabel(unit)
    ax2.set_title('PA='+str(pa)+'\u00B0 Mean S/C Dist ='+str(round(dsun_obs, 3))+' '+titleu)
    #return the image    
    return jimg
        