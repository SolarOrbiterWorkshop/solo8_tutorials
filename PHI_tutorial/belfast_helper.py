import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time
import sunpy.map
from astropy import units as u

def load_field_stop(path = None):
    """load hrt field stop

    Parameters
    ----------
    path: str
        location of the field stop file (optional)

    Returns
    -------
    field_stop: numpy ndarray
        the field stop of the HRT telescope.
    """
    if path is None:
        path = "/scratch/slam/sinjan/solo_attic_fits/demod_mats_field_stops/HRT_field_stop.fits"
    
    hdu_list_tmp = fits.open(path)
    field_stop = np.asarray(hdu_list_tmp[0].data, dtype=np.float32)
    field_stop = np.where(field_stop > 0,1,0)
    
    return field_stop

def plot_fdt_phys_obs(inver_data, suptitle = None): 
    """plot the physical observables from fdt
    
    Parameters
    ----------
    inver_data: numpy array
        FDT physcial observables in format: np.asarray([fdt_icnt, fdt_bmag, fdt_binc, fdt_bazi, fdt_vlos, fdt_blos])
        or '...rte_data_products.fits' file
    suptitle: str
        Name for the plot.
    
    Returns
    -------
    None
    """
    icnt = inver_data[0,:,:]
    bmag = inver_data[1,:,:]
    binc = inver_data[2,:,:]
    bazi = inver_data[3,:,:]
    blos = inver_data[5,:,:]
    
    fig, (ax1, ax2, ax3) = plt.subplots(3,2, figsize = (15,18))

    im4 = ax1[0].imshow(icnt, cmap = "gist_heat", origin="lower") #continuum
    im1 = ax1[1].imshow(bmag, cmap = "plasma", origin="lower") #field strength
    im2 = ax2[0].imshow(binc, cmap = "RdGy", origin="lower") #field inclination
    im3 = ax2[1].imshow(bazi, cmap = "viridis", origin="lower") #field azimuth
    
    #mean correct vlos
    vlos2 = inver_data[4,:,:]
    vlos3 = vlos2 - np.mean(inver_data[4,512:1535,512:1535])
        
    seis = plt.cm.seismic
    norm = plt.Normalize(-2, 2, clip = True)
    rgba_4 = seis(norm(vlos3))
        
    gray = plt.cm.gray
    norm = plt.Normalize(-100, 100, clip = True)
    rgba_5 = gray(norm(blos))
    
    im5 = ax3[0].imshow(rgba_4, cmap = seis, origin="lower") 
    im6 = ax3[1].imshow(rgba_5, cmap = gray, origin="lower")

    fig.colorbar(im4, ax=ax1[0],fraction=0.046, pad=0.04)
    fig.colorbar(im1, ax=ax1[1],fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=ax2[0],fraction=0.046, pad=0.04)
    fig.colorbar(im3, ax=ax2[1],fraction=0.046, pad=0.04)
    fig.colorbar(im5, ax=ax3[0],fraction=0.046, pad=0.04)
    fig.colorbar(im6, ax=ax3[1],fraction=0.046, pad=0.04)

    im1.set_clim(0, 1000)
    im2.set_clim(0,180)
    im3.set_clim(0,180)
    im4.set_clim(0,1.2)
    im5.set_clim(-2,2)
    im6.set_clim(-100,100)

    ax1[1].set_title(r'Magnetic Field Strength [Gauss]')
    ax2[0].set_title(f'Inclination [Degrees]')
    ax2[1].set_title(r'Azimuth [Degrees]')
    ax1[0].set_title("Continuum Intensity")#f'LOS Magnetic Field (Gauss)')
    ax3[0].set_title(r'Vlos [km/s]')
    ax3[1].set_title(r'Blos [Gauss]')
    
    ax1[0].text(35,40, '(a)', color = "white", size = 'x-large')
    ax1[1].text(35,40, '(b)', color = "white", size = 'x-large')
    ax2[0].text(35,40, '(c)', color = "white", size = 'x-large')
    ax2[1].text(35,40, '(d)', color = "white", size = 'x-large')
    ax3[0].text(35,40, '(e)', color = "white", size = 'x-large')
    ax3[1].text(35,40, '(f)', color = "white", size = 'x-large')
    
    if suptitle is not None:
        fig.suptitle(suptitle)
        
    plt.tight_layout()
    plt.show()


def plot_hrt_phys_obs(inver_data, suptitle = None, field_stop = None): 
    """plot the physical observables from hrt
    
    Parameters
    ----------
    inver_data: numpy array
        HRT physcial observables in format: np.asarray([hrt_icnt, hrt_bmag, hrt_binc, hrt_bazi, hrt_vlos, hrt_blos])
        or '...rte_data_products.fits' file
    suptitle: str
        Name for the plot
    field_stop: numpy ndarray
        HRT field stop
    
    Returns
    -------
    None
    """
    if field_stop is None:
        field_stop = load_field_stop()[:,::-1]
        
    inver_data *= field_stop[np.newaxis,:,:]
        
    fs_idx = np.where(field_stop < 1)
    
    #create custom colormaps that are black in the field stop region
    gist = plt.cm.gist_heat
    norm = plt.Normalize(0, 1.2, clip = True)
    rgba_0 = gist(norm(inver_data[0,:,:]))
    rgba_0[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    plasma = plt.cm.plasma
    norm = plt.Normalize(0, 1000, clip = True)
    rgba_1 = plasma(norm(inver_data[1,:,:]))
    rgba_1[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    rdgy = plt.cm.RdGy
    norm = plt.Normalize(0, 180)
    rgba_2 = rdgy(norm(inver_data[2,:,:]))
    rgba_2[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    viridis = plt.cm.viridis
    norm = plt.Normalize(0, 180)
    rgba_3 = viridis(norm(inver_data[3,:,:]))
    rgba_3[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    fig, (ax1, ax2, ax3) = plt.subplots(3,2, figsize = (15,18))

    im4 = ax1[0].imshow(rgba_0, cmap = gist, origin="lower") #continuum
    im1 = ax1[1].imshow(rgba_1, cmap = plasma, origin="lower") #field strength
    im2 = ax2[0].imshow(rgba_2, cmap = rdgy, origin="lower") #field inclination
    im3 = ax2[1].imshow(rgba_3, cmap = viridis, origin="lower") #field azimuth
    
    #mean correct vlos
    vlos2 = inver_data[4,:,:]
    vlos3 = vlos2 - np.mean(inver_data[4,512:1535,512:1535])
    blos = inver_data[1,:,:] * np.cos(inver_data[2,:,:]/180*np.pi) #vlos
    
    if field_stop is not None:
        blos *= field_stop
        vlos3 *= field_stop
        
    seis = plt.cm.seismic
    norm = plt.Normalize(-2, 2, clip = True)
    rgba_4 = seis(norm(vlos3))
    rgba_4[fs_idx[0],fs_idx[1], :3] = 0,0,0
        
    gray = plt.cm.gray
    norm = plt.Normalize(-100, 100, clip = True)
    rgba_5 = gray(norm(blos))
    rgba_5[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    im5 = ax3[0].imshow(rgba_4, cmap = seis, origin="lower") 
    im6 = ax3[1].imshow(rgba_5, cmap = gray, origin="lower")

    fig.colorbar(im4, ax=ax1[0],fraction=0.046, pad=0.04)
    fig.colorbar(im1, ax=ax1[1],fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=ax2[0],fraction=0.046, pad=0.04)
    fig.colorbar(im3, ax=ax2[1],fraction=0.046, pad=0.04)
    fig.colorbar(im5, ax=ax3[0],fraction=0.046, pad=0.04)
    fig.colorbar(im6, ax=ax3[1],fraction=0.046, pad=0.04)

    im1.set_clim(0, 1000)
    im2.set_clim(0,180)
    im3.set_clim(0,180)
    im4.set_clim(0,1.2)
    im5.set_clim(-2,2)
    im6.set_clim(-100,100)

    ax1[1].set_title(r'Magnetic Field Strength [Gauss]')
    ax2[0].set_title(f'Inclination [Degrees]')
    ax2[1].set_title(r'Azimuth [Degrees]')
    ax1[0].set_title("Continuum Intensity")#f'LOS Magnetic Field (Gauss)')
    ax3[0].set_title(r'Vlos [km/s]')
    ax3[1].set_title(r'Blos [Gauss]')
    
    ax1[0].text(35,40, '(a)', color = "white", size = 'x-large')
    ax1[1].text(35,40, '(b)', color = "white", size = 'x-large')
    ax2[0].text(35,40, '(c)', color = "white", size = 'x-large')
    ax2[1].text(35,40, '(d)', color = "white", size = 'x-large')
    ax3[0].text(35,40, '(e)', color = "white", size = 'x-large')
    ax3[1].text(35,40, '(f)', color = "white", size = 'x-large')
    
    if suptitle is not None:
        fig.suptitle(suptitle)
        
    plt.tight_layout()
    plt.show()
    
    
def get_wv_arr_and_ic_wv(file,num_wl = 6,verbose = False):
    '''calculate the wavelength sampling and continuum position from the voltage data in the fits file
    
    Parameters
    ----------
    file: str
        File path to the PHI stokes dataset
    num_wl: int
        Number of wavelengths in the dataset
    verbose: bool
        Print the continuum position to console if True
    
    Returns
    -------
    wave_axis: numpy ndarray
        Wavelength positions in Angstrom
    voltagesData: numpy ndarray
        Voltages of the 6 wavelength positions
    cpos: int
        Index for the Continuum Wavelength position
    '''
    fg_head = 3

    with fits.open(file) as hdu_list:
        header = hdu_list[fg_head].data
        tunning_constant = float(header[0][4])/1e9
        ref_wavelength = float(header[0][5])/1e3
        
        voltagesData = np.zeros(num_wl)
        hi = np.histogram(header['PHI_FG_voltage'],bins=7)
        yi = hi[0]; xi = hi[1]
        j = 0
        for i in range(num_wl + 1):
            if yi[i] != 0 :
                if i < num_wl:
                    idx = np.logical_and(header['PHI_FG_voltage']>=xi[i],header['PHI_FG_voltage']<xi[i+1])
                else:
                    idx = np.logical_and(header['PHI_FG_voltage']>=xi[i],header['PHI_FG_voltage']<=xi[i+1])
                voltagesData[j] = int(np.median(header['PHI_FG_voltage'][idx]))
                j += 1
    
    d1 = voltagesData[0] - voltagesData[1]
    d2 = voltagesData[num_wl-2] - voltagesData[num_wl-1]
    if np.abs(d1) > np.abs(d2):
        cpos = 0
    else:
        cpos = num_wl-1
    if verbose:
        print('Continuum position at wave: ', cpos)
    wave_axis = voltagesData*tunning_constant + ref_wavelength  #6173.3356
    return wave_axis,voltagesData,cpos
    
    
def plot_fdt_stokes(stokes_arr, wv, subsec = None, title = None):
    """plot fdt stokes maps at one wavelength

    Parameters
    ----------
    stokes_arr : numpy ndarray
        Full FDT Stokes Array.
    wv : int
        Index for the desired wavelength position.
    subsec: numpy ndarray
        Region of interest to be plotted [start_x,end_x,start_y,end_y]
    title: str
        Title of figure
        
    Returns
    -------
    None
    """
    fig, (ax1, ax2) = plt.subplots(2,2, figsize = (15,12))

    
    gist = plt.cm.gist_heat
    norm = plt.Normalize(-0.01, 0.01, clip = True)
    rgba_0 = gist(norm(stokes_arr[wv,1,:,:]))
    
    rgba_1 = gist(norm(stokes_arr[wv,2,:,:]))
    
    rgba_2 = gist(norm(stokes_arr[wv,3,:,:]))
    
    if subsec is not None:
        start_row, end_row = subsec[2:4]
        start_col, end_col = subsec[:2]
        assert len(subsec) == 4
        assert start_row >= 0 and start_row < 2048
        assert end_row >= 0 and end_row < 2048
        assert start_col >= 0 and start_col < 2048
        assert end_col >= 0 and end_col < 2048
        
    else:
        start_row, start_col = 0,0
        end_row, end_col = stokes_arr.shape[2]-1,stokes_arr.shape[3]-1
        
    
    im1 = ax1[0].imshow(stokes_arr[wv,0,start_row:end_row,start_col:end_col], cmap = "gist_heat", origin="lower") 
    im2 = ax1[1].imshow(rgba_0[start_row:end_row,start_col:end_col], cmap = gist, origin="lower")
    im3 = ax2[0].imshow(rgba_1[start_row:end_row,start_col:end_col], cmap = gist, origin="lower") 
    im4 = ax2[1].imshow(rgba_2[start_row:end_row,start_col:end_col], cmap = gist, origin="lower")

    fig.colorbar(im1, ax=ax1[0],fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=ax1[1],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    fig.colorbar(im3, ax=ax2[0],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    fig.colorbar(im4, ax=ax2[1],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    
    clim = 0.01

    im1.set_clim(0, 1.2)
    im2.set_clim(-clim, clim)
    im3.set_clim(-clim, clim)
    im4.set_clim(-clim, clim)
    
    ax1[0].set_title(r'I/<I_c>')
    ax1[1].set_title(f'Q/<I_c>')
    ax2[0].set_title(r'U/<I_c>')
    ax2[1].set_title(f'V/<I_c>')

    ax1[0].text(35,40, '(a)', color = "white", size = 'x-large')
    ax1[1].text(35,40, '(b)', color = "white", size = 'x-large')
    ax2[0].text(35,40, '(c)', color = "white", size = 'x-large')
    ax2[1].text(35,40, '(d)', color = "white", size = 'x-large')
    
    if isinstance(title,str):
         plt.suptitle(title)
    else:
        plt.suptitle(f"SO/PHI-FDT Stokes at Wavelength Index: {wv}")
    plt.tight_layout()
    plt.show()
    
    
def plot_hrt_stokes(stokes_arr, wv, subsec = None, title = None, field_stop = None):
    """plot hrt stokes maps at one wavelength

    Parameters
    ----------
    stokes_arr : numpy ndarray
        Full HRT Stokes Array.
    wv : int
        Index for the desired wavelength position.
    subsec: numpy ndarray
        Region of interest to be plotted [start_x,end_x,start_y,end_y]
    title: str
        Title of figure
        
    Returns
    -------
    None
    """
    fig, (ax1, ax2) = plt.subplots(2,2, figsize = (15,12))

    if field_stop is None:
        fs = load_field_stop()[:,::-1]
    else:
        fs = field_stop
    fs_idx = np.where(fs < 1)
    
    gist = plt.cm.gist_heat
    norm = plt.Normalize(-0.01, 0.01, clip = True)
    rgba_0 = gist(norm(stokes_arr[:,:,1,wv]))
    rgba_0[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    rgba_1 = gist(norm(stokes_arr[:,:,2,wv]))
    rgba_1[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    rgba_2 = gist(norm(stokes_arr[:,:,3,wv]))
    rgba_2[fs_idx[0],fs_idx[1], :3] = 0,0,0
    
    if subsec is not None:
        start_row, end_row = subsec[2:4]
        start_col, end_col = subsec[:2]
        assert len(subsec) == 4
        assert start_row >= 0 and start_row < 2048
        assert end_row >= 0 and end_row < 2048
        assert start_col >= 0 and start_col < 2048
        assert end_col >= 0 and end_col < 2048
        
    else:
        start_row, start_col = 0,0
        end_row, end_col = stokes_arr.shape[0]-1,stokes_arr.shape[0]-1
        
    
    im1 = ax1[0].imshow(stokes_arr[start_row:end_row,start_col:end_col,0,wv], cmap = "gist_heat", origin="lower") 
    im2 = ax1[1].imshow(rgba_0[start_row:end_row,start_col:end_col], cmap = gist, origin="lower")
    im3 = ax2[0].imshow(rgba_1[start_row:end_row,start_col:end_col], cmap = gist, origin="lower") 
    im4 = ax2[1].imshow(rgba_2[start_row:end_row,start_col:end_col], cmap = gist, origin="lower")

    fig.colorbar(im1, ax=ax1[0],fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=ax1[1],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    fig.colorbar(im3, ax=ax2[0],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    fig.colorbar(im4, ax=ax2[1],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    
    clim = 0.01

    im1.set_clim(0, 1.2)
    im2.set_clim(-clim, clim)
    im3.set_clim(-clim, clim)
    im4.set_clim(-clim, clim)
    
    ax1[0].set_title(r'I/<I_c>')
    ax1[1].set_title(f'Q/<I_c>')
    ax2[0].set_title(r'U/<I_c>')
    ax2[1].set_title(f'V/<I_c>')

    ax1[0].text(35,40, '(a)', color = "white", size = 'x-large')
    ax1[1].text(35,40, '(b)', color = "white", size = 'x-large')
    ax2[0].text(35,40, '(c)', color = "white", size = 'x-large')
    ax2[1].text(35,40, '(d)', color = "white", size = 'x-large')
    
    if isinstance(title,str):
         plt.suptitle(title)
    else:
        plt.suptitle(f"SO/PHI-HRT Stokes at Wavelength Index: {wv}")
    plt.tight_layout()
    plt.show()
    
    
def plot_hrt_noise_both(stokes_V, blos, field_stop = None):
    """plot histograms with Gaussian fit of Stokes V and Blos of the HRT files

    Parameters
    ----------
    stokes_V : numpy.ndarray
        Stokes V at Ic.
    blos : numpy.ndarray
        Corresponding Blos map.
        
    Returns
    -------
    None
    """
    if field_stop is None:
        field_stop = load_field_stop()
        field_stop = field_stop[:,::-1]
 
    idx = np.where(field_stop == 1)

    def gaussian(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2*sigma**2))

    
    fig, ax = plt.subplots(1, 2, figsize = (14,6))
    
    #stokes V
    stokes_v = stokes_V[idx[0], idx[1]].flatten()
    
    hist, bin_edges = np.histogram(stokes_v, bins = 2000)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

    ax[0].hist(stokes_v, bins = 2000)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    
    p0 = [4e5,0,1e-3]
    coeff, var_matrix = curve_fit(gaussian, bin_centres, hist, p0 = p0)
    y = gaussian(bin_centres, *coeff)
    print(f"Stokes V noise is: {coeff[2]:.4e}")                  
    ax[0].plot(bin_centres, y, linestyle = "--", color= "red", lw = 2, label = f"Mean: {coeff[1]:.2g} Std: {coeff[2]:.1e}")    
    ax[0].set_xlabel("V/<I_c>")
    ax[0].legend(loc = "upper right", prop = {'size':14})
    ax[0].set_ylabel("Count")
    ax[0].set_xlim(-0.01,0.01)
    ylim = max(y)+0.5e4
    ax[0].set_ylim(0,ylim)
    ax[0].set_xticks((-0.01,-0.005,0,0.005,0.01))
    ax[0].set_xticklabels(("-0.01","-0.005","0","0.005","0.01"))
    
    #blos

    blos = blos[idx]
    blos_new = blos[np.where(abs(blos) <= 200)]
    blos_new = blos_new[np.where(abs(blos_new) >= 0.015)]

    hist, bin_edges = np.histogram(blos_new.flatten(), bins = 2000)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    bin_centres = bin_centres[np.where(hist <= 50000)]
    hist = hist[np.where(hist <= 50000)]
    
    p0 = [4e5,0,6]
    coeff, var_matrix = curve_fit(gaussian, bin_centres, hist, p0 = p0)
    y = gaussian(bin_centres, *coeff)
    ax[1].plot(bin_centres, y, linestyle = "--", color= "red", lw = 2, label = f"Mean: {coeff[1]:.2g} Std: {coeff[2]:.2g}")
    ax[1].hist(blos_new.flatten(), bins = 2000, range = (-200,200))[:2]
    ax[1].set_xlabel("Blos [Gauss]")
    ax[1].set_xlim(-50,50)
    ax[1].legend(loc = "upper right", prop = {'size':14})
    ax[1].set_ylabel("Count")
    ylim = max(y)+1e4
    ax[1].set_ylim(0,ylim)
    ax[1].set_xticks((-50,-25,0,25,50))
    ax[1].set_xticklabels(("-50", "-25", "0", "25", "50"))
   
    plt.tight_layout()
    plt.show()