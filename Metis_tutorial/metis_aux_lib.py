import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
from scipy import ndimage
import astropy.units as u
from astropy.coordinates.representation import (
    CartesianRepresentation, 
    CylindricalRepresentation,
)



def get_linspace_arr(val1, val2, dval):
    """
    Return evenly spaced numbers over a specified interval and step.
    This function is aimed to overcome rounding problems when creating
    arrays with ``numpy.linspace`` or ``numpy.arange`` funcitons.
    
    Parameters
    ----------
    val1 : `float`
        The starting value of the sequence.
    val2 : `float`
        The end value of the sequence.
    dval : `float`
        Spacing between values.
    
    Returns
    -------
    polar_img_arr : `numpy.ndarray`
        Array of evenly spaced values.
    
    """
    num = (val2 - val1) / dval
    num = np.around(num, decimals=0)
    arr = np.linspace(val1, val2, int(num)+1)
    return arr


def cart_to_polar_v1(img_arr, r_arr, phi_arr, xc, yc, rsun_pix, rot_angle=None, cval=0.0):
    """
    Convert image array from Cartesian (`x`,`y`) to polar (`Rsun`, `phi`) 
    coordinates.
    
    Parameters
    ----------
    img_arr : `numpy.ndarray`
        The input array with initial image in Cartesian coordinates.
    r_arr : `numpy.ndarray`
        The array of values of polar coordinate `Rsun`, at which the output 
        polar image array will be calculated.
    phi_arr : `numpy.ndarray`
        The array of values of polar coordinate `phi`, at which the output 
        polar image array will be calculated.
    xc : `float`
        The first Cartesian coordinate (along `x`-axis) of the pole of 
        the polar coordinate system.
    yc : `float`
        The second Cartesian coordinate (along `y`-axis) of the pole of 
        the polar coordinate system.
    rsun_pix : `float`
        The radius of the Sun in pixels of ``img_arr``.
    rot_angle : `float`, optional
        If not None, the rotation angle in degrees, applied counter clockwise. 
        Default is None.
    cval : `float`, optional
        The value to fill past edges of input. Default is 0.0.

    Returns
    -------
    polar_img_arr : `numpy.ndarray`
        The result of transforming the input ``img_arr`` to polar coordinates. 
    
    """
    r_pix_arr = r_arr * rsun_pix
    if rot_angle is not None:
        phi_arr = phi_arr - rot_angle
    polar_inds_x = []
    polar_inds_y = []
    for r in r_pix_arr:
        polar_inds_x.append(
            xc + (r) / (1.0+np.tan(np.deg2rad(phi_arr/2.0))**2.0) * (1.0-np.tan(np.deg2rad(phi_arr/2.0))**2.0)
        )
        polar_inds_y.append(
            yc + (r) / (1.0+np.tan(np.deg2rad(phi_arr/2.0))**2.0) * 2.0*np.tan(np.deg2rad(phi_arr/2.0))
        )
    polar_inds = np.array([polar_inds_y, polar_inds_x])  #  NB: index sequence not like in idl: img_arr[x,y]
    polar_img_arr = ndimage.map_coordinates(img_arr, polar_inds, order=1, cval=cval)
    return polar_img_arr


def cart_to_polar(img_arr, r_arr, phi_arr, xc, yc, rsun_pix, rot_angle=None, cval=0.0):
    """
    Convert image array from Cartesian (`x`,`y`) to polar (`Rsun`, `phi`) 
    coordinates by means of CylindricalRepresentation object of 
    astropy.coordinates.representation library.
    
    Parameters
    ----------
    img_arr : `numpy.ndarray`
        The input array with initial image in Cartesian coordinates.
    r_arr : `numpy.ndarray`
        The array of values of polar coordinate `Rsun`, at which the output 
        polar image array will be calculated.
    phi_arr : `numpy.ndarray`
        The array of values of polar coordinate `phi`, at which the output 
        polar image array will be calculated.
    xc : `float`
        The first Cartesian coordinate (along `x`-axis) of the pole of 
        the polar coordinate system.
    yc : `float`
        The second Cartesian coordinate (along `y`-axis) of the pole of 
        the polar coordinate system.
    rsun_pix : `float`
        The radius of the Sun in pixels of ``img_arr``.
    rot_angle : `float`, optional
        If not None, the rotation angle in degrees, applied counter clockwise. 
        Default is None.
    cval : `float`, optional
        The value to fill past edges of input. Default is 0.0.

    Returns
    -------
    polar_img_arr : `numpy.ndarray`
        The result of transforming the input ``img_arr`` to polar coordinates. 
    
    """
    r_pix_arr = r_arr * rsun_pix
    if rot_angle is not None:
        phi_arr = phi_arr - rot_angle
    phi_matrix, r_matrix = np.meshgrid(phi_arr, r_pix_arr)
    polar_repr = CylindricalRepresentation(r_matrix, phi_matrix*u.deg, np.zeros(r_matrix.shape))
    cart_repr = polar_repr.to_cartesian()
    polar_inds = np.array([cart_repr.y + yc, cart_repr.x + xc])  #  NB: index sequence not like in IDL: img_arr[x,y]
    polar_img_arr = ndimage.map_coordinates(img_arr, polar_inds, order=1, cval=cval)
    return polar_img_arr


def cut_img_arr_fov(img_arr, xc, yc, rsun_pix, r1, r2, cval=np.nan):
    """
    Mask (set to ``cval``) the values of the input array outside the 
    field of view defined by ``r1`` and ``r2``.
    
    Parameters
    ----------
    img_arr : `numpy.ndarray`
        The input array.
    xc : `float`
        The first Cartesian coordinate (along `x`-axis) of the Sun center.
    yc : `float`
        The second Cartesian coordinate (along `x`-axis) of the Sun center.
    rsun_pix : `float`
        The radius of the Sun in pixels of ``img_arr``.
    r1, r2 : `float`
        The lower and upper limits (in Rsun) of the field of view covered by 
        the output image array.
    cval : `float`, optional
        The values of maked pixels (outside the field of view). Default is 
        ``np.nan``.

    Returns
    -------
    img_arr : `numpy.ndarray`
        The result of masking the input array.
    
    """
    fov1 = r1 * rsun_pix
    fov2 = r2 * rsun_pix
    x = np.arange(0, img_arr.shape[1], 1)
    y = np.arange(0, img_arr.shape[0], 1)
    xx, yy = np.meshgrid(x, y, sparse=True)
    dist_cen = np.sqrt((xx-xc)**2 + (yy-yc)**2)
    img_arr[dist_cen < fov1] = cval
    img_arr[dist_cen > fov2] = cval
    return img_arr


def polar_to_cart_v1(img_arr, cart_img_shape, xc, yc, rsun_pix, r1, r2, dr, 
                  rot_angle=None, cut_fov=True):
    """
    Convert image array from polar (`Rsun`, `phi`) to Cartesian (`x`,`y`)
    coordinates.
    
    Parameters
    ----------
    img_arr : `numpy.ndarray`
        The input array with initial image in polar coordinates.
    cart_img_shape : `array_like`
        The array, list or tuple which define a shape of the output 
        Cartesian image array.
    xc : `float`
        The first Cartesian coordinate (along `x`-axis) of the Sun center.
    yc : `float`
        The second Cartesian coordinate (along `x`-axis) of the Sun center.
    rsun_pix : `float`
        The radius of the Sun in pixels of ``img_arr``.
    r1, r2 : `float`
        The lower and upper limits (in Rsun) of the field of view covered by 
        the output image array. It corresponds to the range of heliocentric 
        distances covered by polar ``img_arr``.
    dr : `float`
        Spacing between values `Rsun` values at which polar ``img_arr`` was 
        calculated.
    rot_angle : `float`, optional
        If not None, the rotation angle in degrees, applied counter clockwise. 
        Default is None.
    cut_fov : `bool`, optional
        If True (default), the values of the output image array outside the field 
        view defined by ``r1`` and ``r2`` are masked (set to ``np.nan``).

    Returns
    -------
    cart_img_arr : `numpy.ndarray`
        The result of transforming the input ``img_arr`` to Cartesian coordinates. 
    
    """
    n_xpix, n_ypix = cart_img_shape
    x_arr = np.arange(int(n_xpix))
    y_arr = np.arange(int(n_ypix))
    polar_inds_r = []
    polar_inds_phi = []
    for y in y_arr:
        r = np.sqrt((x_arr-xc)**2 + (y-yc)**2)
        if y == yc:
            phi_vals = np.where(x_arr >= xc, 0.0, 180.0)
        else:
            phi_vals = np.rad2deg( 2.0 * np.arctan((y-yc)/(r+x_arr-xc)) )
        if rot_angle is not None:
            phi_vals -= rot_angle
        phi_vals[phi_vals < 0.0] += 360
        polar_inds_r.append(
            (r/rsun_pix - r1)/dr
        )
        polar_inds_phi.append(
            phi_vals
        )
    cart_inds = np.array([polar_inds_r, polar_inds_phi])
    cart_img_arr = ndimage.map_coordinates(img_arr, cart_inds, mode='grid-wrap', 
                                           order=1)
    if cut_fov:
        cart_img_arr = cut_img_arr_fov(cart_img_arr, xc, yc, rsun_pix, r1, r2)
    return cart_img_arr


def polar_to_cart(img_arr, cart_img_shape, xc, yc, rsun_pix, r1, r2, dr, 
                  rot_angle=None, cut_fov=True):
    """
    Convert image array from polar (`Rsun`, `phi`) to Cartesian (`x`,`y`)
    coordinates by means of CartesianRepresentation object and 
    CylindricalRepresentation object of astropy.coordinates.representation 
    library.
    
    Parameters
    ----------
    img_arr : `numpy.ndarray`
        The input array with initial image in polar coordinates.
    cart_img_shape : `array_like`
        The array, list or tuple which define a shape of the output 
        Cartesian image array.
    xc : `float`
        The first Cartesian coordinate (along `x`-axis) of the Sun center.
    yc : `float`
        The second Cartesian coordinate (along `x`-axis) of the Sun center.
    rsun_pix : `float`
        The radius of the Sun in pixels of ``img_arr``.
    r1, r2 : `float`
        The lower and upper limits (in Rsun) of the field of view covered by 
        the output image array. It corresponds to the range of heliocentric 
        distances covered by polar ``img_arr``.
    dr : `float`
        Spacing between values `Rsun` values at which polar ``img_arr`` was 
        calculated.
    rot_angle : `float`, optional
        If not None, the rotation angle in degrees, applied counter clockwise. 
        Default is None.
    cut_fov : `bool`, optional
        If True (default), the values of the output image array outside the field 
        view defined by ``r1`` and ``r2`` are masked (set to ``np.nan``).

    Returns
    -------
    cart_img_arr : `numpy.ndarray`
        The result of transforming the input ``img_arr`` to Cartesian coordinates. 
    
    """
    n_xpix, n_ypix = cart_img_shape
    x_arr = np.arange(int(n_xpix))
    y_arr = np.arange(int(n_ypix))
    x_matrix, y_matrix = np.meshgrid(x_arr, y_arr)
    cart_repr = CartesianRepresentation(
        x_matrix-xc, y_matrix-yc, np.zeros(x_matrix.shape)
    )
    cylin_repr = CylindricalRepresentation(1, 1*u.deg, 0)  # defining an object of the class CylindricalRepresentation
    polar_repr = cylin_repr.from_cartesian(cart_repr)  # converting CartesianRepresentation 'cart_repr' to the class of object 'cylin_repr': i.e. CylindricalRepresentation
    polar_repr_r = (polar_repr.rho/rsun_pix - r1)/dr
    polar_repr_phi = polar_repr.phi.to_value('deg')
    if rot_angle is not None:
        polar_repr_phi -= rot_angle
    polar_repr_phi[polar_repr_phi < 0.0] += 360
    cart_inds = np.array([polar_repr_r, polar_repr_phi])
    cart_img_arr = ndimage.map_coordinates(img_arr, cart_inds, mode='grid-wrap', 
                                           order=1)
    if cut_fov:
        cart_img_arr = cut_img_arr_fov(cart_img_arr, xc, yc, rsun_pix, r1, r2)
    return cart_img_arr


def poly_fun(x, *aks):
    """
    Retrurns a value of multi-dimentional polinomial function. See e.g. 
    equation 6 from Hayes et al. 2001 ApJ, 548, 1081.
    
    Parameters
    ----------
    x : `numpy.ndarray`
        The array of values of variables `x_i`.
    aks : `tuple`
        The tuple with values of `ak` coefficients.

    Returns
    -------
    f : `float`
        The value of polinomial function.

    """
    f = 0
    for i in range(len(aks)):
        f += x[:,i] * aks[i]
    return f


def get_Ar_Br(r, q=0.75):
    """
    Retrun geometric factors `A(r)` and `B(r)`. See equations 2 and 3 from 
    Hayes et al. 2001 ApJ, 548, 1081.
    
    Parameters
    ----------
    r : `float`
        The value of radial height from Sun center.
    q : `float`, optional
        The coefficient of the limb darkening. Default is 0.75 (as in 
        van de Hulst 1950, Bull. Astron. Inst. Netherlands, 11, 135).

    Returns
    -------
    A : `float`
        The value of `A(r)` at ``r``.
    B : `float`
        The value of `B(r)` at ``r``.
    
    """
    sin_g = 1.0/r
    cos_g = np.sqrt(1-sin_g**2)
    var_2A_plus_B = (1-q)/(1-q/3.0)*(2*(1-cos_g)) + q/(1-q/3.0)*(1 - cos_g**2/sin_g * np.log((1+sin_g)/cos_g))
    var_2A_minus_B = (1-q)/(1-q/3.0)*(2./3.*(1-cos_g**3)) + q/(1-q/3.0)*(1./4. + sin_g**2/4.0 - cos_g**4/(4*sin_g) * np.log((1+sin_g)/cos_g))
    A = (var_2A_plus_B + var_2A_minus_B)/4.
    B = (var_2A_plus_B - var_2A_minus_B)/2.0
    return A, B


def fun_to_int_Gkx(r, x, k, q):
    """
    Retrun a value of a function which is integrated for obtaining 
    `G_k(x)`. See equation 5 from Hayes et al. 2001 ApJ, 548, 1081.
    
    Parameters
    ----------
    r : `float`
        The value of radial height from Sun center. ``fun`` is integrated with 
        respect to ``r`` variable.
    x : `float`
        The value of impact distance.
    k : `float`
        The value of exponent.
    q : `float`
        The coefficient of the limb darkening.

    Returns
    -------
    fun : `float`
        The value of function ``fun`` (`G_k(x)`) at ``r``.
    
    """
    A, B = get_Ar_Br(r, q=q)
    fun = A * 2*r**(-k+1)/np.sqrt(r**2-x**2)  +  \
          (B-A) * r**(-k-1) * x**2/np.sqrt(r**2-x**2)
    return fun


def fun_to_int_pGkx(r, x, k, q):
    """
    Retrun a value of a function which is integrated for obtaining 
    `pG_k(x)`. See equation 5 from Hayes et al. 2001 ApJ, 548, 1081.
    
    Parameters
    ----------
    r : `float`
        The value of radial height from Sun center. ``fun`` is integrated with 
        respect to ``r`` variable.
    x : `float`
        The value of impact distance.
    k : `float`
        The value of exponent.
    q : `float`
        The coefficient of the limb darkening.

    Returns
    -------
    fun : `float`
        The value of function ``fun`` (`pG_k(x)`) at ``r``.
    
    """
    A, B = get_Ar_Br(r, q=q)
    fun = (A-B) * r**(-k-1) * x**2/np.sqrt(r**2-x**2)
    return fun


def calc_Gkx_pGkx(x_arr, k_arr, q):
    """
    Calculate integrals `G_k(x)` and `pG_k(x)` from Hayes et al. 2001 ApJ, 548, 1081
    (eqs. 5 and 3).
    
    Parameters
    ----------
    x_arr : `numpy.ndarray`
        The array of values of radial height from Sun center.
    k_arr : `array_like`
        The array or list of values of exponent k used to express electron 
        density along a radial trace in the polynomial form: `sum_k(alpha_k*r**(-k))`.
    q : `float`
        The coefficient of the limb darkening.

    Returns
    -------
    Gkx : `numpy.ndarray`
        The resulting array with `G_k(x)` values, which has a shape of 
        ``len(x_arr)`` x ``len(k_arr)``.
    Gkx_err : `numpy.ndarray`
        The resulting array with `G_k(x)` integration errors, which has 
        a shape of ``len(x_arr)`` x ``len(k_arr)``.
    pGkx : `numpy.ndarray`
        The resulting array with `pG_k(x)` values, which has a shape of 
        ``len(x_arr)`` x ``len(k_arr)``.
    pGkx_err : `numpy.ndarray`
        The resulting array with `pG_k(x)` integration errors, which has 
        a shape of ``len(x_arr)`` x ``len(k_arr)``.
        
    """
    print()
    print('calc_binvertdata:')
    print(' - x = ({x1}-{x2}, dx={dx:0.3f})'.format(x1=x_arr[0], x2=x_arr[-1], dx=x_arr[1]-x_arr[0]))
    print(' - k = {k_arr}'.format(k_arr=k_arr))
    print(' - q = {q} '.format(q=q))
    print('.', end='')
    C = 3/4 * 696340e5 * 6.6524587158e-29*1e4  # adopted from eq. 21 of van de Hulst 1950, Bull. Astron. Inst. Netherlands, 11, 135

    Gkx = np.zeros((len(x_arr), len(k_arr)))
    Gkx_err = np.zeros((len(x_arr), len(k_arr)))
    pGkx = np.zeros((len(x_arr), len(k_arr)))
    pGkx_err = np.zeros((len(x_arr), len(k_arr)))
    for i in range(len(x_arr)):
        print('.', end='', flush=True)
        for j in range(len(k_arr)):
            x = x_arr[i]
            k = k_arr[j]

            Gkx_int_out =  integrate.quad(fun_to_int_Gkx, x, np.inf, args=(x, k, q)) # Gk(x) from eq. 5 of Hayes et al. 2001 ApJ, 548, 1081
            Gkx[i,j] = C*Gkx_int_out[0]
            Gkx_err[i,j] = C*Gkx_int_out[1]

            pGkx_int_out =  integrate.quad(fun_to_int_pGkx, x, np.inf, args=(x, k, q)) # pGk(x) from eq. 3 of Hayes et al. 2001 ApJ, 548, 1081
            pGkx[i,j] = C*pGkx_int_out[0]
            pGkx_err[i,j] = C*pGkx_int_out[1]
    print('Done')
    return Gkx, Gkx_err, pGkx, pGkx_err


def calc_el_dens_and_K(pB_polar, r_arr, k_arr=[1, 2, 3, 4], q=0.63,
                       negative_val=None):
    """
    Calculate the electron density profiles and K-corona maps.
    
    Parameters
    ----------
    pB_polar : `numpy.ndarray`
        The input array with the polarized brightness image in polar coordinates.
    r_arr : `numpy.ndarray`
        The array of values of polar coordinate `Rsun`, at which the output 
        polar image array will be calculated.
    k_arr : `array_like`, optional
        The array or list of values of exponent k used to express electron 
        density along a radial trace in the polynomial form: `sum_k(alpha_k*r**(-k))`.
        See Hayes et al. 2001 ApJ, 548, 1081. Default is [1, 2, 3, 4]
    q : `float`, optional
        The coefficient of the limb darkening. Default is 0.63 (as e.g. in 
        Antonucci et al. 2020, A&A, 642, A10).
    negative_val : `float`, optional
        If not None, the value to fill negative pixels of the K-corona map.
        Default is None.

    Returns
    -------
    dens_polar : `numpy.ndarray`
        The resulting image array of the electron density profiles map.
    K_polar : `numpy.ndarray`
        The resulting image array of the K-corona map.
        
    """
    Gkx, Gkx_err, pGkx, pGkx_err = calc_Gkx_pGkx(r_arr, k_arr=k_arr, q=q)

    ### Calculating electron density ###
    ak_arr = []
    dens_polar = []
    for phi_ind in range(pB_polar.shape[1]):
        ak_p0 = [1., 2., 3., 4.]
        pB_polar_ith = pB_polar[:,phi_ind]
        pB_norm = 10**np.around(np.log10(np.nanmean(pB_polar_ith)))
        x_fit = pGkx / pB_norm
        y_fit = pB_polar_ith / pB_norm
        mask = y_fit > 0
        x_fit = x_fit[mask]
        y_fit = y_fit[mask]
        ak_popt, ak_pcov = curve_fit(poly_fun, x_fit, y_fit, p0=ak_p0)
        ak_perr = np.sqrt(np.diag(ak_pcov))
        ak_arr.append(ak_popt)
        n_r = [ak_popt[i]*r_arr**(-k_arr[i]) for i in range(len(k_arr))]
        n_r = np.sum(np.array(n_r), axis=0)
        dens_polar.append(n_r)
    ak_arr = np.array(ak_arr)
    dens_polar = np.array(dens_polar).transpose()
    
    ### Calculating K-corona ###
    K_polar = np.matmul(Gkx, ak_arr.transpose())
    if negative_val is not None:
        K_polar[K_polar<0.0] = negative_val

    return dens_polar, K_polar


def set_rsun_axes(ax, img_arr, rsun_pix, labelfontsize=18, labelsize=15):
    """
    Convert both axes of input ``ax`` from pixels to Rsun.
    
    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The input axes.
    img_arr : `numpy.ndarray`
        The Cartesian image array plotted in ``ax``.
    rsun_pix : `float`
        The radius of the Sun in pixels of ``img_arr``.
    labelfontsize : `float`, optional
        The fontsize of the axes labels. Default is 18.
    labelsize : `float`, optional
        The fontsize of the tick labels. Default is 15.
    
    """
    cart_img_size1, cart_img_size2 = img_arr.shape
    img_width_rsun_x = cart_img_size1/rsun_pix
    img_width_rsun_y = cart_img_size2/rsun_pix
    ax.set_xticks([])
    secax_x = ax.secondary_xaxis(
        'bottom', 
        functions=(
            lambda x: img_width_rsun_x * (x/(cart_img_size1-1) - 0.5),
            lambda x: (cart_img_size1-1) * (x/img_width_rsun_x + 0.5)
        ),
    )
    secax_x.set_xticks(
        np.linspace(-0.5*img_width_rsun_x, 0.5*img_width_rsun_x, 7)
    )
    secax_x.tick_params(axis='x', which='major', labelsize=labelsize)

    ax.set_yticks([])
    secax_y = ax.secondary_yaxis(
        'left', 
        functions=(
            lambda y: img_width_rsun_y * (y/(cart_img_size2-1) - 0.5),
            lambda y: (cart_img_size2-1) * (y/img_width_rsun_y + 0.5)
        ),
    )
    secax_y.set_yticks(
        np.linspace(-0.5*img_width_rsun_y, 0.5*img_width_rsun_y, 7)
    )
    secax_y.tick_params(axis='y', which='major', labelsize=labelsize)
    
    secax_x.set_xlabel('$R_\odot$', fontsize=labelfontsize)
    secax_y.set_ylabel('$R_\odot$', fontsize=labelfontsize)

    
def set_rsun_phi_axes(ax, img_arr, r_arr, phi_arr, labelfontsize=18, 
                      labelsize=15):
    """
    Convert `x`- and `y`- axes of input ``ax`` from pixels to polar angle and 
    Rsun, respectively.
    
    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The input axes.
    img_arr : `numpy.ndarray`
        The polar image array plotted in ``ax``.
    r_arr : `numpy.ndarray`
        The array of values of polar coordinate `Rsun`, at which the input 
        polar image array ``img_arr`` was calculated.
    phi_arr : `numpy.ndarray`
        The array of values of polar coordinate `phi`, at which the input 
        polar image array ``img_arr`` was calculated.
    labelfontsize : `float`, optional
        The fontsize of the axes labels. Default is 18.
    labelsize : `float`, optional
        The fontsize of the tick labels. Default is 15.
    
    """
    r1 = r_arr[0]
    r2 = r_arr[-1]
    dr = np.diff(r_arr)
    if np.all(np.isclose(dr, dr[0])):
        dr = dr[0]
    else:
        raise ValueError('Error. Provided r_arr is not an evenly spaced array.')
        
    phi1 = phi_arr[0]
    phi2 = phi_arr[-1]
    dphi = np.diff(phi_arr)
    if np.all(np.isclose(dphi, dphi[0])):
        dphi = dphi[0]
    else:
        raise ValueError('Error. Provided phi_arr is not an evenly spaced array.')

    ax.set_xticks([])
    secax_x = ax.secondary_xaxis(
        'bottom', 
        functions=(lambda x: phi1 + x*dphi, lambda x: (x-phi1)/dphi)
    )
    #secax_x.set_xticks(np.arange(phi1, phi2, 90))
    secax_x.set_xticks([phi1, 90, 180, 270, phi2])
    secax_x.tick_params(axis='x', which='major', labelsize=labelsize)

    ax.set_yticks([])
    secax_y = ax.secondary_yaxis(
        'left', 
        functions=(lambda y: r1 + y*dr, lambda y: (y-r1)/dr)
    )
    secax_y.set_yticks(np.linspace(r1, r2, 5))
    secax_y.tick_params(axis='y', which='major', labelsize=labelsize)

    secax_x.set_xlabel('Polar angle [deg]', fontsize=labelfontsize)
    secax_y.set_ylabel('$R_\odot$', fontsize=labelfontsize)

    
def plot_contours(ax, img_arr, vmin=None,  vmax=None):
    """
    Plot contours of ``img_arr`` plotted in ``ax``.
    
    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The input axes.
    img_arr : `numpy.ndarray`
        The image array plotted in ``ax``.
    vmin, vmax : `float`, optional
        If not None, the value range of counters. Default is None.
    
    """
    if vmin is None:
        vmin = np.nanmin(img_arr)
    if vmax is None:
        vmax = np.nanmax(img_arr)
    levels = np.linspace(vmin, 0.5*(vmax+vmin), 5)[1:]
    counters = ax.contour(img_arr, levels=levels, alpha=0.5, colors='k', 
                          linewidths=1, linestyles='dotted')
    ax.clabel(counters, inline=True, fontsize=10, fmt='%0.1e')
    
    
    
