import astropy.units as u
import numpy as np
import cdflib
import pandas as pd

def cdf2df(path):

    with cdflib.cdfread.CDF(path) as file:

        # Extract epoch times
        epoch = file.varget(variable='EPOCH', expand=False, to_np=True)

        # Extract epoch times
        CDF_epoch_class = cdflib.epochs.CDFepoch()
        time = CDF_epoch_class.to_datetime(epoch, to_np=True)

        # Extract B vectors times
        B = file.varget(variable='B_RTN', expand=False)
        norm = np.linalg.norm(B, axis=1)

        # Get data attributes
        attributes = file.globalattsget(expand=True)
        
        df = pd.DataFrame({'BR': B.T[0],
                        'BT': B.T[1],
                        'BN': B.T[2],
                        '|B|': norm}, index = time)
        
        return df

def PS_angle(distance, speed):
    """returns the parker spiral angle. This will be a negative number,
    and is only for positive polarity.

    Parameters
    ----------
    distance : u.Quantity
        Distance of the spacecraft from the Sun
    speed : u.Quantity
        Speed of the solar wind parcel

    Returns
    -------
    float
        Angle in degrees
    """
    # Solar rate of rotation
    omega = (14.713 * np.pi / 180.0) / (24 * 3600 * u.s)
    return (
        np.arctan((-1.0 * distance * omega) / (speed)).to(u.deg).value + 360
    )

def array_dot(r1, t1, r2, t2):
    return r1 * r2 + t1 * t2

def angle_between_components(r1, t1, r2, t2):
    return np.arccos(
        array_dot(r1, t1, r2, t2)
        / (np.sqrt(r1 * r1 + t1 * t1) * np.sqrt(r2 * r2 + t2 * t2))
    )

def polarity_at_sc(r, t, distance, speed, tolerance=45):
    """Works out if magnetic field falls within the Parker Spiral given some tolerance

    Parameters
    ----------
    r : array
        R component
    t : array
        T component
    distance : array with units
        Distance from Sun
    speed : array with units
        Speed of Solar Wind
    tolerance : float, optional
        In degrees, by default 45

    Returns
    -------
    array
        returns an array with -1 for negative polarity and +1 for positive polarity
    """

    # angle of the outwards (positive) Parker spiral
    ps_pos_angle = PS_angle(distance, speed)
    # angle of inwards (negative) polarity Parker spiral
    ps_neg_angle = ps_pos_angle - 180

    # e.g. ps_pos_angle will be around 330

    ps_pos_r = np.cos(ps_pos_angle * np.pi / 180)
    ps_pos_t = np.sin(ps_pos_angle * np.pi / 180)
    ps_neg_r = np.cos(ps_neg_angle * np.pi / 180)
    ps_neg_t = np.sin(ps_neg_angle * np.pi / 180)

    # angle of magnetic field to positive Parker spiral
    angle2pos = angle_between_components(r, t, ps_pos_r, ps_pos_t)
    # angle of magnetic field to negative Parker spiral
    angle2neg = angle_between_components(r, t, ps_neg_r, ps_neg_t)

    is_pos = (abs(angle2pos) * 180 / np.pi < tolerance).astype(int)
    is_neg = (abs(angle2neg) * 180 / np.pi < tolerance).astype(int) * -1
    return is_pos + is_neg

def add_polarity2df(df, ds_period="12H", tolerance=45):
    df_downsampled = df.resample(ds_period).mean()

    df_downsampled["polarity"] = polarity_at_sc(
        df_downsampled["BR"].values,
        df_downsampled["BT"].values,
        df_downsampled["Radius"].values * u.au,
        df_downsampled["V"].values * u.km / u.s,
        tolerance=tolerance,
    )

    # I want to forward fill, so I need to shift the indices back half a window length
    df_downsampled = df_downsampled.shift(periods=int(ds_period[:2])/2, freq='H')

    # now upsample back to df
    df_upsampled = df.combine_first(df_downsampled)
    # # fill in nans
    df_upsampled.fillna(method="ffill", inplace=True)
    df["polarity"] = df_upsampled["polarity"]
    return df

def unwrap_lons(arr, threshold=0):
    try:
        idx = np.argwhere(np.diff(arr) > threshold).flatten()[0]
    except IndexError:
        return arr
    if isinstance(arr[0], u.Quantity):
        val_to_subtract = 360 * u.deg
    else:
        val_to_subtract = 360
    arr[idx + 1 :] -= val_to_subtract
    return arr