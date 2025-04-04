#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Common functions to prepare TPM and to handle/plot/utilize the results.

Read document carefully!
JPL/HORIZONS: https://ssd.jpl.nasa.gov/horizons/manual.html
astroquery  : https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html

- About location
location is like '(MPCcode)' (e.g., 381), 
                 '(MPCcode)@(399, which means Earth)' (e.g., 381@391)
The latter part is to specify major body 
    (Earth '391', Sun '10', Moon '301', solar system barycenter '0' or 'ssb' etc.)
The former part is to specify point on the major body (body center 500 etc.)

J.B. notes that we have to set '@10', or 'None' (not '@0' or 'ssb')
since  the center of Ecliptic coordinate system, 
which is used in Marco's TPM code, is the Sun (not solar system barycenter).


location refers to the coordinate center for the ephemeris, which has 
slightly different physical interpretations depending on the query type: 
  - observer (ephemerides) queries
  - observer location vectors queries
  - (coordinate origin for vectors elements queries)
  - (relative body for orbital elements)

The default value is None for all queries, which corresponds to
  - Earth body center for observer (ephemerides) queries 
    (i.e., location=500, 500@391)
  - Sun body center for orbital elements and vectors queries 
    (i.e., location=@10, 500@10)


2024-05-15 J.B. confirmed that the generated ephem file with @10 
           matched that in DAFEED better than that with @0. 
           J.B. confirmed that the generated obs file with light time correction
           matched that in DAFEED with an accuracy of 0.01 s.
"""
from astroquery.jplhorizons import Horizons
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# To be removed
from scipy import interpolate
from astropy.constants import c, au

# Move to different module ? ==================================================
# e.g., visualization.py
# Color
mycolor = ["#AD002D", "#1e50a2", "#006e54", "#ffd900", 
           "#EFAEA1", "#69821b", "#ec6800", "#afafb0", "#0095b9", "#89c3eb"] 
mycolor = mycolor*100
# Linestyle
myls = ["solid", "dashed", "dashdot", "dotted", (0, (5, 3, 1, 3, 1, 3)), 
        (0, (4,2,1,2,1,2,1,2))]
myls = myls*100
# Marker
mymark = ["o", "^", "x", "D", "+", "v", "<", ">", "h", "H"]
mymark = mymark*100
# Move to different module ? ==================================================


def make_ephemfile(asteroid, df, out, warmuptime_day=30):
    """
    Make ephem file for TPM.
    Light time correction is unnecessary.

    Parameters
    ----------
    asteroid : str
        name of asteroid
    df : pandas.DataFrame
        jd, wavelength, flux, fluxerr, code, and cflag
    out : str
        output ephem filename
    warmuptime_day : float, optional
        time for warming up in day
    """
    # Light-time correction is needless! (private com. with Marco DELBO, July 19 2024)

    # Sort by jd
    jds = sorted(set(df["jd"]))

    # Make chunks. A chunk is composed of obs. with shorter arc
    # This is just to speed up the process.
    # th_arc (50 days) 
    #   i.e., if observation gap is larger than 70 days, 
    #         they are regarded as different chunks
    # Note: th_arc should be larger than margin0 + margin1, otherwise
    #       ephmeris file could be not always ascending order
    th_arc = 50
    jds_chunk = []
    
    # First day
    d0, d1 = jds[0], jds[0]
    for idx_d, d in enumerate(jds):
        if d - d0 < th_arc:
            d1 = d
            # For final day
            if idx_d == len(jds)-1:
                # Save first and last days as a chunk
                ch = (d0, d1)
                jds_chunk.append(ch)
        else:
            # Save first and last days as a chunk
            ch = (d0, d1)
            jds_chunk.append(ch)

            # Update first day
            d0, d1 = d, d

            # For final day
            if idx_d == len(jds)-1:
                # Save first and last days as an independent chunk
                ch = (d, d)
                jds_chunk.append(ch)


    # Make ephem file widh jds_chunk ==========================================
    # Use 1 day time step. Each obs. is linked by interpolation.
    # Just to speed up,  all jds are not used.
    # t0 of chunk 1: (minimum jd of chunk 1) - 50, -50 is margin to reach thermal equilibrium.
    #                Note that this may be enough for slow rotators.
    # t1 of chunk 1: (maximum jd of chunk 1) + 1, 1 is just in case
    #                And add 1 flag to move to next chunk (i.e., t0 of chunk 2)
    # t0 of chunk 2: (minimum jd of chunk 2) - 50, -50 is margin to reach thermal equilibrium.
    # t1 of chunk 2: (maximum jd of chunk 2) + 1, 1 is just in case
    #                And add 1 flag to move to next chunk (i.e., t0 of chunk 3)
    # ......
    # 
    # t0 of final chunk : (minimum jd of final chunk) - 50, -50 is margin to reach thermal equilibrium.
    # t1 of final chunk : (maximum jd of final chunk) + 1, 1 is necessary to run a tpm I don't know why.

    # Note: margin0 of 30 corresopnds to 60 rotations for object wiht rotation period of 12 hr.
    margin0, margin1 = warmuptime_day, 1
    with open(out, "w") as f_eph:
        for idx_ch, ch in enumerate(jds_chunk):
            d0, d1 = ch
            # Add margins
            t0 = d0 - margin0
            t1 = d1 + margin1

            d_list = np.arange(t0, t1, 1)
            for idx, d in enumerate(d_list):
                # Location of @10 means Sun body center (=None) for vectors queries
                # (not solar system barycenter, @0, @ssb).
                S = Horizons(location="@10", id=asteroid, epochs=d)

                # 2024-10-21 for test
                #S = Horizons(location="@0", id=asteroid, epochs=d)

                vec = S.vectors(refplane="ecliptic")
                # Vector from the Sun to asteroid (Sun -> ast)
                x_S, y_S, z_S = vec["x"][0], vec["y"][0], vec["z"][0]
                # Vector from asteroid to the Sun (ast -> Sun)
                x_S, y_S, z_S = -x_S, -y_S, -z_S
                # Last raw should end with 1 to skip (or end) the calculations.
                if (idx == len(d_list)-1):
                    f_eph.write(f"{d} {x_S} {y_S} {z_S} 1\n")
                else:
                    f_eph.write(f"{d} {x_S} {y_S} {z_S}\n")
    # Make ephem file widh jds_chunk ==========================================
     

def make_obsfile(asteroid, df, out, lccor=False, rmnegativeflux=False):
    """
    Parameters
    ----------
    asteroid : str
        name of asteroid
    df : pandas.DataFrame
        jd, wavelength, flux, fluxerr, code, and cflag
    out : str
        output obs filename
    lccor : bool, optional
        wheather perform light-time correction
    """
    # Number of data block or epoch
    N_epoch = len(set(df["jd"]))

    # Sort by jd
    jds = sorted(set(df["jd"]))
    # 0,   1,  2,  3,       4?,     5?,      6?,       7?,       8,        9
    # W1, W2, W3, W4, IRAS 12, IRAS 25, IRAS 60, IRAS 100, Akari 9, Akari 18
    cflag_registered = ["0", "1", "2", "3", "8", "9", ]

    # Make obs file ===========================================================
    with open(out, "w") as f_obs:
        f_obs.write(f"{N_epoch}\n")
        f_obs.write(f"\n")
        for jd in jds:
            # Extract 1 data point. (i.e., merge data points if it is spectrum)
            df_samejd = df[df["jd"]==jd]
            df_samejd = df_samejd.reset_index(drop=True)
            # Number of observations
            N_data = len(df_samejd)

            # Common parameters for all data (i.e., spectrum)
            code, cflag   = df_samejd["code"][0], df_samejd["cflag"][0]
            if str(cflag) in cflag_registered:
                pass
            else:
                cflag = ""

            
            # Location of @10 means Sun body center (=None) for vectors queries
            # (not solar system barycenter, @0, @ssb).
            S = Horizons(location="@10", id=asteroid, epochs=jd)
            # 2024-10-21 for test
            #S = Horizons(location="@0", id=asteroid, epochs=jd)

            vec = S.vectors(refplane="ecliptic")
            # Vector from the Sun to asteroid (Sun -> ast)
            x_S, y_S, z_S = vec["x"][0], vec["y"][0], vec["z"][0]
            # Vector from asteroid to the Sun (ast -> Sun)
            x_S, y_S, z_S = -x_S, -y_S, -z_S

            # Location of code means observer center
            E = Horizons(location=code, id=asteroid, epochs=jd)
            vec = E.vectors(refplane="ecliptic")
            eph = E.ephemerides()

            # Vector from the Earth to asteroid (Earth -> ast)
            x_E, y_E, z_E = vec["x"][0], vec["y"][0], vec["z"][0]
            # Vector from asteroid to the Earth (ast -> Earth)
            x_E, y_E, z_E = -x_E, -y_E, -z_E
            
            if lccor:
                # Lighttime correction
                # Geocentric distance in au
                delta = eph["delta"][0]
                # c: speed of light in m/s
                # au: astronomical unit in m
                c_au_s = c.value/au.value
                c_au_day = c_au_s*24.*3600.
                ltcor = delta/c_au_day
                jd_ltcor = jd - ltcor
            else:
                print("    Do not perform light time correction!")
                jd_ltcor = jd
            
            # Count number of negative flux
            if rmnegativeflux:
                N_nega = 0
                for idx, row in df_samejd.iterrows():
                    flux = row["flux"]
                    if flux < 0:
                        N_nega += 1
                N_data -= N_nega

            # Write 
            #   jd, Ndata
            #   Sun coordinate
            #   Earth coordinate
            f_obs.write(f"{jd_ltcor} {N_data}\n")
            f_obs.write(f"{x_S} {y_S} {z_S}\n")
            f_obs.write(f"{x_E} {y_E} {z_E}\n")
            
            for idx, row in df_samejd.iterrows():
                w             = row["wavelength"]
                flux, fluxerr = row["flux"], row["fluxerr"]
                # Do not use negative flux
                if (rmnegativeflux) and (flux < 0):
                    continue
                f_obs.write(f" {w} {flux} {fluxerr} {cflag}\n")

            f_obs.write("\n")
    # Make obs file ===========================================================


def extract_flux(f0, fixscale=False):
    """
    Extract thermal flux from output of TPM.

    Parameters
    ----------
    f : str
        output file of TPM
    fixscale : bool, optional
        whether fix the scale factor (i.e., trust shape model)

    Return
    ------
    df : pandas.DataFrame
        dataframe with extracted fluxes etc.
    """

    epoch_list, jd_list, w_list = [], [], []
    f_obs_list, ferr_obs_list, f_model_list = [], [], []
    with open (f0, "r") as f:
        f = f.readlines()
        for l in f:
            # l is as below:
            # f> 052     2454199.2663881136   000 18.000        1.698843   0.113268    1.441020    0.000000      1.1789
            #    epoch   JD                   n   wavelength    f_obs      ferr_obs    f_model     ferr_model    f_obs/f_model
            if l[0:2] == "f>":
                # Extract info
                l = l.split()
                epoch, jd, n       = l[1], l[2], l[3]
                w, f_obs, ferr_obs = l[4], l[5], l[6]
                f_model            = l[7]
                epoch_list.append(float(epoch))
                jd_list.append(float(jd))
                w_list.append(float(w))
                f_obs_list.append(float(f_obs))
                ferr_obs_list.append(float(ferr_obs))
                f_model_list.append(float(f_model))
            elif l[0:2] == "r>":
                # Extract scale factor
                # r>     100.0   0.0   0.0  0.00  0.12 1.15125549     19.364     21.576  1.000   2259.587    4.18547
                l = l.split()
                # Fix scale factor to 1
                if fixscale:
                    scalefactor = 1
                # Use scale factor in TPM output
                else:
                    scalefactor = float(l[6])

                # This is useful to check intputs of the TPM calculation
                # Haple angle in deg (t.Bar in the output) 
                TI = float(l[1])
                Htheta = float(l[2])
                A = float(l[5])

            else:
                continue

        # DataFrame
        df = pd.DataFrame({
            "epoch": epoch_list,
            "jd": jd_list,
            "w": w_list,
            "f_obs": f_obs_list,
            "ferr_obs": ferr_obs_list,
            "f_model": f_model_list,
            })
        df["scalefactor"] = scalefactor
        df["TI"] = TI
        df["Htheta"] = Htheta
        df["A"] = A
        #print(f"Scafe factor = {scalefactor}")
    return df


def introduce_var_scalefactor(df, key_t="jd", sf0=0.80, sf1=1.2, sfstep=0.01):
    """
    Introduce variable scale factors per observation.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    key_t : str
        keyword for observation time
    sf0 : float
        minimum scale factor        
    sf1 : float
        maximum scale factor        
    sfstep : float
        step of scale factor

    Return
    ------
    df1 : pandas.DataFrame
        output dataframe with new scale factors
    """

    # Updated dataframe
    df1_list = []

    # Extract observation time and N
    t_list = list(set(df[key_t]))
    # Time in which scale factor to be introduced
    t_cor_list = []
    for t in t_list:
        df_t = df[df[key_t] == t]
        N_t = len(df_t)
        #print(f"{key_t}={t} N={N_t}")
        # Do not introduce scale factor when N = 1 (i.e., photometry)
        if N_t > 1:
            t_cor_list.append(t)
        else:
            # Add into list of updated dataframe
            df1_list.append(df_t)

    N_sf = len(t_cor_list)
    print(f"  Number of scale factors = {N_sf}")
     
    # Determine scale factor
    # List of scale factor to be searched
    sf_list = np.arange(sf0, sf1+sfstep, sfstep)
    for idx_t, t_cor in enumerate(t_cor_list):
        df_t = df[df[key_t] == t_cor]

        # Minimize chi2
        for idx_sf, sf in enumerate(sf_list):
            df_t["scalefactor"] = sf
            if idx_sf == 0:
                # Calculate chi-squared
                # These are bench marks
                chi2 = calc_chi2(df_t)
                scalefactor = sf
            else:
                # Calculate chi2
                chi2_sf = calc_chi2(df_t)
                if chi2_sf < chi2:
                    # Update
                    chi2 = chi2_sf
                    scalefactor = sf
                else:
                    pass
            #print(f"chi2, sf = {chi2:.2f}, {scalefactor}")
        #print(f"{idx_t+1}-th scale factor {scalefactor}")
        df_t["scalefactor"] = scalefactor
        # Do not save for the consistency 
        # (N=1 data has no such column)
        df_t = df_t.drop(["diff"], axis=1)

        df1_list.append(df_t)

    # Make a new dataframe with updated scale factors
    df1 = pd.concat(df1_list)
    # Sort by time
    df1 = df1.sort_values(by=[key_t])
    return df1


    # 2. This is faster =======================================================
    df_blend = df1.copy()
    df_blend["f_model"] = (
        alpha*df1["scalefactor"]**2*df1["f_model"] + 
        (1-alpha)*df2["scalefactor"]**2*df2["f_model"])
    # This is a dummy
    df_blend["scalefactor"] = 1
    # 2. This is faster. =======================================================

    return df


def plot_flux(dfs, fixscale=False, y1range=None, out=None):
    """
    Plot obs and model fluxes.
    Note: model fluxes of input files should be the same.

    Parameters
    ----------
    dfs : array-like
        preprocessed dataframes with fluxes
    out : str, optional
        output filename
    """
    fig = plt.figure(figsize=(24, 16))
    # idx vs. flux (obs. and model, normalized by obs.)
    ax1 = fig.add_axes([0.05, 0.70, 0.90, 0.25])
    # residual
    ax2 = fig.add_axes([0.05, 0.45, 0.90, 0.25])
    # chi2 component 
    ax3 = fig.add_axes([0.05, 0.20, 0.90, 0.25])

    ax3.set_xlabel("Index")
    ax1.set_ylabel("Normalized flux")
    if fixscale:
        ax2.set_ylabel("(f_obs - f_model)/ferr_obs")
    else:
        ax2.set_ylabel("(f_obs - sf**2*f_model)/ferr_obs")
    ax3.set_ylabel("Chi2 component")

    for idx_df, df in enumerate(dfs):
        col = mycolor[idx_df]
        idx_epoch = 0
        chi2 = 0

        for idx_row, row in df.iterrows():

            epoch    = row["epoch"]
            f_obs    = float(row["f_obs"])
            ferr_obs = float(row["ferr_obs"])
            f_model  = float(row["f_model"])
            s        = float(row["scalefactor"])

            # Fit with the size
            if fixscale:
                pass
            else:
                f_model = f_model*s**2

            # Set marker
            if idx_row == 0:
                pass
            else:
                # Change marker for different epoch
                if epoch != epoch0:
                    idx_epoch += 1

            col = mycolor[idx_epoch]
            mark = mymark[idx_epoch]

            # f_obs and s**2 f_model
            # ax.set_ylabel("Flux [Jy]")
            #ax.scatter(idx_row, f_obs, color=col, marker=mark)
            #ax.errorbar(idx_row, f_obs, ferr_obs, color=col)
            #ax.scatter(idx_row, s**2*f_model, color=col, marker=mark)
            
            # Scale factor is already taken into account in f_model above
            # f_obs - s**2 f_model/ferr_obs
            if f_obs > 0:
                res = (f_obs - f_model)/ferr_obs
                chi2 += res**2
            print(f"  f_obs, f_model, res = {f_obs}+-{ferr_obs}, {f_model}, {res}")

            if idx_row == (len(df)-1):
                label = f"Idx {idx_df}: chi2 = {chi2:.3f}"
                label_obs = f"Observations (normalized by themselves, i.e., 1)"
                label_model = f"Models (normalized by observations)"
            else:
                label = None
                label_obs = None
                label_model = None

            # Observation
            ax1.errorbar(idx_row, f_obs/f_obs, ferr_obs/f_obs, color=col, ms=0, marker=mark, label=label_obs)
            # Model
            ax1.scatter(idx_row, f_model/f_obs, color="black", s=5, marker="o", label=label_model)

            # Residual
            ax2.scatter(idx_row, res, color=col, marker=mark, label=label)
            # Chi2 component
            ax3.scatter(idx_row, res**2, color=col, marker=mark, label=label)

            # Update previous epoch
            epoch0 = epoch
    
    if y1range:
        ax1.set_ylim([y1range[0], y1range[1]])
    ymin, ymax = ax2.get_ylim()
    if abs(ymin) > abs(ymax):
        ax2.set_ylim([ymin, -ymin])
    else:
        ax2.set_ylim([-ymax, ymax])

    ax1.legend(ncol=3)
    ax2.legend(ncol=3)
    if out:
        plt.savefig(out)


def crater2Htheta(c_angle, c_density):
    """
    Calculate Haple roughness parameter with 
    semiaperture angle of craters and crater surface density.
    See Hapke 1984, Icarus.
    # Note: This function is no longer useful since
    #       I realised that Htheta was written in the output of TPM!

    Parameters
    ----------
    c_angle : float
        semiaperture angle of craters in degree
    c_density : float
        crater surface density (0 to 1)

    Return
    ------
    Htheta : float
        Haple roughness parameter in degree
        
    """
    # TODO: Calculate theta_bar by myself
    # How to calculate theta_bar?
    Htheta = c_density**0.5*np.tan(np.radians(c_angle))

    # This is temporary one.
    # From Hung+2022, PSJ, Table 2
    # Convert c_angle to theta_bar
    # theta_bar is 12.0 for (gamma, rhoC) = (50.0, 0.50), not 22.0
    # (i.e., Hung+2022 is correct.)
    d = {0.0:0.0, 30.0:3.9, 40.0:12.6, 41.0:16.5, 50.0:12.0, 
         60.0:26.7, 70.0:27.3, 88.0:35.8, 89.0:46.8, 90.0:55.4}
    Htheta = d[c_angle]

    return Htheta
 

def extract_bestparam(df, key_chi2, params):
    """
    Search and return best fit parameters.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    key_chi2 : str
        key of chi-squared
    params : array-like
        key of parameters of interest 

    Returns
    -------
    chi2_min : float
        minimum chi-squared
    best_params : array-like
        best fit parameters and chi2 minimum
    """
    # Extract minimum chi2 and its index
    idx_min = df[key_chi2].idxmin()
    chi2_min = df.loc[idx_min, key_chi2]

    best_params = []
    for key in params:
        best_param = df.loc[idx_min, key]
        best_params.append(best_param)
    return chi2_min, best_params



def blend_flux(df1, df2, alpha):
    """
    Blend fluxes if df1 and df2 with alpha.
      F_blended = alpha*F1*s1^2 + (1 - alpha)*F2*s2^2,
    where s1 and s2 are scale factors.

    Parameters
    ----------
    df1 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    df2 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    alpha : float
        regolith abundance

    Return
    ------
    df_blend : pandas.DataFrame
        dataframe with blended results
    """
    # Sanity check
    assert len(df1) == len(df2), "Check if input dfs are the same dimension!"
    
    # 1. This is safe but slow, because of the for loop. ======================
    # Number of data points
    #N_data = len(df1)
    #epoch_list, jd_list, w_list = [], [], []
    #f_obs_list, ferr_obs_list, f_model_list = [], [], []
    #
    #for n in range(N_data):
    #    epoch1    = df1.loc[n, "epoch"]
    #    jd1       = df1.loc[n, "jd"]
    #    w1        = df1.loc[n, "w"]
    #    f_obs1    = df1.loc[n, "f_obs"]
    #    ferr_obs1 = df1.loc[n, "ferr_obs"]
    #    f_model1  = df1.loc[n, "f_model"]
    #    s1        = df1.loc[n, "scalefactor"]

    #    epoch2    = df2.loc[n, "epoch"]
    #    jd2       = df2.loc[n, "jd"]
    #    w2        = df2.loc[n, "w"]
    #    f_obs2    = df2.loc[n, "f_obs"]
    #    ferr_obs2 = df2.loc[n, "ferr_obs"]
    #    f_model2  = df2.loc[n, "f_model"]
    #    s2        = df2.loc[n, "scalefactor"]

    #    # Sanity checks
    #    assert jd1 == jd2, "Check if input dfs are made from the TPMs with same obs file."
    #    assert w1 == w2, "Check if input dfs are made from the TPMs with same obs file."
    #    assert f_obs1 == f_obs2, "Check if input dfs are made from the TPMs with same obs file."
    #    assert ferr_obs1 == ferr_obs2, "Check if input dfs are made from the TPMs with same obs file."

    #    # Calculate blended flux
    #    f_blend = alpha*s1**2*f_model1 + (1 - alpha)*s2**2*f_model2

    #    # Save info.
    #    epoch_list.append(epoch1)
    #    jd_list.append(epoch1)
    #    w_list.append(w1)
    #    f_obs_list.append(f_obs1)
    #    ferr_obs_list.append(ferr_obs1)
    #    f_model_list.append(f_blend)
    #    
    ## DataFrame
    #df_blend = pd.DataFrame({
    #    "epoch": epoch_list, 
    #    "jd": jd_list, 
    #    "w": w_list,
    #    "f_obs": f_obs_list, 
    #    "ferr_obs": ferr_obs_list, 
    #    "f_model": f_model_list, 
    #    })
    ## This is a dummy
    #df_blend["scalefactor"] = 1
    # 1. This is safe but slow, because of the for loop. ======================

    # 2. This is faster =======================================================
    df_blend = df1.copy()
    df_blend["f_model"] = (
        alpha*df1["scalefactor"]**2*df1["f_model"] + 
        (1-alpha)*df2["scalefactor"]**2*df2["f_model"])
    # This is a dummy
    df_blend["scalefactor"] = 1
    # 2. This is faster. =======================================================

    return df_blend


def calc_chi2(df):
    """
    Calculate chi2 with F_obs, Ferr_obs, F_model, and scale factor.
    Note: 
    The output chi2 could be different from that in the output of TPM.
    This is because runtpm ignores negative fluxes. J.B. thinks negative 
    fluxes are still informative if they are with large errorbars.

    Parameter
    ---------
    df : pandas.DataFrame
        dataframe with fluxes

    Return
    ------
    chi2 : float
        calculated chi-square
    """
    df["diff"] = (df["f_obs"] - df["scalefactor"]**2*df["f_model"])**2/df["ferr_obs"]**2
    chi2 = np.sum(df["diff"])
    return chi2


def search_regolith_abundance(df1, df2, alpha_list, chi2_min=10000, minonly=False):
    """
    Search regolith abundance alpha which minimize chi2.

    Parameters
    ----------
    df1 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    df2 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    alpha_list : array-like
        list of regolith abundance
    chi2_min : float
        initial chi2 minimum
    minonly : bool
        return minimum chi2 and corresponding alpha (i.e., fit by alpha)

    Returns
    -------
    alpha_arr : float
        array of regolith abundance 
    chi2_arr : float
        array of chi2
    """
    alpha_arr, chi2_arr = [], []
    for a in alpha_list:
        # Blend flux as 
        #   F = alpha*F_regolith*s1^2 + (1-alpha)*F_rock*s2^2,
        # where s1 and s2 are scale factors.
        df_blend = blend_flux(df1, df2, a)
        # Calculate chi2 of blended flux
        chi2 = calc_chi2(df_blend)

        if minonly:
            if chi2 < chi2_min:
                chi2_arr = [chi2]
                alpha_arr = [a]
                chi2_min = chi2
            else:
                pass
        else:
            alpha_arr.append(a)
            chi2_arr.append(chi2)
    return alpha_arr, chi2_arr


def calc_C_coord(phi):
    """
    Calculate coordination number (see Sakatani+2018)

    Parameter
    ---------
    phi : float
        macroscopic porosity
        
    Return
    ------
    C_coord : float
        coordination bumber C 
    """
    f = 0.07318 + 2.193 * phi
    C_coord = 2.812 * (1 - phi)**(-1./3.) / (f**2 * (1 + f**2))
    return C_coord


def kappa_Sakatani2018(kappa, D_p, phi, r_c, xi):
    """
    Calculate bulk thermal conductivity 
    with model in Sakatani+2018, Icarus, 309, 13.

    Parameters
    ----------
    kappa : float
        thermal conductivity of solid material
    D_p : float
        particle diameter
    phi : float
        macroscopic porosity
    r_c : float
        radius of the contact area between the spheres
    xi : float
        degree of reduction of the thermal conductance at the contacts 
        owring to the microscopic roughness of the particle surfaces

    Return
    ------
    kappa_bulk : float
        bulk thermal conductivity
    """
    # Particle radius
    R_p = D_p/2.

    # Ratio of 
    # [the effective distance of radiative heat transfer in the voids between particles] 
    # to [the geometric size of the voids]
    # used in Cambioni+2021
    ## 0.68 + 7.6e-5 / D_p
    ## ???
    ## How to calculate r_c?? Assume Young's modulus, Poisson's ratio?
    r_c = xxxx

    # Calculate Coordination number C with phi
    C_coord = calc_C_coord(phi)

    # Equation (8) in Sakatani+2018
    kappa_bulk = 4 / np.pi**2 * kappa * (1 - phi) * C_coord * xi * r_c/R_p
    return kappa_bulk


def calc_Phi(TI_rock, c_p, rho_s):
    """
    Calculate microscopic porosity Phi.

    Parameters
    ----------
    TI_rock : float
        thermal inertia of rock 
    c_p : float
        heat capacity
    rho_s : float
        grain density

    Return
    ------
    Phi : float
        microscopic porosity
    """
    # TI_rock**2/(c_p rho_s (1-Phi)) = 0.11(1-Phi)/Phi
    # Define C = (TI_rock**2)/(0.11 c_p rho_s) + 2 and solve the equation
    # (See note on Cambioni+2021 by JB)
    C = TI_rock**2/(0.11*c_p*rho_s) + 2
    # Solution satisfying 0 < Phi < 1
    Phi = (C - (C**2 - 4)**0.5) / 2
    return Phi
    

def calc_TI_th(TI_rock):
    """
    Calculate threshold of thermal inertia of regolith and rock
    (gamma_c in Cambioni+2021).
    gamma_c is defined as thermal inertia when D_p = l_s,
    where D_p is particle diameter and l_s is thermal skin depth.

    Parameters
    ----------
    TI_rock : float
        thermal inertia of rock

    Return
    ------
    TI_th : float
        threshold of thermal inertia of regolith and rock
    """

    # Fixed parameters in Cambioni+2021 =======================================
    # Macroporosity of Bennu 
    phi = 0.4
    # Grain density of CM meteorites in kg/m^3
    rho_s = 2920
    # Heat capacity for meteorite CM2 Cold Bokkeveld at OTES spot's 
    # mean diurnal temperature (from Figure 3 in Opeil+2020, MPS)
    # TODO: what is mean diurnal temperature? Read by eye?
    c_p = 999
    # Degree of reduction of the thermal conductance at the contacts 
    # owring to the microscopic roughness of the particle surfaces
    # (from Sakatani+2018, Icarus ?)
    xi = 0.12
    # Fixed parameters in Cambioni+2021 =======================================


    # 1. Derive microscipic porosity (microporosity) Phi
    #    with TI_rock, c_p (heat capacity), and rho_s (grain density)
    Phi = calc_Phi(TI_rock, c_p, rho_s)

    # 2. Derive thermal conductivity kappa
    kappa = TI_rock**2/(c_p * rho_s * (1 - Phi))

    # Density considering macroporosity phi and microscopic porosity Phi
    rho = rho_s * (1 - Phi) * (1 - phi)

    # 3.Derive regolith bulk thermal conductivity, kappa_p,
    # with Standard Regolith Model in Sakatani+2018
    # TODO: This kappa is dependent on D_p?
    kappa_p_list = kappa_Sakatani2018(kappa)
    #       Find k_p where D_p = l_s?
    kappa_p = 999


    # 4. Finally derive threshold of thermal inertia
    TI_th = (kappa_p * c_p * rho)**0.5
    return TI_th


def calc_confidence_chi2(paper, chi2_min, dof, n, reduce):
    """
    Calculate condifence levels.
    There are several ways to estimate confidence levels,
    and errors of physical properties using the confidence levels.

    1. P14
       Polishook 2014 Icarus, 241, 79, Cambioni+2021, Nature, 
       chi2 = chi2_min + n * (2/nu)**0.5 (for n-sigma)

    2. V17
       Vokrouhlicky+2017, AJ, Hanus+2018, A&A, Durech+2018 A&A, etc.
       chi2 = chi2_min * (1 + n*(2/nu)**0.5)
            = chi2_min + chi2_min*n*(2/nu)**0.5 (for n-sigma)
       Note: 
           Cambioni+2019 as well with n=1, Equatioin (4) is a typo,
           private com. with Saverio on 2025-01-09

    Since normally chi2_min > 1, 1. gives smaller confidence levels and smaller 
    uncertainties of physical proproperties compared to 2.

    Parameters
    ----------
    paper : str
        P14 or V17
    chi2_min : float
        minimum chi-squared
    dof : int
        degrees of freedom
    n : int
        return n-sigma interval
    reduce : bool
        whether reduced chi-squared or not

    Return
    ------
    chi2_sigma : 
        chi-squared boundary with n-sigma confidence level
    """

    if reduce:
        chi2_sigma = n*np.sqrt(2*dof)/dof
    else:
        chi2_sigma = n*np.sqrt(2*dof)

    # 1. Polishook 2014, Cambioni+2021, 
    if paper == "P14":
        pass
    # 2. Vokrouhlicky+2017, AJ, Hanus+2018, A&A, Durech+2018 A&A, etc.
    elif paper == "V17":
        chi2_sigma = chi2_sigma*chi2_min

    return chi2_sigma
