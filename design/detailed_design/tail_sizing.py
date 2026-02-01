""" tail sizing 
cd {localWorkSpace}
pip install -e ./tools/aero_techt.py
python3 design/tail_sizing.py
"""
import os
import warnings
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from py_planes.avl import Aircraft as AvlAircraft
from py_planes.avl import Avl, get_coef
from mpl_toolkits.mplot3d import Axes3D

# TODO: #9 use classes so that we don't have issues with global and local variables
# TODO: #7 add pickle file implementation so we dont have to run avl over and over again
# we can pickle the returns/ the data that gets plotted
# add if os.pathexists(pickle_path) to the plotting
# else: write to pickle

def plt_3d(x,y,z,title,xlabel,ylabel,zlabel,
           c=None, clabel=None):
    """ plot 3d points """
    plt.figure()
    ax=plt.gcf().add_subplot(111,projection='3d')
    ax.set_title(title)
    if c is None:
        sc=ax.scatter(x,y,z)
    else:
        sc=ax.scatter(x,y,z,c=c,cmap='magma')
        cb=plt.gcf().colorbar(sc, ax=ax, pad=0.1)
        if clabel is None:
            # set to zlabel if not given
            cb.set_label(zlabel)
        else:
            cb.set_label(clabel)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    return sc

def ht_sizing_vht_lt_analysis(min_vht, max_vht, nvht, arht, tail_arm, wing, tail, mass,
                        plot=None):
    """ size horizontal tail based off tail volume and tail arm
    Returns:
        vht_arr: array of horizontal tail volumes (np.array)

        xnps_arr: array of neutral point locations (np.array)
    """
    vht_arr=np.linspace(min_vht, max_vht, nvht)

    xnps=[]
    for vht in vht_arr:
        # iterate through each tail volume in range
        # get tails dimensions and overwrite given dictionary
        sht=vht*wing['chord']*wing['span']/tail_arm
        bht=np.sqrt(arht*sht)
        cht = bht / arht
        tail['h_span']=bht
        tail['h_chord']=cht
        tail['l_tail']=tail_arm
        # create a plane from dimensions and an avl file
        plane=AvlAircraft(wing, tail, mass)
        avl=Avl(f'gen_vht_{vht:.4f}_lt_{tail_arm}.avl',
                plane,
                'Avl',
                'runs/tail_sizing/runs')
        avl.create_avl_file()
        # run the steps from the template associated
        # TODO: change from tail_sizing to horizontal_tail_sizing
        # remove file that may be there already
        file_out_og=f'Avl/runs/tail_sizing/stab/{avl.name}_st.txt'
        if os.path.exists(file_out_og):
            # remove file
            os.remove(file_out_og)
            print(f'removing file: {file_out_og}')
        
        avl.run_avl('Avl/automation/steps_template_tail_sizing.txt')
        # rename the generically named file from the template to something
        # a little more descriptive of our current use cases
        file_out_name=f'Avl/runs/tail_sizing/stab/{avl.name}_{vht:.4f}_st.txt'
        os.rename(
            file_out_og,
            file_out_name)
        # open file, and get its TextIO object to get coefficient from
        with open(file_out_name, 'r') as stab_file:
            file_out = stab_file
            xnp=get_coef(file_out, 'Xnp')
        xnps.append(xnp)
    
    if not plot:
        # return what we need
        return np.asarray(vht_arr), np.asarray(xnps)
    elif plot is True:
        f=plt.figure()
        a=f.add_subplot(111)
        a.plot(vht_arr, xnps, marker='o', linestyle='-')
        a.set_xlabel('horizontal tail volumes')
        a.set_ylabel('neutral point locations (aft leading edge), m')
        a.grid()


def ht_sizing_lt_vht_analysis(ltmin,ltmax,nlt,vhtmin,vhtmax,nvht,arht,wing,tail,mass,
                              plot=None):
    """ horizontal tail, tail arm and tail volume analysis
    Returns:
        lt2plot: 1xn array of tail arms (np.array)

        vht_arr: 1xn array of horizontal tail volumes (np.array)

        xnp_arr: 1xn array of neutral point locations (np.array)
    """
    lt_array=np.linspace(ltmin, ltmax, nlt)
    lt2plot=[]
    vht_arr=[]
    xnp_arr=[]
    for lt in lt_array:
        vhts, xnps=ht_sizing_vht_lt_analysis(vhtmin, vhtmax, nvht, arht, lt, wing, tail, mass)
        # for xnp in xnps:
        #     xnps_arr.append(xnp)
        for vht, xnp in zip(vhts, xnps):
            lt2plot.append(lt)
            vht_arr.append(vht)
            xnp_arr.append(xnp)
    if not plot:
        # return values we need from here
        return np.asarray(lt2plot), np.asarray(vht_arr), np.asarray(xnp_arr)
    elif plot is True:
        # TODO: #8 create a useful data object that we can write avl data to, instead of using get_coef for certain use cases (also should be pickle-able)
        # pickle_path='/Avl/runs/tail_sizing/data/ht_sizing_lt_vht_analysis.p'
        # if os.path.exists(pickle_path):
        #     # read from the pickle file
        # else:
        #     # create the pickle file
        #     with open(pickle_path, 'w') as pfile:
        #         pickle.dump()
        # plot everything
        x=np.asarray(vht_arr)
        y=np.asarray(lt2plot)
        z=np.asarray(xnp_arr)
        plt_3d(x,y,z,
            'horizontal tail analysis',
            'vht', 'tail arm (m)', 'aft le Xnp (m)')

def xnp_from_tail_volume(vht, arht, tail_arm, wing, tail, mass):
    """ find the neutral point for given tail volume """
    sht=vht*wing['chord']*wing['span']/tail_arm
    bht=np.sqrt(arht*sht)
    cht = bht / arht
    tail['h_span']=bht
    tail['h_chord']=cht
    tail['l_tail']=tail_arm
    # create a plane from dimensions and create avl file
    plane=AvlAircraft(wing, tail, mass)
    avl=Avl(f'gen_vht_{vht:.4f}_lt_{tail_arm}',
            plane,
            'Avl',
            'runs/tail_sizing/runs')
    avl.create_avl_file()
    # run the steps from the template associated
    # TODO: change from tail_sizing to horizontal_tail_sizing
    # remove file that may be there already
    file_out_og=f'Avl/runs/tail_sizing/stab/{avl.name}_st.txt'
    if os.path.exists(file_out_og):
        # remove file
        os.remove(file_out_og)
        print(f'>>> AVL removing file: {file_out_og} \n')
    
    avl.run_avl('Avl/automation/steps_template_tail_sizing.txt')
    # rename the generically named file from the template to something
    # a little more descriptive of our current use cases
    file_out_name=f'Avl/runs/tail_sizing/stab/{avl.name}_{vht:.4f}_st.txt'
    os.rename(
        file_out_og,
        file_out_name)
    # open file, and get its TextIO object to get coefficient from
    with open(file_out_name, 'r') as stab_file:
        file_out = stab_file
        xnp=get_coef(file_out, 'Xnp')

    print(f'>> AVL tail volume: {vht} \n     xnp: {xnp} \n')
    return xnp

def xnp_target_from_tail_volume(vht, arht, banner_loc, tail_arm, wing):
    """ find the xnp target from a tail volume """
    mac = wing['chord']
    area_ref = wing['chord']*wing['span']
    l_banner = 9.5/5
    sm_banner = 0.1
    cht=np.sqrt(vht*mac*area_ref/tail_arm/arht)
    x_banner = tail_arm + cht - l_banner + l_banner/2 - banner_loc
    xnp_target=sm_banner*mac + x_banner
    print(f'tail volume: {vht} \nxnp_target: {xnp_target} \ncht: {cht}\n')
    return xnp_target

def tail_volume_from_banner_loc(a, b, arht, banner_loc, tail_arm, wing, tail, mass,
                                tol=1e-4, max_iter=50):
    """ find the horizontal tail volume that meets the xnp target value
    a: minimum tail volume
    b: maximum tail volume
    banner_loc: how far is the aft end of the banner from the end of the tail? (m)
    """

    fa = xnp_from_tail_volume(
        a, arht, tail_arm, wing, tail, mass 
    ) - xnp_target_from_tail_volume(a, arht, banner_loc, tail_arm, wing)
    fb = xnp_from_tail_volume(
        b, arht, tail_arm, wing, tail, mass
    ) - xnp_target_from_tail_volume(b, arht, banner_loc, tail_arm, wing)

    if fa * fb > 0:
        warnings.warn('ROOT NOT BRACKETED', UserWarning, stacklevel=2)

    for i in range(max_iter):
        c = 0.5 * (a + b)
        fc = xnp_from_tail_volume(
            c, arht, tail_arm, wing, tail, mass
        ) - xnp_target_from_tail_volume(c, arht, banner_loc, tail_arm, wing)

        if abs(fc) < tol or abs(b - a) < tol:
            return c

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

    raise RuntimeError("Bisection did not converge")

def tail_volume_from_fuse(x_nose, tail_arm,
                          wing, tail, mass, arht=3.32,
                          vht_min=0.1, vht_max=1.0,sm_banner=0.05,
                          max_iter=50, tol=1e-3):
    """ find the tail volume from the fuselage length
    Returns:
        vht: converged tail volume
        x_np: neutral point, m
        x_banner_excess: how far back banner is aft of the plane
    """
    # the idea is to basically place the LE of the banner 0.05 m from the
    # nose point
    x_banner = 1.9/2 + x_nose + 0.05
    xnp_target = sm_banner*wing['chord'] + x_banner

    # try to find the tail volume that can meet this value
    # use the bisection method
    # try with the minimum value
    a=vht_min
    b=vht_max
    xnp_a=xnp_from_tail_volume(a, arht, tail_arm, wing, tail, mass) - xnp_target
    xnp_b=xnp_from_tail_volume(b, arht, tail_arm, wing, tail, mass) - xnp_target
    
    print(f'>>>> BISECTION: x_nose: {x_nose}, tailarm: {tail_arm} \n'
          f'     xnp_a: {xnp_a+xnp_target}, xnp_b: {xnp_b+xnp_target}, xnp_target: {xnp_target} \n')

    if xnp_a * xnp_b > 0:
        print('>>> BISECTION: no root found \n')
        # warnings.warn('ROOT NOT BRACKETED', UserWarning, stacklevel=2)
        return None, None, None

    for _ in range(max_iter):
        c=0.5 * (a+b)
        x_np_c=xnp_from_tail_volume(c, arht, tail_arm, wing, tail, mass) - xnp_target
        if abs(x_np_c) < tol or abs(b-a) < tol:
            # calculate distance between aft point of banner and the aft point of the tail
            sht=c*wing['chord']*wing['span']/tail_arm
            bht=np.sqrt(arht*sht)
            cht = bht / arht
            x_banner_excess = (tail_arm+cht) - (x_banner + 1.9/2)
            # print to terminal
            print(f'>>> BISECTION: root found... \ntail volume: {c}\nneutral_point: {x_np_c+xnp_target} \n'
                  f'aft_excess: {x_banner_excess} \n'
                  f'x_nose: {x_nose} \ntail_arm: {tail_arm} \nstatic margin: {sm_banner} \n')
            return c, (x_np_c+xnp_target), x_banner_excess
        if xnp_a * x_np_c < 0:
            b=c
            xnp_b = x_np_c
        else:
            a=c
            xnp_a=x_np_c
    warnings.warn('Solution not converged', RuntimeError, stacklevel=2)
    return None, None, None

def tail_volume_from_cog_target(cog_target, sm_target, wing, tail, mass,
                                arht=1.44608, vht_min=0.1, vht_max=1.0,
                                tail_arm=1.5, max_iter=50, tol=1e-3):
    """ tail volume from the cog target 
    Returns:
        vht: converged tail volume, unitless
        xnp: converged neutral point, m
    """
    xnp_target=sm_target*wing['chord'] + cog_target

    a=vht_min
    b=vht_max

    xnp_a=xnp_from_tail_volume(a, arht, tail_arm, wing, tail, mass) - xnp_target
    xnp_b=xnp_from_tail_volume(b, arht, tail_arm, wing, tail, mass) - xnp_target
    
    if xnp_a * xnp_b > 0:
        print('>>> BISECTION: no root found \n')
        # warnings.warn('ROOT NOT BRACKETED', UserWarning, stacklevel=2)
        return None, None, None

    for _ in range(max_iter):
        c=0.5 * (a+b)
        x_np_c=xnp_from_tail_volume(c, arht, tail_arm, wing, tail, mass) - xnp_target
        if abs(x_np_c) < tol or abs(b-a) < tol:
            # calculate distance between aft point of banner and the aft point of the tail
            sht=c*wing['chord']*wing['span']/tail_arm
            bht=np.sqrt(arht*sht)
            cht = bht / arht
            # print to terminal
            print(f'>>> BISECTION: root found... \ntail volume: {c}\nneutral_point: {x_np_c+xnp_target} \n'
                  f'target cog: {cog_target} \ntarget sm: {sm_target}')
            return c, (x_np_c+xnp_target)
        if xnp_a * x_np_c < 0:
            b=c
            xnp_b = x_np_c
        else:
            a=c
            xnp_a=x_np_c
    warnings.warn('Solution not converged', RuntimeError, stacklevel=2)
    return None, None
    


wing={
    'span': 1.524,
    'chord': 0.416,
    'airfoil': 'clark_y.dat',
    'aileron_span': 0.381,   # span of EACH aileron
    'aileron_chord': 0.104,
    'incidence': 2.0
}
tail={
    'h_span': None,                # overwriting
    'h_chord': None,               # overwriting
    'l_tail': None,                # overwriting
    'h_foil': 'naca_0012.dat',
    'h_vertical_offset': 0.1127,
    'elev_chord': 0.08,
    'v_span': 0.2254,
    'v_chord': 0.1543,
    'v_foil': 'naca_0012.dat',
    'v_vertical_offset': None,  # if the tail is offset (skipping for now)
    'rudd_chord': 0.08,
    'incidence': -2.0 # i think idk
}
mass={  # cog of our plane, well have multiple cases, so we can make case yamls
    'x_ref': 0.1543,
    'y_ref': 0.0,
    'z_ref': 0.0,
}

    # ==================== #
    # |     ANALYSIS     | #
    # ==================== #
# when we start focusing on stability, we'll need to make use of mass

"""
# graveyard of failed attempts here...
### original analysis 
### ==================
# ht_sizing_vht_lt_analysis(0.1, 1.0, 12, 1.44608, 1.5, wing, tail, mass, plot=True)


### tail arm + tail volume analysis
### =================================
# ht_sizing_lt_vht_analysis(0.8, 2.0, 6, 0.1, 1.0, 6, 1.44608, wing, tail, mass, plot=True)


### tail arm + tail volume + tail incidence analysis
### =================================================
# dainc_arr=np.linspace(5, -5, 11)
# dainc2plot=[]
# lt2plot=[]  # different than the lt2plot given in the ht_sizing function
# vht2plot=[]
# xnp2plot=[]
# for dainc in dainc_arr:
#     tail['incidence']=dainc
#     lt2plot_old, vhtarr, xnp_arr = ht_sizing_lt_vht_analysis(0.8,2.0,6,0.1,1.0,6,1.4460,wing,tail,mass)
#     for lt, vht, xnp in zip(lt2plot_old, vhtarr, xnp_arr):
#         lt2plot.append(lt)
#         vht2plot.append(vht)
#         xnp2plot.append(xnp)
#         dainc2plot.append(dainc)    # here to ensure homogeneous lengths of arrays
    
# x=np.asarray(vht2plot)
# y=np.asarray(lt2plot)
# z=np.asarray(xnp2plot)
# c=np.asarray(dainc2plot)
# plt_3d(x,y,z,
#        'horizontal tail analysis',
#        'vht', 'tail arm (m)', 'aft_le Xnp (m)',
#        c, 'incidence (deg)')


# the incidence is kind of pointless tbh
# we can assume that the boluem affects the tail mass linearly
# we can try looking into trim cases now I think
# a combination of tail volume and control surface area


# COG sizing
# the idea here is that we want to size the tail and tail arm such that the banner cog
# is at about 10% static margin (let that be our minimum tail volume)
# first well look at the problem with a fixed tail arm and change the htail volume around
tail_arm = 1.5
arht=1.44608
# vht=tail_volume_from_banner_loc(0.1, 1.0, arht, 0.5,
#                                 tail_arm, wing, tail, mass)

# print(vht)
"""


# proof of concept here
# tail_volume_from_fuse(-0.70, 1.5, wing, tail, mass, sm_banner=0.05)

# move around the x nose and tail arm

# # n_iter=5
# lt_arr=np.linspace(0.83, 1.5, 6)
# # lt_arr=np.ones(4,1)*1.5
# xnose_arr=np.linspace(-0.67, -0.7, 6)

# vht2plot=[]
# lt2plot=[]
# xnose2plot=[]
# xnp2plot=[]
# xnose2plot=[]
# xexc2plot=[]
# sht2plot=[]
# bht2plot=[]
# cht2plot=[]
# # just go through each of them
# for lt in lt_arr:
#     for xnose in xnose_arr:
#         print(f'> SOLVER \n tail arm: {lt}')
#         print(f' xnose: {xnose} \n')
#         vht, xnp, xexc = tail_volume_from_fuse(xnose, lt, wing, tail, mass, sm_banner=0.05)
#         if(vht is None or
#            xnp is None or
#            xexc is None):
#             continue
#         else:
#             sht2plot.append(vht*wing['chord']*wing['span']/lt)
#             lt2plot.append(lt)
#             vht2plot.append(vht)
#             xnose2plot.append(xnose)
#             xnp2plot.append(xnp)
#             xexc2plot.append(xexc)

# plt_3d(lt2plot, vht2plot, xnose2plot,
#        'tail volumes from tail arm and xnose',
#        'tail arm, m', 'vht', 'xnose')

# plt_3d(lt2plot, vht2plot, xnose2plot,
#        'test_plot', 'tail arm, m', 'vht', 'xnose',
#        c=xnp2plot, clabel='neutral point location')

# plt_3d(lt2plot, vht2plot, xnose2plot,
#        'tail volume from tail arn and xnose, with excess length',
#        'tail arm, m', 'vht', 'xnose',
#        c=xexc2plot, clabel='excess length')

# plt_3d(lt2plot, sht2plot, xnose2plot,
#        'tail area from tail arm and xnose, with excess length',
#        'tail arm, m', 'tail area, m^2', 'xnose',
#        c=xexc2plot, clabel='excess lengths')


# plt.show()


# lets just pick and choose a tail volume and a tail arm that makes sense
# We'll stick with the 5% static
# it looks like shortening our tail arm may be best
# we can then take one and size the elevators manually I think
# automate tail incidence by solving for a range of CMA, with our cg envelope being say 5-15% static margin
# So maybe at 10% static margin we get this

# for control surface I need to sit down and think about what are the condition I need to meet
# Like what are my trimming conditions, and mess around with my area, and see what my deflections and moments look like

# Another analysis
# ==================

vht = 0.3953
tail_arm = 1.5
xnose = -0.7
sm=0.1
arht=3.32
sht=vht*wing['chord']*wing['span']/tail_arm
bht=np.sqrt(arht*sht)
cht = bht / arht
tail['h_span']=bht
tail['h_chord']=cht
tail['elev_chord']=0.25*cht

# get the vertical tail dimensions
cvt=cht
# AR = b/c ---> b=AR*c
arvt=1.46
bvt=arvt*cvt
tail['v_chord'] = cvt
tail['v_span'] = bvt
tail['h_vertical_offset']

tail['l_tail']=tail_arm
tail['h_vertical_offset'] = bvt / 2
tail['incidence']=-2.0

# create a plane from dimensions and create avl file
plane=AvlAircraft(wing, tail, mass)
avl=Avl(f'gen_vht_{vht:.4f}_lt_{tail_arm}',
        plane,
        'Avl',
        'runs/tail_sizing/runs')
avl.create_avl_file()
print(avl.name)
print(f' \nruns/tail_sizing/runs/{avl.name} \n')
print(f'dimensions \ncht: {cht} \nbht: {bht}, sht: {sht}')

# cog_target=0.2525
sm_target=0.05
vht,xnp = tail_volume_from_cog_target(cog_target, sm_target, wing, tail, mass, arht=3.32)

print(vht, xnp)
