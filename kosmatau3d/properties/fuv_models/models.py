import numpy as np
import astropy.units as u
import astropy.constants as con
import matplotlib.pyplot as plt

from astropy.modeling.physical_models import BlackBody


def u_habing(l):
    '''Spectrum used in Habing (1968), as derived by Bertoldi and Draine (1996).'''
    return (-25/6*(l/1e3)**3 + 25/2*(l/1e3)**2 - 13/3*l/1e3)*1e-14

def u_draine(l, chi=1.71, test=False):
    '''Spectrum used in Draine (1978).'''
    if test:
        return 4e-14 * (53.05461*(l/1e3)**-3 - 85.37824*(l/1e3)**-4 + 34.03412*(l/1e3)**-5)
    else:
        return 4e-14 * chi * (31.02609*(l/1e3)**-3 - 49.92879*(l/1e3)**-4 + 19.90300*(l/1e3)**-5)

def u_mezger(l, floattype=False):
    '''Spectrum used in Mezger, Mathis, and Panagia (1982).'''
    if type(l) == int or type(l) == float:
        floattype = True
        l = np.asarray([l])
    ul = np.zeros(*l.shape)
    bb = BlackBody()
    i = (l>2460)
    ul[i] = l[i]*(bb.evaluate(l[i]*u.AA, scale=1e-18*u.erg/(u.s*u.cm**2*u.AA*u.sr), 
                              temperature=7500*u.K) 
                  + bb.evaluate(l[i]*u.AA, scale=1.65e-17*u.erg/(u.cm**2*u.s*u.AA*u.sr), 
                                temperature=4000*u.K)) * 4*np.pi / 2.9979e10
    i = (l>1340) & (l<=2460)
    ul[i] = 2.373e-14*4*np.pi * (l[i]/1e4)**-0.6678
    i = (l>1100) & (l<=1340)
    ul[i] = 6.821e-13*4*np.pi * (l[i]/1e4)
    i = (l>=912) & (l<=1100)
    ul[i] = 1.287e-9*4*np.pi * (l[i]/1e4)**4.4172
    if floattype:
        return ul[0]
    else:
        return ul

def u_zucconi(l, temp=22000, scale=5.7e-17):
    '''Spectrum used in Zucconi, Walmsley, and Galli (2003).'''
    ul = (2*con.h*con.c/l**3/u.AA**3).to(u.g/u.s**2).value * scale / (np.exp(con.h*con.c/l/u.AA/con.k_B/temp/u.K).value - 1)
    return ul

def u_kosma(l, floattype=False, **kwargs):
    ''' Spectrum used in KOSMA-tau (Röllig, Szczerba, Ossenkopf, and Glück  2013).'''
    if type(l) == int or type(l) == float:
        floattype = True
        l = np.asarray([l])
    ul = u_draine(l, **kwargs)
    i = l >= 2000
    ul[i] = 1.39148e-5/con.c.to(u.cm/u.s).value * l[i]**0.7
    if floattype:
        return ul[0]
    else:
        return ul

def compare_models(l_range=(912, 2066)):
    '''
    Print spectrum information to screen.
    '''
    l = np.linspace(l_range[0], l_range[1], num=1000)
    t = 29000
    w = 1.5e-11

    habing = u_habing(l)
    draine = u_draine(l)
    mezger = u_mezger(l)
    zucconi = u_zucconi(l, temp=t, scale=w)
    kosma = u_kosma(l)

    uint_habing = np.trapz(habing/l, l)
    uint_draine = np.trapz(draine/l, l)
    uint_mezger = np.trapz(mezger/l, l)
    uint_zucconi = np.trapz(zucconi/l, l)
    uint_kosma = np.trapz(kosma/l, l)

    print('integrated Habing spectrum: {:.3e} erg/cm^3\n'.format(uint_habing))

    # chi
    print('Chi')
    print('  Draine: {:.3f}'.format(u_draine(1000)/4e-14))
    print('  Mezger: {:.3f}'.format(u_mezger(1000)/4e-14))
    print('  Zucconi: {:.3f}'.format(u_zucconi(1000, temp=t, scale=w)/4e-14))
    print('  KOSMA-tau: {:.3f}'.format(u_kosma(1000)/4e-14))

    # integrated energy density
    print('Total energy density ratio')
    print('  Draine: {:.3f}'.format(uint_draine/uint_habing))
    print('  Mezger: {:.3f}'.format(uint_mezger/uint_habing))
    print('  Zucconi: {:.3f}'.format(uint_zucconi/uint_habing))
    print('  KOSMA-tau: {:.3f}'.format(uint_kosma/uint_habing))

    return


def plot_comparison(l_range=(912, 2066), num=10000, ticklabelsize=16, labelsize=24):
    '''
    Plot spectra together.
    '''

    l = np.linspace(l_range[0], l_range[1], num=num)
    t = 29000
    w = 1.5e-11

    habing = u_habing(l)
    draine = u_draine(l)
    mezger = u_mezger(l)
    zucconi = u_zucconi(l, temp=t, scale=w)
    kosma = u_kosma(l)

    fig, ax = plt.subplots(1, 1, figsize=(14, 10))

    ax.plot(l, habing/1e-13, lw=3, label='Habing (1968)')
    ax.plot(l, draine/1e-13, lw=3, label='Draine (1978)')
    ax.plot(l, mezger/1e-13, ls=':', lw=3, label='Mezger et al. (1982)')
    ax.plot(l, zucconi/1e-13, ls='-.', lw=3, label='Zucconi et al. (2003)')
    ax.plot(l, kosma/1e-13, ls='--', lw=3, label='Röllig et al. (2013)')

    ax.set_xlabel(r'$\lambda$ $(\AA)$', fontsize=labelsize)
    ax.set_ylabel(r'$\lambda u_\lambda$ $\left( 10^{-13}~\mathrm{erg}\,\mathrm{cm^3} \right)$', fontsize=labelsize)
    ax.tick_params(axis='both', labelsize=ticklabelsize)
    ax.ticklabel_format(axis='both', style='sci', useLocale=False, useOffset=False, useMathText=True)
    ax.legend(fontsize=labelsize)
    ax.set_xlim(l.min(), l.max())

    plt.show()
    plt.savefig('/home/yanitski/projects/pdr/FUV_comparison.png')

    return

