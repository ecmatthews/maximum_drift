import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord

class binary():
    '''class that will handle the calculation of secondary mass and RV drift.'''

    def __init__(self, csv):
        '''requires a csv file from Gaia containing information about nearby stars'''

        ## separation in arcsec, plx in mas
        self.csv = csv
        
        df = pd.read_csv(self.csv)
        ## HACK: this is an imperfect, but usually effective, way of identifying host
        ## >> brightest star in the search!
        ix_host = np.nanargmin(np.array(df['phot_g_mean_mag']))
        c1 = SkyCoord(ra=df['ra'][ix_host]*u.degree, dec=df['dec'][ix_host]*u.degree, frame='icrs')
        self.g1 = df['phot_g_mean_mag'][ix_host]

        self.plx = df['parallax'][ix_host]
        
        ## lists of separation and gmag and 'type'
        ## b=bound, v=vc/unknown, n=not a companion based on color
        self.s_list = []
        self.g_list = []
        self.bprp_list = []
        self.t_list = []
        
        for ix, row in df.iterrows():
                        
            ## get those without parallax - these are *unknown* VCs
            if np.isnan(row['parallax']):
                c2 = SkyCoord(ra=row['ra']*u.degree, dec=row['dec']*u.degree, frame='icrs')
                self.s_list.append(c1.separation(c2).arcsecond)
                self.g_list.append(row['phot_g_mean_mag'])  
                self.bprp_list.append(row['bp_rp'])  
                self.t_list.append('v')
                
            ## get those with parallax matching foreground star (with huge error tolerance)
            ## there are *bound* companions!        
            ## HACK the tolerance is hardcoded here
            if 1.1*self.plx > row['parallax'] > 0.9*self.plx:
                if not row['parallax'] == self.plx:
                    c2 = SkyCoord(ra=row['ra']*u.degree, dec=row['dec']*u.degree, frame='icrs')
                    self.s_list.append(c1.separation(c2).arcsecond)
                    self.g_list.append(row['phot_g_mean_mag'])  
                    self.bprp_list.append(row['bp_rp'])  
                    self.t_list.append('b')    
                                    
    def calc_m2(self):
        '''function to calculate the secondary mass using inputs from __init__'''

        ## data-table
        table = './data/mamajek_pecaut.dat'
        df = pd.read_csv(table,skiprows=22,nrows=118,delim_whitespace=True).replace('....',np.nan).replace('...',np.nan)

        self.m2list = []
        self.mclist = []
        
        for ix, gmag in enumerate(self.g_list):
            ## abs_mag
            abs_mag = gmag - 2.5*np.log10(((1000/self.plx)/10)**2)
            
            ## faint objects abs_mag > 17.3
            ## these are in the degenerate BD region - set to largest BD mass
            if abs_mag > 17.3:
                m2 = 0.075
                
            ## otherwise we can convert apparent mag to a mass (**IN G BAND ONLY**)
            else:
                mags = np.array(df['M_G']).astype(float)
                f = interpolate.interp1d(mags,df['Msun'])
                print('abs_mag > ', abs_mag)
                m2 = f(abs_mag)            
                
            self.m2list.append(m2 * u.Msun)  
                    
            ## also look at colors - if it's too blue, it has to be a background
            mags = np.array(df['Bp-Rp']).astype(float)
            f = interpolate.interp1d(mags,df['Msun'])
            ## HACK hardcoded bprp error 
            ## super generous estimate - mass from color min and max
            e_flux = 0.4
            #mc = f((self.bprp_list[ix])*u.Msun)
            if self.bprp_list[ix] < -0.12+e_flux:
                mc_min = f((self.bprp_list[ix]+e_flux)*u.Msun)
                mc_max = 10
            elif self.bprp_list[ix] > 4.86-e_flux:
                mc_min = 0.01
                mc_max = f((self.bprp_list[ix]-e_flux)*u.Msun)
            else:
                mc_min = f((self.bprp_list[ix]+e_flux)*u.Msun)
                mc_max = f((self.bprp_list[ix]-e_flux)*u.Msun)
            self.mclist.append([mc_min,mc_max])
            if (m2 > mc_max) or (m2 < mc_min):
                self.t_list[ix] = 'n'
            
            
    def manual_mag2(self, mag2):
        # ELISABETH TODO - delete this?
        
        for ix, v in enumerate(self.g_list):
            if np.isnan(v):
                self.g_list[ix] = mag2
        
            
    def plot_drift(self, obs=-100):
        '''
        plot the maximum drift as a function of projected separation of the companion.

        inputs:
         - obs: the observed drift, as a vector [value, uncertainty].
                default value is -100 -> in this case, no drift is plotted

        for "unknown" companions, this method assumes the companion is at the host star distance.
        '''
                
        ## define y-axis unit
        msy = u.m / u.s / u.yr
        
        for ix, _ in enumerate(self.s_list):
        
            ## axes
            xmin, xmax = 10, 5000
            
            ## companion properties
            sep_au = self.s_list[ix]*1000/self.plx*u.au
            print('on-sky separation: {0:.2f}"'.format(self.s_list[ix]))
            print('projected separation: {0:.0f}'.format(sep_au))
            print('mass at host distance:', self.m2list[ix])
            print('color-mass range: {0:.2f}-{1:.2f}'.format(self.mclist[ix][0],self.mclist[ix][1]))
            
            ## color-code
            if self.t_list[ix] == 'b':
                col = 'limegreen'
                print('this companion is bound! (matching plx)')
            elif self.t_list[ix] == 'v':
                col = 'dodgerblue'
                print('companion status unknown')
            elif self.t_list[ix] == 'n':
                col = 'red'
                print('companion *not* bound based on BP RP')
            else:
                col = 'black'
                print('companion type missing!')
            
            ## calculate linear drift for sep > minimum separation
            ## this formula is from Montagnier+2008
            rho = np.linspace(sep_au.value,xmax,num=50000) * sep_au.unit
            acc = (const.G*self.m2list[ix]/rho**2 * np.cos(np.arcsin(sep_au/rho))).to(msy)
            
            ## plot
            plt.plot(rho,acc,color='k')
            plt.xscale('log')
            plt.xlabel('rho [AU]')
            plt.xlim(xmin,xmax)
            plt.ylabel('linear drift [m/s/yr]')

            ## add the observed drift
            if not obs == -100:
                print('observed drift: {}+/-{}'.format(obs[0],obs[1]))
                plt.fill_between([xmin,xmax],[obs[0]-obs[1],obs[0]-obs[1]],[obs[0]+obs[1],obs[0]+obs[1]],
                                 color=col,alpha=0.3)
                plt.plot([xmin,xmax],[obs[0],obs[0]],color=col)

            plt.show()
                     