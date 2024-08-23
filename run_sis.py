"""
Create figure comparing different SIS worlds
"""

import sciris as sc
import numpy as np
import starsim as ss
import matplotlib.pyplot as plt

t_inf = 40

class one_more(ss.Intervention):
    
    def __init__(self, ti=t_inf, n=2):
        super().__init__()
        self.ti = ti
        self.n = n
    
    def apply(self, sim):
        if sim.ti == self.ti:
            sis = sim.diseases.sis
            uninfected = sis.susceptible.uids[:self.n]
            sim.diseases.sis.set_prognoses(uninfected)
            
            
def window_corr(tvec, x, y, wind=40):
    tout = []
    cout = []
    for i in range(len(tvec)-wind):
        twind = tvec[i:i+wind]
        xwind = x[i:i+wind]
        ywind = y[i:i+wind]
        this_t = twind.max()
        this_c = np.corrcoef(xwind, ywind)[0,1]**2
        tout.append(this_t)
        cout.append(this_c)
    if t_inf <= wind:
        tout = list(range(t_inf)) + tout
        cout = [1]*t_inf + cout
    return tout, cout
    

n_agents = 50

pars = dict(
    n_agents = n_agents,
    diseases = 'sis',
    networks = 'static',
    n_years = 200,
    rand_seed = 8,
)

ss.options(_centralized=False)

s1 = ss.Sim(pars)
s2 = ss.Sim(pars, interventions=one_more())
s1, s2 = ss.parallel(s1, s2).sims

ss.options(_centralized=True)

s3 = ss.Sim(pars)
s4 = ss.Sim(pars, interventions=one_more())
s3, s4 = ss.parallel(s3, s4).sims

r1 = s1.results.sis.n_infected
r2 = s2.results.sis.n_infected
r3 = s3.results.sis.n_infected
r4 = s4.results.sis.n_infected
tvec = np.arange(len(r1))

t_crn, c_crn = window_corr(tvec, r1, r2)
t_cen, c_cen = window_corr(tvec, r3, r4)

sc.options(dpi=200)
plt.figure(figsize=(14,6))
kw = dict(alpha=0.7, lw=2)
lkw = dict(frameon=False, loc='upper right', ncols=2)
ytop = max(r1.max(), r3.max())*1.3
vkw = dict(linestyle=':', c='k', **kw)


ax1 = plt.subplot(3,1,1)
ax1.plot(tvec, r2, label='Added infection', c='cornflowerblue', **kw)
ax1.plot(tvec, r1, label='Baseline', c='darkblue', **kw)
ax1.set_title('Default (centralized)', fontweight='bold') # WARNING: these seem to be reversed!!
plt.ylim(top=ytop)
plt.ylabel('Number infected')

ax2 = plt.subplot(3,1,2)
ax2.plot(tvec, r4, label='Added infection', c='darkorange', **kw)
ax2.plot(tvec, r3, label='Baseline', c='darkred', **kw)
ax2.set_title('Common random numbers', fontweight='bold')
plt.ylim(top=ytop)
plt.ylabel('Number infected')

ax3 = plt.subplot(3,1,3)
ax3.plot(t_crn, c_crn, label='Centralized', c='darkblue')
ax3.plot(t_cen, c_cen, label='Common', c='darkorange') # CK: WARNING: also swapped!
ax3.set_title('Correlation', fontweight='bold')
plt.xlabel('Time (days)')
plt.ylabel('Correlation')

for ax in [ax1, ax2, ax3]:
    ax.legend(**lkw)
    ax.axvline(t_inf, **vkw)
    ax.set_xlim(left=tvec.min(), right=tvec.max())
    sc.boxoff(ax=ax)
    

sc.figlayout()
plt.show()
sc.savefig('time_series_v2.png')