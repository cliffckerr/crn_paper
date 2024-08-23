"""
Create figure comparing different SIS worlds
"""

import sciris as sc
import numpy as np
import starsim as ss
import matplotlib.pyplot as plt

class one_more(ss.Intervention):
    
    def __init__(self, ti=30, n=2):
        super().__init__()
        self.ti = ti
        self.n = n
    
    def apply(self, sim):
        if sim.ti == self.ti:
            sis = sim.diseases.sis
            uninfected = sis.susceptible.uids[:self.n]
            sim.diseases.sis.set_prognoses(uninfected)
            
            
def window_corr(tvec, x, y, wind=50):
    tout = []
    cout = []
    for i in range(len(tvec)-wind):
        twind = tvec[i:i+wind]
        xwind = x[i:i+wind]
        ywind = y[i:i+wind]
        this_t = twind.mean()
        this_c = np.corrcoef(xwind, ywind)[0,1]**2
        tout.append(this_t)
        cout.append(this_c)
    return tout, cout
    

n_agents = 100

pars = dict(
    n_agents = n_agents,
    diseases = 'sis',
    networks = 'static',
    n_years = 500,
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
plt.figure(figsize=(12,8))
kw = dict(alpha=0.7, lw=2)

plt.subplot(3,1,1)
plt.plot(tvec, r1, label='Baseline', **kw)
plt.plot(tvec, r2, label='Added infection', **kw)
plt.title('Common random numbers')
plt.legend()

plt.subplot(3,1,2)
plt.plot(tvec, r3, label='Baseline', **kw)
plt.plot(tvec, r4, label='Added infection', **kw)
plt.title('Centralized random numbers')
plt.legend()

plt.subplot(3,1,3)
plt.plot(t_crn, c_crn, label='Common')
plt.plot(t_cen, c_cen, label='Centralized')
plt.legend()
plt.title('Correlation')

sc.figlayout()
plt.show()
sc.savefig('time_series.png')