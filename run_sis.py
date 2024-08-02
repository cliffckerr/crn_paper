"""
Create figure comparing different SIS worlds
"""

# import sciris as sc
import starsim as ss
import matplotlib.pyplot as plt

class one_more(ss.Intervention):
    
    def __init__(self, ti=30):
        super().__init__()
        self.ti = ti
    
    def apply(self, sim):
        if sim.ti == self.ti:
            sis = sim.diseases.sis
            first_uninfected = ss.uids([sis.susceptible.uids[0]])
            sim.diseases.sis.set_prognoses(first_uninfected)

n_agents = 100

pars = dict(
    n_agents = n_agents,
    diseases = 'sis',
    networks = 'static',
    n_years = 500,
    rand_seed = 2,
)

s1 = ss.Sim(pars)
s2 = ss.Sim(pars, interventions = one_more())
s1, s2 = ss.parallel(s1, s2).sims

r1 = s1.results.sis.n_infected
r2 = s2.results.sis.n_infected
plt.plot(r1)
plt.plot(r2)
# sim.plot()
# sim.diseases.sis.plot()

plt.show()