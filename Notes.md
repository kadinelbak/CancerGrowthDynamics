Cut off after xyz day.
Make sure to keep n constant on non IC50 curves 
Optimize the treatment with calculate interaction term
Optimize the treatment without interaction term
Do simulation on actual data to match 



Untreated Models
u0: {20k: 67, 30k: 100}

Theta Logistic 
r, K, v
Logistic
r, K
Gompertz
r, K
Logistic_death
r, K, d

r: [0, 2]
K: [500, 10,000}]
v: [.5, 5]
d: [0, 1]

Treated Models

Data: {IC25: 1.47 uM, IC50: 1.00 uM, IC75: 0.62 uM} at 20k; add Day 0 = 67 cells if missing and cut off at Day 12. Baseline r, K, v pulled from untreated monoculture logisitic fits.

Hill_ramp
Emax, n, lambda, IC50 (bounds: [0, 1.5], [0.1, 5], [0.01, 10], [1e-3, 100])

Hill_ramp_tonset
Emax, n, lambda, IC50, t_onset (t_onset: [0, 5])

two_pop_static_res
f_sens0, E_sens, n, IC50, lambda (f_sens0: [0, 1]; E_sens: [0, 2])

PARPi model1
alpha, beta, d ([0, 1], [0, 0.5], [0, 4])

PARPi model2 (n damage compartments collapsed, n=3)
alpha, beta, d ([0, 1], [0, 0.5], [0, 4])

PARPi model3
alpha, beta, gamma, d ([0, 1], [0, 0.5], [0, 4], [0, 4])

PARPi model4
alpha, phi, d ([0, 1], [0, 2], [0, 4])

two_pop_ramp_tonset_static
f_s0, E_sens, n, IC50, lambda, t_onset, delta_s (t_onset: [0, 5]; delta_s: [0, 1])

two_pop_ramp_tonset_switch
f_s0, E_sens, E_res, n, IC50, lambda, t_onset, delta_s, delta_r, q_sr (E_res: [0, 1]; delta_s/delta_r/q_sr: [0, 1]; t_onset: [0, 5])

