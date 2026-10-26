# ==============================================================================
# Program D : P2_ET1.ipynb
# ==============================================================================

# ------------------------------ Cell 1 ------------------------------
'''>>>>> Extra Task ET1 - fractional sensitivity of the heat flux predicted by Eq. (1.5)
   q" = n1 (Tw-Tsat)^n2 (g + n4 gen gamma)^n3 P^n5
   with the best-fit constants of Task 2.3 (read from results/task23.json written by P2_Task2.3).'''
from __future__ import print_function
import math, json, numpy
%matplotlib inline
import matplotlib.pyplot as plt

gen_en = 9.8
best = json.load(open('results/task23.json'))['best']
n1, n2, n3, n4, n5 = best['n1'], best['n2'], best['n3'], best['n4'], best['n5']
print('Task 2.3 best-fit constants: n1 = %.5e  n2 = %.4f  n3 = %.4f  n4 = %.4f  n5 = %.4f' % (n1, n2, n3, n4, n5))

def qpp(dT, g, gam, P):
    return n1 * dT**n2 * (g + n4 * gen_en * gam)**n3 * P**n5

# typical sets of inputs [(Tw-Tsat)1, g1, gamma1, P1]
STATES = {'S1 (earth gravity)': (40.0, 9.8, 1.79, 5.5),
          'S2 (microgravity)':  (40.0, 0.098, 1.79, 5.5),
          'S3 (2 gen, 9.5 kPa)': (35.0, 19.6, 1.79, 9.5)}
DX_SUPERHEAT = 2.0      # deg C  (step 2 of the statement)
FRAC = 0.05             # 5 % step for g, gamma and P (same relative size as 2 C / 40 C)

et1 = {}
for sname, (dT1, g1, gam1, P1) in STATES.items():
    q1 = qpp(dT1, g1, gam1, P1)                                    # step 1
    out = {'q1': q1}
    for var, x1 in [('dT', dT1), ('g', g1), ('gamma', gam1), ('P', P1)]:
        dx = DX_SUPERHEAT if var == 'dT' else FRAC * x1
        args = dict(dT=dT1, g=g1, gam=gam1, P=P1)
        args['gam' if var == 'gamma' else var] = x1 + dx              # step 2: change only this input
        q2 = qpp(**args)
        S = (q2 - q1) * x1 / (dx * q1)                                 # step 3: fractional change ratio
        out[var] = dict(x1=x1, dx=dx, q2=q2, S=S)
    # analytical logarithmic derivatives (limit dx -> 0), for comparison
    s = g1 + n4 * gen_en * gam1
    out['analytic'] = dict(dT=n2, g=n3 * g1 / s, gamma=n3 * n4 * gen_en * gam1 / s, P=n5)
    et1[sname] = out
    print('\n%s: Tw-Tsat = %.1f C, g = %.3f m/s2, gamma = %.2f, P = %.1f kPa  ->  q"1 = %.2f W/cm2' % (sname, dT1, g1, gam1, P1, q1))
    for var in ['dT', 'g', 'gamma', 'P']:
        o = out[var]
        print('   %-6s x1 = %7.3f  dx = %6.3f  q"2 = %7.2f  dq"/q" / (dx/x) = %7.4f   (analytic %.4f)'
              % (var, o['x1'], o['dx'], o['q2'], o['S'], out['analytic'][var]))

# ------------------------------ Cell 2 ------------------------------
# identical to Program A, Cell 2 (listed above)

# ------------------------------ Cell 3 ------------------------------
'''ET1 - limiting fractional inaccuracy of the inputs for |dq"/q"| < 0.15
   To first order  dq"/q" = S_dT dT/dT + S_g dg/g + S_gamma dgamma/gamma + S_P dP/P, so:
   (i)   one input inaccurate at a time   : eps_x < 0.15 / |S_x|
   (ii)  all inputs with the same eps, errors adding in the worst direction : eps < 0.15 / sum|S_x|
   (iii) same eps, independent random errors (root-sum-square)            : eps < 0.15 / sqrt(sum S_x^2)'''
QTOL = 0.15
VARS = ['dT', 'g', 'gamma', 'P']
lim = {}
for sname, out in et1.items():
    S = numpy.array([abs(out[v]['S']) for v in VARS])
    single = QTOL / S
    worst = QTOL / S.sum()
    rss = QTOL / math.sqrt((S**2).sum())
    lim[sname] = dict(single=dict(zip(VARS, single.tolist())), worst=worst, rss=rss)
    # nonlinear check: all inputs off by the worst-case eps in the direction that increases q"
    dT1, g1, gam1, P1 = [out[v]['x1'] for v in VARS]
    sgn = numpy.sign([out[v]['S'] for v in VARS])
    qw = qpp(dT1 * (1 + sgn[0] * worst), g1 * (1 + sgn[1] * worst), gam1 * (1 + sgn[2] * worst), P1 * (1 + sgn[3] * worst))
    lim[sname]['check_worst'] = qw / out['q1'] - 1
    print('%s:' % sname)
    print('   (i)   single input : ' + ', '.join('%s < %.3f' % (v, e) for v, e in zip(VARS, single)) +
          '   -> superheat within +/- %.2f C' % (single[0] * dT1))
    print('   (ii)  worst case   : eps < %.4f  (superheat +/- %.2f C, g +/- %.3f m/s2, gamma +/- %.3f, P +/- %.3f kPa);'
          ' nonlinear check dq"/q" = %.4f' % (worst, worst * dT1, worst * g1, worst * gam1, worst * P1, lim[sname]['check_worst']))
    print('   (iii) RSS          : eps < %.4f  (superheat +/- %.2f C)' % (rss, rss * dT1))

fig, ax = plt.subplots(figsize=(6.6, 3.6))
w = 0.26
for c, (sname, out) in enumerate(et1.items()):
    ax.bar(numpy.arange(4) + (c - 1) * w, [out[v]['S'] for v in VARS], width=w - 0.03, color=C[c], label=sname)
ax.set_xticks(range(4))
ax.set_xticklabels([r'superheat $T_w-T_{sat}$', r'gravity $g$', r'surface tension $\gamma$', r'pressure $P$'])
ax.set_ylabel(r'$(\Delta q^{\prime\prime}/q^{\prime\prime})/(\Delta x/x)$')
ax.axhline(0, color='#52514e', lw=0.8)
ax.legend(loc='upper right')
ax.grid(axis='x', visible=False)
savefig('et1_sensitivity')
json.dump(dict(constants=best, states=STATES, sens=et1, limits=lim), open(RESDIR + 'et1.json', 'w'), indent=1)
