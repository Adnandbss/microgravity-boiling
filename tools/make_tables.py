'''Generate LaTeX tables and number macros for the report from code/results/*.json'''
import json, math, os
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent          # repository root
R = str(ROOT / 'code' / 'results') + '/'
T = str(ROOT / 'report' / 'tables') + '/'
os.makedirs(T, exist_ok=True)
macros = []

def mac(name, val):
    macros.append('\\newcommand{\\%s}{%s}' % (name, val))

def sci(x, d=3):
    if x == 0:
        return '0'
    e = int(math.floor(math.log10(abs(x))))
    if -2 <= e <= 3:
        return ('%.' + str(max(d - 1 - e, 0)) + 'f') % x if e >= 0 else ('%.' + str(d - 1 - e) + 'f') % x
    return '%.*f\\times10^{%d}' % (d - 1, x / 10**e, e)

def f(x, d=4):
    return ('%.' + str(d) + 'f') % x

def write(name, s):
    open(T + name + '.tex', 'w').write(s)

def tup(t):
    return '(' + ',\\,'.join(sci(v, 3) if abs(v) < 0.01 else ('%g' % v) for v in t) + ')'

# ------------------------------------------------------------------ Task 2.1
if os.path.exists(R + 'task21.json'):
    t21 = json.load(open(R + 'task21.json'))
    mac('TwoOneMlow', f(t21['0.098']['m'], 2)); mac('TwoOneClow', sci(t21['0.098']['C']))
    mac('TwoOneMhigh', f(t21['9.8']['m'], 2)); mac('TwoOneChigh', sci(t21['9.8']['C']))

# ------------------------------------------------------------------ Task 2.2
if os.path.exists(R + 'task22.json'):
    t = json.load(open(R + 'task22.json'))
    b, ref, base, st = t['best'], t['ref'], t['base'], t['stats']
    mac('TwoTwoNone', sci(b['n1'])); mac('TwoTwoNtwo', f(b['n2'])); mac('TwoTwoNthree', f(b['n3']))
    mac('TwoTwoF', f(b['F'], 5)); mac('TwoTwoRefF', f(ref['F'], 5))
    mac('TwoTwoRefNone', sci(ref['n1'])); mac('TwoTwoRefNtwo', f(ref['n2'])); mac('TwoTwoRefNthree', f(ref['n3']))
    mac('TwoTwoBaseF', f(base['F'], 5)); mac('TwoTwoBaseKconv', '%d' % base['k035']); mac('TwoTwoBaseKmin', '%d' % base['kmin'])
    mac('TwoTwoRmsAbs', f(st['rms_abs'], 1)); mac('TwoTwoRmsRel', f(100 * st['rms_rel'], 1)); mac('TwoTwoMeanAbsRel', f(100 * st['mean_abs_rel'], 1))
    mac('TwoTwoMaxRel', f(100 * st['max_abs_rel'], 0)); mac('TwoTwoWithin', '%d' % round(100 * st['within10'])); mac('TwoTwoMeanQ', f(st['mean_q'], 0))
    for gk, nm in [('0.098', 'Low'), ('9.8', 'Mid'), ('19.6', 'High')]:
        mac('TwoTwoBias' + nm, '%+.1f' % (100 * t['bias'][gk]['mean_dev'])); mac('TwoTwoRms' + nm, f(100 * t['bias'][gk]['rms_dev'], 1))
    mac('TwoTwoBestSet', b['set']); mac('TwoTwoBestSeed', '%d' % b['seed'])
    # constants table
    s = '\\begin{tabular}{@{}l c c c c@{}}\n\\toprule\n & $n_1$ & $n_2$ & $n_3$ & $\\Ferrmean$\\\\\n\\midrule\n'
    s += 'Skeleton guesses (Set 1, MUT $=0.085$): MINIMUM & $%s$ & %s & %s & %s\\\\\n' % (sci(base['n1']), f(base['n2']), f(base['n3']), f(base['F'], 5))
    s += '\\phantom{Skeleton guesses (Set 1, MUT $=0.085$):} ENDING & $%s$ & %s & %s & %s\\\\\n' % (sci(base['n1end']), f(base['n2end']), f(base['n3end']), f(base['Fend'], 5))
    s += '\\phantom{Skeleton guesses (Set 1, MUT $=0.085$):} TIME AVG & $%s$ & %s & %s & %s$^\\dagger$\\\\\n' % (sci(base['n1ta']), f(base['n2ta']), f(base['n3ta']), f(base['Fta'], 5))
    s += '\\midrule\n\\textbf{Best GA fit} (%s, two-stage, seed %d): MINIMUM & $\\mathbf{%s}$ & \\textbf{%s} & \\textbf{%s} & \\textbf{%s}\\\\\n' % (b['set'], b['seed'], sci(b['n1']), f(b['n2']), f(b['n3']), f(b['F'], 5))
    s += '\\midrule\nDeterministic reference minimum (check) & $%s$ & %s & %s & %s\\\\\n' % (sci(ref['n1']), f(ref['n2']), f(ref['n3']), f(ref['F'], 5))
    s += '\\bottomrule\n\\end{tabular}\n'
    write('t22_constants', s)
    # study table
    s = '\\begin{tabular}{@{}l l l l c c c c r c@{}}\n\\toprule\n'
    s += 'Set & Initial guesses $(n_1,n_2,n_3)$ & MUT & NGEN & best & median & worst & $<0.035$ & $k_{0.035}$ & $n_2$ (min--max)\\\\\n'
    s += '& & & & \\multicolumn{3}{c}{$\\Ferrmean$ min.\\ over the run} & runs & & \\\\\n\\midrule\n'
    rows = t['table']
    for i, r in enumerate(rows):
        first = i % 5 == 0
        mut = r['MUT'] if isinstance(r['MUT'], str) else '%.3g' % r['MUT']
        ng = r['NGEN'] if isinstance(r['NGEN'], str) else '%d' % r['NGEN']
        k = '%d' % r['k035med'] if r['k035med'] >= 0 else 'never'
        col = '' if r['nsucc'] == r['nruns'] else '\\color{bad}'
        s += '%s & %s & %s & %s & %s & %s & %s & %s%d/%d & %s & %s--%s\\\\\n' % (
            r['set'] if first else '', ('$' + tup(r['init']) + '$') if first else '', mut, ng.replace('+', '+'), f(r['Fbest'], 5), f(r['Fmed'], 5), f(r['Fworst'], 5),
            col, r['nsucc'], r['nruns'], k, f(r['n2min'], 2), f(r['n2max'], 2))
        if i % 5 == 4 and i < len(rows) - 1:
            s += '\\midrule\n'
    s += '\\bottomrule\n\\end{tabular}\n'
    write('t22_study', s)


# ------------------------------------------------------------------ Task 2.3
if os.path.exists(R + 'task23.json'):
    t = json.load(open(R + 'task23.json'))
    b, ref, base, st = t['best'], t['ref'], t['base'], t['stats']
    K = ['n1', 'n2', 'n3', 'n4', 'n5']
    for k, nm in zip(K, ['None', 'Ntwo', 'Nthree', 'Nfour', 'Nfive']):
        mac('TwoThree' + nm, sci(b[k]) if k == 'n1' else f(b[k]))
        mac('TwoThreeRef' + nm, sci(ref[k]) if k == 'n1' else f(ref[k]))
    mac('TwoThreeF', f(b['F'], 5)); mac('TwoThreeFval', f(b['Fval'], 5)); mac('TwoThreeRefF', f(ref['F'], 5))
    mac('TwoThreeFinit', f(t['F_init_vii'], 5)); mac('TwoThreeBaseF', f(base['F'], 5)); mac('TwoThreeBaseKmin', '%d' % base['kmin'])
    mac('TwoThreeBaseKconv', '%d' % base['k030'])
    mac('TwoThreeRmsTrain', f(100 * st['train']['rms_rel'], 1)); mac('TwoThreeRmsVal', f(100 * st['val']['rms_rel'], 1))
    mac('TwoThreeMadTrain', f(100 * st['train']['mean_abs_rel'], 1)); mac('TwoThreeMadVal', f(100 * st['val']['mean_abs_rel'], 1))
    mac('TwoThreeWithinTrain', '%d' % round(100 * st['train']['within10'])); mac('TwoThreeWithinVal', '%d' % round(100 * st['val']['within10']))
    mac('TwoThreeRmsAbsTrain', f(st['train']['rms_abs'], 1)); mac('TwoThreeRmsAbsVal', f(st['val']['rms_abs'], 1))
    mac('TwoThreeBestSet', t['bestseq']['set']); mac('TwoThreeBestSeed', '%d' % t['bestseq']['seed'])
    geq = b['n4'] * 9.8 * 1.79
    mac('TwoThreeGeq', f(geq, 1)); mac('TwoThreeGeqRatio', f(geq / 9.8, 1))
    s = '\\begin{tabular}{@{}l c c c c c c c@{}}\n\\toprule\n & $n_1$ & $n_2$ & $n_3$ & $n_4$ & $n_5$ & $\\Ferrmean$ & $\\Ferrmean$\\\\\n & & & & & & (77 training) & (15 validation)\\\\\n\\midrule\n'
    s += 'Initial guesses (vii) & $4.76\\times10^{-4}$ & 3.028 & 0.2249 & 1.054 & 0.217 & %s & --\\\\\n' % f(t['F_init_vii'], 5)
    s += 'Run with guesses (vii): MINIMUM & $%s$ & %s & %s & %s & %s & %s & --\\\\\n' % (sci(base['n1']), f(base['n2']), f(base['n3']), f(base['n4']), f(base['n5']), f(base['F'], 5))
    s += '\\phantom{Run with guesses (vii):} ENDING & $%s$ & %s & %s & %s & %s & %s & --\\\\\n' % (sci(base['n1end']), f(base['n2end']), f(base['n3end']), f(base['n4end']), f(base['n5end']), f(base['Fend'], 5))
    s += '\\phantom{Run with guesses (vii):} TIME AVG & $%s$ & %s & %s & %s & %s & %s$^\\dagger$ & --\\\\\n' % (sci(base['n1ta']), f(base['n2ta']), f(base['n3ta']), f(base['n4ta']), f(base['n5ta']), f(base['Fta'], 5))
    s += '\\midrule\n\\textbf{Best GA fit} (set %s, refined) & $\\mathbf{%s}$ & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s}\\\\\n' % (
        t['bestseq']['set'], sci(b['n1']), f(b['n2']), f(b['n3']), f(b['n4']), f(b['n5']), f(b['F'], 5), f(b['Fval'], 5))
    s += '\\midrule\nDeterministic reference minimum & $%s$ & %s & %s & %s & %s & %s & --\\\\\n' % (sci(ref['n1']), f(ref['n2']), f(ref['n3']), f(ref['n4']), f(ref['n5']), f(ref['F'], 5))
    s += '\\bottomrule\n\\end{tabular}\n'
    write('t23_constants', s)
    # study table
    s = '\\begin{tabular}{@{}l l c c c c c c c c c@{}}\n\\toprule\n'
    s += 'Set & Initial guesses $(n_1,\\dots,n_5)$ & MUT & MFRAC & best & median & worst & $<0.03$ & $k_{0.03}$ & frac.\\ & refined\\\\\n'
    s += '& & & & \\multicolumn{3}{c}{single run, $\\Ferrmean$} & runs & & nkeep$=0$ & best--worst\\\\\n\\midrule\n'
    rows = t['table']; prev = None
    for r in rows:
        first = r['set'] != prev
        if first and prev is not None:
            s += '\\midrule\n'
        ref_runs = [d['F'] for d in t['refine'] if d['set'] == r['set']]
        k = '%d' % r['k030med'] if r['k030med'] >= 0 else 'never'
        col = '' if r['nsucc'] == r['nruns'] else '\\color{bad}'
        s += '%s & %s & %.2f & %.1f & %s & %s & %s & %s%d/%d & %s & %.2f & %s\\\\\n' % (
            r['set'] if first else '', ('$' + tup(r['init']) + '$') if first else '', r['MUT'], r['MFRAC'], f(r['Fbest'], 5), f(r['Fmed'], 5), f(r['Fworst'], 5),
            col, r['nsucc'], r['nruns'], k, r['nkeep0'], ('%s--%s' % (f(min(ref_runs), 5), f(max(ref_runs), 5))) if first else '')
        prev = r['set']
    s += '\\bottomrule\n\\end{tabular}\n'
    write('t23_study', s)

    if 'bias' in t:
        t22b = json.load(open(R + 'task22.json'))['bias'] if os.path.exists(R + 'task22.json') else {}
        s = '\\begin{tabular}{@{}l c c c c@{}}\n\\toprule\nTest condition & points & \\multicolumn{2}{c}{Eq.~(1.5), five constants} & Eq.~(1.1), three constants\\\\\n'
        s += '\\cmidrule(lr){3-4}\\cmidrule(l){5-5}\n & & mean deviation & rms deviation & mean deviation\\\\\n\\midrule\n'
        for r in t['bias']:
            k22 = ('%g' % r['g']) if (r['gamma'] == 1.79 and r['P'] == 5.5) else None
            old = ('$%+.1f\\,\\%%$' % (100 * t22b[k22]['mean_dev'])) if k22 in t22b else '--'
            flag = abs(r['mean_dev']) > 0.10
            mean = ('$\\mathbf{%+.1f\\,\\%%}$' if flag else '$%+.1f\\,\\%%$') % (100 * r['mean_dev'])
            s += '$g=%g$~m/s$^2$, $\\gamma=%g$, $P=%g$~kPa & %d & %s & $%.1f\\,\\%%$ & %s\\\\\n' % (
                r['g'], r['gamma'], r['P'], r['n'], mean, 100 * r['rms_dev'], old)
        s += '\\bottomrule\n\\end{tabular}\n'
        write('t23_bias', s)
        bmap = {(r['g'], r['gamma'], r['P']): r for r in t['bias']}
        mac('TwoThreeBiasOneG', '%+.0f' % (100 * bmap[(9.8, 1.79, 5.5)]['mean_dev']))
        mac('TwoThreeBiasTwoG', '%+.0f' % (100 * bmap[(19.6, 1.79, 5.5)]['mean_dev']))
        others = [r['mean_dev'] for k, r in bmap.items() if k not in [(9.8, 1.79, 5.5), (19.6, 1.79, 5.5)]]
        mac('TwoThreeBiasOthersMin', '%+.0f' % (100 * min(others))); mac('TwoThreeBiasOthersMax', '%+.0f' % (100 * max(others)))

# ------------------------------------------------------------------ ET1
if os.path.exists(R + 'et1.json'):
    t = json.load(open(R + 'et1.json'))
    V = ['dT', 'g', 'gamma', 'P']
    s = '\\begin{tabular}{@{}l c c c c c c c c c@{}}\n\\toprule\n'
    s += 'State & $(T_w-T_{sat})_1$ & $g_1$ & $\\gamma_1$ & $P_1$ & $\\qpp_1$ & \\multicolumn{4}{c}{fractional change ratio $(\\Delta\\qpp/\\qpp_1)/(\\Delta x/x_1)$}\\\\\n'
    s += '\\cmidrule(l){7-10}\n & \\textdegree C & m/s$^2$ & & kPa & W/cm$^2$ & $x=T_w-T_{sat}$ & $x=g$ & $x=\\gamma$ & $x=P$\\\\\n\\midrule\n'
    for sn, o in t['sens'].items():
        x = [o[v]['x1'] for v in V]
        s += '%s & %g & %g & %g & %g & %.1f & %s & %s & %s & %s\\\\\n' % (sn.replace('gen', '$g_{en}$'), x[0], x[1], x[2], x[3], o['q1'],
             *['%.3f \\scriptsize{(%.3f)}' % (o[v]['S'], o['analytic'][v]) for v in V])
    s += '\\bottomrule\n\\end{tabular}\n'
    write('et1_sens', s)
    s = '\\begin{tabular}{@{}l c c c c c c c@{}}\n\\toprule\n'
    s += 'State & \\multicolumn{4}{c}{(i) one input at a time: $0.15/|S_x|$} & (ii) worst case & (iii) RSS & superheat\\\\\n'
    s += '\\cmidrule(lr){2-5}\n & $T_w-T_{sat}$ & $g$ & $\\gamma$ & $P$ & $0.15/\\sum|S_x|$ & $0.15/\\sqrt{\\sum S_x^2}$ & limit (ii)\\\\\n\\midrule\n'
    for sn, l in t['limits'].items():
        o = t['sens'][sn]
        sing = ['%.3f' % l['single'][v] if l['single'][v] < 10 else '$>10$' for v in V]
        s += '%s & %s & %s & %s & %s & %.4f & %.4f & $\\pm$%.2f\\,\\textdegree C\\\\\n' % (sn.replace('gen', '$g_{en}$'), *sing, l['worst'], l['rss'], l['worst'] * o['dT']['x1'])
    s += '\\bottomrule\n\\end{tabular}\n'
    write('et1_limits', s)
    s1 = t['sens']['S1 (earth gravity)']; l1 = t['limits']['S1 (earth gravity)']
    mac('EtOneSdT', f(s1['dT']['S'], 2)); mac('EtOneSg', f(s1['g']['S'], 3)); mac('EtOneSgam', f(s1['gamma']['S'], 3)); mac('EtOneSP', f(s1['P']['S'], 3))
    mac('EtOneWorst', f(100 * l1['worst'], 1)); mac('EtOneRss', f(100 * l1['rss'], 1)); mac('EtOneSingle', f(100 * l1['single']['dT'], 1))
    mac('EtOneWorstC', f(l1['worst'] * 40, 1)); mac('EtOneSingleC', f(l1['single']['dT'] * 40, 1)); mac('EtOneCheck', f(100 * l1['check_worst'], 1))
    mac('EtOneSgLow', f(t['sens']['S2 (microgravity)']['g']['S'], 4)); mac('EtOneSgamLow', f(t['sens']['S2 (microgravity)']['gamma']['S'], 3))

# ------------------------------------------------------------------ ET2
if os.path.exists(R + 'et2.json'):
    t = json.load(open(R + 'et2.json'))
    s = '\\begin{tabular}{@{}l c c c c c c@{}}\n\\toprule\n'
    s += ' & \\multicolumn{3}{c}{genetic algorithm} & \\multicolumn{3}{c}{exact minima (check)}\\\\\n\\cmidrule(lr){2-4}\\cmidrule(l){5-7}\n'
    s += ' & GA fit of (ET2.1) & from Task 2.3 GA $n_i$ & difference & min.\\ of (ET2.1) & from min.\\ of (1.5) & difference\\\\\n\\midrule\n'
    rel = {'n1': '$(n_1)^{-1/n_2}$', 'n2': '$1/n_2$', 'n3': '$-n_3/n_2$', 'n4': '$n_4$', 'n5': '$-n_5/n_2$'}
    for r in t['table']:
        s += '$\\hat n_%s$ = %s & %s & %s & $%+.1f\\,\\%%$ & %s & %s & $%+.1f\\,\\%%$\\\\\n' % (r['k'][1], rel[r['k']], ('$%.4f$' % r['ga']), ('$%.4f$' % r['imp']), r['diff'],
                                                                     '$%.4f$' % r['inv_ref'], '$%.4f$' % r['imp_ref'], r['diff_ref'])
    s += '\\midrule\n$\\Ferrmean$ of (ET2.1) & %s & %s & & %s & & \\\\\n' % (f(t['ga_F'], 5), f(t['F_implied'], 5), f(t['inv_ref_F'], 5))
    s += '\\bottomrule\n\\end{tabular}\n'
    write('et2_table', s)
    s = '\\begin{tabular}{@{}l c c c c@{}}\n\\toprule\n & measured data & (a) noise-free & (b) scatter on $\\qpp$ & (c) scatter on $T_w-T_{sat}$\\\\\n'
    s += ' & (exact minima) & synthetic, GA & synthetic, 30 realisations & synthetic, 30 realisations\\\\\n\\midrule\n'
    for r in t['noise_table']:
        nf = 100 * (t['noisefree'][r['k']] - t['implied'][r['k']]) / abs(t['implied'][r['k']])
        s += '$\\hat n_%s$ & $%+.1f\\,\\%%$ & $%+.2f\\,\\%%$ & $%+.1f\\pm%.1f\\,\\%%$ & $%+.1f\\pm%.1f\\,\\%%$\\\\\n' % (r['k'][1], r['data'], nf, r['b_mean'], r['b_std'], r['c_mean'], r['c_std'])
    s += '\\bottomrule\n\\end{tabular}\n'
    write('et2_noise', s)
    mac('EtTwoSigma', f(t['sigma'], 3)); mac('EtTwoGaF', f(t['ga_F'], 5)); mac('EtTwoImpF', f(t['F_implied'], 5)); mac('EtTwoNoiseFreeF', sci(t['noisefree']['F'], 2))
    for r in t['table']:
        mac('EtTwoDiff' + {'n1': 'One', 'n2': 'Two', 'n3': 'Three', 'n4': 'Four', 'n5': 'Five'}[r['k']], '%+.1f' % r['diff'])
        mac('EtTwoDiffRef' + {'n1': 'One', 'n2': 'Two', 'n3': 'Three', 'n4': 'Four', 'n5': 'Five'}[r['k']], '%+.1f' % r['diff_ref'])

# ------------------------------------------------------------------ Task 2.5
if os.path.exists(R + 'task25.json'):
    t = json.load(open(R + 'task25.json'))
    b = t['best']
    mac('TwoFiveNone', f(b['n1'], 3)); mac('TwoFiveNtwo', f(b['n2'], 3)); mac('TwoFiveNthree', f(b['n3'], 3))
    mac('TwoFiveFtrain', f(b['Ftrain'], 4)); mac('TwoFiveFval', f(b['Fval'], 4)); mac('TwoFiveFall', f(b['Fall'], 4))
    mac('TwoFiveUA', f(b['UAeff'], 3)); mac('TwoFiveHr', f(b['hr_mean'], 2)); mac('TwoFiveKmin', '%d' % b['kmin'])
    mac('TwoFiveBestSet', b['set']); mac('TwoFiveBestNgen', '%d' % b['NGEN']); mac('TwoFiveBestSeed', '%d' % b['seed'])
    if 'ref' in t:
        r = t['ref']; mac('TwoFiveRefF', f(r['F'], 5)); mac('TwoFiveRefFval', f(r['Fval'], 4))
        mac('TwoFiveRefNone', f(r['n1'], 2)); mac('TwoFiveRefNtwo', f(r['n2'], 3)); mac('TwoFiveRefNthree', f(r['n3'], 3))
        c = t['conv_range']
        mac('TwoFiveNconv', '%d' % c['n']); mac('TwoFiveConvNoneMin', f(c['n1'][0], 2)); mac('TwoFiveConvNoneMax', f(c['n1'][1], 2))
        mac('TwoFiveConvNtwoMin', f(c['n2'][0], 2)); mac('TwoFiveConvNtwoMax', f(c['n2'][1], 2))
        mac('TwoFiveConvNthreeMin', f(c['n3'][0], 2)); mac('TwoFiveConvNthreeMax', f(c['n3'][1], 2))
        mac('TwoFiveConvUAMin', f(c['UAeff'][0], 3)); mac('TwoFiveConvUAMax', f(c['UAeff'][1], 3))
        pf = t['profile']; mac('TwoFiveProfMaxF', f(max(p[3] for p in pf), 6)); mac('TwoFiveProfMinF', f(min(p[3] for p in pf), 6))
    s = '\\begin{tabular}{@{}l l c c c c c c c c c c c@{}}\n\\toprule\n'
    s += 'Set & Initial guesses & MUT & MFRAC & \\multicolumn{3}{c}{training MAE (3 seeds)} & valid. & $<0.045$ & $k_{0.045}$ & frac. & \\multicolumn{2}{c}{best run}\\\\\n'
    s += '\\cmidrule(lr){5-7}\\cmidrule(l){12-13}\n & $(n_1,n_2,n_3)$ & & & best & median & worst & MAE & runs & & nkeep$=0$ & $n_1$, $n_2$ & $n_3$\\\\\n\\midrule\n'
    prev = None
    for r in t['table']:
        first = r['set'] != prev
        if first and prev is not None:
            s += '\\midrule\n'
        k = '%d' % r['k045med'] if r['k045med'] >= 0 else 'never'
        col = '' if r['nsucc'] == r['nruns'] else '\\color{bad}'
        s += '%s & %s & %.2f & %.1f & %s & %s & %s & %s & %s%d/%d & %s & %.2f & %.2f, %.2f & %.2f\\\\\n' % (
            r['set'] if first else '', ('$' + tup(r['init']) + '$') if first else '', r['MUT'], r['MFRAC'], f(r['Fbest'], 4), f(r['Fmed'], 4), f(r['Fworst'], 4),
            f(r['Fval'], 4), col, r['nsucc'], r['nruns'], k, r['nkeep0'], r['n1'], r['n2'], r['n3'])
        prev = r['set']
    s += '\\bottomrule\n\\end{tabular}\n'
    write('t25_study', s)
    s = '\\begin{tabular}{@{}r c c c c@{}}\n\\toprule\nNGEN & \\multicolumn{2}{c}{Set 3 $(5,5,1)$, MFRAC $=0.5$} & \\multicolumn{2}{c}{Set 1 $(1,1,1)$, MFRAC $=1.0$}\\\\\n'
    s += '\\cmidrule(lr){2-3}\\cmidrule(l){4-5}\n & training MAE & validation MAE & training MAE & validation MAE\\\\\n\\midrule\n'
    for NG in [100, 500, 2000, 6000, 20000]:
        cells = []
        for nm in ['Set 3', 'Set 1']:
            runs = [d for d in t['ngen'] if d['set'] == nm and d['NGEN'] == NG]
            cells.append('%s--%s' % (f(min(d['F'] for d in runs), 4), f(max(d['F'] for d in runs), 4)))
            cells.append('%s--%s' % (f(min(d['Fval'] for d in runs), 4), f(max(d['Fval'] for d in runs), 4)))
        s += '%d & %s & %s & %s & %s\\\\\n' % (NG, *cells)
    s += '\\bottomrule\n\\end{tabular}\n'
    write('t25_ngen', s)

open(str(ROOT / 'report' / 'numbers.tex'), 'w').write('% generated from code/results/*.json\n' + '\n'.join(macros) + '\n')
print('%d macros, tables written to %s' % (len(macros), T))
