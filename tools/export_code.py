'''Export the code cells of the executed notebooks to report/code/*.py for the LaTeX appendix.
   A cell identical to a cell already printed in an earlier program is replaced by a reference.'''
import json, os
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent          # repository root
CODE = str(ROOT / 'code') + '/'
OUT = str(ROOT / 'report' / 'code') + '/'
os.makedirs(OUT, exist_ok=True)
ORDER = [('A', 'P2_Task2.1'), ('B', 'P2_Task2.2'), ('C', 'P2_Task2.3'), ('D', 'P2_ET1'), ('E', 'P2_ET2'), ('F', 'P2_Task2.4-2.5')]
seen = {}
nonascii = set()
for letter, name in ORDER:
    nb = json.load(open(CODE + name + '.ipynb'))
    parts = ['# ' + '=' * 78, '# Program %s : %s.ipynb' % (letter, name), '# ' + '=' * 78]
    k = 0
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        k += 1
        src = ''.join(c['source']).rstrip()
        parts.append('')
        parts.append('# ' + '-' * 30 + ' Cell %d ' % k + '-' * 30)
        if src in seen:
            parts.append('# identical to Program %s, Cell %d (listed above)' % seen[src])
        else:
            seen[src] = (letter, k)
            parts.append(src)
            nonascii |= {ch for ch in src if ord(ch) > 127}
    open(OUT + name.replace('.', '_').replace('-', '_') + '.py', 'w').write('\n'.join(parts) + '\n')
    print(name, '->', OUT + name.replace('.', '_').replace('-', '_') + '.py', sum(1 for _ in open(OUT + name.replace('.', '_').replace('-', '_') + '.py')), 'lines')
print('non-ASCII characters in code:', nonascii)
