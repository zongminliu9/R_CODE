"""PyMOL side of the lab-meeting movies and stills (run inside PyMOL: pymol -cq a33_pymol_scene.py -- spec.json).

spec: {"panels": [{"name", "pdb", "dcd", "style", "alpha", "partner", "iso", "units"?}], "mode": "frames"|"still"|"pse",
       "states": [...] (1-based, per panel list or one list), "out": path or dir, "width", "height", "view_sel", "buffer",
       "turn": [[axis, deg], ...], "smooth": window (frames mode only; 0 = raw)}
Styles: loop (32-45 highlighted, side chains 39/42/43/45), groove (TBCB groove + tail aromatic), pocket (residue 114
and its neighbours), unfold (alpha coloured by folding unit), overview (whole alpha with the sites labelled).
"""
import sys, json, os
from pymol import cmd

spec = json.load(open(sys.argv[1]))
HEX = {'tuba1a': '2a78d6', 'tuba8': 'eb6834', 'aqua': '1baf7a', 'yel': 'eda100', 'mag': 'e87ba4', 'grn': '008300', 'vio': '6250d6',
       'red': 'e34948', 'ink': '3a3a38', 'alpha': 'd9d8d2', 'partner': 'b8c4cf'}
for k, v in HEX.items():
    cmd.set_color(k, [int(v[i:i + 2], 16) / 255 for i in (0, 2, 4)])
cmd.set('bg_rgb', 'white'); cmd.set('ray_opaque_background', 1); cmd.set('antialias', 2); cmd.set('ray_shadows', 0)
cmd.set('cartoon_fancy_helices', 1); cmd.set('cartoon_side_chain_helper', 1); cmd.set('stick_radius', 0.22)
cmd.set('specular', 0.2); cmd.set('ambient', 0.45); cmd.set('depth_cue', 1); cmd.set('fog_start', 0.55)
cmd.set('label_size', 12); cmd.set('label_position', (1.2, 1.8, 4.0)); cmd.set('label_color', 'ink'); cmd.set('label_outline_color', 'white'); cmd.set('label_font_id', 7)
cmd.set('dash_color', 'ink'); cmd.set('dash_width', 3); cmd.set('dash_gap', 0.3); cmd.feedback('disable', 'all', 'everything')

UNIT_COL = ['tuba1a', 'tuba8', 'aqua', 'yel', 'mag', 'grn']


def style(p):
    o, A, B, iso = p['name'], p['alpha'], p.get('partner'), p['iso']
    al = f'{o} and chain {A}'
    cmd.hide('everything', o)
    cmd.show('cartoon', f'{o} and polymer'); cmd.color('alpha', f'{o} and polymer')
    if B:
        cmd.color('partner', f'{o} and chain {B}')
    cmd.show('sticks', f'{o} and resn GTP'); cmd.color('grey50', f'{o} and resn GTP and elem C')
    cmd.show('spheres', f'{o} and resn MG'); cmd.set('sphere_scale', 0.35, f'{o} and resn MG'); cmd.color('grn', f'{o} and resn MG')
    s = p['style']
    if s in ('loop', 'overview'):
        cmd.color(iso, f'{al} and resi 32-45')
        cmd.set('cartoon_loop_radius', 0.45, f'{al} and resi 32-45')
        cmd.show('sticks', f'{al} and resi 39+42+43+45 and not name N+C+O'); cmd.set('stick_radius', 0.32, f'{al} and resi 39+42+43+45')
        cmd.color(iso, f'{al} and resi 39+42+43+45 and elem C'); cmd.util.cnc(f'{al} and resi 39+42+43+45')
    if s in ('pocket', 'overview'):
        cmd.show('sticks', f'{al} and resi 114 and not name N+C+O'); cmd.color(iso, f'{al} and resi 114 and elem C'); cmd.util.cnc(f'{al} and resi 114')
    if s == 'overview':
        cmd.color('aqua', f'{al} and resi 438-451'); cmd.set('cartoon_loop_radius', 0.45, f'{al} and resi 438-451')
        cmd.show('spheres', f'{al} and resi 114 and not name N+C+O+CA'); cmd.set('sphere_scale', 0.7, f'{al} and resi 114')
        cmd.set('label_size', 10)
    if s == 'pocket':
        cmd.set('cartoon_transparency', 0.55, o)
        cmd.select('pk', f'byres ({o} and polymer and not hydro) within 4.5 of ({al} and resi 114 and not name N+C+O+CA)')
        cmd.show('sticks', f'pk and not name N+C+O'); cmd.color('grey60', 'pk and elem C'); cmd.util.cnc('pk')
        cmd.color(iso, f'{al} and resi 114 and elem C'); cmd.show('sticks', f'{al} and resi 110 and name C+O+CA')
        if cmd.count_atoms(f'{al} and resi 114 and name OG'):
            cmd.distance(f'{o}_hb', f'{al} and resi 114 and name OG', f'{al} and resi 110 and name O', mode=0)
            cmd.hide('labels', f'{o}_hb')
    if s == 'groove':
        cmd.color(iso, f'{al} and resi 438-451')
        tail = f'{al} and (resi 449 and resn PHE or resi 451 and resn TYR)'
        cmd.show('sticks', f'{tail} and not name N+C+O'); cmd.color(iso, f'{tail} and elem C'); cmd.util.cnc(tail)
        cmd.set('cartoon_transparency', 0.35, f'{o} and chain {B}')
        gr = f'{o} and chain {B} and resi 187+188+192+203+223'
        cmd.show('sticks', f'{gr} and not name N+C+O'); cmd.color('grey40', f'{gr} and elem C'); cmd.util.cnc(gr)
    if s == 'unfold':
        for u, res in enumerate(p['units']):
            cmd.color(UNIT_COL[u % 6], f'{al} and resi {"+".join(str(r) for r in res)}')
        cmd.color('vio', f'{al} and resi 32-45'); cmd.color('red', f'{al} and resi 110-120')


def show_only(name):
    for q in cmd.get_names('all'):
        cmd.disable(q)
    cmd.enable(name)
    if f'{name}_hb' in cmd.get_names('all'):
        cmd.enable(f'{name}_hb')


panels = spec['panels']
for p in panels:
    cmd.load(p['pdb'], p['name'])
    if p.get('dcd'):
        cmd.load_traj(p['dcd'], p['name'], state=1)
    style(p)
def look(spec, obj):
    # camera on the 'from' side of 'target', looking at it; 'up' = [low_sel, high_sel] sets the vertical
    # PyMOL view[0:9] holds the model->camera rotation column-major: row i of R is (v[i], v[i+3], v[i+6])
    def c(sel):
        m = cmd.get_model(f'{obj} and ({sel})', state=1).atom
        assert m, sel
        return [sum(at.coord[k] for at in m) / len(m) for k in range(3)]
    sub = lambda p, q: [p[k] - q[k] for k in range(3)]
    dot = lambda p, q: sum(p[k] * q[k] for k in range(3))
    nrm = lambda p: [x / dot(p, p) ** 0.5 for x in p]
    cross = lambda p, q: [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]]
    tgt, frm = spec['look']
    z = nrm(sub(c(frm), c(tgt)))
    up = sub(c(spec['up'][1]), c(spec['up'][0])) if spec.get('up') else ([0.0, 1.0, 0.0] if abs(z[1]) < 0.9 else [1.0, 0.0, 0.0])
    up = [up[k] - z[k] * dot(up, z) for k in range(3)]; y = nrm(up); x = cross(y, z)
    v = list(cmd.get_view())
    for i, row in enumerate((x, y, z)):
        v[i], v[i + 3], v[i + 6] = row
    cmd.set_view(v)


first = panels[0]


def make_view(obj, alpha):
    cmd.orient(f'{obj} and chain {alpha} and name CA')
    for ax_, deg in spec.get('turn', []):
        cmd.turn(ax_, deg)
    if spec.get('look'):
        look(spec, obj)
    if spec.get('view_sel'):
        cmd.zoom(f'{obj} and ({spec["view_sel"]})', spec.get('buffer', 8), state=1)
    if spec.get('move_z'):
        cmd.move('z', spec['move_z'])
    return cmd.get_view()


VIEW = make_view(first['name'], first['alpha'])
PVIEW = {p['name']: (make_view(p['name'], p['alpha']) if spec.get('per_panel_view') else VIEW) for p in panels}
W, H = spec.get('width', 1280), spec.get('height', 720)

if spec['mode'] == 'pse':
    for p in panels[1:]:
        cmd.disable(p['name']) if spec.get('only_first') else None
    cmd.set_view(VIEW); cmd.frame(1)
    cmd.save(spec['out'])
elif spec['mode'] == 'still':
    for i, p in enumerate(panels):
        show_only(p['name'])
        st = spec['states'][i] if isinstance(spec['states'][0], list) else spec['states']
        for k, s_ in enumerate(st):
            cmd.set_view(PVIEW[p['name']]); cmd.set('state', s_)
            if spec.get('labels') and k == 0:
                for sel, txt in spec['labels']:
                    cmd.label(f'{p["name"]} and {sel}', f'"{txt}"')
            cmd.png(os.path.join(spec['out'], f'{p["name"]}_s{s_:04d}.png'), width=W, height=H, ray=1, quiet=1)
else:                                                  # frames for a movie, one panel at a time
    cmd.set('antialias', spec.get('antialias', 1))
    for p in panels:
        if spec.get('smooth', 0):
            cmd.smooth(f'{p["name"]} and polymer', 1, spec['smooth'])
    for i, p in enumerate(panels):
        show_only(p['name'])
        n = cmd.count_states(p['name'])
        d = os.path.join(spec['out'], p['name']); os.makedirs(d, exist_ok=True)
        for s_ in range(1, n + 1):
            cmd.set_view(VIEW); cmd.set('state', s_)
            cmd.png(os.path.join(d, f'f{s_:04d}.png'), width=W, height=H, ray=1, quiet=1)
        print('rendered', p['name'], n, flush=True)
