from pathlib import Path
import html
import numpy as np

OUT = Path(__file__).resolve().parent / 'ergebnisse'

def vandermonde(x, y, raster):
    """Monombasis in aufsteigender Reihenfolge; Auswertung mit Horner."""
    V = np.vander(x, N=len(x), increasing=True)
    a = np.linalg.solve(V, y)
    p = np.zeros_like(raster)
    for coefficient in a[::-1]:
        p = p * raster + coefficient
    return p, a, np.linalg.cond(V, np.inf)


def baryzentrische_basis(x, raster):
    """Gewichte einmal berechnen, dann normierte Lagrange-Basis auswerten.

    Gemeinsame Skalierung der Gewichte verändert die Formel nicht.
    Der Sonderfall x == Stützstelle wird zusätzlich korrekt behandelt.
    """
    d = x[:, None] - x[None, :]
    np.fill_diagonal(d, 1.)
    w = 1. / np.prod(d, axis=1)
    w /= np.max(np.abs(w))
    delta = raster[:, None] - x[None, :]
    hit = delta == 0
    safe = np.where(hit, 1., delta)
    terms = w / safe
    L = terms / terms.sum(axis=1, keepdims=True)
    for row in np.flatnonzero(hit.any(axis=1)):
        L[row] = hit[row].astype(float)
    return L


def plot(name, title, x, curves, ylabel='Fehler / °C', log=True, right=None):
    """Abhängigkeitenarme SVG-Diagramme; optional zweite logarithmische Achse."""
    colors = ['#1261a0', '#d34a2c', '#25834b', '#8944a0']
    left, top, width, height = 92, 68, 720, 355
    def transform(v):
        return np.log10(np.maximum(v, 1e-18)) if log else v
    values = np.concatenate([transform(np.asarray(y)) for _, y in curves])
    lo, hi = np.nanmin(values), np.nanmax(values)
    if log:
        lo, hi = np.floor(lo), np.ceil(hi)
    pad = max((hi-lo)*.05, .01)
    lo, hi = lo-pad, hi+pad
    xmin, xmax = min(x), max(x)
    px = left + (np.asarray(x)-xmin)/(xmax-xmin)*width
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="960" height="560" viewBox="0 0 960 560"><rect width="960" height="560" fill="white"/><g font-family="Arial" font-size="13"><text x="40" y="30" font-size="19">{html.escape(title)}</text>']
    def py(y, low=lo, high=hi):
        return top+height-(np.asarray(y)-low)/(high-low)*height
    for tick in np.linspace(xmin, xmax, 9):
        xx=left+(tick-xmin)/(xmax-xmin)*width
        parts.append(f'<path d="M{xx},{top}v{height}" stroke="#eee"/><text x="{xx}" y="445" text-anchor="middle">{tick:.3g}</text>')
    ticks = np.arange(np.ceil(lo), hi, max(1, np.ceil((hi-lo)/7))) if log else np.linspace(lo, hi, 7)
    for tick in ticks:
        yy=py(tick)
        label=f'10^{int(tick)}' if log else f'{tick:.3g}'
        parts.append(f'<path d="M{left},{yy}h{width}" stroke="#ddd"/><text x="82" y="{yy+4}" text-anchor="end">{label}</text>')
    def line(label, y, color, index, low=lo, high=hi):
        points=' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in zip(px, py(transform(np.asarray(y)), low, high)))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="1.8"/><text x="95" y="{492+index*19}" fill="{color}">{html.escape(label)}</text>')
    for j,(label,y) in enumerate(curves):
        line(label,y,colors[j],j)
    if right:
        label,y=right
        ry=transform(np.asarray(y)); rlo=np.floor(min(ry)); rhi=np.ceil(max(ry))
        if rhi==rlo: rhi+=1
        line(label,y,'#8944a0',len(curves),rlo,rhi)
        for tick in np.linspace(rlo,rhi,6):
            parts.append(f'<text x="825" y="{py(tick,rlo,rhi)+4}" fill="#8944a0">10^{tick:.0f}</text>')
        parts.append('<text x="890" y="56" fill="#8944a0" text-anchor="middle">κ∞ (rechts)</text>')
    parts.append(f'<path d="M{left},{top}v{height}h{width}" fill="none" stroke="black"/><text x="450" y="469" text-anchor="middle">{"Polynomgrad n" if xmax==40 else "Thermospannung E / mV"}</text><text x="18" y="255" transform="rotate(-90 18 255)" text-anchor="middle">{html.escape(ylabel)}</text></g></svg>')
    (OUT/name).write_text(''.join(parts),encoding='utf-8')


