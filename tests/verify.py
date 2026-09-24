#!/usr/bin/env python3
"""cosmos-cinema dogrulama kapsulu: Kepler, enerji, momentum, determinizm."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import cosmos

OK = 0


def kontrol(ad, kosul, detay=""):
    global OK
    if kosul:
        OK += 1
        print(f"  [ok] {ad}")
    else:
        print(f"  [FAIL] {ad} {detay}")
        sys.exit(1)


print("1) Kepler: iki cisim dairesel yorunge")
pos = np.array([[0.0, 0, 0], [2.0, 0, 0]])
vel = np.array([[0.0, 0.5, 0], [0.0, -0.5, 0]])  # dairesel: v_rel = sqrt(GM/d) = 1
mass = np.array([1.0, 1.0])
k0, p0 = cosmos.enerji(pos, vel, mass, 0.01)
r_maks = 0.0
for _ in range(4000):
    cosmos.adimat(pos, vel, mass, 0.001, 0.01)
    r_maks = max(r_maks, float(np.linalg.norm(pos[0] - pos[1])))
k1, p1 = cosmos.enerji(pos, vel, mass, 0.01)
sapma = abs((k1 + p1) - (k0 + p0)) / abs(k0 + p0)
kontrol("yorunge bagli (r < 1.3x r0)", r_maks < 2.6, f"r_maks={r_maks:.2f}")
kontrol("enerji sapmasi < %0.5", sapma < 0.005, f"%{100 * sapma:.4f}")

print("2) momentum korunumu (mutlak)")
print(f"  p = {(mass @ vel).round(9).tolist()}")
kontrol("momentum == [0,0,0]", bool(np.allclose(mass @ vel, 0, atol=1e-9)))

print("3) kuvvet simetrisi (Newton 3. yasasi)")
pp = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]])
mm = np.array([2.0, 3.0])
a = cosmos.ivme(pp, mm, 0.01)
# a_i * m_i = -a_j * m_j olmali
sol = a[0] * mm[0]
sag = a[1] * mm[1]
kontrol("a1*m1 == -a2*m2", bool(np.allclose(sol, -sag)), f"{sol} vs {-sag}")

print("4) determinizm")
cfg_test = {"galaxies": [{"center": [0, 0, 0], "velocity": [0.2, 0, 0], "count": 200,
                          "radius": 5, "mass": 1.0}],
            "clusters": [], "singles": [],
            "sim": {"dt": 0.01, "steps": 50, "softening": 0.15}}


def kos():
    rng_state = np.random.default_rng(2026)
    p, v, m = cosmos.sahne_uret(cfg_test)
    # sahne_uret rng'yi icte uretir — determinizmi seed uzerinden test etmek icin
    # ayni cfg ayni sonucu vermeli
    return p, v


# sahne_uret her cagride yeni rng uretir (seed cfg'den) — ayni cfg ayni sonuc:
def uret():
    import json
    return cosmos.sahne_uret(json.loads(json.dumps(cfg_test)))


p1, v1, _ = uret() if False else (None, None, None)
import json as _json
p1, v1, m1 = cosmos.sahne_uret(_json.loads(_json.dumps(cfg_test)))
p2, v2, m2 = cosmos.sahne_uret(_json.loads(_json.dumps(cfg_test)))
kontrol("ayni cfg -> ayni baslangic", bool(np.array_equal(p1, p2) and np.array_equal(v1, v2)))

print(f"\nTUM KAPILAR GECTI ({OK} kontrol)")
