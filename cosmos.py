#!/usr/bin/env python3
"""
cosmos-cinema — JSON tanımlı yerçekimi N-kütle evreni.

Doğrudan özetleme (direct-summation) kütleçekimi + leapfrog (kick-drift-kick)
entegratörü. Milyonlarca kütleçekim çifti, bloklu NumPy ile. Çıktı: hız-ağırlıklı
renk haritasıyla filmlenmiş GIF.

    python3 cosmos.py --scene scenes/collision.json
    python3 cosmos.py --scene scenes/colossus.json --preview
"""

import argparse
import json
import math
import os
import shutil
import struct
import subprocess
import zlib

import numpy as np

G = 1.0
BLOK = 2048


def png_yaz(path, img):
    h, w, _ = img.shape
    raw = bytearray()
    for y in range(h):
        raw.append(0)
        raw += img[y].tobytes()

    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr)
                + chunk(b"IDAT", zlib.compress(bytes(raw), 6)) + chunk(b"IEND", b""))


# ---------------------------------------------------------------- fizik

def ivme(pos, mass, soft):
    """Bloklu doğrudan özetleme: a_i = -G Σ m_j r_ij / (|r_ij|²+eps²)^1.5."""
    n = pos.shape[0]
    a = np.zeros((n, 3))
    for i0 in range(0, n, BLOK):
        i1 = min(i0 + BLOK, n)
        dx = pos[None, :, :] - pos[i0:i1, None, :]          # (b, n, 3)
        r2 = np.einsum("ijk,ijk->ij", dx, dx)
        np.fill_diagonal(r2[:, i0:i1], np.inf)              # kendini atla
        r2 += soft * soft
        inv_r3 = 1.0 / (r2 * np.sqrt(r2))
        a[i0:i1] = G * np.einsum("ij,ijk->ik", (mass * inv_r3), dx)
    return a


def enerji(pos, vel, mass, soft):
    kin = 0.5 * float(np.sum(mass * np.einsum("ij,ij->i", vel, vel)))
    pot = 0.0
    n = pos.shape[0]
    for i0 in range(0, n, BLOK):
        i1 = min(i0 + BLOK, n)
        dx = pos[None, :, :] - pos[i0:i1, None, :]
        r2 = np.einsum("ijk,ijk->ij", dx, dx)
        np.fill_diagonal(r2[:, i0:i1], np.inf)
        r = np.sqrt(r2 + soft * soft)
        pot += G * float(np.sum(mass[None, :] * mass[i0:i1, None] / r)) * 0.5
    return kin, pot


def adimat(pos, vel, mass, dt, soft):
    """Leapfrog KDK."""
    a = ivme(pos, mass, soft)
    vel += 0.5 * dt * a
    pos += dt * vel
    a = ivme(pos, mass, soft)
    vel += 0.5 * dt * a


# ---------------------------------------------------------------- sahne üretimi

def disk_galaksi(merkez, hiz, n, yarıcap, rng, kalınlık=0.1, core_kutle=3000.0):
    """Merkez kittle + donen yildiz diski: v = sqrt(G*M_core/r)."""
    r = yarıcap * np.sqrt(rng.random(n))
    r = np.clip(r, yarıcap * 0.15, yarıcap)
    theta = 2 * np.pi * rng.random(n)
    x = merkez[0] + r * np.cos(theta)
    y = merkez[1] + rng.normal(0, yarıcap * kalınlık, n)
    z = merkez[2] + r * np.sin(theta)
    v = np.sqrt(G * core_kutle / r)
    vx = hiz[0] - v * np.sin(theta)
    vy = hiz[1] + rng.normal(0, v * 0.04, n)
    vz = hiz[2] + v * np.cos(theta)
    pos = np.stack([x, y, z], axis=1)
    vel = np.stack([vx, vy, vz], axis=1)
    # merkez cekirdek (kara delik rolunde)
    pos_core = np.array([merkez])
    vel_core = np.array([hiz])
    return np.vstack([pos_core, pos]), np.vstack([vel_core, vel])


def plummer_kume(merkez, n, yarıcap, rng):
    """Plummer küresi örnekleme (küresel yıldız kümesi)."""
    u = rng.random(n)
    r = yarıcap / np.sqrt(u ** (-2.0 / 3.0) - 1.0)
    r = np.minimum(r, yarıcap * 8)
    cos_t = 2 * rng.random(n) - 1
    sin_t = np.sqrt(1 - cos_t ** 2)
    phi = 2 * np.pi * rng.random(n)
    x = merkez[0] + r * sin_t * np.cos(phi)
    y = merkez[1] + r * cos_t
    z = merkez[2] + r * sin_t * np.sin(phi)
    # kolay baslangic: kismen kirisik (cokus izlemek icin)
    vel = rng.normal(0, 0.12, (n, 3))
    pos = np.stack([x, y, z], axis=1)
    return pos, vel


def sahne_uret(cfg):
    rng = np.random.default_rng(cfg.get("seed", 1))
    pos_l, vel_l, mass_l = [], [], []
    for gal in cfg["galaxies"]:
        p, v = disk_galaksi(
            np.array(gal["center"], float),
            np.array(gal["velocity"], float),
            gal["count"], gal["radius"], rng,
            gal.get("thickness", 0.1), gal.get("core_mass", 3000.0))
        pos_l.append(p)
        vel_l.append(v)
        mass_l.append(np.concatenate([[gal.get("core_mass", 3000.0)],
                                      np.full(gal["count"], gal.get("mass", 1.0))]))
    for k in cfg.get("clusters", []):
        p, v = plummer_kume(np.array(k["center"], float), k["count"],
                            k["radius"], rng)
        pos_l.append(p)
        vel_l.append(v)
        mass_l.append(np.full(k["count"], k.get("mass", 1.0)))
    for y in cfg.get("singles", []):
        pos_l.append(np.array([y["pos"]], float))
        vel_l.append(np.array([y["vel"]], float))
        mass_l.append(np.array([y.get("mass", 400.0)]))
    return np.vstack(pos_l), np.vstack(vel_l), np.concatenate(mass_l)


# ---------------------------------------------------------------- render

def renk_haritasi(hiz_norm):
    """hiz 0..1 -> mavi/beyaz/sari/kirmizi termal palet."""
    h = np.clip(hiz_norm, 0, 1)
    r = np.clip(2 * h, 0, 1) * 255
    g_ = np.clip(1 - 2 * np.abs(h - 0.5), 0, 1) * 255
    b = np.clip(1 - 2 * h, 0, 1) * 255
    return np.stack([r, g_, b], axis=-1)


def kare_ciz(pos, vel, cfg):
    W, H = cfg["render"]["width"], cfg["render"]["height"]
    merkez = np.array(cfg["render"].get("center", [0, 0, 0]))
    genis = cfg["render"].get("span", 24)
    olcek = W / genis
    x = (pos[:, 0] - merkez[0]) * olcek + W / 2
    y = ((pos[:, 2] - merkez[2]) * olcek + H / 2)
    hiz = np.linalg.norm(vel, axis=1)
    hn = hiz / (cfg["render"].get("vmax", 3.0))

    kütük_yoğunluk = np.zeros((H, W))
    hız_toplam = np.zeros((H, W))
    xi = np.clip(x.astype(int), 0, W - 1)
    yi = np.clip(y.astype(int), 0, H - 1)
    np.add.at(kütük_yoğunluk, (yi, xi), 1.0)
    np.add.at(hız_toplam, (yi, xi), hn)
    dolu = kütük_yoğunluk > 0

    yog = np.log1p(kütük_yoğunluk * 3.0)
    yog /= max(np.percentile(yog, 99.7), 1e-9)
    yog = np.clip(yog, 0, 1)
    ort_hiz = np.where(dolu, hız_toplam / np.maximum(kütük_yoğunluk, 1), 0)

    parlak = yog ** 0.65
    renk = renk_haritasi(ort_hiz) * parlak[..., None]
    art = np.clip(renk, 0, 255).astype(np.uint8)
    return art


# ---------------------------------------------------------------- ana

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene", required=True)
    ap.add_argument("--preview", action="store_true")
    args = ap.parse_args()

    cfg = json.load(open(args.scene, encoding="utf-8"))
    pos, vel, mass = sahne_uret(cfg)
    n = pos.shape[0]
    dt = cfg["sim"]["dt"]
    steps = cfg["sim"]["steps"]
    soft = cfg["sim"].get("softening", 0.15)
    render_araligi = cfg["render"].get("every", max(1, steps // 240))
    if args.preview:
        steps = min(steps, 40)
        render_araligi = max(1, steps // 8)

    out = cfg["output"]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    tmp = out + ".frames"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)

    k0, p0 = enerji(pos, vel, mass, soft)
    e0 = k0 + p0
    p_vec0 = mass @ vel
    print(f"  parcacik: {n}, adim: {steps}, dt: {dt}")
    print(f"  E0 = {e0:.2f} (kin {k0:.2f} + pot {p0:.2f})")

    kare_i = 0
    kacis = cfg["render"].get("capture_at", [])
    for adim in range(steps):
        adimat(pos, vel, mass, dt, soft)
        if (adim + 1) % render_araligi == 0 or adim == steps - 1:
            img = kare_ciz(pos, vel, cfg)
            png_yaz(f"{tmp}/f{kare_i:05d}.png", img)
            kare_i += 1
        if adim % 50 == 0:
            k, p = enerji(pos, vel, mass, soft)
            print(f"  adim {adim}: E sapmasi %{100 * abs((k + p) - e0) / abs(e0):.3f}")

    shutil.copy(f"{tmp}/f{kare_i - 1:05d}.png", out.replace(".gif", ".png"))

    fps = cfg["render"].get("fps", 20)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", f"{tmp}/f%05d.png",
                    "-vf", "palettegen=max_colors=256:stats_mode=diff", f"{tmp}/pal.png"],
                   check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
                    "-i", f"{tmp}/f%05d.png", "-i", f"{tmp}/pal.png",
                    "-lavfi", "paletteuse=dither=bayer:bayer_scale=3:diff_mode=rectangle",
                    "-loop", "0", out], check=True)
    shutil.rmtree(tmp, ignore_errors=True)

    k, p = enerji(pos, vel, mass, soft)
    sapma = 100 * abs((k + p) - e0) / abs(e0)
    p_vec1 = mass @ vel
    print(f"  enerji sapmasi: %{sapma:.4f} | momentum suruklendi: "
          f"{np.abs(p_vec1 - p_vec0).max():.2e}")
    print(f"== BİTTİ -> {out}")


if __name__ == "__main__":
    main()
