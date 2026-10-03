# ============================================================
# analizar_erd_ers.py — Paso 3: ERD/ERS (versión para el curso)
# ============================================================
#
# Entrada : filtrado/*_filtrado.csv   (salida de filter_eeg.py)
# Salida  : resultados/...
#
# Calcula el mismo análisis del proyecto completo:
#   A. ERD/ERS en dominio tiempo (% de cambio de potencia vs baseline)
#   B. Mapas tiempo-frecuencia (wavelets de Morlet) + FDR
#   C. PSD de Welch (CROSS vs INTENTION)
#   D. Potencia absoluta (uV²)
#
# en TRES niveles de la jerarquía:
#
#   POR TRIAL   -> un trial suelto, sin promediar nada
#   POR CORRIDA -> promedio de los trials de UN archivo
#   GLOBAL      -> todas las corridas juntas (estadística con n = n_corridas)
#
# Cada carpeta contiene reporte_erd_ers.csv + 7 gráficas PNG:
#   erd_ers_timecourse.png, erd_ers_barras.png,
#   potencia_timecourse_completa.png, potencia_timecourse_bandas.png,
#   tfr_ersp.png, tfr_ersp_potencia_absoluta.png, psd_baseline_vs_tarea.png
#
# Resultados en:
#   resultados/global/
#   resultados/por_corrida/<archivo>/
#   resultados/por_trial/<archivo>/trial_<n>/
#
# Protocolo de cada trial (eventos en la columna event_id):
#   100 CROSS (0-3 s, baseline) | 200 INTENTION (3-6 s)
#   300 MOTION (6-9 s)          | 400 RESET | 999 (se descarta)
#
# Uso:   python analizar_erd_ers.py
#        python analizar_erd_ers.py --sin-trial   (omite el nivel por trial: más rápido)
# ============================================================
import argparse
import textwrap
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                      # sin ventana: solo guarda PNGs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.signal import butter, sosfiltfilt, hilbert, fftconvolve, welch
from scipy.stats import ttest_1samp

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "filtrado"
OUT_ROOT = ROOT / "resultados"
OUT_GLOBAL = OUT_ROOT / "global"
OUT_POR_CORRIDA = OUT_ROOT / "por_corrida"
OUT_POR_TRIAL = OUT_ROOT / "por_trial"

FS = 500.0
EVENTS = {"CROSS": 100, "INTENTION": 200, "MOTION": 300, "RESET": 400}
META_COLS = {"datetime", "timestamp", "event_id", "event_label", "trial"}

# Canales EEG: se detectan automáticamente al leer el primer CSV.
CHANNELS = []

# ── Duraciones nominales del protocolo (eje de tiempo real, segundos) ──
T_CROSS_END = 3.0       # fin de CROSS / inicio de INTENTION
T_INTENTION_END = 6.0   # fin de INTENTION / inicio de MOTION
T_MOTION_END = 9.0      # fin de MOTION

# ── Bandas ──────────────────────────────────────────────────────────────
BANDS = {                                   # timecourse de %ΔP
    "Mu (8-13 Hz)": (8.0, 13.0),
    "Beta (13-30 Hz)": (13.0, 30.0),
}
BAND_GAMMA = {"Gamma (30-45 Hz)": (30.0, 45.0)}
BANDS_BARRAS = {**BANDS, **BAND_GAMMA}      # reporte y barras (3 bandas)

BAND_ZONES = {                              # potencia multi-banda y franjas del PSD
    "Mu (8-13 Hz)":         (8.0, 13.0),
    "Beta Low (13-21 Hz)":  (13.0, 21.0),
    "Beta High (21-30 Hz)": (21.0, 30.0),
    "Gamma (30-45 Hz)":     (30.0, 45.0),
}
BAND_ZONE_COLORS = {
    "Mu (8-13 Hz)":         "#8ecae6",
    "Beta Low (13-21 Hz)":  "#ffd166",
    "Beta High (21-30 Hz)": "#e63946",
    "Gamma (30-45 Hz)":     "#2ca02c",
}
POWER_BAND_FULL = {"Banda completa (8-45 Hz)": (8.0, 45.0)}
POWER_BANDS_ALL = {**POWER_BAND_FULL, **BAND_ZONES}
ALL_BANDS = {**BANDS_BARRAS, **POWER_BANDS_ALL}   # unión: se filtra una sola vez

N_POINTS_TIMECOURSE = 200
TFR_FREQS = np.arange(8.0, 46.0, 1.0)
TFR_MIN_CYCLES = 3.0
TFR_CYCLES_PER_HZ = 0.5
N_POINTS_TFR = 60
FDR_ALPHA = 0.05
PSD_YLIM = (1e-1, 1e1)                      # eje Y (log, uV²/Hz) del PSD


# ── Utilidad de títulos ───────────────────────────────────────────────────

def set_wrapped_title(fig, title: str, fontsize: int = 12, width: int = 95,
                      extra_bottom_in: float = 0.0) -> float:
    """Título envuelto en varias líneas, reservando el espacio vertical en pulgadas reales."""
    wrapped = "\n".join(textwrap.wrap(title, width=width)) or title
    n_lines = wrapped.count("\n") + 1

    fig_h_in = fig.get_size_inches()[1]
    line_h_in = (fontsize / 72.0) * 1.35
    title_h_in = n_lines * line_h_in
    pad_in = 0.05

    title_top_frac = 1.0 - (pad_in * 0.5) / fig_h_in
    fig.suptitle(wrapped, fontsize=fontsize, y=title_top_frac, va="top")

    title_bottom_frac = title_top_frac - title_h_in / fig_h_in
    top = title_bottom_frac - (pad_in * 0.5 + extra_bottom_in) / fig_h_in
    top = max(top, 0.55)

    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.subplots_adjust(top=top)
    return title_bottom_frac


# ── Señal ──────────────────────────────────────────────────────────────────

def bandpass(signal, low, high, fs, order=4):
    sos = butter(order, [low, high], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, signal)


def band_power(signal, low, high, fs):
    """Potencia instantánea = envolvente de Hilbert al cuadrado."""
    return np.abs(hilbert(bandpass(signal, low, high, fs))) ** 2


def morlet_wavelet(freq, n_cycles, fs):
    sigma_t = n_cycles / (2 * np.pi * freq)
    t = np.arange(-4 * sigma_t, 4 * sigma_t, 1.0 / fs)
    w = np.exp(2j * np.pi * freq * t) * np.exp(-(t ** 2) / (2 * sigma_t ** 2))
    w /= np.sqrt(0.5 * np.sum(np.abs(w) ** 2))
    return w


def tfr_power(signal, freqs, fs):
    tfr = np.zeros((len(freqs), len(signal)))
    for i, f in enumerate(freqs):
        n_cycles = max(TFR_MIN_CYCLES, f * TFR_CYCLES_PER_HZ)
        w = morlet_wavelet(f, n_cycles, fs)
        tfr[i] = np.abs(fftconvolve(signal, w, mode="same")) ** 2
    return tfr


def resample(arr, n_out):
    if len(arr) < 2:
        return np.full(n_out, arr[0] if len(arr) else np.nan)
    return np.interp(np.linspace(0, 1, n_out), np.linspace(0, 1, len(arr)), arr)


def eje_tiempo_real(n_out_por_fase: int) -> np.ndarray:
    """Eje X en segundos (0-9 s) para 3 fases concatenadas: CROSS, INTENTION, MOTION."""
    x_cross = np.linspace(0.0, T_CROSS_END, n_out_por_fase, endpoint=False)
    x_int = np.linspace(T_CROSS_END, T_INTENTION_END, n_out_por_fase, endpoint=False)
    x_mot = np.linspace(T_INTENTION_END, T_MOTION_END, n_out_por_fase)
    return np.concatenate([x_cross, x_int, x_mot])


XLABEL_FASES = (f"Tiempo (s)  →  CROSS 0–{T_CROSS_END:.0f}s | INTENTION "
                f"{T_CROSS_END:.0f}–{T_INTENTION_END:.0f}s (línea sólida en "
                f"t={T_CROSS_END:.0f}) | MOTION {T_INTENTION_END:.0f}–"
                f"{T_MOTION_END:.0f}s (línea punteada en t={T_INTENTION_END:.0f})")


# ── Carga y segmentación ──────────────────────────────────────────────────

def detect_channels(df):
    return [c for c in df.select_dtypes(include=[np.number]).columns
            if c not in META_COLS]


def load_data(path):
    df = pd.read_csv(path)
    df = df[df["event_id"] != 999]
    return df.reset_index(drop=True)


def get_trial_segments(df):
    """dict trial -> índices de CROSS / INTENTION / MOTION (solo trials completos)."""
    segs, omitidos = {}, []
    for t in sorted(df["trial"].dropna().unique()):
        sub = df[df["trial"] == t]
        cross = sub.index[sub["event_id"] == EVENTS["CROSS"]].to_numpy()
        inten = sub.index[sub["event_id"] == EVENTS["INTENTION"]].to_numpy()
        motion = sub.index[sub["event_id"] == EVENTS["MOTION"]].to_numpy()
        if len(cross) and len(inten) and len(motion):
            segs[int(t)] = {"cross": cross, "intention": inten, "motion": motion}
        else:
            omitidos.append(int(t))
    return segs, omitidos


def prepare_corrida(path: Path) -> dict:
    """Carga un CSV filtrado y calcula UNA vez los pasos costosos (Hilbert y wavelets)."""
    df = load_data(path)
    if not CHANNELS:
        CHANNELS.extend(detect_channels(df))
    segs, omitidos = get_trial_segments(df)

    power = {}
    for ch in CHANNELS:
        for band_name, (low, high) in ALL_BANDS.items():
            power[(ch, band_name)] = band_power(df[ch].values, low, high, FS)
    tfr = {ch: tfr_power(df[ch].values, TFR_FREQS, FS) for ch in CHANNELS}

    return {"name": path.stem.replace("_filtrado", ""), "df": df, "segs": segs,
            "omitidos": omitidos, "pow": power, "tfr": tfr}


# ── SECCIÓN A: ERD/ERS dominio tiempo ─────────────────────────────────────

def erd_ers_for_corrida(c, trials):
    resumen = {(ch, b): [] for ch in CHANNELS for b in BANDS_BARRAS}
    timecourses = {(ch, b): [] for ch in CHANNELS for b in BANDS_BARRAS}
    for t in trials:
        seg = c["segs"][t]
        for ch in CHANNELS:
            for band_name in BANDS_BARRAS:
                power = c["pow"][(ch, band_name)]
                p_base = power[seg["cross"]].mean()
                if p_base <= 0:
                    continue
                p_task = power[seg["intention"]].mean()
                resumen[(ch, band_name)].append((p_task - p_base) / p_base * 100.0)

                cross_r = resample((power[seg["cross"]] - p_base) / p_base * 100.0, N_POINTS_TIMECOURSE)
                int_r = resample((power[seg["intention"]] - p_base) / p_base * 100.0, N_POINTS_TIMECOURSE)
                mot_r = resample((power[seg["motion"]] - p_base) / p_base * 100.0, N_POINTS_TIMECOURSE)
                timecourses[(ch, band_name)].append(np.concatenate([cross_r, int_r, mot_r]))
    return resumen, timecourses


def fdr_bh(pvals, alpha=0.05):
    pvals = np.asarray(pvals)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    thresh = alpha * (np.arange(1, n + 1) / n)
    below = ranked <= thresh
    if not np.any(below):
        return np.zeros(n, dtype=bool)
    cutoff = ranked[np.max(np.where(below))]
    return pvals <= cutoff


def build_report(resumen, n_label="n_trials"):
    """
    t-test (vs 0) por canal/banda + corrección FDR (Benjamini-Hochberg).
    Con n=1 (nivel por trial) no hay prueba de hipótesis posible: p=1.0 y
    el tipo se reporta como "(ERD/ERS no signif.)" según el signo.
    """
    rows, pvals = [], []
    for (ch, band_name), values in resumen.items():
        values = np.array(values)
        if len(values) == 0:
            continue
        mean_erd, std_erd = values.mean(), values.std()
        if len(values) > 1:
            _, pval = ttest_1samp(values, 0.0)
        else:
            pval = np.nan
        pval = pval if pval == pval else 1.0
        pvals.append(pval)
        rows.append({
            "canal": ch, "banda": band_name, n_label: len(values),
            "%cambio_promedio": round(mean_erd, 2), "std": round(std_erd, 2),
            "p_value": round(pval, 4),
        })

    cols = ["canal", "banda", n_label, "%cambio_promedio", "std", "tipo",
            "p_value", "significativo_FDR<0.05"]
    if not rows:
        return pd.DataFrame(columns=cols)

    sig_mask = fdr_bh(np.array(pvals), alpha=FDR_ALPHA)
    for row, sig in zip(rows, sig_mask):
        tipo = "ERD" if row["%cambio_promedio"] < 0 else "ERS"
        row["tipo"] = tipo if sig else f"({tipo} no signif.)"
        row["significativo_FDR<0.05"] = "sí" if sig else "no"
    return pd.DataFrame(rows)[cols].sort_values(["canal", "banda"]).reset_index(drop=True)


def plot_timecourses(timecourses, out_path, titulo_extra, bands_subset=None):
    bands_subset = bands_subset or BANDS
    fig, axes = plt.subplots(len(CHANNELS), 1, figsize=(9, 2.3 * len(CHANNELS)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, ch in zip(axes, CHANNELS):
        for band_name in bands_subset:
            arrs = timecourses.get((ch, band_name), [])
            if not arrs:
                continue
            mean_curve = np.nanmean(np.vstack(arrs), axis=0)
            ax.plot(eje_tiempo_real(len(mean_curve) // 3), mean_curve, label=band_name)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.axvline(T_CROSS_END, color="black", linewidth=0.8)
        ax.axvline(T_INTENTION_END, color="black", linewidth=0.8, linestyle=":")
        ax.set_xlim(0, T_MOTION_END)
        ax.set_ylabel(f"{ch}\n% vs baseline")
        ax.legend(fontsize=7, loc="upper right")
    axes[-1].set_xlabel(XLABEL_FASES)
    set_wrapped_title(fig, f"ERD (-) / ERS (+) promedio entre trials — {titulo_extra}")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_bars(report, out_path, titulo_extra):
    fig, ax = plt.subplots(figsize=(10, 4.5))
    bands = list(BANDS_BARRAS.keys())
    n_bands = len(bands)
    x = np.arange(len(CHANNELS))
    width = 0.8 / n_bands
    for i, band_name in enumerate(bands):
        sub = report[report["banda"] == band_name].set_index("canal").reindex(CHANNELS)
        vals = sub["%cambio_promedio"].values.astype(float)
        sig = sub["significativo_FDR<0.05"].values
        ax.bar(x + i * width, vals, width, label=band_name)
        for xi, v, s in zip(x + i * width, vals, sig):
            if s == "sí":
                ax.text(xi, v, "*", ha="center", va="bottom" if v >= 0 else "top", fontsize=14)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xticks(x + width * (n_bands - 1) / 2)
    ax.set_xticklabels(CHANNELS)
    ax.set_ylabel("% cambio (INTENTION vs baseline)")
    ax.legend()
    set_wrapped_title(fig, f"ERD/ERS promedio por canal y banda (* = p<0.05 FDR) — {titulo_extra}")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ── SECCIÓN D: potencia absoluta ──────────────────────────────────────────

def potencia_for_corrida(c, trials, bands_dict):
    timecourses = {(ch, b): [] for ch in CHANNELS for b in bands_dict}
    for t in trials:
        seg = c["segs"][t]
        for ch in CHANNELS:
            for band_name in bands_dict:
                power = c["pow"][(ch, band_name)]
                cross_r = resample(power[seg["cross"]], N_POINTS_TIMECOURSE)
                int_r = resample(power[seg["intention"]], N_POINTS_TIMECOURSE)
                mot_r = resample(power[seg["motion"]], N_POINTS_TIMECOURSE)
                timecourses[(ch, band_name)].append(np.concatenate([cross_r, int_r, mot_r]))
    return timecourses


def plot_potencia(timecourses, out_path, titulo_extra, bands_subset, colors=None):
    fig, axes = plt.subplots(len(CHANNELS), 1, figsize=(9, 2.3 * len(CHANNELS)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, ch in zip(axes, CHANNELS):
        for band_name in bands_subset:
            arrs = timecourses.get((ch, band_name), [])
            if not arrs:
                continue
            mean_curve = np.nanmean(np.vstack(arrs), axis=0)
            color = colors.get(band_name) if colors else None
            ax.plot(eje_tiempo_real(len(mean_curve) // 3), mean_curve, label=band_name, color=color)
        ax.axvline(T_CROSS_END, color="black", linewidth=0.8)
        ax.axvline(T_INTENTION_END, color="black", linewidth=0.8, linestyle=":")
        ax.set_xlim(0, T_MOTION_END)
        ax.set_ylabel(f"{ch}\nPotencia (uV²)")
        ax.legend(fontsize=7, loc="upper right")
    axes[-1].set_xlabel(XLABEL_FASES)
    set_wrapped_title(fig, f"Potencia absoluta promedio entre trials — {titulo_extra}")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ── SECCIÓN B: tiempo-frecuencia ──────────────────────────────────────────

def tfr_stack_for_corrida(c, trials):
    """Por canal, apila (n_trials, n_freqs, 3*N_POINTS_TFR): %ΔP vs CROSS y potencia absoluta."""
    stacks_pct, stacks_abs = {}, {}
    for ch in CHANNELS:
        tfr_full = c["tfr"][ch]
        filas_pct, filas_abs = [], []
        for t in trials:
            seg = c["segs"][t]
            pct_rows, abs_rows = [], []
            for fi in range(len(TFR_FREQS)):
                pot = tfr_full[fi]
                p_base = pot[seg["cross"]].mean()

                abs_rows.append(np.concatenate([
                    resample(pot[seg["cross"]], N_POINTS_TFR),
                    resample(pot[seg["intention"]], N_POINTS_TFR),
                    resample(pot[seg["motion"]], N_POINTS_TFR)]))

                if p_base <= 0:
                    pct_rows.append(np.full(3 * N_POINTS_TFR, np.nan))
                    continue
                pct_rows.append(np.concatenate([
                    resample((pot[seg["cross"]] - p_base) / p_base * 100.0, N_POINTS_TFR),
                    resample((pot[seg["intention"]] - p_base) / p_base * 100.0, N_POINTS_TFR),
                    resample((pot[seg["motion"]] - p_base) / p_base * 100.0, N_POINTS_TFR)]))
            filas_pct.append(np.vstack(pct_rows))
            filas_abs.append(np.vstack(abs_rows))
        if filas_pct:
            stacks_pct[ch] = np.stack(filas_pct, axis=0)
            stacks_abs[ch] = np.stack(filas_abs, axis=0)
    return stacks_pct, stacks_abs


def combinar_tfr(stacks_pct_list, stacks_abs_list):
    """Mapas promedio (%ΔP y potencia absoluta) + máscara de significancia FDR."""
    resultados = {}
    for ch in CHANNELS:
        piezas_pct = [s[ch] for s in stacks_pct_list if ch in s]
        piezas_abs = [s[ch] for s in stacks_abs_list if ch in s]
        if not piezas_pct:
            continue
        stack_pct = np.concatenate(piezas_pct, axis=0)
        stack_abs = np.concatenate(piezas_abs, axis=0)

        mapa_pct = np.nanmean(stack_pct, axis=0)
        mapa_potencia = np.nanmean(stack_abs, axis=0)

        n_trials, n_freqs, n_time = stack_pct.shape
        pvals = np.ones((n_freqs, n_time))
        if n_trials > 1:                       # con 1 solo trial no hay prueba t
            _, p = ttest_1samp(stack_pct, 0.0, axis=0, nan_policy="omit")
            pvals = np.where(np.isnan(p), 1.0, p)
        sig = fdr_bh(pvals.ravel(), alpha=FDR_ALPHA).reshape(n_freqs, n_time)
        resultados[ch] = {"mapa_pct": mapa_pct, "mapa_potencia": mapa_potencia,
                          "sig": sig, "n_trials": n_trials}
    return resultados


def plot_tfr_pct(resultados, out_path, titulo_extra):
    n_ch = len(resultados)
    fig, axes = plt.subplots(n_ch, 1, figsize=(9, 2.6 * n_ch), sharex=True)
    axes = np.atleast_1d(axes)
    x = eje_tiempo_real(N_POINTS_TFR)
    vmax = max(np.nanmax(np.abs(r["mapa_pct"])) for r in resultados.values())
    vmax = vmax if vmax > 0 else 1.0
    for ax, (ch, r) in zip(axes, resultados.items()):
        im = ax.pcolormesh(x, TFR_FREQS, r["mapa_pct"], cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, shading="auto")
        if r["sig"].any():
            ax.contour(x, TFR_FREQS, r["sig"].astype(float), levels=[0.5],
                       colors="black", linewidths=1.2)
        ax.axvline(T_CROSS_END, color="black", linewidth=0.8)
        ax.axvline(T_INTENTION_END, color="black", linewidth=0.8, linestyle=":")
        ax.set_xlim(0, T_MOTION_END)
        ax.set_ylabel(f"{ch}\nFrecuencia (Hz)")
        fig.colorbar(im, ax=ax, label="% vs baseline")
    axes[-1].set_xlabel(XLABEL_FASES)
    set_wrapped_title(fig, f"Mapa tiempo-frecuencia (ERSP), contorno=p<0.05 FDR — {titulo_extra}")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tfr_abs(resultados, out_path, titulo_extra):
    n_ch = len(resultados)
    fig, axes = plt.subplots(n_ch, 1, figsize=(9, 2.6 * n_ch), sharex=True)
    axes = np.atleast_1d(axes)
    x = eje_tiempo_real(N_POINTS_TFR)
    vmax = max(np.nanmax(r["mapa_potencia"]) for r in resultados.values())
    vmax = vmax if vmax > 0 else 1.0
    for ax, (ch, r) in zip(axes, resultados.items()):
        im = ax.pcolormesh(x, TFR_FREQS, r["mapa_potencia"], cmap="RdBu_r",
                           vmin=0, vmax=vmax, shading="auto")
        ax.axvline(T_CROSS_END, color="black", linewidth=0.8)
        ax.axvline(T_INTENTION_END, color="black", linewidth=0.8, linestyle=":")
        ax.set_xlim(0, T_MOTION_END)
        ax.set_ylabel(f"{ch}\nFrecuencia (Hz)")
        fig.colorbar(im, ax=ax, label="Potencia (uV²)")
    axes[-1].set_xlabel(XLABEL_FASES)
    set_wrapped_title(fig, f"Mapa tiempo-frecuencia (potencia absoluta) — {titulo_extra}")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ── SECCIÓN C: PSD Welch ──────────────────────────────────────────────────

def psd_segments_for_corrida(c, trials):
    segs = {ch: {"cross": [], "intention": []} for ch in CHANNELS}
    for t in trials:
        seg = c["segs"][t]
        for ch in CHANNELS:
            segs[ch]["cross"].append(c["df"][ch].values[seg["cross"]])
            segs[ch]["intention"].append(c["df"][ch].values[seg["intention"]])
    return {ch: (np.concatenate(segs[ch]["cross"]), np.concatenate(segs[ch]["intention"]))
            for ch in CHANNELS if segs[ch]["cross"]}


def combinar_psd(segmentos_list):
    resultados = {}
    for ch in CHANNELS:
        cross_piezas = [s[ch][0] for s in segmentos_list if ch in s]
        int_piezas = [s[ch][1] for s in segmentos_list if ch in s]
        if not cross_piezas:
            continue
        cross_sig = np.concatenate(cross_piezas)
        int_sig = np.concatenate(int_piezas)
        nperseg = max(int(min(FS * 1.0, len(cross_sig), len(int_sig))), 32)
        f_c, psd_c = welch(cross_sig, fs=FS, nperseg=nperseg)
        _, psd_i = welch(int_sig, fs=FS, nperseg=nperseg)
        resultados[ch] = {"freqs": f_c, "psd_cross": psd_c, "psd_intention": psd_i}
    return resultados


def plot_psd(resultados, out_path, titulo_extra):
    fig, axes = plt.subplots(len(resultados), 1, figsize=(8, 2.6 * len(resultados)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, (ch, r) in zip(axes, resultados.items()):
        for band_name, (lo, hi) in BAND_ZONES.items():
            ax.axvspan(lo, hi, color=BAND_ZONE_COLORS[band_name], alpha=0.15)
        ax.semilogy(r["freqs"], r["psd_cross"], label="CROSS (baseline)", color="tab:gray")
        ax.semilogy(r["freqs"], r["psd_intention"], label="INTENTION (tarea)", color="tab:blue")
        ax.set_xlim(6, 46)
        ax.set_ylim(*PSD_YLIM)
        ax.set_ylabel(f"{ch}\nPSD (uV²/Hz)")
        ax.legend(fontsize=7, loc="upper right")
    axes[-1].set_xlabel("Frecuencia (Hz)")

    legend_h_in = 0.30
    title_bottom = set_wrapped_title(
        fig, f"PSD Welch: CROSS vs INTENTION — {titulo_extra}",
        extra_bottom_in=legend_h_in + 0.08,
    )
    band_patches = [mpatches.Patch(facecolor=BAND_ZONE_COLORS[b], alpha=0.3, label=b)
                    for b in BAND_ZONES]
    fig.legend(handles=band_patches, loc="upper center", ncol=4, fontsize=7,
               bbox_to_anchor=(0.5, title_bottom - 0.01))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ── Orquestador ───────────────────────────────────────────────────────────

def analizar(items, out_dir: Path, titulo: str, modo: str = "trials", imprimir: bool = True):
    """
    items: lista de (corrida, [trials]).
    modo="trials"   -> estadística sobre los trials (POR CORRIDA y POR TRIAL).
    modo="corridas" -> las curvas juntan todos los trials, pero la estadística
                       usa el promedio de ERD% por corrida (n = n_corridas). GLOBAL.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    resumen_tot = {(ch, b): [] for ch in CHANNELS for b in BANDS_BARRAS}
    tc_tot = {(ch, b): [] for ch in CHANNELS for b in BANDS_BARRAS}
    pot_tot = {(ch, b): [] for ch in CHANNELS for b in POWER_BANDS_ALL}
    corrida_means = {(ch, b): [] for ch in CHANNELS for b in BANDS_BARRAS}
    filas_por_corrida = []
    stacks_pct, stacks_abs, segs_psd = [], [], []

    for c, trials in items:
        r, tc = erd_ers_for_corrida(c, trials)
        for k in r:
            resumen_tot[k].extend(r[k])
            tc_tot[k].extend(tc[k])
            vals = np.array(r[k])
            if len(vals) > 0:
                corrida_means[k].append(vals.mean())
                filas_por_corrida.append({
                    "corrida": c["name"], "canal": k[0], "banda": k[1],
                    "n_trials": len(vals), "%cambio_promedio": round(vals.mean(), 2),
                })
        pot = potencia_for_corrida(c, trials, POWER_BANDS_ALL)
        for k in pot:
            pot_tot[k].extend(pot[k])
        st_pct, st_abs = tfr_stack_for_corrida(c, trials)
        stacks_pct.append(st_pct)
        stacks_abs.append(st_abs)
        segs_psd.append(psd_segments_for_corrida(c, trials))

    if modo == "corridas":
        report = build_report(corrida_means, n_label="n_corridas")
    else:
        report = build_report(resumen_tot, n_label="n_trials")
    report.to_csv(out_dir / "reporte_erd_ers.csv", index=False)
    if imprimir:
        print(f"\n=== REPORTE ERD/ERS — {titulo} ===")
        print(report.to_string(index=False))

    if modo == "corridas" and filas_por_corrida:
        (pd.DataFrame(filas_por_corrida)
           .sort_values(["corrida", "canal", "banda"]).reset_index(drop=True)
           .to_csv(out_dir / "reporte_por_corrida.csv", index=False))

    # A: %ΔP
    plot_timecourses(tc_tot, out_dir / "erd_ers_timecourse.png", titulo, bands_subset=BANDS)
    plot_bars(report, out_dir / "erd_ers_barras.png", titulo)
    # D: potencia absoluta
    plot_potencia(pot_tot, out_dir / "potencia_timecourse_completa.png", titulo,
                  bands_subset=POWER_BAND_FULL)
    plot_potencia(pot_tot, out_dir / "potencia_timecourse_bandas.png", titulo,
                  bands_subset=BAND_ZONES, colors=BAND_ZONE_COLORS)
    # B: tiempo-frecuencia
    tfr_res = combinar_tfr(stacks_pct, stacks_abs)
    plot_tfr_pct(tfr_res, out_dir / "tfr_ersp.png", titulo)
    plot_tfr_abs(tfr_res, out_dir / "tfr_ersp_potencia_absoluta.png", titulo)
    # C: PSD
    plot_psd(combinar_psd(segs_psd), out_dir / "psd_baseline_vs_tarea.png", titulo)

    if imprimir:
        print(f"Guardado en: {out_dir}")


def main():
    ap = argparse.ArgumentParser(description="ERD/ERS: global, por corrida y por trial.")
    ap.add_argument("--sin-trial", action="store_true",
                    help="omite el nivel por trial (genera muchas gráficas; esto lo acelera)")
    args = ap.parse_args()

    files = sorted(INPUT_DIR.glob("*_filtrado.csv"))
    if not files:
        print(f"No hay archivos *_filtrado.csv en:\n  {INPUT_DIR}\n"
              f"Corre primero:  python filter_eeg.py")
        return

    print(f"\nCargando {len(files)} archivo(s) filtrado(s) y calculando Hilbert + wavelets...\n")
    corridas = []
    for f in files:
        c = prepare_corrida(f)
        if not c["segs"]:
            print(f"  ⚠  {f.name}: sin trials completos (CROSS+INTENTION+MOTION), se omite.")
            continue
        msg = f"  ✓  {c['name']}: {len(c['segs'])} trials completos"
        if c["omitidos"]:
            msg += f" (omitidos por eventos incompletos: {c['omitidos']})"
        print(msg)
        corridas.append(c)
    if not corridas:
        return
    print(f"\nCanales detectados: {', '.join(CHANNELS)}")

    # ---------- POR CORRIDA ----------
    print(f"\n{'=' * 70}\nANÁLISIS POR CORRIDA\n{'=' * 70}")
    for c in corridas:
        analizar([(c, sorted(c["segs"]))], OUT_POR_CORRIDA / c["name"],
                 f"POR CORRIDA: {c['name']}", modo="trials")

    # ---------- POR TRIAL ----------
    if not args.sin_trial:
        print(f"\n{'=' * 70}\nANÁLISIS POR TRIAL\n{'=' * 70}")
        n = 0
        for c in corridas:
            for t in sorted(c["segs"]):
                analizar([(c, [t])], OUT_POR_TRIAL / c["name"] / f"trial_{t}",
                         f"TRIAL {t} — {c['name']}", modo="trials", imprimir=False)
                n += 1
            print(f"  ✓  {c['name']}: {len(c['segs'])} trials")
        print(f"Por trial: {n} carpeta(s) generadas.")

    # ---------- GLOBAL ----------
    print(f"\n{'=' * 70}\nANÁLISIS GLOBAL\n{'=' * 70}")
    analizar([(c, sorted(c["segs"])) for c in corridas], OUT_GLOBAL,
             f"GLOBAL ({len(corridas)} corridas)", modo="corridas")

    print(f"\nListo. Revisa las carpetas en:\n  {OUT_ROOT}")


if __name__ == "__main__":
    main()
