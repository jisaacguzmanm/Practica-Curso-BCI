# ============================================================
# filter_eeg.py — Paso 2: Preprocesamiento (filtrado) de la señal EEG
# ============================================================
#
# Entrada : raw/*.csv           (señal cruda en uV, 4 corridas)
# Salida  : filtrado/*_filtrado.csv   (señal filtrada en uV)
#
# Pasos:
#   1. Lectura de los datos crudos.
#   2. Filtro Notch 60 Hz (+ armónicos hasta Nyquist) -> interferencia de red.
#   3. Filtro pasa-banda 6-47 Hz FIR fase cero -> banda de interés motor.
#   4. Manejo de bordes: reflect-padding para absorber el transitorio FIR.
#
# Las columnas de metadatos (datetime, timestamp, event_id, event_label,
# trial) se copian sin modificación. Los canales EEG se detectan
# automáticamente (todas las columnas numéricas que no son metadatos).
#
# Uso:   python filter_eeg.py
# ============================================================
import numpy as np
import pandas as pd
import mne
from pathlib import Path

mne.set_log_level("ERROR")

ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "raw"
OUTPUT_DIR = ROOT / "filtrado"

# ── Parámetros de filtrado ──────────────────────────────────────────────
SFREQ = 500.0       # frecuencia de muestreo (Hz)
LINE_FREQ = 60.0    # frecuencia de la red eléctrica (Hz)
L_FREQ = 6.0        # corte inferior del pasa-banda (Hz)
H_FREQ = 47.0       # corte superior del pasa-banda (Hz)
PAD_SAMPLES = 500   # muestras de reflect-padding en cada borde

META_COLS = {"datetime", "timestamp", "event_id", "event_label", "trial"}

# ─────────────────────────────────────────────────────────────────────────


def detect_channels(df: pd.DataFrame) -> list:
    """Canales EEG = columnas numéricas que no son metadatos."""
    return [c for c in df.select_dtypes(include=[np.number]).columns
            if c not in META_COLS]


def filter_signal(df: pd.DataFrame, channels: list) -> pd.DataFrame:
    """Aplica Notch + pasa-banda con reflect-padding. Devuelve un df nuevo (misma forma)."""
    data_V = df[channels].values.T * 1e-6          # uV -> V (MNE trabaja en V)

    # Reflect-padding al inicio para absorber el transitorio del filtro FIR
    pad = data_V[:, 1:PAD_SAMPLES + 1][:, ::-1]
    data_padded = np.concatenate([pad, data_V], axis=1)

    info = mne.create_info(ch_names=channels, sfreq=SFREQ, ch_types="eeg")
    raw = mne.io.RawArray(data_padded, info, verbose=False)

    raw.notch_filter(
        freqs=np.arange(LINE_FREQ, SFREQ / 2, LINE_FREQ),   # 60, 120, 180, 240 Hz
        fir_design="firwin", phase="zero-double", verbose=False,
    )
    raw.filter(
        l_freq=L_FREQ, h_freq=H_FREQ,
        fir_design="firwin", phase="zero-double", verbose=False,
    )

    data_filt = raw.get_data()[:, PAD_SAMPLES:]    # se descarta el padding
    df_out = df.copy()
    df_out[channels] = (data_filt * 1e6).T         # V -> uV
    return df_out


def filter_file(src: Path, dst: Path) -> None:
    df_raw = pd.read_csv(src)
    channels = detect_channels(df_raw)
    df_out = filter_signal(df_raw, channels)

    dst.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(dst, index=False)

    n_trials = df_out["trial"].nunique() if "trial" in df_out.columns else 0
    print(f"  ✓  {src.name}  →  {dst.name}   ({len(df_out)} muestras, {n_trials} valores de trial)")
    print(f"     Canales: {', '.join(channels)}")
    for ch in channels:
        v = df_out[ch].values
        print(f"       {ch}: min={v.min():.1f}  max={v.max():.1f}  std={v.std():.2f}  (uV, post-filtro)")


def main() -> None:
    csv_files = sorted(INPUT_DIR.glob("*.csv"))
    if not csv_files:
        print(f"No se encontraron archivos .csv en:\n  {INPUT_DIR}")
        return

    print(f"\nFiltrando {len(csv_files)} archivo(s): Notch {LINE_FREQ:.0f} Hz + "
          f"Pasa-banda {L_FREQ:.0f}-{H_FREQ:.0f} Hz  |  Padding: {PAD_SAMPLES} muestras\n")

    for src in csv_files:
        dst = OUTPUT_DIR / f"{src.stem}_filtrado.csv"
        filter_file(src, dst)
        print()

    print(f"Listo. Archivos filtrados en:\n  {OUTPUT_DIR}")
    print("\nSiguiente paso: python analizar_erd_ers.py")


if __name__ == "__main__":
    main()
