# Guía de ejecución — ERD/ERS con EEG (Curso BCI + EEG)

En esta práctica vas a recorrer el flujo completo de análisis de imaginería motora:
**señal cruda → filtrado → ERD/ERS**, revisando tú mismo cada resultado.

## Contenido del repositorio

```
├── raw/                    ← 4 corridas de una misma persona (CSV), 10 trials cada una
├── filter_eeg.py           ← Paso 3: filtrado de la señal
├── analizar_erd_ers.py     ← Paso 4: análisis ERD/ERS y gráficas
├── requirements.txt        ← librerías necesarias
└── README.md               ← esta guía
```

## Requisitos previos

- Descarga el repositorio (botón verde **Code → Download ZIP**, y descomprime) o clónalo con `git clone`.
- Abre una terminal **dentro de la carpeta del repositorio** (donde está `requirements.txt`). En Windows, usa **PowerShell**.
- Comprueba tu versión de **Python**:

  ```
  python --version
  ```

  > En Mac/Linux, si `python` no funciona, usa `python3` en todos los comandos.
  > En Windows, si no funciona, prueba con `py`.

### ¿Qué versión de Python tengo?

**Si tu versión es 3.10, 3.11 o 3.12** → todo bien, **sáltate esta sección** y continúa directo con el **Paso 1**.

**Si tu versión NO es 3.10, 3.11 ni 3.12** (por ejemplo 3.9, 3.13 o mayor, o no tienes Python) → sigue estos pasos en Windows:

**A. Instalar uv (una sola vez) y reiniciar PowerShell**

```
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Cierra PowerShell por completo y vuelve a abrirlo dentro de la carpeta del repositorio.

**B. Crear y activar un entorno con Python 3.11**

```
uv venv --python 3.11 --seed
.venv\Scripts\activate
```

Comprueba que quedó bien: `python --version` debe mostrar `Python 3.11.x`. Luego continúa con el **Paso 1**.

> - Cada vez que abras una terminal nueva, repite solo `.venv\Scripts\activate` antes de ejecutar los scripts.
> - Si al activar aparece un error de "ejecución de scripts deshabilitada", ejecuta
>   `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass` y vuelve a activar.

---

## Paso 1 — Observa la señal cruda (manual)

Abre la carpeta `raw/` aquí mismo en GitHub y mira los 4 archivos `.csv`.

| Columna | Qué es |
|---|---|
| `datetime`, `timestamp` | marcas de tiempo |
| `event_id` / `event_label` | fase del trial: `100` CROSS, `200` INTENTION, `300` MOTION, `400` RESET |
| `trial` | número de trial (1 a 10) |
| resto de columnas | canales EEG en microvoltios (uV), muestreados a 500 Hz |

Cada trial dura: **CROSS 3 s** (cruz de fijación, es el *baseline*) → **INTENTION 3 s** (imaginas el movimiento) → **MOTION 3 s** → **RESET**.

**Qué observar:** la señal cruda es "ruidosa": tiene un desfase (offset) y mucha interferencia de la red eléctrica (60 Hz).

---

## Paso 2 — Instalar las librerías

Copia y pega en la terminal:

```
python -m pip install -r requirements.txt
```

## Paso 3 — Filtrar la señal

Copia y pega en la terminal:

```
python filter_eeg.py
```

Esto aplica:

- **Notch a 60 Hz** (y armónicos: 120, 180 y 240 Hz) → elimina la interferencia de la red eléctrica.
- **Pasa-banda 6–47 Hz** (FIR, fase cero) → conserva la banda de interés para imaginería motora (Mu, Beta y parte de Gamma).

Se crea la carpeta **`filtrado/`** con un archivo `*_filtrado.csv` por cada corrida.

### Revisión manual (Paso 3)

Abre un archivo de `filtrado/` junto a su equivalente de `raw/` y compara:

- ¿Qué pasó con el offset de la señal (valores promedio)?
- ¿Cambió el rango de amplitud (min/max)?
- ¿Se conservan las columnas de eventos y trials igual que antes?

---

## Paso 4 — Calcular ERD/ERS

Copia y pega en la terminal:

```
python analizar_erd_ers.py
```

Tarda unos minutos. Se crea la carpeta **`resultados/`** con:

```
resultados/
├── global/                       ← las 4 corridas juntas (n = 4 corridas)
├── por_corrida/<corrida>/        ← una carpeta por corrida (promedio de sus 10 trials)
└── por_trial/<corrida>/trial_N/  ← un trial suelto, sin promediar
```

> Si tu computador es lento y solo quieres ver los niveles global y por corrida:
> `python analizar_erd_ers.py --sin-trial`

### Qué hay en cada carpeta

Todas las carpetas contienen el mismo conjunto de archivos:

| Archivo | Qué muestra |
|---|---|
| `reporte_erd_ers.csv` | tabla con el % de cambio de potencia por canal y banda, tipo (ERD/ERS), p-value y significancia (FDR) |
| `erd_ers_timecourse.png` | curva de % de cambio de potencia (Mu y Beta) durante los 9 s del trial |
| `erd_ers_barras.png` | barras de % de cambio por canal y banda (Mu, Beta, Gamma); `*` = significativo |
| `potencia_timecourse_completa.png` | potencia absoluta (uV²) en la banda completa 8–45 Hz |
| `potencia_timecourse_bandas.png` | potencia absoluta separada por bandas (Mu, Beta Low, Beta High, Gamma) |
| `tfr_ersp.png` | mapa tiempo-frecuencia en % vs baseline; el contorno negro marca zonas significativas |
| `tfr_ersp_potencia_absoluta.png` | mapa tiempo-frecuencia en potencia absoluta (uV²) |
| `psd_baseline_vs_tarea.png` | espectro de potencia: CROSS (baseline) vs INTENTION (tarea) |

`resultados/global/` incluye además `reporte_por_corrida.csv`, que permite ver si las 4 corridas son consistentes entre sí.

---

## Paso 5 — Revisión manual de resultados

Recorre las gráficas, empezando por `resultados/global/`:

1. **`erd_ers_barras.png`**: ¿en qué canales y bandas aparece ERD (barra negativa) o ERS (barra positiva)? ¿cuáles tienen `*`?
2. **`erd_ers_timecourse.png`**: ¿en qué momento exacto cae la potencia? ¿ocurre al empezar INTENTION (t = 3 s)?
3. **`tfr_ersp.png`**: ¿en qué frecuencias aparece el cambio (azul = ERD, rojo = ERS)?
4. **`psd_baseline_vs_tarea.png`**: ¿las curvas de CROSS e INTENTION se separan en la banda Mu (8–13 Hz)?
5. Compara `global/` con `por_corrida/`: ¿el patrón se repite en las 4 corridas?
6. Abre algunos `por_trial/`: ¿por qué un trial suelto es mucho más ruidoso que el promedio?

### Cómo interpretar

- **ERD (Event-Related Desynchronization)**: la potencia **baja** respecto al baseline (% negativo). Típico de Mu/Beta al imaginar movimiento.
- **ERS (Event-Related Synchronization)**: la potencia **sube** respecto al baseline (% positivo).
- `%cambio = (P_INTENTION − P_CROSS) / P_CROSS × 100`.
- En el nivel **por trial** no hay prueba estadística (solo hay un dato): `p_value = 1.0` y el tipo aparece como `(ERD no signif.)` o `(ERS no signif.)`. Son gráficas descriptivas.

---

## Problemas frecuentes

| Problema | Solución |
|---|---|
| `python: command not found` | Prueba `python3` (Mac/Linux) o `py` (Windows), o sigue la sección "¿Qué versión de Python tengo?" |
| `No module named ...` | Repite el Paso 2 (y si usas el entorno de uv, actívalo antes con `.venv\Scripts\activate`) |
| `'uv' no se reconoce como un comando` | Cierra PowerShell por completo y ábrelo de nuevo |
| Error al instalar `numpy`, `scipy` o `mne` | Tu Python probablemente es muy nuevo o muy viejo: sigue la sección "¿Qué versión de Python tengo?" |
| `No se encontraron archivos .csv` | Ejecuta los comandos desde la carpeta raíz del repositorio |
| `No hay archivos *_filtrado.csv` | Ejecuta primero el Paso 3 |
