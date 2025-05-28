This repository provides a collection of MATLAB scripts for analyzing **action potential (AP) features** from slice recordings, and **fiber photometry** and **ΔF/F₀ signal dynamics** from dopamine sensor imaging. 

| Script Name                    | Description                                                                 |
|--------------------------------|-----------------------------------------------------------------------------|
| `AP_Statistic.m`              | Extracts waveform parameters from action potential traces (e.g., spike amplitude, half-width, AHP). |
| `Code_FIcurve_extract_HERE.m` | Extracts firing properties across current injection steps (F-I curves).    |
| `Code_FIcurve_statistic_HERE.m` | Aggregates and visualizes F-I data; exports to Excel.                     |
| `eventedgefinder.m`           | Identifies onset/offset times of signal events using threshold detection.   |
| `Spontaneous_Freq.m`          | Measures spontaneous firing frequency using peak detection.                 |
| `fiberphotometryIp_AUC.m`     | Computes AUC for fiber photometry traces around L-DOPA injection; includes detrending. |
| `Sig_process_auto_deltaFF0.m` | Batch-processing script to compute ΔF/F₀ across trials using baseline fitting. |
| `SlowMode_AUC.m`              | Fits and subtracts slow exponential baseline changes; outputs AUC. |
