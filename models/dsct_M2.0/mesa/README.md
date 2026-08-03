# dsct_M2.0

2 M☉ δ Scuti model for the nonlinear mode-coupling pipeline. Target `log g = 3.900`.

```sh
mesasdk_init
./mk
./rn
python3 select_model.py
```

`select_model.py` picks the profile closest in `log g`, checks `R`, `ν_dyn` and
`log L`, and prints the path of the paired `.GYRE` file. Non-zero exit means the
model is wrong — don't run GYRE on it.

Current result: `LOGS/profile331.data.GYRE` (model 988), `log g = 3.8999`,
`R = 2.6278 R☉`, `Teff = 7586 K`, `√(GM/R³) = 2.866 c/d`.

## Extra profile columns

| column | |
|---|---|
| `dGamma1_dlnRho_s` | in neither GYRE's output nor MESA's defaults |
| `dGamma1_dlnRho_T`, `dGamma1_dlnT_Rho` | the two pieces of it |
| `t_thermal` | flags where adiabatic eigenfunctions stop being valid |
| `omega_conv` | for turbulent damping later |

Profile headers also carry `E_star_erg` and `dyn_freq_rad_per_s`.
