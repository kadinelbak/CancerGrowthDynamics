| Variant | Loss | 2b structure | 3a model (k) | β_C (3a) | ε_N (4) | ΣBIC 2b+3a+4 | Held-out Ce1 RMSE |
|---|---|---|---|---|---|---|---|
| V0_current (floor 50% (current)) | raw | 3 compartments → Hill kill | cc_K_alpha_deathN_beta (6) | 5.77 | 7.58 | 3443.4 | 711 |
| V1_floor0 (no floor) | raw | 3 compartments → Hill kill | cc_K_alpha_deathN_beta (6) | 9.36e-13 | 7.73 | 3419.5 | 711 |
| V2_floor25 (floor 25%) | raw | 3 compartments → Hill kill | cc_K_alpha_deathN_beta (6) | 2.89 | 7.66 | 3430.0 | 711 |
| V3_floor75 (floor 75%) | raw | 3 compartments → Hill kill | cc_K_alpha_deathN_beta (6) | 8.66 | 7.48 | 3460.0 | 711 |
| V4_ns (Norton-Simon kill (2b, 3a, 4)) | raw | 3 compartments → linear kill | cc_K_alpha_deathN_beta (6) | 0.202 | 8.25 | 3437.6 | 775 |
| V6_sd (current, SD-weighted loss) | sd | 3 compartments → Hill kill | cc_K_alpha_deathN_beta (6) | 4.22 | 9.35 | 626.0 | 681 |
| V5_ns_sd (Norton-Simon, SD-weighted loss) | sd | 3 compartments → linear kill | cc_K_alpha_beta (5) | 0.236 | 9.64 | 600.9 | 572 |
| V8_sd_floor0 (current, SD-weighted, no floor) | sd | 3 compartments → Hill kill | cc_full (7) | 1.53e-13 | 10.1 | 606.5 | — |
| V7_ns_sd_floor0 (Norton-Simon, SD-weighted, no floor) | sd | 3 compartments → linear kill | cc_full (7) | 5.52e-16 | 11 | 583.7 | — |
