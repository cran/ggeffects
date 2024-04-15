# test_predictions, engine emmeans, 3-way interaction

    Code
      print(test_predictions(m, terms = c("x1", "x2", "x3"), engine = "emmeans"))
    Output
      # Pairwise comparisons
      
      x1  |  x2 |  x3 | Contrast |       95% CI |     p
      -------------------------------------------------
      0-1 | 1-1 | b-b |     1.60 | -0.14,  3.34 | 0.071
      0-0 | 1-2 | b-b |     0.15 | -0.87,  1.17 | 0.774
      0-1 | 1-2 | b-b |    -0.19 | -1.50,  1.12 | 0.773
      0-0 | 1-3 | b-b |     0.03 | -0.82,  0.87 | 0.952
      0-1 | 1-3 | b-b |     0.80 | -0.32,  1.93 | 0.160
      0-0 | 1-1 | b-a |     0.17 | -0.70,  1.04 | 0.701
      0-1 | 1-1 | b-a |     1.05 |  0.03,  2.08 | 0.043
      0-0 | 1-2 | b-a |     0.61 | -0.26,  1.48 | 0.169
      0-1 | 1-2 | b-a |     1.44 | -0.31,  3.18 | 0.105
      0-0 | 1-3 | b-a |    -0.16 | -1.00,  0.69 | 0.712
      0-1 | 1-3 | b-a |     1.13 | -0.62,  2.87 | 0.202
      0-0 | 1-1 | b-d |    -0.60 | -1.56,  0.35 | 0.214
      0-1 | 1-1 | b-d |    -0.88 | -2.00,  0.25 | 0.124
      0-0 | 1-2 | b-d |    -0.05 | -0.87,  0.78 | 0.912
      0-1 | 1-2 | b-d |    -0.23 | -1.53,  1.08 | 0.733
      0-0 | 1-3 | b-d |    -0.56 | -1.59,  0.46 | 0.275
      0-1 | 1-3 | b-d |     0.65 | -1.10,  2.39 | 0.463
      0-0 | 1-1 | b-c |     0.30 | -0.52,  1.12 | 0.469
      0-1 | 1-1 | b-c |     1.17 | -0.57,  2.92 | 0.185
      0-0 | 1-2 | b-c |     0.48 | -0.48,  1.43 | 0.325
      0-1 | 1-2 | b-c |     0.12 | -1.63,  1.86 | 0.894
      0-0 | 1-3 | b-c |     0.32 | -0.58,  1.23 | 0.479
      0-1 | 1-3 | b-c |    -0.39 | -2.13,  1.36 | 0.660
      1-0 | 1-2 | b-b |    -1.45 | -3.28,  0.37 | 0.116
      1-1 | 1-2 | b-b |    -1.79 | -3.79,  0.21 | 0.078
      1-0 | 1-3 | b-b |    -1.58 | -3.31,  0.15 | 0.074
      1-1 | 1-3 | b-b |    -0.80 | -2.68,  1.08 | 0.400
      1-0 | 1-1 | b-a |    -1.43 | -3.18,  0.31 | 0.106
      1-1 | 1-1 | b-a |    -0.55 | -2.37,  1.28 | 0.552
      1-0 | 1-2 | b-a |    -0.99 | -2.74,  0.75 | 0.260
      1-1 | 1-2 | b-a |    -0.17 | -2.47,  2.14 | 0.886
      1-0 | 1-3 | b-a |    -1.76 | -3.49, -0.03 | 0.046
      1-1 | 1-3 | b-a |    -0.48 | -2.78,  1.83 | 0.682
      1-0 | 1-1 | b-d |    -2.20 | -3.99, -0.42 | 0.016
      1-1 | 1-1 | b-d |    -2.48 | -4.36, -0.60 | 0.011
      1-0 | 1-2 | b-d |    -1.65 | -3.37,  0.07 | 0.060
      1-1 | 1-2 | b-d |    -1.83 | -3.82,  0.17 | 0.072
      1-0 | 1-3 | b-d |    -2.17 | -3.99, -0.34 | 0.021
      1-1 | 1-3 | b-d |    -0.96 | -3.26,  1.35 | 0.412
      1-0 | 1-1 | b-c |    -1.30 | -3.02,  0.42 | 0.136
      1-1 | 1-1 | b-c |    -0.43 | -2.74,  1.88 | 0.712
      1-0 | 1-2 | b-c |    -1.13 | -2.91,  0.66 | 0.213
      1-1 | 1-2 | b-c |    -1.48 | -3.79,  0.82 | 0.204
      1-0 | 1-3 | b-c |    -1.28 | -3.04,  0.48 | 0.153
      1-1 | 1-3 | b-c |    -1.99 | -4.29,  0.32 | 0.090
      0-1 | 2-2 | b-b |    -0.34 | -1.75,  1.07 | 0.636
      0-0 | 2-3 | b-b |    -0.12 | -1.12,  0.88 | 0.809
      0-1 | 2-3 | b-b |     0.65 | -0.59,  1.90 | 0.299
      0-0 | 2-1 | b-a |     0.02 | -1.00,  1.04 | 0.967
      0-1 | 2-1 | b-a |     0.91 | -0.25,  2.06 | 0.121
      0-0 | 2-2 | b-a |     0.46 | -0.56,  1.48 | 0.372
      0-1 | 2-2 | b-a |     1.29 | -0.54,  3.11 | 0.164
      0-0 | 2-3 | b-a |    -0.30 | -1.30,  0.69 | 0.545
      0-1 | 2-3 | b-a |     0.98 | -0.85,  2.80 | 0.289
      0-0 | 2-1 | b-d |    -0.75 | -1.84,  0.35 | 0.177
      0-1 | 2-1 | b-d |    -1.03 | -2.27,  0.22 | 0.105
      0-0 | 2-2 | b-d |    -0.19 | -1.17,  0.79 | 0.695
      0-1 | 2-2 | b-d |    -0.37 | -1.79,  1.04 | 0.601
      0-0 | 2-3 | b-d |    -0.71 | -1.87,  0.44 | 0.223
      0-1 | 2-3 | b-d |     0.50 | -1.32,  2.32 | 0.588
      0-0 | 2-1 | b-c |     0.15 | -0.83,  1.13 | 0.757
      0-1 | 2-1 | b-c |     1.02 | -0.80,  2.85 | 0.267
      0-0 | 2-2 | b-c |     0.33 | -0.77,  1.42 | 0.553
      0-1 | 2-2 | b-c |    -0.03 | -1.85,  1.79 | 0.973
      0-0 | 2-3 | b-c |     0.18 | -0.88,  1.23 | 0.740
      0-1 | 2-3 | b-c |    -0.53 | -2.36,  1.29 | 0.562
      1-0 | 2-3 | b-b |     0.22 | -1.07,  1.50 | 0.740
      1-1 | 2-3 | b-b |     0.99 | -0.50,  2.48 | 0.189
      1-0 | 2-1 | b-a |     0.36 | -0.95,  1.67 | 0.586
      1-1 | 2-1 | b-a |     1.24 | -0.17,  2.66 | 0.083
      1-0 | 2-2 | b-a |     0.80 | -0.51,  2.11 | 0.228
      1-1 | 2-2 | b-a |     1.62 | -0.37,  3.62 | 0.109
      1-0 | 2-3 | b-a |     0.03 | -1.26,  1.32 | 0.960
      1-1 | 2-3 | b-a |     1.32 | -0.68,  3.31 | 0.194
      1-0 | 2-1 | b-d |    -0.41 | -1.77,  0.95 | 0.551
      1-1 | 2-1 | b-d |    -0.69 | -2.18,  0.80 | 0.360
      1-0 | 2-2 | b-d |     0.14 | -1.13,  1.42 | 0.823
      1-1 | 2-2 | b-d |    -0.04 | -1.67,  1.60 | 0.966
      1-0 | 2-3 | b-d |    -0.37 | -1.79,  1.04 | 0.599
      1-1 | 2-3 | b-d |     0.84 | -1.16,  2.83 | 0.407
      1-0 | 2-1 | b-c |     0.49 | -0.78,  1.77 | 0.446
      1-1 | 2-1 | b-c |     1.36 | -0.64,  3.36 | 0.179
      1-0 | 2-2 | b-c |     0.67 | -0.70,  2.03 | 0.335
      1-1 | 2-2 | b-c |     0.31 | -1.69,  2.30 | 0.761
      1-0 | 2-3 | b-c |     0.51 | -0.82,  1.85 | 0.445
      1-1 | 2-3 | b-c |    -0.20 | -2.19,  1.80 | 0.845
      0-1 | 3-3 | b-b |     0.78 | -0.33,  1.88 | 0.165
      0-0 | 3-1 | b-a |     0.14 | -0.70,  0.99 | 0.736
      0-1 | 3-1 | b-a |     1.03 |  0.03,  2.03 | 0.044
      0-0 | 3-2 | b-a |     0.58 | -0.26,  1.43 | 0.173
      0-1 | 3-2 | b-a |     1.41 | -0.32,  3.14 | 0.109
      0-0 | 3-3 | b-a |    -0.18 | -1.00,  0.63 | 0.657
      0-1 | 3-3 | b-a |     1.10 | -0.63,  2.83 | 0.209
      0-0 | 3-1 | b-d |    -0.63 | -1.56,  0.30 | 0.184
      0-1 | 3-1 | b-d |    -0.90 | -2.01,  0.20 | 0.107
      0-0 | 3-2 | b-d |    -0.07 | -0.86,  0.72 | 0.857
      0-1 | 3-2 | b-d |    -0.25 | -1.54,  1.04 | 0.699
      0-0 | 3-3 | b-d |    -0.59 | -1.59,  0.41 | 0.243
      0-1 | 3-3 | b-d |     0.62 | -1.11,  2.35 | 0.477
      0-0 | 3-1 | b-c |     0.27 | -0.52,  1.07 | 0.492
      0-1 | 3-1 | b-c |     1.15 | -0.58,  2.88 | 0.191
      0-0 | 3-2 | b-c |     0.45 | -0.48,  1.38 | 0.339
      0-1 | 3-2 | b-c |     0.09 | -1.64,  1.82 | 0.917
      0-0 | 3-3 | b-c |     0.30 | -0.58,  1.18 | 0.502
      0-1 | 3-3 | b-c |    -0.41 | -2.14,  1.32 | 0.637
      1-0 | 3-1 | b-a |    -0.63 | -1.76,  0.49 | 0.266
      1-1 | 3-1 | b-a |     0.25 | -0.99,  1.50 | 0.687
      1-0 | 3-2 | b-a |    -0.19 | -1.32,  0.93 | 0.733
      1-1 | 3-2 | b-a |     0.63 | -1.25,  2.52 | 0.505
      1-0 | 3-3 | b-a |    -0.96 | -2.06,  0.15 | 0.088
      1-1 | 3-3 | b-a |     0.32 | -1.56,  2.21 | 0.733
      1-0 | 3-1 | b-d |    -1.40 | -2.59, -0.21 | 0.022
      1-1 | 3-1 | b-d |    -1.68 | -3.01, -0.35 | 0.014
      1-0 | 3-2 | b-d |    -0.85 | -1.94,  0.24 | 0.124
      1-1 | 3-2 | b-d |    -1.03 | -2.52,  0.46 | 0.173
      1-0 | 3-3 | b-d |    -1.37 | -2.61, -0.12 | 0.032
      1-1 | 3-3 | b-d |    -0.16 | -2.04,  1.73 | 0.870
      1-0 | 3-1 | b-c |    -0.50 | -1.59,  0.59 | 0.361
      1-1 | 3-1 | b-c |     0.37 | -1.51,  2.25 | 0.697
      1-0 | 3-2 | b-c |    -0.33 | -1.52,  0.86 | 0.586
      1-1 | 3-2 | b-c |    -0.69 | -2.57,  1.20 | 0.471
      1-0 | 3-3 | b-c |    -0.48 | -1.63,  0.67 | 0.411
      1-1 | 3-3 | b-c |    -1.19 | -3.07,  0.69 | 0.213
      0-1 | 1-1 | a-a |     0.89 | -0.14,  1.91 | 0.088
      0-0 | 1-2 | a-a |     0.44 | -0.43,  1.31 | 0.319
      0-1 | 1-2 | a-a |     1.27 | -0.48,  3.01 | 0.152
      0-0 | 1-3 | a-a |    -0.33 | -1.17,  0.52 | 0.444
      0-1 | 1-3 | a-a |     0.96 | -0.79,  2.70 | 0.278
      0-0 | 1-1 | a-d |    -0.77 | -1.72,  0.19 | 0.113
      0-1 | 1-1 | a-d |    -1.05 | -2.17,  0.08 | 0.068
      0-0 | 1-2 | a-d |    -0.21 | -1.04,  0.61 | 0.604
      0-1 | 1-2 | a-d |    -0.39 | -1.70,  0.91 | 0.550
      0-0 | 1-3 | a-d |    -0.73 | -1.76,  0.29 | 0.157
      0-1 | 1-3 | a-d |     0.48 | -1.27,  2.22 | 0.587
      0-0 | 1-1 | a-c |     0.13 | -0.69,  0.95 | 0.751
      0-1 | 1-1 | a-c |     1.00 | -0.74,  2.75 | 0.256
      0-0 | 1-2 | a-c |     0.31 | -0.65,  1.26 | 0.525
      0-1 | 1-2 | a-c |    -0.05 | -1.80,  1.69 | 0.952
      0-0 | 1-3 | a-c |     0.15 | -0.75,  1.06 | 0.735
      0-1 | 1-3 | a-c |    -0.56 | -2.30,  1.19 | 0.528
      1-0 | 1-2 | a-a |    -0.45 | -1.47,  0.58 | 0.387
      1-1 | 1-2 | a-a |     0.38 | -1.44,  2.20 | 0.679
      1-0 | 1-3 | a-a |    -1.21 | -2.21, -0.21 | 0.018
      1-1 | 1-3 | a-a |     0.07 | -1.75,  1.89 | 0.939
      1-0 | 1-1 | a-d |    -1.65 | -2.75, -0.56 | 0.004
      1-1 | 1-1 | a-d |    -1.93 | -3.18, -0.69 | 0.003
      1-0 | 1-2 | a-d |    -1.10 | -2.08, -0.12 | 0.028
      1-1 | 1-2 | a-d |    -1.28 | -2.69,  0.13 | 0.075
      1-0 | 1-3 | a-d |    -1.62 | -2.77, -0.47 | 0.007
      1-1 | 1-3 | a-d |    -0.41 | -2.23,  1.41 | 0.657
      1-0 | 1-1 | a-c |    -0.75 | -1.73,  0.23 | 0.130
      1-1 | 1-1 | a-c |     0.12 | -1.71,  1.94 | 0.899
      1-0 | 1-2 | a-c |    -0.58 | -1.67,  0.51 | 0.295
      1-1 | 1-2 | a-c |    -0.94 | -2.76,  0.89 | 0.309
      1-0 | 1-3 | a-c |    -0.73 | -1.78,  0.32 | 0.171
      1-1 | 1-3 | a-c |    -1.44 | -3.26,  0.38 | 0.120
      0-1 | 2-2 | a-a |     0.83 | -0.92,  2.57 | 0.348
      0-0 | 2-3 | a-a |    -0.77 | -1.61,  0.08 | 0.075
      0-1 | 2-3 | a-a |     0.52 | -1.23,  2.26 | 0.556
      0-0 | 2-1 | a-d |    -1.21 | -2.16, -0.25 | 0.014
      0-1 | 2-1 | a-d |    -1.49 | -2.61, -0.36 | 0.010
      0-0 | 2-2 | a-d |    -0.65 | -1.48,  0.17 | 0.117
      0-1 | 2-2 | a-d |    -0.83 | -2.14,  0.47 | 0.208
      0-0 | 2-3 | a-d |    -1.17 | -2.19, -0.15 | 0.025
      0-1 | 2-3 | a-d |     0.04 | -1.71,  1.78 | 0.965
      0-0 | 2-1 | a-c |    -0.31 | -1.13,  0.51 | 0.458
      0-1 | 2-1 | a-c |     0.56 | -1.18,  2.31 | 0.522
      0-0 | 2-2 | a-c |    -0.13 | -1.09,  0.82 | 0.782
      0-1 | 2-2 | a-c |    -0.49 | -2.23,  1.25 | 0.576
      0-0 | 2-3 | a-c |    -0.28 | -1.19,  0.62 | 0.534
      0-1 | 2-3 | a-c |    -0.99 | -2.74,  0.75 | 0.260
      1-0 | 2-3 | a-a |    -1.59 | -3.32,  0.14 | 0.071
      1-1 | 2-3 | a-a |    -0.31 | -2.62,  2.00 | 0.790
      1-0 | 2-1 | a-d |    -2.04 | -3.82, -0.25 | 0.026
      1-1 | 2-1 | a-d |    -2.31 | -4.20, -0.43 | 0.017
      1-0 | 2-2 | a-d |    -1.48 | -3.20,  0.24 | 0.090
      1-1 | 2-2 | a-d |    -1.66 | -3.66,  0.34 | 0.102
      1-0 | 2-3 | a-d |    -2.00 | -3.82, -0.18 | 0.032
      1-1 | 2-3 | a-d |    -0.79 | -3.10,  1.52 | 0.498
      1-0 | 2-1 | a-c |    -1.13 | -2.85,  0.58 | 0.193
      1-1 | 2-1 | a-c |    -0.26 | -2.57,  2.04 | 0.821
      1-0 | 2-2 | a-c |    -0.96 | -2.75,  0.83 | 0.288
      1-1 | 2-2 | a-c |    -1.32 | -3.62,  0.99 | 0.258
      1-0 | 2-3 | a-c |    -1.11 | -2.87,  0.65 | 0.213
      1-1 | 2-3 | a-c |    -1.82 | -4.13,  0.49 | 0.120
      0-1 | 3-3 | a-a |     1.28 | -0.45,  3.01 | 0.144
      0-0 | 3-1 | a-d |    -0.44 | -1.37,  0.49 | 0.345
      0-1 | 3-1 | a-d |    -0.72 | -1.82,  0.38 | 0.197
      0-0 | 3-2 | a-d |     0.11 | -0.68,  0.90 | 0.781
      0-1 | 3-2 | a-d |    -0.07 | -1.36,  1.22 | 0.916
      0-0 | 3-3 | a-d |    -0.41 | -1.41,  0.59 | 0.419
      0-1 | 3-3 | a-d |     0.80 | -0.93,  2.53 | 0.358
      0-0 | 3-1 | a-c |     0.46 | -0.33,  1.25 | 0.254
      0-1 | 3-1 | a-c |     1.33 | -0.40,  3.06 | 0.130
      0-0 | 3-2 | a-c |     0.63 | -0.30,  1.56 | 0.180
      0-1 | 3-2 | a-c |     0.27 | -1.46,  2.00 | 0.754
      0-0 | 3-3 | a-c |     0.48 | -0.40,  1.36 | 0.280
      0-1 | 3-3 | a-c |    -0.23 | -1.96,  1.50 | 0.793
      1-0 | 3-1 | a-d |    -1.73 | -3.51,  0.06 | 0.058
      1-1 | 3-1 | a-d |    -2.00 | -3.89, -0.12 | 0.037
      1-0 | 3-2 | a-d |    -1.17 | -2.89,  0.55 | 0.179
      1-1 | 3-2 | a-d |    -1.35 | -3.35,  0.65 | 0.182
      1-0 | 3-3 | a-d |    -1.69 | -3.51,  0.13 | 0.069
      1-1 | 3-3 | a-d |    -0.48 | -2.79,  1.83 | 0.680
      1-0 | 3-1 | a-c |    -0.82 | -2.54,  0.89 | 0.342
      1-1 | 3-1 | a-c |     0.05 | -2.26,  2.35 | 0.968
      1-0 | 3-2 | a-c |    -0.65 | -2.44,  1.14 | 0.471
      1-1 | 3-2 | a-c |    -1.01 | -3.31,  1.30 | 0.386
      1-0 | 3-3 | a-c |    -0.80 | -2.56,  0.96 | 0.368
      1-1 | 3-3 | a-c |    -1.51 | -3.82,  0.79 | 0.196
      0-1 | 1-1 | d-d |    -0.28 | -1.47,  0.91 | 0.644
      0-0 | 1-2 | d-d |     0.55 | -0.36,  1.46 | 0.229
      0-1 | 1-2 | d-d |     0.38 | -0.99,  1.74 | 0.585
      0-0 | 1-3 | d-d |     0.04 | -1.06,  1.13 | 0.948
      0-1 | 1-3 | d-d |     1.25 | -0.54,  3.03 | 0.169
      0-0 | 1-1 | d-c |     0.90 | -0.01,  1.81 | 0.052
      0-1 | 1-1 | d-c |     1.77 | -0.01,  3.56 | 0.052
      0-0 | 1-2 | d-c |     1.08 |  0.04,  2.11 | 0.041
      0-1 | 1-2 | d-c |     0.72 | -1.07,  2.50 | 0.427
      0-0 | 1-3 | d-c |     0.92 | -0.06,  1.91 | 0.066
      0-1 | 1-3 | d-c |     0.21 | -1.57,  2.00 | 0.812
      1-0 | 1-2 | d-d |     0.83 | -0.26,  1.92 | 0.132
      1-1 | 1-2 | d-d |     0.65 | -0.84,  2.14 | 0.385
      1-0 | 1-3 | d-d |     0.31 | -0.93,  1.56 | 0.618
      1-1 | 1-3 | d-d |     1.52 | -0.36,  3.41 | 0.111
      1-0 | 1-1 | d-c |     1.18 |  0.09,  2.27 | 0.034
      1-1 | 1-1 | d-c |     2.05 |  0.17,  3.93 | 0.033
      1-0 | 1-2 | d-c |     1.35 |  0.16,  2.54 | 0.027
      1-1 | 1-2 | d-c |     0.99 | -0.89,  2.88 | 0.296
      1-0 | 1-3 | d-c |     1.20 |  0.05,  2.35 | 0.041
      1-1 | 1-3 | d-c |     0.49 | -1.39,  2.37 | 0.605
      0-1 | 2-2 | d-d |    -0.18 | -1.45,  1.10 | 0.780
      0-0 | 2-3 | d-d |    -0.52 | -1.50,  0.46 | 0.296
      0-1 | 2-3 | d-d |     0.69 | -1.03,  2.41 | 0.425
      0-0 | 2-1 | d-c |     0.35 | -0.42,  1.12 | 0.372
      0-1 | 2-1 | d-c |     1.22 | -0.50,  2.94 | 0.162
      0-0 | 2-2 | d-c |     0.52 | -0.39,  1.43 | 0.257
      0-1 | 2-2 | d-c |     0.16 | -1.56,  1.88 | 0.851
      0-0 | 2-3 | d-c |     0.37 | -0.49,  1.23 | 0.394
      0-1 | 2-3 | d-c |    -0.34 | -2.06,  1.38 | 0.695
      1-0 | 2-3 | d-d |    -0.34 | -1.75,  1.07 | 0.634
      1-1 | 2-3 | d-d |     0.87 | -1.13,  2.87 | 0.388
      1-0 | 2-1 | d-c |     0.53 | -0.75,  1.80 | 0.414
      1-1 | 2-1 | d-c |     1.40 | -0.60,  3.39 | 0.168
      1-0 | 2-2 | d-c |     0.70 | -0.66,  2.06 | 0.310
      1-1 | 2-2 | d-c |     0.34 | -1.66,  2.34 | 0.734
      1-0 | 2-3 | d-c |     0.55 | -0.78,  1.88 | 0.414
      1-1 | 2-3 | d-c |    -0.16 | -2.16,  1.84 | 0.873
      0-1 | 3-3 | d-d |     1.21 | -0.61,  3.03 | 0.190
      0-0 | 3-1 | d-c |     0.86 | -0.12,  1.84 | 0.083
      0-1 | 3-1 | d-c |     1.74 | -0.09,  3.56 | 0.062
      0-0 | 3-2 | d-c |     1.04 | -0.05,  2.13 | 0.062
      0-1 | 3-2 | d-c |     0.68 | -1.14,  2.50 | 0.459
      0-0 | 3-3 | d-c |     0.89 | -0.16,  1.94 | 0.097
      0-1 | 3-3 | d-c |     0.18 | -1.65,  2.00 | 0.846
      1-0 | 3-1 | d-c |    -0.35 | -2.06,  1.37 | 0.690
      1-1 | 3-1 | d-c |     0.53 | -1.78,  2.83 | 0.651
      1-0 | 3-2 | d-c |    -0.17 | -1.96,  1.62 | 0.849
      1-1 | 3-2 | d-c |    -0.53 | -2.84,  1.78 | 0.649
      1-0 | 3-3 | d-c |    -0.32 | -2.08,  1.44 | 0.716
      1-1 | 3-3 | d-c |    -1.03 | -3.34,  1.27 | 0.375
      0-1 | 1-1 | c-c |     0.87 | -0.85,  2.59 | 0.316
      0-0 | 1-2 | c-c |     0.17 | -0.73,  1.08 | 0.703
      0-1 | 1-2 | c-c |    -0.18 | -1.90,  1.54 | 0.832
      0-0 | 1-3 | c-c |     0.02 | -0.84,  0.88 | 0.957
      0-1 | 1-3 | c-c |    -0.69 | -2.41,  1.03 | 0.429
      1-0 | 1-2 | c-c |    -0.70 | -2.48,  1.09 | 0.440
      1-1 | 1-2 | c-c |    -1.06 | -3.36,  1.25 | 0.365
      1-0 | 1-3 | c-c |    -0.85 | -2.61,  0.91 | 0.341
      1-1 | 1-3 | c-c |    -1.56 | -3.86,  0.75 | 0.183
      0-1 | 2-2 | c-c |    -0.36 | -2.15,  1.43 | 0.690
      0-0 | 2-3 | c-c |    -0.15 | -1.14,  0.84 | 0.761
      0-1 | 2-3 | c-c |    -0.86 | -2.65,  0.93 | 0.340
      1-0 | 2-3 | c-c |     0.21 | -1.55,  1.97 | 0.815
      1-1 | 2-3 | c-c |    -0.50 | -2.81,  1.80 | 0.665
      0-1 | 3-3 | c-c |    -0.71 | -2.47,  1.05 | 0.425

---

    Code
      print(test_predictions(m, terms = c("x1", "x2", "x3"), engine = "emmeans",
      test = "interaction"))
    Output
      # Interaction contrasts
      
      x1  |      x2 |      x3 | Contrast |      95% CI |     p
      --------------------------------------------------------
      0-1 | 1 and 2 | b and a |     1.88 | -1.14, 4.90 | 0.219
      0-1 | 1 and 3 | b and a |     1.22 | -1.66, 4.10 | 0.401
      0-1 | 2 and 3 | b and a |    -0.66 | -3.70, 2.38 | 0.667
      0-1 | 1 and 2 | b and d |     2.04 | -0.80, 4.88 | 0.158
      0-1 | 1 and 3 | b and d |     2.31 | -0.69, 5.31 | 0.129
      0-1 | 2 and 3 | b and d |     0.28 | -2.58, 3.13 | 0.848
      0-1 | 1 and 2 | b and c |     0.71 | -2.63, 4.05 | 0.674
      0-1 | 1 and 3 | b and c |    -0.76 | -3.97, 2.46 | 0.641
      0-1 | 2 and 3 | b and c |    -1.46 | -4.55, 1.62 | 0.347
      0-1 | 1 and 2 | a and d |     0.16 | -2.51, 2.83 | 0.907
      0-1 | 1 and 3 | a and d |     1.09 | -1.87, 4.05 | 0.465
      0-1 | 2 and 3 | a and d |     0.93 | -2.38, 4.25 | 0.576
      0-1 | 1 and 2 | a and c |    -1.17 | -4.37, 2.03 | 0.468
      0-1 | 1 and 3 | a and c |    -1.98 | -5.15, 1.20 | 0.219
      0-1 | 2 and 3 | a and c |    -0.81 | -4.32, 2.70 | 0.649
      0-1 | 1 and 2 | d and c |    -1.33 | -4.36, 1.70 | 0.386
      0-1 | 1 and 3 | d and c |    -3.07 | -6.36, 0.22 | 0.067
      0-1 | 2 and 3 | d and c |    -1.74 | -5.09, 1.61 | 0.304

# test_predictions, engine emmeans, by and variable name = level value

    Code
      print(test_predictions(m, c("time", "coffee"), engine = "emmeans"))
    Output
      # Pairwise comparisons
      
      time                |          coffee | Contrast |       95% CI |      p
      ------------------------------------------------------------------------
      morning-noon        |   coffee-coffee |     1.93 | -2.02,  5.88 | 0.336 
      morning-afternoon   |   coffee-coffee |    -1.93 | -5.88,  2.02 | 0.336 
      morning-morning     |  coffee-control |     5.78 |  1.83,  9.73 | 0.004 
      morning-noon        |  coffee-control |     0.00 | -3.95,  3.95 | > .999
      morning-afternoon   |  coffee-control |     0.00 | -3.95,  3.95 | > .999
      noon-afternoon      |   coffee-coffee |    -3.86 | -7.81,  0.09 | 0.056 
      noon-morning        |  coffee-control |     3.86 | -0.09,  7.81 | 0.056 
      noon-noon           |  coffee-control |    -1.93 | -5.88,  2.02 | 0.336 
      noon-afternoon      |  coffee-control |    -1.93 | -5.88,  2.02 | 0.336 
      afternoon-morning   |  coffee-control |     7.71 |  3.76, 11.66 | < .001
      afternoon-noon      |  coffee-control |     1.93 | -2.02,  5.88 | 0.336 
      afternoon-afternoon |  coffee-control |     1.93 | -2.02,  5.88 | 0.336 
      morning-noon        | control-control |    -5.78 | -9.73, -1.83 | 0.004 
      morning-afternoon   | control-control |    -5.78 | -9.73, -1.83 | 0.004 
      noon-afternoon      | control-control |     0.00 | -3.95,  3.95 | > .999

---

    Code
      print(test_predictions(m, c("time", "coffee"), by = "sex", engine = "emmeans"))
    Output
      # Pairwise comparisons
      
      sex = female
      
      time                |          coffee | Contrast |       95% CI |      p
      ------------------------------------------------------------------------
      morning-noon        |   coffee-coffee |     1.93 | -2.02,  5.88 | 0.336 
      morning-afternoon   |   coffee-coffee |    -1.93 | -5.88,  2.02 | 0.336 
      morning-morning     |  coffee-control |     5.78 |  1.83,  9.73 | 0.004 
      morning-noon        |  coffee-control |     0.00 | -3.95,  3.95 | > .999
      morning-afternoon   |  coffee-control |     0.00 | -3.95,  3.95 | > .999
      noon-afternoon      |   coffee-coffee |    -3.86 | -7.81,  0.09 | 0.056 
      noon-morning        |  coffee-control |     3.86 | -0.09,  7.81 | 0.056 
      noon-noon           |  coffee-control |    -1.93 | -5.88,  2.02 | 0.336 
      noon-afternoon      |  coffee-control |    -1.93 | -5.88,  2.02 | 0.336 
      afternoon-morning   |  coffee-control |     7.71 |  3.76, 11.66 | < .001
      afternoon-noon      |  coffee-control |     1.93 | -2.02,  5.88 | 0.336 
      afternoon-afternoon |  coffee-control |     1.93 | -2.02,  5.88 | 0.336 
      morning-noon        | control-control |    -5.78 | -9.73, -1.83 | 0.004 
      morning-afternoon   | control-control |    -5.78 | -9.73, -1.83 | 0.004 
      noon-afternoon      | control-control |     0.00 | -3.95,  3.95 | > .999
      
      sex = male
      
      time                |          coffee | Contrast |       95% CI |      p
      ------------------------------------------------------------------------
      morning-noon        |   coffee-coffee |     1.93 | -2.02,  5.88 | 0.336 
      morning-afternoon   |   coffee-coffee |    -1.93 | -5.88,  2.02 | 0.336 
      morning-morning     |  coffee-control |     5.78 |  1.83,  9.73 | 0.004 
      morning-noon        |  coffee-control |     0.00 | -3.95,  3.95 | > .999
      morning-afternoon   |  coffee-control |     0.00 | -3.95,  3.95 | > .999
      noon-afternoon      |   coffee-coffee |    -3.86 | -7.81,  0.09 | 0.056 
      noon-morning        |  coffee-control |     3.86 | -0.09,  7.81 | 0.056 
      noon-noon           |  coffee-control |    -1.93 | -5.88,  2.02 | 0.336 
      noon-afternoon      |  coffee-control |    -1.93 | -5.88,  2.02 | 0.336 
      afternoon-morning   |  coffee-control |     7.71 |  3.76, 11.66 | < .001
      afternoon-noon      |  coffee-control |     1.93 | -2.02,  5.88 | 0.336 
      afternoon-afternoon |  coffee-control |     1.93 | -2.02,  5.88 | 0.336 
      morning-noon        | control-control |    -5.78 | -9.73, -1.83 | 0.004 
      morning-afternoon   | control-control |    -5.78 | -9.73, -1.83 | 0.004 
      noon-afternoon      | control-control |     0.00 | -3.95,  3.95 | > .999
