import logging
import math
import os.path
import pprint
import statistics
from pathlib import Path
from typing import List, Union

from util import PDBtool

"""
Taken from https://www.ebi.ac.uk/thornton-srv/databases/sidechains/index.html on 2023-12-01.
"""
web_atlas = {
    "total": 482555,
    # Alanine
    "numALA": 40341,
    "ALA-ALA": 6850,
    "ALA-ARG": 3967,  # ARG-ALA: 3955
    "ALA-ASN": 2304,  # ASN-ALA: 2302
    "ALA-ASP": 2500,  # ASP-ALA: 2495
    "ALA-CYS": 1353,  # CYS-ALA: 1353
    "ALA-GLN": 2310,  # GLN-ALA: 2307
    "ALA-GLU": 2825,  # GLU-ALA: 2808
    "ALA-GLY": 3678,  # GLY-ALA: 3678
    "ALA-HIS": 1812,  # HIS-ALA: 1811
    "ALA-ILE": 9519,  # ILE-ALA: 9518
    "ALA-LEU": 15282,  # LEU-ALA: 15280
    "ALA-LYS": 2180,  # LYS-ALA: 2140
    "ALA-MET": 2857,  # MET-ALA: 2855
    "ALA-PHE": 7233,  # PHE-ALA: 7232
    "ALA-PRO": 2975,  # PRO-ALA: 2975
    "ALA-SER": 2934,  # SER-ALA: 2931
    "ALA-THR": 4081,  # THR-ALA: 4080
    "ALA-TRP": 2680,  # TRP-ALA: 2680
    "ALA-TYR": 5144,  # TYR-ALA: 5144
    "ALA-VAL": 10941,  # VAL-ALA: 10941
    # Arginine
    "numARG": 24542,
    "ARG-ALA": 3955,
    "ARG-ARG": 4441,
    "ARG-ASN": 2986,  # ASN-ARG: 3003
    "ARG-ASP": 7391,  # ASP-ARG: 7416
    "ARG-CYS": 856,  # CYS-ARG: 860
    "ARG-GLN": 3185,  # GLN-ARG: 3189
    "ARG-GLU": 9671,  # GLU-ARG: 9678
    "ARG-GLY": 3814,  # GLY-ARG: 3831
    "ARG-HIS": 2092,  # HIS-ARG: 2096
    "ARG-ILE": 4608,  # ILE-ARG: 4636
    "ARG-LEU": 8397,  # LEU-ARG: 8423
    "ARG-LYS": 2282,  # LYS-ARG: 2232
    "ARG-MET": 1620,  # MET-ARG: 1628
    "ARG-PHE": 4114,  # PHE-ARG: 4135
    "ARG-PRO": 2942,  # PRO-ARG: 2955
    "ARG-SER": 3230,  # SER-ARG: 3239
    "ARG-THR": 3862,  # THR-ARG: 3886
    "ARG-TRP": 2291,  # TRP-ARG: 2303
    "ARG-TYR": 4587,  # TYR-ARG: 4613
    "ARG-VAL": 5178,  # VAL-ARG: 5190
    # Asparagine
    "numASN": 20874,
    "ASN-ALA": 2302,
    "ASN-ARG": 3003,
    "ASN-ASN": 2821,
    "ASN-ASP": 2844,  # ASP-ASN: 2845
    "ASN-CYS": 603,  # CYS-ASN: 604
    "ASN-GLN": 2227,  # GLN-ASN: 2217
    "ASN-GLU": 2861,  # GLU-ASN: 2840
    "ASN-GLY": 2753,  # GLY-ASN: 2755
    "ASN-HIS": 1332,  # HIS-ASN: 1336
    "ASN-ILE": 2740,  # ILE-ASN: 2744
    "ASN-LEU": 3996,  # LEU-ASN: 3999
    "ASN-LYS": 2435,  # LYS-ASN: 2395
    "ASN-MET": 1052,  # MET-ASN: 1053
    "ASN-PHE": 2470,  # PHE-ASN: 2468
    "ASN-PRO": 1695,  # PRO-ASN: 1695
    "ASN-SER": 2488,  # SER-ASN: 2487
    "ASN-THR": 2991,  # THR-ASN: 2993
    "ASN-TRP": 1317,  # TRP-ASN: 1317
    "ASN-TYR": 2735,  # TYR-ASN: 2736
    "ASN-VAL": 3001,  # VAL-ASN: 3003
    # Aspartate
    "numASP": 27965,
    "ASP-ALA": 2495,
    "ASP-ARG": 7416,
    "ASP-ASN": 2845,
    "ASP-ASP": 2058,
    "ASP-CYS": 560,  # CYS-ASP: 560
    "ASP-GLN": 2382,  # GLN-ASP: 2375
    "ASP-GLU": 1897,  # GLU-ASP: 1891
    "ASP-GLY": 3133,  # GLY-ASP: 3136
    "ASP-HIS": 2383,  # HIS-ASP: 2384
    "ASP-ILE": 2567,  # ILE-ASP: 2586
    "ASP-LEU": 4050,  # LEU-ASP: 4056
    "ASP-LYS": 5989,  # LYS-ASP: 5883
    "ASP-MET": 946,  # MET-ASP: 948
    "ASP-PHE": 2290,  # PHE-ASP: 2291
    "ASP-PRO": 1753,  # PRO-ASP: 1754
    "ASP-SER": 3276,  # SER-ASP: 3279
    "ASP-THR": 3435,  # THR-ASP: 3435
    "ASP-TRP": 1302,  # TRP-ASP 1302
    "ASP-TYR": 3533,  # TYR-ASP: 3531
    "ASP-VAL": 3001,  # VAL-ASP: 3004
    # CYSTEINE
    "numCYS": 6510,
    "CYS-ALA": 1353,
    "CYS-ARG": 860,
    "CYS-ASN": 604,
    "CYS-ASP": 560,
    "CYS-CYS": 3490,
    "CYS-GLN": 606,  # GLN-CYS: 603
    "CYS-GLU": 600,  # GLU-CYS: 596
    "CYS-GLY": 999,  # GLY-CYS: 999
    "CYS-HIS": 639,  # HIS-CYS: 639
    "CYS-ILE": 2072,  # ILE-CYS: 2072
    "CYS-LEU": 3296,  # LEU-CYS: 3295
    "CYS-LYS": 597,  # LYS-CYS: 587
    "CYS-MET": 641,  # MET-CYS: 641
    "CYS-PHE": 1987,  # PHE-CYS: 1987
    "CYS-PRO": 924,  # PRO-CYS: 924
    "CYS-SER": 815,  # SER-CYS: 814
    "CYS-THR": 967,  # THR-CYS: 967
    "CYS-TRP": 709,  # TRP-CYS: 709
    "CYS-TYR": 1287,  # TYR-CYS: 1287
    "CYS-VAL": 2341,  # VAL-CYS: 2340
    # GLUTAMINE
    "numGLN": 18401,
    "GLN-ALA": 2307,
    "GLN-ARG": 3189,
    "GLN-ASN": 2217,
    "GLN-ASP": 2375,
    "GLN-CYS": 603,
    "GLN-GLN": 2172,
    "GLN-GLU": 2591,  # GLU-GLN: 2594
    "GLN-GLY": 2087,  # GLY-GLN: 2091
    "GLN-HIS": 1224,  # HIS-GLN: 1230
    "GLN-ILE": 3105,  # ILE-GLN: 3115
    "GLN-LEU": 5433,  # LEU-GLN: 5448
    "GLN-LYS": 2473,  # LYS-GLN: 2427
    "GLN-MET": 1154,  # MET-GLN: 1154
    "GLN-PHE": 2668,  # PHE-GLN: 2678
    "GLN-PRO": 1802,  # PRO-GLN: 1805
    "GLN-SER": 2234,  # SER-GLN: 2242
    "GLN-THR": 2785,  # THR-GLN: 2791
    "GLN-TRP": 1435,  # TRP-GLN: 1437
    "GLN-TYR": 2737,  # TYR-GLN: 2745
    "GLN-VAL": 3359,  # VAL-GLN: 3364
    # GLUTAMATE
    "numGLU": 32564,
    "GLU-ALA": 2808,
    "GLU-ARG": 9678,
    "GLU-ASN": 2840,
    "GLU-ASP": 1891,
    "GLU-CYS": 596,
    "GLU-GLN": 2594,
    "GLU-GLU": 2475,
    "GLU-GLY": 2607,  # GLY-GLU: 2613
    "GLU-HIS": 2576,  # HIS-GLU: 2583
    "GLU-ILE": 3842,  # ILE-GLU: 3853
    "GLU-LEU": 6285,  # LEU-GLU: 6314
    "GLU-LYS": 7861,  # LYS-GLU: 7755
    "GLU-MET": 1324,  # MET-GLU: 1330
    "GLU-PHE": 3082,  # PHE-GLU: 3094
    "GLU-PRO": 2369,  # PRO-GLU: 2382
    "GLU-SER": 3260,  # SER-GLU: 3277
    "GLU-THR": 3623,  # THR-GLU: 3637
    "GLU-TRP": 1803,  # TRP-GLU: 1808
    "GLU-TYR": 3967,  # TYR-GLU: 3980
    "GLU-VAL": 4154,  # VAL-GLU: 4172
    # GLYCINE
    "numGLY": 35841,
    "GLY-ALA": 3678,
    "GLY-ARG": 3831,
    "GLY-ASN": 2755,
    "GLY-ASP": 3136,
    "GLY-CYS": 999,
    "GLY-GLN": 2091,
    "GLY-GLU": 2613,
    "GLY-GLY": 2816,
    "GLY-HIS": 1537,  # HIS-GLY: 1537
    "GLY-ILE": 4182,  # ILE-GLY: 4180
    "GLY-LEU": 6369,  # LEU-GLY: 6368
    "GLY-LYS": 2460,  # LYS-GLY: 2427
    "GLY-MET": 1624,  # MET-GLY: 1623
    "GLY-PHE": 3871,  # PHE-GLY: 3870
    "GLY-PRO": 2211,  # PRO-GLY: 2211
    "GLY-SER": 2859,  # SER-GLY: 2859
    "GLY-THR": 3481,  # THR-GLY: 3481
    "GLY-TRP": 1860,  # TRP-GLY: 1860
    "GLY-TYR": 3617,  # TYR-GLY: 3615
    "GLY-VAL": 4697,  # VAL-GLY: 4696
    # HISTIDINE
    "numHIS": 11167,
    "HIS-ALA": 1811,
    "HIS-ARG": 2096,
    "HIS-ASN": 1336,
    "HIS-ASP": 2384,
    "HIS-CYS": 639,
    "HIS-GLN": 1230,
    "HIS-GLU": 2583,
    "HIS-GLY": 1537,
    "HIS-HIS": 1750,
    "HIS-ILE": 2329,  # ILE-HIS: 2328
    "HIS-LEU": 3845,  # LEU-HIS: 3845
    "HIS-LYS": 1193,  # LYS-HIS: 1171
    "HIS-MET": 978,  # MET-HIS: 979
    "HIS-PHE": 2292,  # PHE-HIS: 2291
    "HIS-PRO": 1318,  # PRO-HIS: 1218
    "HIS-SER": 1785,  # SER-HIS: 1789
    "HIS-THR": 1989,  # THR-HIS: 1992
    "HIS-TRP": 1179,  # TRP-HIS: 1181
    "HIS-TYR": 2384,  # TYR-HIS: 2384
    "HIS-VAL": 2635,  # VAL-HIS: 2637
    # ISOLEUCINE
    "numILE": 27100,
    "ILE-ALA": 9518,
    "ILE-ARG": 4636,
    "ILE-ASN": 2744,
    "ILE-ASP": 2568,
    "ILE-CYS": 2072,
    "ILE-GLN": 3115,
    "ILE-GLU": 3853,
    "ILE-GLY": 4180,
    "ILE-HIS": 2328,
    "ILE-ILE": 18624,
    "ILE-LEU": 26652,  # LEU-ILE: 26651
    "ILE-LYS": 3816,  # LYS-ILE: 3758
    "ILE-MET": 4498,  # MET-ILE: 4493
    "ILE-PHE": 11586,  # PHE-ILE: 11585
    "ILE-PRO": 3726,  # PRO-ILE: 3727
    "ILE-SER": 3514,  # SER-ILE: 3513
    "ILE-THR": 5706,  # THR-ILE: 5704
    "ILE-TRP": 3672,  # TRP-ILE: 3672
    "ILE-TYR": 7310,  # TYR-ILE: 7311
    "ILE-VAL": 17789,  # VAL-ILE: 17790
    # LEUCINE
    "numLEU": 44478,
    "LEU-ALA": 15280,
    "LEU-ARG": 8423,
    "LEU-ASN": 3999,
    "LEU-ASP": 4056,
    "LEU-CYS": 3295,
    "LEU-GLN": 5448,
    "LEU-GLU": 6314,
    "LEU-GLY": 6368,
    "LEU-HIS": 3845,
    "LEU-ILE": 26651,
    "LEU-LEU": 47638,
    "LEU-LYS": 5696,  # LYS-LEU: 5592
    "LEU-MET": 6993,  # MET-LEU: 6993
    "LEU-PHE": 19454,  # PHE-LEU: 19454
    "LEU-PRO": 6433,  # PRO-LEU: 6438
    "LEU-SER": 5619,  # SER-LEU: 5614
    "LEU-THR": 8233,  # THR-LEU: 8236
    "LEU-TRP": 6487,  # TRP-LEU: 6485
    "LEU-TYR": 12030,  # TYR-LEU: 12039
    "LEU-VAL": 27217,  # VAL-LEU: 27218
    # LYSINE
    "numLYS": 27812,
    "LYS-ALA": 2140,
    "LYS-ARG": 2232,
    "LYS-ASN": 2395,
    "LYS-ASP": 5883,
    "LYS-CYS": 587,
    "LYS-GLN": 2428,
    "LYS-GLU": 7755,
    "LYS-GLY": 2427,
    "LYS-HIS": 1171,
    "LYS-ILE": 3758,
    "LYS-LEU": 5592,
    "LYS-LYS": 1616,
    "LYS-MET": 1120,  # MET-LYS: 1140
    "LYS-PHE": 2890,  # PHE-LYS: 2935
    "LYS-PRO": 1523,  # PRO-LYS: 1568
    "LYS-SER": 2365,  # SER-LYS: 2409
    "LYS-THR": 2834,  # THR-LYS: 2884
    "LYS-TRP": 1407,  # TRP-LYS: 1425
    "LYS-TYR": 3656,  # TYR-LYS: 3724
    "LYS-VAL": 3894,  # VAL-LYS: 3963
    # METHIONINE
    "numMET": 8855,
    "MET-ALA": 2855,
    "MET-ARG": 1628,
    "MET-ASN": 1053,
    "MET-ASP": 948,
    "MET-CYS": 641,
    "MET-GLN": 1154,
    "MET-GLU": 1330,
    "MET-GLY": 1623,
    "MET-HIS": 979,
    "MET-ILE": 4493,
    "MET-LEU": 6993,
    "MET-LYS": 1140,
    "MET-MET": 1973,
    "MET-PHE": 3859,  # PHE-MET: 3865
    "MET-PRO": 1485,  # PRO-MET: 1486
    "MET-SER": 1243,  # SER-MET: 1245
    "MET-THR": 1811,  # THR-MET: 1813
    "MET-TRP": 1411,  # TRP-MET: 1412
    "MET-TYR": 2613,  # TYR-MET: 2613
    "MET-VAL": 4576,  # VAL-MET: 4572
    # PHENYLALANINE
    "numPHE": 19578,
    "PHE-ALA": 7232,
    "PHE-ARG": 4135,
    "PHE-ASN": 2468,
    "PHE-ASP": 2291,
    "PHE-CYS": 1987,
    "PHE-GLN": 2678,
    "PHE-GLU": 3094,
    "PHE-GLY": 3870,
    "PHE-HIS": 2291,
    "PHE-ILE": 11585,
    "PHE-LEU": 19454,
    "PHE-LYS": 2935,
    "PHE-MET": 3865,
    "PHE-PHE": 11127,
    "PHE-PRO": 3755,  # PRO-PHE: 3756
    "PHE-SER": 3307,  # SER-PHE: 3306
    "PHE-THR": 4283,  # THR-PHE: 4283
    "PHE-TRP": 3676,  # TRP-PHE: 3676
    "PHE-TYR": 6660,  # TYR-PHE: 6658
    "PHE-VAL": 12109,  # VAL-PHE: 12106
    # PROLINE
    "numPRO": 22507,
    "PRO-ALA": 2975,
    "PRO-ARG": 2955,
    "PRO-ASN": 1695,
    "PRO-ASP": 1754,
    "PRO-CYS": 924,
    "PRO-GLN": 1805,
    "PRO-GLU": 2382,
    "PRO-GLY": 2211,
    "PRO-HIS": 1318,
    "PRO-ILE": 3727,
    "PRO-LEU": 6438,
    "PRO-LYS": 1568,
    "PRO-MET": 1486,
    "PRO-PHE": 3756,
    "PRO-PRO": 2142,
    "PRO-SER": 1990,  # SER-PRO: 1990
    "PRO-THR": 2471,  # THR-PRO: 2470
    "PRO-TRP": 2135,  # TRP-PRO: 2135
    "PRO-TYR": 4149,  # TYR-PRO: 4149
    "PRO-VAL": 4353,  # VAL-PRO: 4352
    # SERINE
    "numSER": 28561,
    "SER-ALA": 2931,
    "SER-ARG": 3239,
    "SER-ASN": 2487,
    "SER-ASP": 3279,
    "SER-CYS": 814,
    "SER-GLN": 2242,
    "SER-GLU": 3277,
    "SER-GLY": 2859,
    "SER-HIS": 1789,
    "SER-ILE": 3513,
    "SER-LEU": 5614,
    "SER-LYS": 2409,
    "SER-MET": 1245,
    "SER-PHE": 3306,
    "SER-PRO": 1990,
    "SER-SER": 2809,
    "SER-THR": 3132,  # THR-SER: 3132
    "SER-TRP": 1571,  # TRP-SER: 1571
    "SER-TYR": 2951,  # TYR-SER: 2952
    "SER-VAL": 4126,  # VAL-SER: 4130
    # THREONINE
    "numTHR": 26716,
    "THR-ALA": 4080,
    "THR-ARG": 3886,
    "THR-ASN": 2993,
    "THR-ASP": 3435,
    "THR-CYS": 967,
    "THR-GLN": 2791,
    "THR-GLU": 3637,
    "THR-GLY": 3481,
    "THR-HIS": 1992,
    "THR-ILE": 5704,
    "THR-LEU": 8236,
    "THR-LYS": 2884,
    "THR-MET": 1813,
    "THR-PHE": 4283,
    "THR-PRO": 2470,
    "THR-SER": 3132,
    "THR-THR": 4262,
    "THR-TRP": 1773,  # TRP-THR: 1773
    "THR-TYR": 3567,  # TYR-THR: 3568
    "THR-VAL": 6411,  # VAL-THR: 6412
    # TRYPTOPHAN
    "numTRP": 7118,
    "TRP-ALA": 2680,
    "TRP-ARG": 2303,
    "TRP-ASN": 1317,
    "TRP-ASP": 1302,
    "TRP-CYS": 709,
    "TRP-GLN": 1437,
    "TRP-GLU": 1808,
    "TRP-GLY": 1860,
    "TRP-HIS": 1181,
    "TRP-ILE": 3672,
    "TRP-LEU": 6485,
    "TRP-LYS": 1425,
    "TRP-MET": 1412,
    "TRP-PHE": 3676,
    "TRP-PRO": 2135,
    "TRP-SER": 1571,
    "TRP-THR": 1773,
    "TRP-TRP": 1800,
    "TRP-TYR": 2737,  # TYR-TRP: 2738
    "TRP-VAL": 3870,  # VAL-TRP: 3869
    # TYROSINE
    "numTYR": 17058,
    "TYR-ALA": 5144,
    "TYR-ARG": 4613,
    "TYR-ASN": 2736,
    "TYR-ASP": 3531,
    "TYR-CYS": 1287,
    "TYR-GLN": 2745,
    "TYR-GLU": 3980,
    "TYR-GLY": 3615,
    "TYR-HIS": 2384,
    "TYR-ILE": 7311,
    "TYR-LEU": 12039,
    "TYR-LYS": 3724,
    "TYR-MET": 2613,
    "TYR-PHE": 6658,
    "TYR-PRO": 4149,
    "TYR-SER": 2952,
    "TYR-THR": 3568,
    "TYR-TRP": 2738,
    "TYR-TYR": 5179,
    "TYR-VAL": 7633,  # VAL-TYR: 7632
    # VALINE
    "numVAL": 34567,
    "VAL-ALA": 10941,
    "VAL-ARG": 5190,
    "VAL-ASN": 3003,
    "VAL-ASP": 3004,
    "VAL-CYS": 2340,
    "VAL-GLN": 3364,
    "VAL-GLU": 4172,
    "VAL-GLY": 4696,
    "VAL-HIS": 2637,
    "VAL-ILE": 17790,
    "VAL-LEU": 27218,
    "VAL-LYS": 3963,
    "VAL-MET": 4572,
    "VAL-PHE": 12106,
    "VAL-PRO": 4352,
    "VAL-SER": 4130,
    "VAL-THR": 6412,
    "VAL-TRP": 3869,
    "VAL-TYR": 7632,
    "VAL-VAL": 19723,
}

"""
Based on:
Singh, J., and Thornton, J. M. (1992) Atlas of Protein Side Chain Interactions, 
IRL Press at Oxford University Press, Oxford.
"""
book_atlas = {
    # ALANINE
    "ALA:ALA": 1.2,
    "ALA:ARG": 0.6,
    "ALA:ASN": 0.9,
    "ALA:ASP": 0.8,
    "ALA:CYS": 0.5,
    "ALA:GLN": 1.0,
    "ALA:GLU": 0.6,
    "ALA:GLY": 1.1,
    "ALA:HIS": 0.8,
    "ALA:ILE": 1.3,
    "ALA:LEU": 1.2,
    "ALA:LYS": 0.8,
    "ALA:MET": 0.8,
    "ALA:PHE": 1.1,
    "ALA:PRO": 0.8,
    "ALA:SER": 0.9,
    "ALA:THR": 1.0,
    "ALA:TRP": 1.1,
    "ALA:TYR": 1.1,
    "ALA:VAL": 1.3,
    # ARGININE
    "ARG:ARG": 0.7,
    "ARG:ASN": 0.9,
    "ARG:ASP": 2.4,
    "ARG:CYS": 0.6,
    "ARG:GLN": 1.5,
    "ARG:GLU": 2.2,
    "ARG:GLY": 0.8,
    "ARG:HIS": 0.9,
    "ARG:ILE": 0.6,
    "ARG:LEU": 0.6,
    "ARG:LYS": 0.4,
    "ARG:MET": 0.6,
    "ARG:PHE": 1.0,
    "ARG:PRO": 0.9,
    "ARG:SER": 1.0,
    "ARG:THR": 0.9,
    "ARG:TRP": 1.1,
    "ARG:TYR": 1.1,
    "ARG:VAL": 0.5,
    # ASPARAGINE
    "ASN:ASN": 1.6,
    "ASN:ASP": 1.6,
    "ASN:CYS": 0.7,
    "ASN:GLN": 1.6,
    "ASN:GLU": 1.4,
    "ASN:GLY": 1.4,
    "ASN:HIS": 1.2,
    "ASN:ILE": 0.4,
    "ASN:LEU": 0.6,
    "ASN:LYS": 1.4,
    "ASN:MET": 0.6,
    "ASN:PHE": 0.5,
    "ASN:PRO": 0.8,
    "ASN:SER": 1.3,
    "ASN:THR": 1.3,
    "ASN:TRP": 0.9,
    "ASN:TYR": 1.1,
    "ASN:VAL": 0.7,
    # ASPARTATE
    "ASP:ASP": 0.9,
    "ASP:CYS": 0.6,
    "ASP:GLN": 0.9,
    "ASP:GLU": 0.9,
    "ASP:GLY": 1.4,
    "ASP:HIS": 1.4,
    "ASP:ILE": 0.5,
    "ASP:LEU": 0.4,
    "ASP:LYS": 2.7,
    "ASP:MET": 0.4,
    "ASP:PHE": 0.5,
    "ASP:PRO": 0.6,
    "ASP:SER": 1.8,
    "ASP:THR": 1.5,
    "ASP:TRP": 0.6,
    "ASP:TYR": 0.9,
    "ASP:VAL": 0.4,
    # CYSTEINE
    "CYS:CYS": 13.6,
    "CYS:GLN": 0.8,
    "CYS:GLU": 0.3,
    "CYS:GLY": 0.9,
    "CYS:HIS": 0.8,
    "CYS:ILE": 1.0,
    "CYS:LEU": 0.8,
    "CYS:LYS": 0.3,
    "CYS:MET": 1.3,
    "CYS:PHE": 0.8,
    "CYS:PRO": 1.0,
    "CYS:SER": 0.8,
    "CYS:THR": 0.6,
    "CYS:TRP": 0.7,
    "CYS:TYR": 0.9,
    "CYS:VAL": 0.8,
    # GLUTAMINE
    "GLN:GLN": 1.1,
    "GLN:GLU": 0.9,
    "GLN:GLY": 1.3,
    "GLN:HIS": 0.7,
    "GLN:ILE": 0.7,
    "GLN:LEU": 0.8,
    "GLN:LYS": 1.0,
    "GLN:MET": 0.9,
    "GLN:PHE": 0.7,
    "GLN:PRO": 1.4,
    "GLN:SER": 1.1,
    "GLN:THR": 1.3,
    "GLN:TRP": 0.9,
    "GLN:TYR": 1.1,
    "GLN:VAL": 0.9,
    # GLUTAMATE
    "GLU:GLU": 1.0,
    "GLU:GLY": 0.6,
    "GLU:HIS": 1.4,
    "GLU:ILE": 0.6,
    "GLU:LEU": 0.6,
    "GLU:LYS": 2.8,
    "GLU:MET": 0.7,
    "GLU:PHE": 0.6,
    "GLU:PRO": 0.7,
    "GLU:SER": 1.5,
    "GLU:THR": 1.3,
    "GLU:TRP": 0.7,
    "GLU:TYR": 0.8,
    "GLU:VAL": 0.7,
    # GLYCINE
    "GLY:GLY": 1.5,
    "GLY:HIS": 0.9,
    "GLY:ILE": 0.8,
    "GLY:LEU": 0.8,
    "GLY:LYS": 0.8,
    "GLY:MET": 0.7,
    "GLY:PHE": 0.6,
    "GLY:PRO": 0.9,
    "GLY:SER": 1.2,
    "GLY:THR": 1.5,
    "GLY:TRP": 1.0,
    "GLY:TYR": 1.3,
    "GLY:VAL": 1.0,
    # HISTIDINE
    "HIS:HIS": 2.3,
    "HIS:ILE": 0.6,
    "HIS:LEU": 0.7,
    "HIS:LYS": 0.7,
    "HIS:MET": 0.8,
    "HIS:PHE": 1.1,
    "HIS:PRO": 0.8,
    "HIS:SER": 0.9,
    "HIS:THR": 1.1,
    "HIS:TRP": 1.1,
    "HIS:TYR": 1.0,
    "HIS:VAL": 0.7,
    # ISOLEUCINE
    "ILE:ILE": 1.7,
    "ILE:LEU": 1.5,
    "ILE:LYS": 0.6,
    "ILE:MET": 1.2,
    "ILE:PHE": 1.2,
    "ILE:PRO": 0.6,
    "ILE:SER": 0.6,
    "ILE:THR": 0.8,
    "ILE:TRP": 1.3,
    "ILE:TYR": 1.0,
    "ILE:VAL": 1.5,
    # LEUCINE
    "LEU:LEU": 1.5,
    "LEU:LYS": 0.5,
    "LEU:MET": 1.4,
    "LEU:PHE": 1.4,
    "LEU:PRO": 0.7,
    "LEU:SER": 0.8,
    "LEU:THR": 0.7,
    "LEU:TRP": 1.3,
    "LEU:TYR": 1.0,
    "LEU:VAL": 1.5,
    # LYSINE
    "LYS:LYS": 0.4,
    "LYS:MET": 0.7,
    "LYS:PHE": 0.8,
    "LYS:PRO": 0.7,
    "LYS:SER": 0.9,
    "LYS:THR": 1.0,
    "LYS:TRP": 0.8,
    "LYS:TYR": 1.3,
    "LYS:VAL": 0.6,
    # METHIONINE
    "MET:MET": 2.4,
    "MET:PHE": 1.5,
    "MET:PRO": 1.0,
    "MET:SER": 0.8,
    "MET:THR": 0.7,
    "MET:TRP": 1.4,
    "MET:TYR": 1.0,
    "MET:VAL": 1.0,
    # PHENYLALANINE
    "PHE:PHE": 1.5,
    "PHE:PRO": 1.0,
    "PHE:SER": 0.7,
    "PHE:THR": 0.7,
    "PHE:TRP": 1.2,
    "PHE:TYR": 0.8,
    "PHE:VAL": 1.3,
    # PROLINE
    "PRO:PRO": 1.1,
    "PRO:SER": 0.9,
    "PRO:THR": 0.9,
    "PRO:TRP": 1.7,
    "PRO:TYR": 1.8,
    "PRO:VAL": 1.1,
    # SERINE
    "SER:SER": 1.6,
    "SER:THR": 1.6,
    "SER:TRP": 0.5,
    "SER:TYR": 0.9,
    "SER:VAL": 0.9,
    # THREONINE
    "THR:THR": 1.6,
    "THR:TRP": 0.7,
    "THR:TYR": 0.7,
    "THR:VAL": 0.8,
    # TRYPTOPHAN
    "TRP:TRP": 0.9,
    "TRP:TYR": 0.9,
    "TRP:VAL": 1.1,
    # TYROSINE
    "TYR:TYR": 1.1,
    "TYR:VAL": 0.9,
    # VALINE
    "VAL:VAL": 1.6,
}

# per residue sum of normalized contact propensities >= 1.0
book_summed_propensities = {'ALA': 11.4,
                            'ARG': 10.3,
                            'ASN': 13.9,
                            'ASP': 12.8,
                            'CYS': 16.9,
                            'GLN': 12.4,
                            'GLU': 11.6,
                            'GLY': 12.7,
                            'HIS': 10.6,
                            'ILE': 11.7,
                            'LEU': 10.8,
                            'LYS': 10.2,
                            'MET': 12.2,
                            'PHE': 12.3,
                            'PRO': 10.1,
                            'SER': 11.1,
                            'THR': 13.2,
                            'TRP': 12.3,
                            'TYR': 13.9,
                            'VAL': 11.4}


def scii_from_seq(seq: str) -> float:
    """
    Calculates the side chain interaction index for the given amino acid sequence.
    Based on:
    Gehenn, K., Pipkorn, R., & Reed, J. (2003). Successful Design and Synthesis of a Polarity-Triggered β → α
    Conformational Switch Using the Side Chain Interaction Index (SCII) as a Measure of Local Stuctural Stability.
    In Biochemistry (Vol. 43, Issue 3, pp. 607–612). American Chemical Society (ACS). https://doi.org/10.1021/bi0301744

    :param seq: A string of amino acids in single letter format. Should not exceed 20 amino acids.
    :return: The side chain interaction index for the given sequence.
    """
    if len(seq) > 20:
        logging.warning(f"The passed sequence {seq} is {len(seq)} amino acids long. This is likely too long to "
                        f"calculate a meaningful SCII value for!")
    aas = [PDBtool.one_to_three(a, "") for a in seq]
    buf = {aa: {"partners": [], "scores": []} for aa in aas}
    for n, aa in enumerate(aas):
        if buf[aa].get("visited"):  # Already have calculated score for this amino acid
            continue
        buf[aa]["visited"] = True
        for i in range(len(aas)):
            if i == n or aas[i] in buf[aa]["partners"]:  # Don't count self, only count unique amino acids
                continue
            prop = 0.
            if _ := book_atlas.get(f"{aa}:{aas[i]}"):
                prop = _
            else:
                prop = book_atlas[f"{aas[i]}:{aa}"]
            if prop < 1.0:  # Only include favorable interactions
                continue
            buf[aa]["partners"].append(aas[i])
            buf[aa]["scores"].append(prop)
        buf[aa]["index_number"] = sum(buf[aa]["scores"]) / book_summed_propensities[aa]
    # print(pprint.pformat(buf))
    inums = [buf[aa]["index_number"] for aa in aas]
    return statistics.mean(inums)


def scii_for_res(res: str, partners: Union[str, List[str]]) -> float:
    """
    Calculates the side chain interaction index for the given amino acid and its given partners.
    Based on:
    Gehenn, K., Pipkorn, R., & Reed, J. (2003). Successful Design and Synthesis of a Polarity-Triggered β → α
    Conformational Switch Using the Side Chain Interaction Index (SCII) as a Measure of Local Stuctural Stability.
    In Biochemistry (Vol. 43, Issue 3, pp. 607–612). American Chemical Society (ACS). https://doi.org/10.1021/bi0301744

    :param res: An amino acid in abbreviated form (one letter) to calculate the index number for
    :param partners: A list of partner amino acids to consider when calculating the index number. Should be either a
        string of one letter abbreviations or a list of abbreviations.
    :return: The index number / side chain interaction index for the given residue and its partners.
    """
    r = PDBtool.one_to_three(res, "")
    if r == "":
        logging.error(f"Could not recognize residue {res}!")
        raise ValueError(f"Could not recognize residue {res}!")
    aas = [PDBtool.one_to_three(aa, "") for aa in partners]
    buf = {aa: {"partners": [], "scores": []} for aa in aas}
    for aa in aas:
        if buf[aa].get("visited"):  # Already have calculated score for this amino acid
            continue
        buf[aa]["visited"] = True
        for i in range(len(aas)):
            if aas[i] in buf[aa]["partners"]:  # Don't count self, only count unique amino acids
                continue
            prop = 0.
            if _ := book_atlas.get(f"{aa}:{aas[i]}"):
                prop = _
            else:
                prop = book_atlas[f"{aas[i]}:{aa}"]
            if prop < 1.0:  # Only include favorable interactions
                continue
            buf[aa]["partners"].append(aas[i])
            buf[aa]["scores"].append(prop)
        buf[aa]["index_number"] = sum(buf[aa]["scores"]) / book_summed_propensities[aa]
    # print(pprint.pformat(buf))
    return statistics.mean([buf[aa]["index_number"] for aa in aas])


def scii_for_pdb(path: Union[str, Path], chain: str = None, radius: float = 7.0, segment_size: int = False,
                 use_centroid: bool = False) -> Union[float, dict]:
    """
    Calculates the side chain interaction index for a given pdb by calculating the average SCII for each residue and
    every partner residue within a certain radius around it. The distance is calculated between the centroid of all
    atoms in the residue.
    Based on:
    Gehenn, K., Pipkorn, R., & Reed, J. (2003). Successful Design and Synthesis of a Polarity-Triggered β → α
    Conformational Switch Using the Side Chain Interaction Index (SCII) as a Measure of Local Stuctural Stability.
    In Biochemistry (Vol. 43, Issue 3, pp. 607–612). American Chemical Society (ACS). https://doi.org/10.1021/bi0301744

    :param path: The path to the PDB file to calculate the SCII for.
    :param chain: Which chain to analyze (optional). If not set, uses the first chain encountered in the PDB.
    :param radius: How many Å away the centroid of a residue may be to still be included in the neighborhood to
        calculate the SCII for.
    :param segment_size: What size of segment to use, if per-segment SCII is to be calculated. Default is false, set to
        any integer to enable segment-wise SCII
    :param use_centroid: Whether to use centroids for the inter-residue distance calculations, or the distance between
        the closest two atoms of each residue
    :return: The SCII for the molecule in the PDB.
    """
    try:
        path = Path(path).resolve(strict=True)
        path = PDBtool.renumber_ascending(os.path.normpath(path))
    except (FileNotFoundError, RuntimeError) as e:
        logging.error(f"Couldn't resolve the path to the PDB file. Stacktrace: {e}")
        raise e
    if chain is None:
        chain = PDBtool.get_chains(path)
        if chain is None:
            logging.error("Couldn't autodiscover which chain to use!")
            raise ValueError("Couldn't autodiscover which chain to use!")
        chain = chain[0]
    if use_centroid:
        aas = PDBtool.get_amino_acids_on_chain(path, chain, with_id=True)
        res_coord = {aa: None for aa in aas}
        for aa, res_id in aas:
            res_coord[(aa, res_id)] = PDBtool.get_centroid(path, res_id)
        # print(pprint.pformat(res_coord))
    else:
        _ = PDBtool.get_atoms_on_chain(path, chain)
        residues = {}
        for atom in _:
            if atom["comp_num"] not in residues:
                residues[atom["comp_num"]] = []
            residues[atom["comp_num"]].append(atom)
        # print(pprint.pformat(residues))
    sciis = {}
    if use_centroid:
        for residue1, centroid1 in res_coord.items():
            partners = []
            for residue2, centroid2 in res_coord.items():
                if residue1[1] == residue2[1]:
                    continue
                if math.sqrt((centroid1[0] - centroid2[0]) ** 2 + (centroid1[1] - centroid2[1]) ** 2 + (
                        centroid1[2] - centroid2[2]) ** 2) <= radius:
                    partners.append(residue2[0])
            # print(f"partners for res {residue1}: " + pprint.pformat(partners))
            sciis[residue1[1]] = scii_for_res(residue1[0], partners)
    else:
        # Maybe faster with KD-Tree?
        def _dist(coord1: List[float], coord2: List[float]) -> float:
            return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2)

        for residue1, atoms1 in residues.items():
            partners = []
            for residue2, atoms2 in residues.items():
                if residue1 == residue2:
                    continue
                _ = False
                for a1 in atoms1:
                    if _:
                        _ = False
                        break
                    for a2 in atoms2:
                        dist = _dist([a1["X"], a1["Y"], a1["Z"]], [a2["X"], a2["Y"], a2["Z"]])
                        if dist >= 22:
                            # tryptophan is ~10 angstrom along longest axis, if an atom is more than 22 angstrom away.
                            # it can't be within the set radius, even if both residues are tryptophans and the sampled
                            # atoms were the worst possible choice. Has an added safety margin of 2 angstrom
                            _ = True
                            break
                        if dist <= radius:
                            _ = True
                            partners.append(a2["atom_comp_id"])
                            break
            sciis[residue1] = scii_for_res(PDBtool.three_to_one(atoms1[0]["atom_comp_id"]), partners)
    if not segment_size:
        return statistics.mean(sciis.values())
    else:
        buf = [[] for _ in range(len(sciis) // segment_size + 1)]
        for res_id, scii in sciis.items():
            buf[res_id // segment_size].append((res_id, scii))
        if len(buf[-1]) < segment_size // 2:
            buf[-2] = buf[-2] + buf[-1]
            buf.pop()
        segment_sciis = {}
        for segment in buf:
            segment_sciis[f"{segment[0][0]}-{segment[-1][0]}"] = statistics.mean([v[1] for v in segment])
        return segment_sciis
