# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                 Institute for Mathematics and Computer-based Simulation
#                 University of the Bundeswehr Munich
#                 https://www.unibw.de/imcs-en
#
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

fiber_ref = {
    "frames": [1, 2, 3],
    "connectivity": [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]],
    "tangents": [
        [
            [1.6666666666666727e-01, 7.4014868308343765e-17, 2.9605947323337506e-16],
            [1.6666666666666755e-01, 7.4014868308343765e-17, 2.9605947323337506e-16],
            [1.6666666666666755e-01, 7.4014868308343765e-17, 2.9605947323337506e-16],
            [1.6666666666666755e-01, 7.4014868308343765e-17, 2.9605947323337506e-16],
            [1.6666666666666785e-01, 7.4014868308343765e-17, 2.9605947323337506e-16],
            [1.6666666666666666e-01, -7.4014868308343765e-17, -2.9605947323337506e-16],
            [-1.1564823173178713e-12, 3.0000000000000021e-01, 1.1572964808692632e-12],
            [-1.1564823173178713e-12, 3.0000000000000021e-01, 1.1572964808692632e-12],
            [-1.1564823173178713e-12, 3.0000000000000154e-01, 1.1572964808692632e-12],
            [-1.1564823173178713e-12, 3.0000000000000027e-01, 1.1572964808692632e-12],
            [-1.1564823173178713e-12, 3.0000000000000043e-01, 1.1572964808692632e-12],
            [-1.1566303470544881e-12, 2.9999999999999921e-01, 1.1567043619227964e-12],
        ],
        [
            [1.6273909132413109e-01, 3.1199812227798291e-02, 1.7929351197317239e-02],
            [1.6212850911143031e-01, 3.3215613382746444e-02, 1.9748124639686965e-02],
            [1.6154438995719422e-01, 3.5037566875858052e-02, 2.1323995252048904e-02],
            [1.6099888393687789e-01, 3.6670196057711268e-02, 2.2662805532595327e-02],
            [1.6050233419613327e-01, 3.8117569552254350e-02, 2.3769921938070038e-02],
            [1.6006341490375275e-01, 3.9383233826328391e-02, 2.4650167105336518e-02],
            [-5.1777607021698924e-02, 2.9282361101573989e-01, -3.9642010116601810e-02],
            [-4.5252402072889386e-02, 2.9453341908381842e-01, -3.4644415267251670e-02],
            [-3.6826604709168285e-02, 2.9639641939042632e-01, -2.8120340767823659e-02],
            [-2.6473203074376690e-02, 2.9814737825513865e-01, -2.0122330863413158e-02],
            [-1.4189645869202830e-02, 2.9946803994657234e-01, -1.0720314595960900e-02],
            [-9.5492977419953107e-06, 2.9999522535170414e-01, -9.5492954284755651e-06],
        ],
        [
            [1.5726141409634017e-01, 4.9745513804322815e-02, 2.3967893847206323e-02],
            [1.5607257201520541e-01, 5.2178742018658873e-02, 2.6438589317587063e-02],
            [1.5498495807745569e-01, 5.4309098924245659e-02, 2.8469018403483853e-02],
            [1.5402211433184290e-01, 5.6151900286578162e-02, 3.0074561124326864e-02],
            [1.5320255241907180e-01, 5.7720130834906493e-02, 3.1268384577404364e-02],
            [1.5254038163770986e-01, 5.9024327199983483e-02, 3.2061245105525403e-02],
            [-8.0087818312943204e-02, 2.8257305236335956e-01, -6.1120467500650401e-02],
            [-7.2139907619065177e-02, 2.8596819312396465e-01, -5.4908934505147812e-02],
            [-6.0279586754371305e-02, 2.9032207440837604e-01, -4.5562497617658661e-02],
            [-4.4333297933022375e-02, 2.9483874179409703e-01, -3.3174150506372214e-02],
            [-2.4217151334562887e-02, 2.9847645291126107e-01, -1.7900851555000024e-02],
            [-1.9098594327434288e-05, 2.9999045070341318e-01, -1.9098592015875937e-05],
        ],
    ],
    "fiber_radii": [0.1, 0.1],
    "coordinates": [
        [
            [0.5, 0.5, 2.0],
            [1.0, 0.5, 2.0],
            [1.5, 0.5, 2.0],
            [2.0, 0.5, 2.0],
            [2.5, 0.5, 2.0],
            [3.0, 0.5, 2.0],
            [0.5, 0.5, 2.0],
            [0.5, 1.4, 2.0],
            [0.5, 2.3, 2.0],
            [0.5, 3.2, 2.0],
            [0.5, 4.1, 2.0],
            [0.5, 5.0, 2.0],
        ],
        [
            [0.948286838466663, 0.5447475509007359, 2.342076817894439],
            [1.435583292108013, 0.641419746127858, 2.398654518353514],
            [1.921084418679117, 0.7438473916302395, 2.460322670974626],
            [2.40488830275625, 0.851455833639283, 2.526361440180104],
            [2.887126732139883, 0.9636833248973671, 2.596067822956782],
            [3.367959988263933, 1.079979556021063, 2.668754114771083],
            [0.948286838466663, 0.5447475509007359, 2.342076817894439],
            [0.8022713320361106, 1.425716080814451, 2.23025989897235],
            [0.6786733606603241, 2.312103054699511, 2.13573703069914],
            [0.583239960963886, 3.203984202470153, 2.063012228932619],
            [0.5217653035337566, 4.100558428259391, 2.016407726903026],
            [0.5, 5.0, 2.0],
        ],
        [
            [1.228097757969286, 0.6166652990143313, 2.549960586810679],
            [1.698076808895335, 0.7696294826804005, 2.625682450542643],
            [2.164634605638421, 0.9294349353296336, 2.708151899466971],
            [2.62811149182423, 1.095196581351745, 2.796071726407185],
            [3.08891071850994, 1.266071830554488, 2.88818760466178],
            [3.547484327176075, 1.441253399164251, 2.983281132442972],
            [1.228097757969286, 0.6166652990143313, 2.549960586810679],
            [0.9988027227855376, 1.469158244893169, 2.375123247716616],
            [0.7991710579441209, 2.333441298989852, 2.223642061066096],
            [0.6412135106078395, 3.211261125172633, 2.104791761444802],
            [0.5373459605299356, 4.10160116972303, 2.027485364727713],
            [0.5, 5.0, 2.0],
        ],
    ],
}
