%% INPUT DATA SAMPLE

input.xpoints = [
    -7.0, 0.0, 0.0;
    -7.0, 0.0, -0.2;
    -7.0, 1.5, 0.0;
    -7.0, 1.5, -0.2;
    -5.0, 0.0, 0.0;
    -3.0, 0.0, 0.0;
    -1.0, 0.0, 0.0;
    1.0, 0.0, 0.0;
    3.0, 0.0, 0.0;
    5.0, 0.0, 0.0;
    -5.0, 1.6, 0.0;
    -3.0, 1.6, 0.0;
    -1.0, 1.6, 0.0;
    1.0, 1.6, 0.0;
    3.0, 1.6, 0.0;
    5.0, 1.6, 0.0;
    -5.0, 0.0, -0.2;
    -3.0, 0.0, -0.2;
    -1.0, 0.0, -0.2;
    1.0, 0.0, -0.2;
    3.0, 0.0, -0.2;
    5.0, 0.0, -0.2;
    -5.0, 1.6, -0.2;
    -3.0, 1.6, -0.2;
    -1.0, 1.6, -0.2;
    1.0, 1.6, -0.2;
    3.0, 1.6, -0.2;
    5.0, 1.6, -0.2;
    -0.5, 2.6, -0.6;
    0.5, 2.6, -0.6;
    -0.5, 2.6, -1.1;
    0.5, 2.6, -1.1;
    0.0, 3.1, -0.8;
    -1.0, 1.5, -1.6;
    1.0, 1.5, -1.6;
    1.0, 0.0, -1.2;
    -1.0, 0.0, -1.2;
    1.5, 1.0, -2.2;
    -1.5, 1.0, -2.2;
    0.8, -2.0, -0.2;
    -0.8, -2.0, -0.2;
    0.8, -2.0, -1.0;
    -0.8, -2.0, -1.0;
    0.4, -3.8, -0.2;
    -0.4, -3.8, -0.2;
    0.4, -3.8, -0.6;
    -0.4, -3.8, -0.6;
    0.0, -5.0, -0.2;
    0.0, -5.0, -0.4;
    7.0, 0.0, 0.0;
    7.0, 0.0, -0.2;
    7.0, 1.5, 0.0;
    7.0, 1.5, -0.2;
];


input.T = [
     1     3
     1     5
     1     2
     2    17
     2     4
     3    11
     4    23
     3     4
     1    11
     2    23
     5     6
     6     7
     7     8
     8     9
     9    10
    10    16
    11    12
    12    13
    13    14
    14    15
    15    16
     5    11
     6    12
     7    13
     8    14
     9    15
     5    17
     6    18
     7    19
     8    20
     9    21
    10    22
    11    23
    12    24
    13    25
    14    26
    15    27
    16    28
    17    18
    18    19
    19    20
    20    21
    21    22
    22    28
    17    23
    23    24
    24    25
    25    26
    26    27
    27    28
    18    24
    19    25
    20    26
    21    27
     6    11
    18    23
    18    25
     6    13
     7    14
    20    25
     9    14
    21    26
     9    16
    21    28
    25    29
    26    30
    29    30
    31    32
    29    31
    30    32
    29    33
    31    33
    30    33
    32    33
    33    31
    31    34
    32    35
    34    37
    35    36
    26    35
    25    34
    37    19
    20    36
    37    39
    39    34
    35    38
    36    38
    36    42
    20    40
    40    42
    19    41
    37    43
    41    43
    40    44
    41    45
    45    47
    44    46
    42    46
    43    47
    47    49
    45    48
    44    48
    46    49
    48    49
    10    50
    22    51
    50    51
    50    52
    51    53
    16    52
    28    53
    52    53
    16    50
    28    51
    38    34
    39    35
    42    43
    36    37
    46    47
     2     3
     2     5
     5    18
    18     7
     8    21
    21    10
    10    51
    11    24
    24    13
    14    27
    27    16
     4    11
    16    53
    52    51
    40    41
    44    45
    36    40
    37    41
    36    43
    43    46
    19    40
    40    45
    40    46
    41    47
    47    48
    46    48
    34    35
    30    33
    29    32
    33    32
    11    17
    12    18
    13    19
    14    20
    15    21
    16    22
     7    20
    14    25
    20    37
    41    42
    26    34
    44    47
    35    37
    34    36
    35    20
    34    19
    26    32
    25    31
     1    23
    23     6
     6    25
    50    28
    28     9
     9    26
    29    32
    26    37
    36    25
    37    40
    41    46
    29    26
    32    34
    20    41
    36    41
    ];
    
% Element properties

input.E = ones(size(input.T,1),1) * 210e9 ; % 210e9
input.A = ones(size(input.T,1),1) * 3*10^-4 ; % m^2 
input.rho = ones(size(input.T,1),1) * 7850 ; % kg/m^3

% BC
input.fixnodes = [29, 1, 0;...
            29, 2, 0;...
            29, 3, 0;...  %Rueda derecha
            30, 1, 0;...
            30, 2, 0;...
            30, 3, 0;...   %Rueda derecha
            45, 1, 0;...
            48-3, 2, 0;...
            48-3, 3, 0;... %Nodo superior cola
];
   
%% MASAS ADICIONALES
input.g = 9.81;

masa_motor = 65;
masa_piloto = 80;
masa_pasajero = 80;

input.mass = [
    33, masa_motor ;
    39-3, masa_piloto ;
    40-3, masa_pasajero ;
    ];

clear masa_motor masa_piloto masa_pasajero;



%% FORCES
Lift = 940*input.g; % L = Weight
lift_parcial = Lift / 14;
    
input.Fext = [1, 3, lift_parcial;...
              5, 3, lift_parcial;...
              6, 3, lift_parcial;...
              7, 3, lift_parcial; 
              8, 3, lift_parcial;...
              9, 3, lift_parcial;...
              10, 3, lift_parcial;...
              50, 3, lift_parcial;...
              
              3, 3, lift_parcial;...
              11, 3, lift_parcial;...
              12, 3, lift_parcial;...
              13, 3, lift_parcial;...
              14, 3, lift_parcial;...
              15, 3, lift_parcial;...
              16, 3, lift_parcial;...
              52, 3, lift_parcial;...
              ];
          
          
%% TERMIC
input.alpha = 12*10^-6;

input.Delta_T = ones(size(input.T,1),1) * -40;

              
