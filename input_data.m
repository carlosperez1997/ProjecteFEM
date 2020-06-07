%% INPUT DATA SAMPLE

input.xpoints = [
  0, 0, 0;           %(1)
  0, 0, -0.2;        %(2)
  0, 1.6, 0;         %(3)
  0, 1.6, -0.2;      %(4)
  2, 0, 0;           %(5)
  4, 0, 0;           %(6)
  6, 0, 0;           %(7)
  8, 0, 0;           %(8)
  10, 0, 0;          %(9)
  12, 0, 0;          %(10)
  2, 1.65, 0;        %(11)
  4, 1.65, 0;        %(12)
  6, 1.65, 0;        %(13)
  8, 1.65, 0;        %(14)
  10, 1.65, 0;        %(15)
  12, 1.65, 0;       %(16)
    2, 0, -0.2;         %(17)
  4, 0, -0.2;           %(18)
  6, 0, -0.2;           %(19)
  8, 0, -0.2;           %(20)
  10, 0, -0.2;           %(21)
  12, 0, -0.2;          %(22)
  2, 1.65, -0.2;        %(23)
  4, 1.65, -0.2;        %(24)
  6, 1.65, -0.2;        %(25)
  8, 1.65, -0.2;        %(26)
  10, 1.65, -0.2;        %(27)
  12, 1.65, -0.2;       %(28)
  6.5, 2.6, -0.6;       %(29)
  7.5, 2.6, -0.6;       %(30)
  6.5, 2.6, -1.1;       %(31)
  7.5, 2.6, -1.1;       %(32)
  7, 3.1, -0.8;       %(33)
  %7.5, 3.1, -0.6;       %(34)
  %6.5, 3.1, -1.1;       %(35)
  %7.5, 3.1, -1.1;       %(36)
  6, 1.6, -1.6;         %(37)
  8, 1.6, -1.6;         %(38)
  8, 0, -1.2;           %(39)
  6, 0, -1.2;           %(40)
  8.5, 1, -2.2;         %(41)
  5.5, 1, -2.2;         %(42)
  7.8,-2, -0.2;       %(43)
  6.2,-2, -0.2;       %(44)
  7.8,-2, -1;         %(45)
  6.2,-2, -1;         %(46)
  7.4,-3.8, -0.2;     %(47)
  6.6,-3.8, -0.2;     %(48)
  7.4,-3.8, -0.6;     %(49)
  6.6,-3.8, -0.6;     %(50)
  7,-5, -0.2;         %(51)
  7,-5, -0.4;         %(52)
  14,0, 0;            %(53)
  14,0, -0.2;         %(54)
  14,1.6, 0;          %(55)
  14,1.6, -0.2;       %(56)
  
  %7-0.8, -4.5, -0.2;          %Añado (57)
  %7+0.8, -4.5 , -0.2;       % Añado (58)
];

for i = 1:size(input.xpoints,1)
    input.xpoints(i,1) = input.xpoints(i,1) - 7;
end

input.T = [
     1     3
     1     5
     1     2
     2     17
     2     4
     3     11
     4     23
     3     4
     1     11
     2     23
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
    31 33 %31    35
    30 33 %30    34
    32 33 %32    36
    %34    36
    33  31
    %33    35
    %35    36
    %33    34
    31    37
    32    38
    37    40
    38    39
    26    38
    25    37
    40    19
    20    39
    40    42
    42    37
    38    41
    39    41
    39    45
    20    43
    43    45
    19    44
    40    46
    44    46
    43    47
    44    48
    48    50
    47    49
    45    49
    46    50
    50    52
    48    51
    47    51
    49    52
    51    52
    10    53
    22    54
    53    54
    53    55
    54    56
    16    55
    28    56
    55    56
    16    53
    28    54
    41    37
    42    38
    45    46
    39    40
    49    50
    2     3
    2     5
    5     18
    18    7
    8     21
    21    10
    10    54     
    11    24
    24    13
    14    27
    27    16
    4     11
    16    56
    55    54
    43    44
    47    48
    39    43
    40    44
    39    46
    46    49
    19    43
    43    48
    43    49
    44    50
    50    51
    49    51
    37    38
    %34    35
    %33    36
    30 33%30    36
    29    35
    % Repe: 33 29 %34    29
    33 32 %35    32
    11    17
    12    18
    13    19
    14    20
    15    21
    16    22
    7     20
    14    25
    20    40
    44    45
    26    37
    47    50
    38    40
    37    39
    38    20
    37    19
    26    32
    25    31
    1     23
    23    6
    6     25
    53    28
    28    9
    9     26
    29    32
    26    40
    39    25
    40    43
    44    49
    29    26
    32    37
    20    44 % Añado para simetria
    39    44 % Añado para simetria
    
    %ANYADO
    %57 51 
    %57 48
    
    %58 51
    %58 47
    ];

for i=1:size(input.T,1)
    if input.T(i,1) > 33
        input.T(i,1) = input.T(i,1) -3;
    end
    
    if input.T(i,2) > 33
        input.T(i,2) = input.T(i,2) -3;
    end 
end

% input.xpoints = [
%   0, 0, 0;           %(1)
%   0, 0, -0.2;        %(2)
%   0, 1.6, 0;         %(3)
%   0, 1.6, -0.2;      %(4)
%   2, 0, 0;           %(5)
%   4, 0, 0;           %(6)
%   6, 0, 0;           %(7)
%   8, 0, 0;           %(8)
%   10, 0, 0;          %(9)
%   12, 0, 0;          %(10)
%   2, 1.65, 0;        %(11)
%   4, 1.65, 0;        %(12)
%   6, 1.65, 0;        %(13)
%   8, 1.65, 0;        %(14)
%   10, 1.65, 0;        %(15)
%   12, 1.65, 0;       %(16)
%     2, 0, -0.2;         %(17)
%   4, 0, -0.2;           %(18)
%   6, 0, -0.2;           %(19)
%   8, 0, -0.2;           %(20)
%   10, 0, -0.2;           %(21)
%   12, 0, -0.2;          %(22)
%   2, 1.65, -0.2;        %(23)
%   4, 1.65, -0.2;        %(24)
%   6, 1.65, -0.2;        %(25)
%   8, 1.65, -0.2;        %(26)
%   10, 1.65, -0.2;        %(27)
%   12, 1.65, -0.2;       %(28)
%   6.5, 2.6, -0.6;       %(29)
%   7.5, 2.6, -0.6;       %(30)
%   6.5, 2.6, -1.1;       %(31)
%   7.5, 2.6, -1.1;       %(32)
%   6.5, 3.1, -0.6;       %(33)
%   7.5, 3.1, -0.6;       %(34)
%   6.5, 3.1, -1.1;       %(35)
%   7.5, 3.1, -1.1;       %(36)
%   6, 1.6, -1.6;         %(37)
%   8, 1.6, -1.6;         %(38)
%   8, 0, -1.2;           %(39)
%   6, 0, -1.2;           %(40)
%   8.5, 1, -2.2;         %(41)
%   5.5, 1, -2.2;         %(42)
%   7.8,-2, -0.2;       %(43)
%   6.2,-2, -0.2;       %(44)
%   7.8,-2, -1;         %(45)
%   6.2,-2, -1;         %(46)
%   7.4,-3.8, -0.2;     %(47)
%   6.6,-3.8, -0.2;     %(48)
%   7.4,-3.8, -0.6;     %(49)
%   6.6,-3.8, -0.6;     %(50)
%   7,-5, -0.2;         %(51)
%   7,-5, -0.4;         %(52)
%   14,0, 0;            %(53)
%   14,0, -0.2;         %(54)
%   14,1.6, 0;          %(55)
%   14,1.6, -0.2;       %(56)
% ];
% 
% input.T = [
%      1     3
%      1     5
%      1     2
%      2     17
%      2     4
%      3     11
%      4     23
%      3     4
%      1     11
%      2     23
%      5     6
%      6     7
%      7     8
%      8     9
%      9    10
%     10    16
%     11    12
%     12    13
%     13    14
%     14    15
%     15    16
%      5    11
%      6    12
%      7    13
%      8    14
%      9    15
%      5    17
%      6    18
%      7    19
%      8    20
%      9    21
%     10    22
%     11    23
%     12    24
%     13    25
%     14    26
%     15    27
%     16    28
%     17    18
%     18    19
%     19    20
%     20    21
%     21    22
%     22    28
%     17    23
%     23    24
%     24    25
%     25    26
%     26    27
%     27    28
%     18    24
%     19    25
%     20    26
%     21    27
%      6    11
%     18    23
%     18    25
%      6    13
%      7    14
%     20    25
%      9    14
%     21    26
%      9    16
%     21    28
%     25    29
%     26    30
%     29    30
%     31    32
%     29    31
%     30    32
%     29    33
%     31    35
%     30    34
%     32    36
%     34    36
%     33    35
%     35    36
%     33    34
%     31    37
%     32    38
%     37    40
%     38    39
%     26    38
%     25    37
%     40    19
%     20    39
%     40    42
%     42    37
%     38    41
%     39    41
%     39    45
%     20    43
%     43    45
%     19    44
%     40    46
%     44    46
%     43    47
%     44    48
%     48    50
%     47    49
%     45    49
%     46    50
%     50    52
%     48    51
%     47    51
%     49    52
%     51    52
%     10    53
%     22    54
%     53    54
%     53    55
%     54    56
%     16    55
%     28    56
%     55    56
%     16    53
%     28    54
%     41    37
%     42    38
%     45    46
%     39    40
%     49    50
%     2     3
%     2     5
%     5     18
%     18    7
%     8     21
%     21    10
%     10    54     
%     11    24
%     24    13
%     14    27
%     27    16
%     4     11
%     16    56
%     55    54
%     43    44
%     47    48
%     39    43
%     40    44
%     39    46
%     46    49
%     19    43
%     43    48
%     43    49
%     44    50
%     50    51
%     49    51
%     37    38
%     34    35
%     33    36
%     30    36
%     29    35
%     34    29
%     35    32
%     11    17
%     12    18
%     13    19
%     14    20
%     15    21
%     16    22
%     7     20
%     14    25
%     20    40
%     44    45
%     26    37
%     47    50
%     38    40
%     37    39
%     38    20
%     37    19
%     26    32
%     25    31
%     1     23
%     23    6
%     6     25
%     53    28
%     28    9
%     9     26
%     29    32
%     26    40
%     39    25
%     40    43
%     44    49
%     29    26
%     32    37
%     ];


clear T
    
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
            48-3, 1, 0;...
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
              3, 3, lift_parcial;...
              5, 3, lift_parcial;...
              6, 3, lift_parcial; 
              7, 3, lift_parcial;...
              8, 3, lift_parcial;...
              9, 3, lift_parcial;...
              10, 3, lift_parcial;...
              11, 3, lift_parcial;...
              12, 3, lift_parcial;...
              13, 3, lift_parcial;...
              14, 3, lift_parcial;...
              15, 3, lift_parcial;...
              16, 3, lift_parcial;...
              53, 3, lift_parcial];
          
          
%% TERMIC
input.alpha = 1;

input.Delta_T = ones(size(input.T,1),1) * -40;

              
