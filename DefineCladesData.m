function clade=DefineClades

% function Clade=DefineClades
%
% Indo-European clades of interest
% GKN 28/10

GlobalSwitches;

CN=1;%*
clade{CN}.language={'Irish_A' 'Irish_B' 'Welsh_N' 'Welsh_C' 'Breton_List' 'Breton_SE' 'Breton_ST'};
clade{CN}.name='Celtic';
clade{CN}.rootrange=[1700 inf];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Welsh_N' 'Welsh_C' 'Breton_List' 'Breton_SE' 'Breton_ST'};
clade{CN}.name='Brythonic';
clade{CN}.rootrange=[1450 1600]; 
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Romanian_List' 'Vlach' 'Italian' 'Ladin' 'Provencal' 'French' 'Walloon'...
      'French_Creole_C' 'French_Creole_D' 'Sardinian_N' 'Sardinian_L' 'Sardinian_C' 'Spanish'...
      'Portuguese_ST' 'Brazilian' 'Catalan'};
clade{CN}.name='Italic';
clade{CN}.rootrange=[1700 1850];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Italian' 'Ladin' 'Provencal' 'French' 'Walloon' 'French_Creole_C' 'French_Creole_D'...
      'Spanish' 'Portuguese_ST' 'Brazilian' 'Catalan'};
clade{CN}.name='IberianFrench';
clade{CN}.rootrange=[1200 1550];
clade{CN}.adamrange=[];

CN=CN+1;
clade{CN}.language={'Irish_A' 'Irish_B' 'Welsh_N' 'Welsh_C' 'Breton_List' 'Breton_SE' 'Breton_ST'...
      'Romanian_List' 'Vlach' 'Italian' 'Ladin' 'Provencal' 'French' 'Walloon' 'French_Creole_C'...
      'French_Creole_D' 'Sardinian_N' 'Sardinian_L' 'Sardinian_C' 'Spanish' 'Portuguese_ST'...
      'Brazilian' 'Catalan'};
clade{CN}.name='ItaloCeltic';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'German_ST' 'Penn_Dutch' 'Dutch_List' 'Afrikaans' 'Flemish' 'Frisian'...
      'Swedish_Up' 'Swedish_VL' 'Swedish_List' 'Danish' 'Riksmal' 'Icelandic_ST' 'Faroese'...
      'English_ST' 'Takitaki'};
clade{CN}.name='Germanic';
clade{CN}.rootrange=[1750 1950];
clade{CN}.adamrange=[];

CN=CN+1;
clade{CN}.language={'German_ST' 'Penn_Dutch' 'Dutch_List' 'Afrikaans' 'Flemish' 'Frisian'...         
      'English_ST' 'Takitaki'};
clade{CN}.name='WestGermanic';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[];

CN=CN+1;
clade{CN}.language={'Swedish_Up' 'Swedish_VL' 'Swedish_List' 'Danish' 'Riksmal' 'Icelandic_ST' 'Faroese'};
clade{CN}.name='NorthGermanic';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Lithuanian_O' 'Lithuanian_ST' 'Latvian' 'Slovenian' 'Lusatian_L' 'Lusatian_U'...
      'Czech' 'Slovak' 'Czech_E' 'Ukrainian' 'Byelorussian' 'Polish' 'Russian' 'Macedonian'...
      'Bulgarian' 'Serbocroatian'};
clade{CN}.name='BaltoSlav';
clade{CN}.rootrange=[1900 3400];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Slovenian' 'Lusatian_L' 'Lusatian_U' 'Czech' 'Slovak' 'Czech_E' 'Ukrainian'...
      'Byelorussian' 'Polish' 'Russian' 'Macedonian' 'Bulgarian' 'Serbocroatian'};
clade{CN}.name='Slav';
clade{CN}.rootrange=[1300 inf];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Gypsy_Gk' 'Singhalese' 'Kashmiri' 'Marathi' 'Gujarati' 'Panjabi_ST' 'Lahnda'...
      'Hindi' 'Bengali' 'Nepali_List' 'Khaskura'};
clade{CN}.name='Indic';
clade{CN}.rootrange=[2200 inf];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Gypsy_Gk' 'Singhalese' 'Kashmiri' 'Marathi' 'Gujarati' 'Panjabi_ST' 'Lahnda'...
      'Hindi' 'Bengali' 'Nepali_List' 'Khaskura' 'Ossetic' 'Afghan' 'Waziri' 'Persian_List'...
      'Tadzik' 'Baluchi' 'Wakhi'};
clade{CN}.name='IndoIranian';
clade{CN}.rootrange=[3000 inf];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Ossetic' 'Afghan' 'Waziri' 'Persian_List' 'Tadzik' 'Baluchi' 'Wakhi'};
clade{CN}.name='Iranian';
clade{CN}.rootrange=[2500 inf];
clade{CN}.adamrange=[];

CN=CN+1;
clade{CN}.language={'Albanian_T' 'Albanian_Top' 'Albanian_G' 'Albanian_K' 'Albanian_C'};
clade{CN}.name='Albanian';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[];

CN=CN+1;%*
clade{CN}.language={'Greek_ML' 'Greek_MD' 'Greek_Mod' 'Greek_D' 'Greek_K'};
clade{CN}.name='Greek';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[3500 inf];

CN=CN+1;
clade{CN}.language={'Armenian_Mod' 'Armenian_List'};
clade{CN}.name='Armenian';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[];

CN=CN+1;
clade{CN}.language={'TOCHARIAN_A' 'TOCHARIAN_B'};
clade{CN}.name='Tocharic';
clade{CN}.rootrange=[1650 inf];
clade{CN}.adamrange=[];

CN=CN+1;
clade{CN}.language={'Irish_A'    'Irish_B'    'Welsh_N'    'Welsh_C'...
    'Breton_List'    'Breton_SE'    'Breton_ST'    'Romanian_List'    'Vlach'    'Italian'    'Ladin'...
    'Provencal'    'French'    'Walloon'    'French_Creole_C'    'French_Creole_D'    'Sardinian_N'...
    'Sardinian_L'    'Sardinian_C'    'Spanish'    'Portuguese_ST'    'Brazilian'    'Catalan'...
    'German_ST'    'Penn_Dutch'    'Dutch_List'    'Afrikaans'    'Flemish'    'Frisian'    'Swedish_Up'...
    'Swedish_VL'    'Swedish_List'    'Danish'    'Riksmal'    'Icelandic_ST'    'Faroese'    'English_ST'...
    'Takitaki'    'Lithuanian_O'    'Lithuanian_ST'    'Latvian'    'Slovenian'    'Lusatian_L'...
    'Lusatian_U'    'Czech'    'Slovak'    'Czech_E'    'Ukrainian'    'Byelorussian'    'Polish'...
    'Russian'    'Macedonian'    'Bulgarian'    'Serbocroatian'    'Gypsy_Gk'    'Singhalese'    'Kashmiri'...
    'Marathi'    'Gujarati'    'Panjabi_ST'    'Lahnda'    'Hindi'    'Bengali'    'Nepali_List'...
    'Khaskura'    'Greek_ML'    'Greek_MD'    'Greek_Mod'    'Greek_D'    'Greek_K'    'Armenian_Mod'...
    'Armenian_List'    'Ossetic'    'Afghan'    'Waziri'    'Persian_List'    'Tadzik'    'Baluchi'...
    'Wakhi'    'Albanian_T'    'Albanian_Top'    'Albanian_G'    'Albanian_K'    'Albanian_C'};
clade{CN}.name='Modern';
clade{CN}.rootrange=[];
clade{CN}.adamrange=[];
