    function targMap = targDataMap(),

    ;%***********************
    ;% Create Parameter Map *
    ;%***********************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 2;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc paramMap
        ;%
        paramMap.nSections           = nTotSects;
        paramMap.sectIdxOffset       = sectIdxOffset;
            paramMap.sections(nTotSects) = dumSection; %prealloc
        paramMap.nTotData            = -1;

        ;%
        ;% Auto data (rtP)
        ;%
            section.nData     = 54;
            section.data(54)  = dumData; %prealloc

                    ;% rtP.C_ts
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtP.J_eq
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% rtP.K_sa
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 2;

                    ;% rtP.K_sia
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 3;

                    ;% rtP.L_d
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 4;

                    ;% rtP.L_ls
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 5;

                    ;% rtP.L_q
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 6;

                    ;% rtP.Pp
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 7;

                    ;% rtP.R_s_ref
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 8;

                    ;% rtP.R_ts
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 9;

                    ;% rtP.T_s_ref
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 10;

                    ;% rtP.alpha_Cu
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 11;

                    ;% rtP.b_a
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 12;

                    ;% rtP.b_eq
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 13;

                    ;% rtP.k_i
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 14;

                    ;% rtP.k_l
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 15;

                    ;% rtP.k_omega2
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 16;

                    ;% rtP.k_theta2
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 17;

                    ;% rtP.lambda_m
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 18;

                    ;% rtP.r
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 19;

                    ;% rtP.fromWS_Signal1_Time0
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 20;

                    ;% rtP.fromWS_Signal1_Data0
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 22;

                    ;% rtP.Integrator_IC
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 24;

                    ;% rtP.Integrator1_IC
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 25;

                    ;% rtP.Integrator_IC_i2wxdirj0b
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 26;

                    ;% rtP.Integrator1_IC_oww15y3mqv
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 27;

                    ;% rtP.Integrator1_IC_gd0uumlv2e
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 28;

                    ;% rtP.Gain4_Gain
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 29;

                    ;% rtP.FromWorkspace_Time0
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 30;

                    ;% rtP.FromWorkspace_Data0
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 171;

                    ;% rtP.FromWorkspace_Time0_o2ayzzhu4m
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 312;

                    ;% rtP.FromWorkspace_Data0_knbhfzr1k4
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 314;

                    ;% rtP.Integrator_IC_gpeemvg0a0
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 316;

                    ;% rtP.Integrator1_IC_epj0vehlxx
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 317;

                    ;% rtP.Integrator2_IC
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 318;

                    ;% rtP.Gain_Gain
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 319;

                    ;% rtP.Gain1_Gain
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 320;

                    ;% rtP.Gain2_Gain
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 321;

                    ;% rtP.Integrator2_IC_cltf442w1a
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 322;

                    ;% rtP.Gain3_Gain
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 323;

                    ;% rtP.Integrator2_IC_lurynlyjlu
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 324;

                    ;% rtP.Integrator_IC_n5atnjvnfo
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 325;

                    ;% rtP.Gain_Gain_bpfw2sk42p
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 326;

                    ;% rtP.Gain1_Gain_o051fegvds
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 327;

                    ;% rtP.Gain2_Gain_fpelshvuq5
                    section.data(45).logicalSrcIdx = 44;
                    section.data(45).dtTransOffset = 328;

                    ;% rtP.FromWorkspace_Time0_bh3fvyhkqh
                    section.data(46).logicalSrcIdx = 45;
                    section.data(46).dtTransOffset = 329;

                    ;% rtP.FromWorkspace_Data0_ahxyd5aapb
                    section.data(47).logicalSrcIdx = 46;
                    section.data(47).dtTransOffset = 331;

                    ;% rtP.Gain1_Gain_bfommv4lak
                    section.data(48).logicalSrcIdx = 47;
                    section.data(48).dtTransOffset = 333;

                    ;% rtP.Gain2_Gain_gt5xi45qct
                    section.data(49).logicalSrcIdx = 48;
                    section.data(49).dtTransOffset = 334;

                    ;% rtP.FromWorkspace_Time0_fivuiqtoah
                    section.data(50).logicalSrcIdx = 49;
                    section.data(50).dtTransOffset = 335;

                    ;% rtP.FromWorkspace_Data0_kevgbgk3mx
                    section.data(51).logicalSrcIdx = 50;
                    section.data(51).dtTransOffset = 337;

                    ;% rtP.Integrator2_IC_dzcpbao03e
                    section.data(52).logicalSrcIdx = 51;
                    section.data(52).dtTransOffset = 339;

                    ;% rtP.Constant_Value
                    section.data(53).logicalSrcIdx = 52;
                    section.data(53).dtTransOffset = 340;

                    ;% rtP.Constant2_Value
                    section.data(54).logicalSrcIdx = 53;
                    section.data(54).dtTransOffset = 341;

            nTotData = nTotData + section.nData;
            paramMap.sections(1) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% rtP.ManualSwitch1_CurrentSetting
                    section.data(1).logicalSrcIdx = 54;
                    section.data(1).dtTransOffset = 0;

                    ;% rtP.ManualSwitch_CurrentSetting
                    section.data(2).logicalSrcIdx = 55;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            paramMap.sections(2) = section;
            clear section


            ;%
            ;% Non-auto Data (parameter)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        paramMap.nTotData = nTotData;



    ;%**************************
    ;% Create Block Output Map *
    ;%**************************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 1;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc sigMap
        ;%
        sigMap.nSections           = nTotSects;
        sigMap.sectIdxOffset       = sectIdxOffset;
            sigMap.sections(nTotSects) = dumSection; %prealloc
        sigMap.nTotData            = -1;

        ;%
        ;% Auto data (rtB)
        ;%
            section.nData     = 123;
            section.data(123)  = dumData; %prealloc

                    ;% rtB.k2sn1moatl
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtB.m4zsvrd15a
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% rtB.jpgl5uqnen
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 2;

                    ;% rtB.hpaqwswri4
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 3;

                    ;% rtB.ojegsw1r4w
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 4;

                    ;% rtB.hpqeqr0o0i
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 5;

                    ;% rtB.efmorfp5hk
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 6;

                    ;% rtB.homtptr1wn
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 7;

                    ;% rtB.dbz2he40vm
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 8;

                    ;% rtB.mbfoowmnkk
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 9;

                    ;% rtB.dlox32o0ua
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 10;

                    ;% rtB.oh5qhuuk4l
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 11;

                    ;% rtB.dx12nyer4s
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 12;

                    ;% rtB.c5fdscj3uw
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 13;

                    ;% rtB.lnhoszrvh2
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 14;

                    ;% rtB.c11k2clsrb
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 15;

                    ;% rtB.mqiqfyapeu
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 16;

                    ;% rtB.kylrj4an0z
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 17;

                    ;% rtB.dshebsf5x3
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 18;

                    ;% rtB.ly2xklr4fz
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 19;

                    ;% rtB.h1lgwz5ovu
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 20;

                    ;% rtB.gifca5zh5d
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 21;

                    ;% rtB.gysjfnsk1k
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 22;

                    ;% rtB.ajowt5z1q5
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 23;

                    ;% rtB.fwvyoitwh3
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 24;

                    ;% rtB.netv40pb4f
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 25;

                    ;% rtB.jo0ai5lhtc
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 26;

                    ;% rtB.k532jfpdyn
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 27;

                    ;% rtB.nxola2u2g3
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 28;

                    ;% rtB.kmls4v4rl0
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 29;

                    ;% rtB.af03zirre2
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 30;

                    ;% rtB.dhufv0kct4
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 31;

                    ;% rtB.fd152jijlt
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 32;

                    ;% rtB.klsfhtexlj
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 33;

                    ;% rtB.p3yvsgsymc
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 34;

                    ;% rtB.elrtlardao
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 35;

                    ;% rtB.pxphckphjm
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 36;

                    ;% rtB.k25kf23xp1
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 37;

                    ;% rtB.g5whoifsex
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 38;

                    ;% rtB.fb0zj32vbu
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 39;

                    ;% rtB.kga0ttwxih
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 40;

                    ;% rtB.kf3sjuqaae
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 41;

                    ;% rtB.fwd4ryyrfq
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 42;

                    ;% rtB.iamwwft4eh
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 45;

                    ;% rtB.edszxpjalm
                    section.data(45).logicalSrcIdx = 44;
                    section.data(45).dtTransOffset = 46;

                    ;% rtB.azzpx0o04y
                    section.data(46).logicalSrcIdx = 45;
                    section.data(46).dtTransOffset = 47;

                    ;% rtB.nd4l40ce0g
                    section.data(47).logicalSrcIdx = 46;
                    section.data(47).dtTransOffset = 48;

                    ;% rtB.pr2t5fqsxs
                    section.data(48).logicalSrcIdx = 47;
                    section.data(48).dtTransOffset = 49;

                    ;% rtB.nly0v00iut
                    section.data(49).logicalSrcIdx = 48;
                    section.data(49).dtTransOffset = 50;

                    ;% rtB.cri0rzhvjd
                    section.data(50).logicalSrcIdx = 49;
                    section.data(50).dtTransOffset = 51;

                    ;% rtB.o1eaam1gvh
                    section.data(51).logicalSrcIdx = 50;
                    section.data(51).dtTransOffset = 52;

                    ;% rtB.i1cztljifv
                    section.data(52).logicalSrcIdx = 51;
                    section.data(52).dtTransOffset = 53;

                    ;% rtB.hst2mc5opt
                    section.data(53).logicalSrcIdx = 52;
                    section.data(53).dtTransOffset = 54;

                    ;% rtB.kpwiznui1d
                    section.data(54).logicalSrcIdx = 53;
                    section.data(54).dtTransOffset = 55;

                    ;% rtB.ollu0ludgl
                    section.data(55).logicalSrcIdx = 54;
                    section.data(55).dtTransOffset = 56;

                    ;% rtB.bpf3qazr44
                    section.data(56).logicalSrcIdx = 55;
                    section.data(56).dtTransOffset = 57;

                    ;% rtB.nc1t1u3huh
                    section.data(57).logicalSrcIdx = 56;
                    section.data(57).dtTransOffset = 58;

                    ;% rtB.dnzl4cfy5s
                    section.data(58).logicalSrcIdx = 57;
                    section.data(58).dtTransOffset = 59;

                    ;% rtB.ougga4l4aa
                    section.data(59).logicalSrcIdx = 58;
                    section.data(59).dtTransOffset = 60;

                    ;% rtB.pxdilwla4w
                    section.data(60).logicalSrcIdx = 59;
                    section.data(60).dtTransOffset = 61;

                    ;% rtB.gl1fdweeu4
                    section.data(61).logicalSrcIdx = 60;
                    section.data(61).dtTransOffset = 62;

                    ;% rtB.fhjbzzgwbv
                    section.data(62).logicalSrcIdx = 61;
                    section.data(62).dtTransOffset = 63;

                    ;% rtB.kl3fuewomp
                    section.data(63).logicalSrcIdx = 62;
                    section.data(63).dtTransOffset = 64;

                    ;% rtB.h1yubcelv5
                    section.data(64).logicalSrcIdx = 63;
                    section.data(64).dtTransOffset = 65;

                    ;% rtB.k1afv5mxmb
                    section.data(65).logicalSrcIdx = 64;
                    section.data(65).dtTransOffset = 66;

                    ;% rtB.lggf4lsbsd
                    section.data(66).logicalSrcIdx = 65;
                    section.data(66).dtTransOffset = 67;

                    ;% rtB.evhxefbhrl
                    section.data(67).logicalSrcIdx = 66;
                    section.data(67).dtTransOffset = 68;

                    ;% rtB.gt3nq2cscb
                    section.data(68).logicalSrcIdx = 67;
                    section.data(68).dtTransOffset = 69;

                    ;% rtB.aqulfpaq14
                    section.data(69).logicalSrcIdx = 68;
                    section.data(69).dtTransOffset = 70;

                    ;% rtB.bd4jar1csl
                    section.data(70).logicalSrcIdx = 69;
                    section.data(70).dtTransOffset = 71;

                    ;% rtB.dvu22fpphc
                    section.data(71).logicalSrcIdx = 70;
                    section.data(71).dtTransOffset = 72;

                    ;% rtB.fnvmvvqodz
                    section.data(72).logicalSrcIdx = 71;
                    section.data(72).dtTransOffset = 73;

                    ;% rtB.dz4hiszfdk
                    section.data(73).logicalSrcIdx = 72;
                    section.data(73).dtTransOffset = 76;

                    ;% rtB.lswxnkpgli
                    section.data(74).logicalSrcIdx = 73;
                    section.data(74).dtTransOffset = 77;

                    ;% rtB.oqeh2ikwco
                    section.data(75).logicalSrcIdx = 74;
                    section.data(75).dtTransOffset = 78;

                    ;% rtB.jxvum3zvoh
                    section.data(76).logicalSrcIdx = 75;
                    section.data(76).dtTransOffset = 79;

                    ;% rtB.ou4cqm5an5
                    section.data(77).logicalSrcIdx = 76;
                    section.data(77).dtTransOffset = 80;

                    ;% rtB.onoapnnemw
                    section.data(78).logicalSrcIdx = 77;
                    section.data(78).dtTransOffset = 81;

                    ;% rtB.iktm1emn42
                    section.data(79).logicalSrcIdx = 78;
                    section.data(79).dtTransOffset = 82;

                    ;% rtB.kwquioc0a2
                    section.data(80).logicalSrcIdx = 79;
                    section.data(80).dtTransOffset = 83;

                    ;% rtB.caar3rhznz
                    section.data(81).logicalSrcIdx = 80;
                    section.data(81).dtTransOffset = 84;

                    ;% rtB.gbixzpfrni
                    section.data(82).logicalSrcIdx = 81;
                    section.data(82).dtTransOffset = 85;

                    ;% rtB.jwtzlohzhs
                    section.data(83).logicalSrcIdx = 82;
                    section.data(83).dtTransOffset = 86;

                    ;% rtB.cbv3kpprwx
                    section.data(84).logicalSrcIdx = 83;
                    section.data(84).dtTransOffset = 87;

                    ;% rtB.cdjjap5tb4
                    section.data(85).logicalSrcIdx = 84;
                    section.data(85).dtTransOffset = 88;

                    ;% rtB.k5iqwa5dub
                    section.data(86).logicalSrcIdx = 85;
                    section.data(86).dtTransOffset = 89;

                    ;% rtB.lbeji2oewc
                    section.data(87).logicalSrcIdx = 86;
                    section.data(87).dtTransOffset = 90;

                    ;% rtB.crxm0t2udf
                    section.data(88).logicalSrcIdx = 87;
                    section.data(88).dtTransOffset = 91;

                    ;% rtB.avres0c4zw
                    section.data(89).logicalSrcIdx = 88;
                    section.data(89).dtTransOffset = 92;

                    ;% rtB.buftnulajc
                    section.data(90).logicalSrcIdx = 89;
                    section.data(90).dtTransOffset = 93;

                    ;% rtB.pwtikpf0j0
                    section.data(91).logicalSrcIdx = 90;
                    section.data(91).dtTransOffset = 94;

                    ;% rtB.marzp3jauv
                    section.data(92).logicalSrcIdx = 91;
                    section.data(92).dtTransOffset = 95;

                    ;% rtB.dymq05epmd
                    section.data(93).logicalSrcIdx = 92;
                    section.data(93).dtTransOffset = 96;

                    ;% rtB.d3woaso1vk
                    section.data(94).logicalSrcIdx = 93;
                    section.data(94).dtTransOffset = 97;

                    ;% rtB.lata252i34
                    section.data(95).logicalSrcIdx = 94;
                    section.data(95).dtTransOffset = 98;

                    ;% rtB.p5w0nljbtt
                    section.data(96).logicalSrcIdx = 95;
                    section.data(96).dtTransOffset = 99;

                    ;% rtB.jkeaxbjqzp
                    section.data(97).logicalSrcIdx = 96;
                    section.data(97).dtTransOffset = 100;

                    ;% rtB.lpe0g0uuxy
                    section.data(98).logicalSrcIdx = 97;
                    section.data(98).dtTransOffset = 101;

                    ;% rtB.kwzodnsn51
                    section.data(99).logicalSrcIdx = 98;
                    section.data(99).dtTransOffset = 102;

                    ;% rtB.kjg2pniu00
                    section.data(100).logicalSrcIdx = 99;
                    section.data(100).dtTransOffset = 103;

                    ;% rtB.c4ogmaamk3
                    section.data(101).logicalSrcIdx = 100;
                    section.data(101).dtTransOffset = 104;

                    ;% rtB.jivucflgxz
                    section.data(102).logicalSrcIdx = 101;
                    section.data(102).dtTransOffset = 105;

                    ;% rtB.b0djz0xuyx
                    section.data(103).logicalSrcIdx = 102;
                    section.data(103).dtTransOffset = 106;

                    ;% rtB.moi1nwsely
                    section.data(104).logicalSrcIdx = 103;
                    section.data(104).dtTransOffset = 107;

                    ;% rtB.apbcudoiru
                    section.data(105).logicalSrcIdx = 104;
                    section.data(105).dtTransOffset = 108;

                    ;% rtB.bptosmp44x
                    section.data(106).logicalSrcIdx = 105;
                    section.data(106).dtTransOffset = 109;

                    ;% rtB.fxltoioz2o
                    section.data(107).logicalSrcIdx = 106;
                    section.data(107).dtTransOffset = 110;

                    ;% rtB.k3waxxng0u
                    section.data(108).logicalSrcIdx = 107;
                    section.data(108).dtTransOffset = 111;

                    ;% rtB.ek2tbc2vjr
                    section.data(109).logicalSrcIdx = 108;
                    section.data(109).dtTransOffset = 112;

                    ;% rtB.j1emwd2qae
                    section.data(110).logicalSrcIdx = 109;
                    section.data(110).dtTransOffset = 113;

                    ;% rtB.lsyqe5i55w
                    section.data(111).logicalSrcIdx = 110;
                    section.data(111).dtTransOffset = 114;

                    ;% rtB.c4n4jb2dgo
                    section.data(112).logicalSrcIdx = 111;
                    section.data(112).dtTransOffset = 115;

                    ;% rtB.ny2zw5nsgp
                    section.data(113).logicalSrcIdx = 112;
                    section.data(113).dtTransOffset = 116;

                    ;% rtB.fzxuttsywh
                    section.data(114).logicalSrcIdx = 113;
                    section.data(114).dtTransOffset = 117;

                    ;% rtB.lb5025wpa3
                    section.data(115).logicalSrcIdx = 114;
                    section.data(115).dtTransOffset = 118;

                    ;% rtB.ib34wfcjm3
                    section.data(116).logicalSrcIdx = 115;
                    section.data(116).dtTransOffset = 119;

                    ;% rtB.nf03pyptgh
                    section.data(117).logicalSrcIdx = 116;
                    section.data(117).dtTransOffset = 120;

                    ;% rtB.iss1dct3mr
                    section.data(118).logicalSrcIdx = 117;
                    section.data(118).dtTransOffset = 121;

                    ;% rtB.hpveetbx1j
                    section.data(119).logicalSrcIdx = 118;
                    section.data(119).dtTransOffset = 122;

                    ;% rtB.eufxrbxe4z
                    section.data(120).logicalSrcIdx = 119;
                    section.data(120).dtTransOffset = 123;

                    ;% rtB.j0wuzgingg
                    section.data(121).logicalSrcIdx = 120;
                    section.data(121).dtTransOffset = 124;

                    ;% rtB.aduksi4izn
                    section.data(122).logicalSrcIdx = 121;
                    section.data(122).dtTransOffset = 125;

                    ;% rtB.pia0trm2jv
                    section.data(123).logicalSrcIdx = 122;
                    section.data(123).dtTransOffset = 126;

            nTotData = nTotData + section.nData;
            sigMap.sections(1) = section;
            clear section


            ;%
            ;% Non-auto Data (signal)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        sigMap.nTotData = nTotData;



    ;%*******************
    ;% Create DWork Map *
    ;%*******************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 6;
        sectIdxOffset = 1;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc dworkMap
        ;%
        dworkMap.nSections           = nTotSects;
        dworkMap.sectIdxOffset       = sectIdxOffset;
            dworkMap.sections(nTotSects) = dumSection; %prealloc
        dworkMap.nTotData            = -1;

        ;%
        ;% Auto data (rtDW)
        ;%
            section.nData     = 4;
            section.data(4)  = dumData; %prealloc

                    ;% rtDW.parg3eaod2
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.na0wt4oepn
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.dejwxn1xy1
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.dcb22robm1
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 3;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 21;
            section.data(21)  = dumData; %prealloc

                    ;% rtDW.jrzquc22ff.TimePtr
                    section.data(1).logicalSrcIdx = 4;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.ml4gzok05f.TimePtr
                    section.data(2).logicalSrcIdx = 5;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.diriss0ycm.TimePtr
                    section.data(3).logicalSrcIdx = 6;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.dsna2h3s4j.LoggedData
                    section.data(4).logicalSrcIdx = 7;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.pv5y2snvjf.LoggedData
                    section.data(5).logicalSrcIdx = 8;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.imrctj1fis.LoggedData
                    section.data(6).logicalSrcIdx = 9;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.ak1r52kzjl.LoggedData
                    section.data(7).logicalSrcIdx = 10;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.o0dx25v4bj.LoggedData
                    section.data(8).logicalSrcIdx = 11;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.opv5nm4cbf.AQHandles
                    section.data(9).logicalSrcIdx = 12;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.cqxqnkszbm.AQHandles
                    section.data(10).logicalSrcIdx = 13;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.dt1kwx4eay.AQHandles
                    section.data(11).logicalSrcIdx = 14;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.dnhkap2j1u.AQHandles
                    section.data(12).logicalSrcIdx = 15;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.put0eavhrg.LoggedData
                    section.data(13).logicalSrcIdx = 16;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.n02okhcqno.LoggedData
                    section.data(14).logicalSrcIdx = 17;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.a0ojjrdjng.LoggedData
                    section.data(15).logicalSrcIdx = 18;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.ebooqd5zoc.LoggedData
                    section.data(16).logicalSrcIdx = 19;
                    section.data(16).dtTransOffset = 17;

                    ;% rtDW.i2e1lohrez.LoggedData
                    section.data(17).logicalSrcIdx = 20;
                    section.data(17).dtTransOffset = 20;

                    ;% rtDW.mvlzl14qts.LoggedData
                    section.data(18).logicalSrcIdx = 21;
                    section.data(18).dtTransOffset = 21;

                    ;% rtDW.ekcesp1k2p.TimePtr
                    section.data(19).logicalSrcIdx = 22;
                    section.data(19).dtTransOffset = 22;

                    ;% rtDW.hrwve4sncv.TimePtr
                    section.data(20).logicalSrcIdx = 23;
                    section.data(20).dtTransOffset = 23;

                    ;% rtDW.pqdji2nimh.AQHandles
                    section.data(21).logicalSrcIdx = 24;
                    section.data(21).dtTransOffset = 24;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section

            section.nData     = 4;
            section.data(4)  = dumData; %prealloc

                    ;% rtDW.asp1sby2rd
                    section.data(1).logicalSrcIdx = 25;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.ekogk2wdou
                    section.data(2).logicalSrcIdx = 26;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.axuxe1t12v
                    section.data(3).logicalSrcIdx = 27;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.fpyg5ifl2j
                    section.data(4).logicalSrcIdx = 28;
                    section.data(4).dtTransOffset = 3;

            nTotData = nTotData + section.nData;
            dworkMap.sections(3) = section;
            clear section

            section.nData     = 5;
            section.data(5)  = dumData; %prealloc

                    ;% rtDW.jptvtcgs5c.PrevIndex
                    section.data(1).logicalSrcIdx = 29;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.l2ca51dh03.PrevIndex
                    section.data(2).logicalSrcIdx = 30;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.h33gajte4x.PrevIndex
                    section.data(3).logicalSrcIdx = 31;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.gjvhwykizl.PrevIndex
                    section.data(4).logicalSrcIdx = 32;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.ojtfq2nxq4.PrevIndex
                    section.data(5).logicalSrcIdx = 33;
                    section.data(5).dtTransOffset = 4;

            nTotData = nTotData + section.nData;
            dworkMap.sections(4) = section;
            clear section

            section.nData     = 4;
            section.data(4)  = dumData; %prealloc

                    ;% rtDW.iqsjj14to5
                    section.data(1).logicalSrcIdx = 34;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.derxydwg12
                    section.data(2).logicalSrcIdx = 35;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.ef11rieqzy
                    section.data(3).logicalSrcIdx = 36;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.aa1tc5fl4n
                    section.data(4).logicalSrcIdx = 37;
                    section.data(4).dtTransOffset = 3;

            nTotData = nTotData + section.nData;
            dworkMap.sections(5) = section;
            clear section

            section.nData     = 4;
            section.data(4)  = dumData; %prealloc

                    ;% rtDW.hpxndlhxow
                    section.data(1).logicalSrcIdx = 38;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.d3jkgjc1oh
                    section.data(2).logicalSrcIdx = 39;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.gh0bm4uzee
                    section.data(3).logicalSrcIdx = 40;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.lljtkeuzy4
                    section.data(4).logicalSrcIdx = 41;
                    section.data(4).dtTransOffset = 3;

            nTotData = nTotData + section.nData;
            dworkMap.sections(6) = section;
            clear section


            ;%
            ;% Non-auto Data (dwork)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        dworkMap.nTotData = nTotData;



    ;%
    ;% Add individual maps to base struct.
    ;%

    targMap.paramMap  = paramMap;
    targMap.signalMap = sigMap;
    targMap.dworkMap  = dworkMap;

    ;%
    ;% Add checksums to base struct.
    ;%


    targMap.checksum0 = 883753185;
    targMap.checksum1 = 1892521921;
    targMap.checksum2 = 2740970809;
    targMap.checksum3 = 2593317770;

