supera::ImageMeta3D meta;
meta.set(-500, -500, -500, 500, 500, 500, 2500, 2500, 2500);

supera::EventInput evInput;
evInput.reserve(16);
supera::ParticleInput particleInput0;
particleInput0.pcloud.reserve(112);
supera::EDep particleInput0_edep0;
particleInput0_edep0.x = 0.6;
particleInput0_edep0.y = 0.6;
particleInput0_edep0.z = -25;
particleInput0_edep0.t = 1.97256;
particleInput0_edep0.e = 0.500036;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep0));
supera::EDep particleInput0_edep1;
particleInput0_edep1.x = 0.6;
particleInput0_edep1.y = 0.6;
particleInput0_edep1.z = -25.4;
particleInput0_edep1.t = 1.98512;
particleInput0_edep1.e = 1.00007;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep1));
supera::EDep particleInput0_edep2;
particleInput0_edep2.x = 0.6;
particleInput0_edep2.y = 0.6;
particleInput0_edep2.z = -25.8;
particleInput0_edep2.t = 1.99768;
particleInput0_edep2.e = 0.921565;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep2));
supera::EDep particleInput0_edep3;
particleInput0_edep3.x = 0.6;
particleInput0_edep3.y = 0.6;
particleInput0_edep3.z = -25.8;
particleInput0_edep3.t = 2.01024;
particleInput0_edep3.e = 0.0570227;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep3));
supera::EDep particleInput0_edep4;
particleInput0_edep4.x = 0.6;
particleInput0_edep4.y = 0.6;
particleInput0_edep4.z = -26.2;
particleInput0_edep4.t = 2.01357;
particleInput0_edep4.e = 0.252422;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep4));
supera::EDep particleInput0_edep5;
particleInput0_edep5.x = 0.6;
particleInput0_edep5.y = 0.6;
particleInput0_edep5.z = -26.2;
particleInput0_edep5.t = 2.0169;
particleInput0_edep5.e = 0.572607;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep5));
supera::EDep particleInput0_edep6;
particleInput0_edep6.x = 0.6;
particleInput0_edep6.y = 0.6;
particleInput0_edep6.z = -26.6;
particleInput0_edep6.t = 2.03086;
particleInput0_edep6.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep6));
supera::EDep particleInput0_edep7;
particleInput0_edep7.x = 0.6;
particleInput0_edep7.y = 0.6;
particleInput0_edep7.z = -27;
particleInput0_edep7.t = 2.04481;
particleInput0_edep7.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep7));
supera::EDep particleInput0_edep8;
particleInput0_edep8.x = 0.6;
particleInput0_edep8.y = 0.6;
particleInput0_edep8.z = -27.4;
particleInput0_edep8.t = 2.05876;
particleInput0_edep8.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep8));
supera::EDep particleInput0_edep9;
particleInput0_edep9.x = 0.6;
particleInput0_edep9.y = 0.6;
particleInput0_edep9.z = -27.8;
particleInput0_edep9.t = 2.07271;
particleInput0_edep9.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep9));
supera::EDep particleInput0_edep10;
particleInput0_edep10.x = 0.6;
particleInput0_edep10.y = 0.6;
particleInput0_edep10.z = -28.2;
particleInput0_edep10.t = 2.08666;
particleInput0_edep10.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep10));
supera::EDep particleInput0_edep11;
particleInput0_edep11.x = 0.6;
particleInput0_edep11.y = 0.6;
particleInput0_edep11.z = -28.6;
particleInput0_edep11.t = 2.10062;
particleInput0_edep11.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep11));
supera::EDep particleInput0_edep12;
particleInput0_edep12.x = 0.6;
particleInput0_edep12.y = 0.6;
particleInput0_edep12.z = -29;
particleInput0_edep12.t = 2.11457;
particleInput0_edep12.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep12));
supera::EDep particleInput0_edep13;
particleInput0_edep13.x = 0.6;
particleInput0_edep13.y = 0.6;
particleInput0_edep13.z = -29.4;
particleInput0_edep13.t = 2.12852;
particleInput0_edep13.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep13));
supera::EDep particleInput0_edep14;
particleInput0_edep14.x = 0.6;
particleInput0_edep14.y = 0.6;
particleInput0_edep14.z = -29.8;
particleInput0_edep14.t = 2.14247;
particleInput0_edep14.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep14));
supera::EDep particleInput0_edep15;
particleInput0_edep15.x = 0.6;
particleInput0_edep15.y = 0.6;
particleInput0_edep15.z = -30.2;
particleInput0_edep15.t = 2.15642;
particleInput0_edep15.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep15));
supera::EDep particleInput0_edep16;
particleInput0_edep16.x = 0.6;
particleInput0_edep16.y = 0.6;
particleInput0_edep16.z = -30.6;
particleInput0_edep16.t = 2.17038;
particleInput0_edep16.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep16));
supera::EDep particleInput0_edep17;
particleInput0_edep17.x = 0.6;
particleInput0_edep17.y = 0.6;
particleInput0_edep17.z = -31;
particleInput0_edep17.t = 2.18433;
particleInput0_edep17.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep17));
supera::EDep particleInput0_edep18;
particleInput0_edep18.x = 0.6;
particleInput0_edep18.y = 0.6;
particleInput0_edep18.z = -31.4;
particleInput0_edep18.t = 2.19828;
particleInput0_edep18.e = 0.283858;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep18));
supera::EDep particleInput0_edep19;
particleInput0_edep19.x = 0.2;
particleInput0_edep19.y = 0.6;
particleInput0_edep19.z = -31.4;
particleInput0_edep19.t = 2.21223;
particleInput0_edep19.e = 0.593697;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep19));
supera::EDep particleInput0_edep20;
particleInput0_edep20.x = 0.2;
particleInput0_edep20.y = 0.6;
particleInput0_edep20.z = -31.8;
particleInput0_edep20.t = 2.22618;
particleInput0_edep20.e = 0.877555;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep20));
supera::EDep particleInput0_edep21;
particleInput0_edep21.x = 0.2;
particleInput0_edep21.y = 0.6;
particleInput0_edep21.z = -32.2;
particleInput0_edep21.t = 2.24014;
particleInput0_edep21.e = 0.414264;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep21));
supera::EDep particleInput0_edep22;
particleInput0_edep22.x = 0.2;
particleInput0_edep22.y = 0.6;
particleInput0_edep22.z = -32.2;
particleInput0_edep22.t = 2.25409;
particleInput0_edep22.e = 0.459219;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep22));
supera::EDep particleInput0_edep23;
particleInput0_edep23.x = 0.2;
particleInput0_edep23.y = 0.6;
particleInput0_edep23.z = -32.6;
particleInput0_edep23.t = 2.25897;
particleInput0_edep23.e = 0.0666798;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep23));
supera::EDep particleInput0_edep24;
particleInput0_edep24.x = 0.2;
particleInput0_edep24.y = 0.6;
particleInput0_edep24.z = -32.6;
particleInput0_edep24.t = 2.26384;
particleInput0_edep24.e = 0.824885;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep24));
supera::EDep particleInput0_edep25;
particleInput0_edep25.x = 0.2;
particleInput0_edep25.y = 0.6;
particleInput0_edep25.z = -33;
particleInput0_edep25.t = 2.27507;
particleInput0_edep25.e = 0.103992;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep25));
supera::EDep particleInput0_edep26;
particleInput0_edep26.x = 0.6;
particleInput0_edep26.y = 0.6;
particleInput0_edep26.z = -33;
particleInput0_edep26.t = 2.28629;
particleInput0_edep26.e = 0.789377;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep26));
supera::EDep particleInput0_edep27;
particleInput0_edep27.x = 0.6;
particleInput0_edep27.y = 0.6;
particleInput0_edep27.z = -33.4;
particleInput0_edep27.t = 2.29751;
particleInput0_edep27.e = 0.893369;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep27));
supera::EDep particleInput0_edep28;
particleInput0_edep28.x = 0.6;
particleInput0_edep28.y = 0.6;
particleInput0_edep28.z = -33.8;
particleInput0_edep28.t = 2.30874;
particleInput0_edep28.e = 0.893369;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep28));
supera::EDep particleInput0_edep29;
particleInput0_edep29.x = 0.6;
particleInput0_edep29.y = 0.6;
particleInput0_edep29.z = -34.2;
particleInput0_edep29.t = 2.31996;
particleInput0_edep29.e = 0.0672888;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep29));
supera::EDep particleInput0_edep30;
particleInput0_edep30.x = 0.6;
particleInput0_edep30.y = 1;
particleInput0_edep30.z = -34.2;
particleInput0_edep30.t = 2.33118;
particleInput0_edep30.e = 0.82608;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep30));
supera::EDep particleInput0_edep31;
particleInput0_edep31.x = 0.6;
particleInput0_edep31.y = 1;
particleInput0_edep31.z = -34.6;
particleInput0_edep31.t = 2.3424;
particleInput0_edep31.e = 0.893369;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep31));
supera::EDep particleInput0_edep32;
particleInput0_edep32.x = 0.6;
particleInput0_edep32.y = 1;
particleInput0_edep32.z = -35;
particleInput0_edep32.t = 2.35363;
particleInput0_edep32.e = 0.893369;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep32));
supera::EDep particleInput0_edep33;
particleInput0_edep33.x = 0.6;
particleInput0_edep33.y = 1;
particleInput0_edep33.z = -35.4;
particleInput0_edep33.t = 2.36485;
particleInput0_edep33.e = 0.00333362;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep33));
supera::EDep particleInput0_edep34;
particleInput0_edep34.x = 0.6;
particleInput0_edep34.y = 1;
particleInput0_edep34.z = -35.4;
particleInput0_edep34.t = 2.37607;
particleInput0_edep34.e = 0.869494;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep34));
supera::EDep particleInput0_edep35;
particleInput0_edep35.x = 0.6;
particleInput0_edep35.y = 1;
particleInput0_edep35.z = -35.4;
particleInput0_edep35.t = 2.38999;
particleInput0_edep35.e = 0.140379;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep35));
supera::EDep particleInput0_edep36;
particleInput0_edep36.x = 0.6;
particleInput0_edep36.y = 1;
particleInput0_edep36.z = -35.8;
particleInput0_edep36.t = 2.40362;
particleInput0_edep36.e = 0.918773;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep36));
supera::EDep particleInput0_edep37;
particleInput0_edep37.x = 0.6;
particleInput0_edep37.y = 1;
particleInput0_edep37.z = -36.2;
particleInput0_edep37.t = 2.41725;
particleInput0_edep37.e = 0.918773;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep37));
supera::EDep particleInput0_edep38;
particleInput0_edep38.x = 0.6;
particleInput0_edep38.y = 1;
particleInput0_edep38.z = -36.6;
particleInput0_edep38.t = 2.43088;
particleInput0_edep38.e = 0.918773;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep38));
supera::EDep particleInput0_edep39;
particleInput0_edep39.x = 0.6;
particleInput0_edep39.y = 1;
particleInput0_edep39.z = -37;
particleInput0_edep39.t = 2.44452;
particleInput0_edep39.e = 0.918773;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep39));
supera::EDep particleInput0_edep40;
particleInput0_edep40.x = 0.6;
particleInput0_edep40.y = 1;
particleInput0_edep40.z = -37.4;
particleInput0_edep40.t = 2.45815;
particleInput0_edep40.e = 0.717649;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep40));
supera::EDep particleInput0_edep41;
particleInput0_edep41.x = 0.6;
particleInput0_edep41.y = 1;
particleInput0_edep41.z = -37.4;
particleInput0_edep41.t = 2.47178;
particleInput0_edep41.e = 0.228553;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep41));
supera::EDep particleInput0_edep42;
particleInput0_edep42.x = 0.6;
particleInput0_edep42.y = 1;
particleInput0_edep42.z = -37.8;
particleInput0_edep42.t = 2.48437;
particleInput0_edep42.e = 1.04407;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep42));
supera::EDep particleInput0_edep43;
particleInput0_edep43.x = 0.6;
particleInput0_edep43.y = 1;
particleInput0_edep43.z = -38.2;
particleInput0_edep43.t = 2.49696;
particleInput0_edep43.e = 1.04407;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep43));
supera::EDep particleInput0_edep44;
particleInput0_edep44.x = 0.6;
particleInput0_edep44.y = 1;
particleInput0_edep44.z = -38.6;
particleInput0_edep44.t = 2.50955;
particleInput0_edep44.e = 1.04407;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep44));
supera::EDep particleInput0_edep45;
particleInput0_edep45.x = 0.6;
particleInput0_edep45.y = 1;
particleInput0_edep45.z = -39;
particleInput0_edep45.t = 2.52214;
particleInput0_edep45.e = 0.544715;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep45));
supera::EDep particleInput0_edep46;
particleInput0_edep46.x = 0.6;
particleInput0_edep46.y = 1;
particleInput0_edep46.z = -39;
particleInput0_edep46.t = 2.53473;
particleInput0_edep46.e = 0.498381;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep46));
supera::EDep particleInput0_edep47;
particleInput0_edep47.x = 0.6;
particleInput0_edep47.y = 1;
particleInput0_edep47.z = -39.4;
particleInput0_edep47.t = 2.55013;
particleInput0_edep47.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep47));
supera::EDep particleInput0_edep48;
particleInput0_edep48.x = 0.6;
particleInput0_edep48.y = 1;
particleInput0_edep48.z = -39.8;
particleInput0_edep48.t = 2.56553;
particleInput0_edep48.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep48));
supera::EDep particleInput0_edep49;
particleInput0_edep49.x = 0.6;
particleInput0_edep49.y = 1;
particleInput0_edep49.z = -40.2;
particleInput0_edep49.t = 2.58093;
particleInput0_edep49.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep49));
supera::EDep particleInput0_edep50;
particleInput0_edep50.x = 0.6;
particleInput0_edep50.y = 1;
particleInput0_edep50.z = -40.6;
particleInput0_edep50.t = 2.59633;
particleInput0_edep50.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep50));
supera::EDep particleInput0_edep51;
particleInput0_edep51.x = 0.6;
particleInput0_edep51.y = 1;
particleInput0_edep51.z = -41;
particleInput0_edep51.t = 2.61173;
particleInput0_edep51.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep51));
supera::EDep particleInput0_edep52;
particleInput0_edep52.x = 0.6;
particleInput0_edep52.y = 1;
particleInput0_edep52.z = -41.4;
particleInput0_edep52.t = 2.62713;
particleInput0_edep52.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep52));
supera::EDep particleInput0_edep53;
particleInput0_edep53.x = 0.6;
particleInput0_edep53.y = 1;
particleInput0_edep53.z = -41.8;
particleInput0_edep53.t = 2.64253;
particleInput0_edep53.e = 1.04203;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep53));
supera::EDep particleInput0_edep54;
particleInput0_edep54.x = 0.6;
particleInput0_edep54.y = 1;
particleInput0_edep54.z = -42.2;
particleInput0_edep54.t = 2.65793;
particleInput0_edep54.e = 0.65184;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep54));
supera::EDep particleInput0_edep55;
particleInput0_edep55.x = 0.6;
particleInput0_edep55.y = 1;
particleInput0_edep55.z = -42.2;
particleInput0_edep55.t = 2.67333;
particleInput0_edep55.e = 0.419798;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep55));
supera::EDep particleInput0_edep56;
particleInput0_edep56.x = 0.6;
particleInput0_edep56.y = 1;
particleInput0_edep56.z = -42.6;
particleInput0_edep56.t = 2.68472;
particleInput0_edep56.e = 1.1211;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep56));
supera::EDep particleInput0_edep57;
particleInput0_edep57.x = 0.6;
particleInput0_edep57.y = 1;
particleInput0_edep57.z = -43;
particleInput0_edep57.t = 2.69612;
particleInput0_edep57.e = 1.1211;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep57));
supera::EDep particleInput0_edep58;
particleInput0_edep58.x = 0.6;
particleInput0_edep58.y = 1;
particleInput0_edep58.z = -43.4;
particleInput0_edep58.t = 2.70751;
particleInput0_edep58.e = 0.218871;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep58));
supera::EDep particleInput0_edep59;
particleInput0_edep59.x = 0.6;
particleInput0_edep59.y = 1;
particleInput0_edep59.z = -43.4;
particleInput0_edep59.t = 2.7189;
particleInput0_edep59.e = 0.903237;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep59));
supera::EDep particleInput0_edep60;
particleInput0_edep60.x = 0.6;
particleInput0_edep60.y = 1;
particleInput0_edep60.z = -43.8;
particleInput0_edep60.t = 2.73211;
particleInput0_edep60.e = 1.12235;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep60));
supera::EDep particleInput0_edep61;
particleInput0_edep61.x = 0.6;
particleInput0_edep61.y = 1;
particleInput0_edep61.z = -44.2;
particleInput0_edep61.t = 2.74531;
particleInput0_edep61.e = 1.12235;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep61));
supera::EDep particleInput0_edep62;
particleInput0_edep62.x = 0.6;
particleInput0_edep62.y = 1;
particleInput0_edep62.z = -44.6;
particleInput0_edep62.t = 2.75851;
particleInput0_edep62.e = 1.12235;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep62));
supera::EDep particleInput0_edep63;
particleInput0_edep63.x = 0.6;
particleInput0_edep63.y = 1;
particleInput0_edep63.z = -45;
particleInput0_edep63.t = 2.77171;
particleInput0_edep63.e = 1.12235;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep63));
supera::EDep particleInput0_edep64;
particleInput0_edep64.x = 0.6;
particleInput0_edep64.y = 1;
particleInput0_edep64.z = -45.4;
particleInput0_edep64.t = 2.78492;
particleInput0_edep64.e = 0.2676;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep64));
supera::EDep particleInput0_edep65;
particleInput0_edep65.x = 0.2;
particleInput0_edep65.y = 1;
particleInput0_edep65.z = -45.4;
particleInput0_edep65.t = 2.79812;
particleInput0_edep65.e = 0.854753;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep65));
supera::EDep particleInput0_edep66;
particleInput0_edep66.x = 0.2;
particleInput0_edep66.y = 1;
particleInput0_edep66.z = -45.8;
particleInput0_edep66.t = 2.81132;
particleInput0_edep66.e = 0.0611669;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep66));
supera::EDep particleInput0_edep67;
particleInput0_edep67.x = 0.2;
particleInput0_edep67.y = 1;
particleInput0_edep67.z = -45.8;
particleInput0_edep67.t = 2.82452;
particleInput0_edep67.e = 1.28088;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep67));
supera::EDep particleInput0_edep68;
particleInput0_edep68.x = 0.2;
particleInput0_edep68.y = 1;
particleInput0_edep68.z = -46.2;
particleInput0_edep68.t = 2.84155;
particleInput0_edep68.e = 1.35471;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep68));
supera::EDep particleInput0_edep69;
particleInput0_edep69.x = 0.2;
particleInput0_edep69.y = 1;
particleInput0_edep69.z = -46.6;
particleInput0_edep69.t = 2.85858;
particleInput0_edep69.e = 1.35471;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep69));
supera::EDep particleInput0_edep70;
particleInput0_edep70.x = 0.2;
particleInput0_edep70.y = 1;
particleInput0_edep70.z = -47;
particleInput0_edep70.t = 2.87561;
particleInput0_edep70.e = 1.35471;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep70));
supera::EDep particleInput0_edep71;
particleInput0_edep71.x = 0.2;
particleInput0_edep71.y = 1;
particleInput0_edep71.z = -47.4;
particleInput0_edep71.t = 2.89264;
particleInput0_edep71.e = 0.810099;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep71));
supera::EDep particleInput0_edep72;
particleInput0_edep72.x = 0.2;
particleInput0_edep72.y = 1;
particleInput0_edep72.z = -47.4;
particleInput0_edep72.t = 2.90967;
particleInput0_edep72.e = 0.508074;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep72));
supera::EDep particleInput0_edep73;
particleInput0_edep73.x = 0.2;
particleInput0_edep73.y = 1;
particleInput0_edep73.z = -47.8;
particleInput0_edep73.t = 2.92397;
particleInput0_edep73.e = 1.26383;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep73));
supera::EDep particleInput0_edep74;
particleInput0_edep74.x = 0.2;
particleInput0_edep74.y = 1;
particleInput0_edep74.z = -48.2;
particleInput0_edep74.t = 2.93827;
particleInput0_edep74.e = 1.26383;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep74));
supera::EDep particleInput0_edep75;
particleInput0_edep75.x = 0.2;
particleInput0_edep75.y = 1;
particleInput0_edep75.z = -48.6;
particleInput0_edep75.t = 2.95257;
particleInput0_edep75.e = 1.26383;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep75));
supera::EDep particleInput0_edep76;
particleInput0_edep76.x = 0.2;
particleInput0_edep76.y = 1;
particleInput0_edep76.z = -49;
particleInput0_edep76.t = 2.96687;
particleInput0_edep76.e = 0.283419;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep76));
supera::EDep particleInput0_edep77;
particleInput0_edep77.x = 0.2;
particleInput0_edep77.y = 1;
particleInput0_edep77.z = -49;
particleInput0_edep77.t = 2.98117;
particleInput0_edep77.e = 1.12244;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep77));
supera::EDep particleInput0_edep78;
particleInput0_edep78.x = 0.2;
particleInput0_edep78.y = 1;
particleInput0_edep78.z = -49.4;
particleInput0_edep78.t = 2.99627;
particleInput0_edep78.e = 1.44692;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep78));
supera::EDep particleInput0_edep79;
particleInput0_edep79.x = 0.2;
particleInput0_edep79.y = 1;
particleInput0_edep79.z = -49.8;
particleInput0_edep79.t = 3.01137;
particleInput0_edep79.e = 1.44692;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep79));
supera::EDep particleInput0_edep80;
particleInput0_edep80.x = 0.2;
particleInput0_edep80.y = 1;
particleInput0_edep80.z = -50.2;
particleInput0_edep80.t = 3.02647;
particleInput0_edep80.e = 0.219572;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep80));
supera::EDep particleInput0_edep81;
particleInput0_edep81.x = 0.2;
particleInput0_edep81.y = 1;
particleInput0_edep81.z = -50.2;
particleInput0_edep81.t = 3.04157;
particleInput0_edep81.e = 1.36006;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep81));
supera::EDep particleInput0_edep82;
particleInput0_edep82.x = 0.2;
particleInput0_edep82.y = 1;
particleInput0_edep82.z = -50.6;
particleInput0_edep82.t = 3.05422;
particleInput0_edep82.e = 1.40495;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep82));
supera::EDep particleInput0_edep83;
particleInput0_edep83.x = 0.2;
particleInput0_edep83.y = 1.4;
particleInput0_edep83.z = -50.6;
particleInput0_edep83.t = 3.06686;
particleInput0_edep83.e = 0.198421;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep83));
supera::EDep particleInput0_edep84;
particleInput0_edep84.x = 0.2;
particleInput0_edep84.y = 1.4;
particleInput0_edep84.z = -51;
particleInput0_edep84.t = 3.07951;
particleInput0_edep84.e = 0.805664;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep84));
supera::EDep particleInput0_edep85;
particleInput0_edep85.x = 0.2;
particleInput0_edep85.y = 1.4;
particleInput0_edep85.z = -51;
particleInput0_edep85.t = 3.09216;
particleInput0_edep85.e = 0.893761;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep85));
supera::EDep particleInput0_edep86;
particleInput0_edep86.x = 0.2;
particleInput0_edep86.y = 1.4;
particleInput0_edep86.z = -51.4;
particleInput0_edep86.t = 3.10625;
particleInput0_edep86.e = 1.79643;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep86));
supera::EDep particleInput0_edep87;
particleInput0_edep87.x = 0.2;
particleInput0_edep87.y = 1.4;
particleInput0_edep87.z = -51.8;
particleInput0_edep87.t = 3.12034;
particleInput0_edep87.e = 0.652922;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep87));
supera::EDep particleInput0_edep88;
particleInput0_edep88.x = 0.2;
particleInput0_edep88.y = 1.4;
particleInput0_edep88.z = -51.8;
particleInput0_edep88.t = 3.13443;
particleInput0_edep88.e = 0.323813;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep88));
supera::EDep particleInput0_edep89;
particleInput0_edep89.x = 0.2;
particleInput0_edep89.y = 1.4;
particleInput0_edep89.z = -51.8;
particleInput0_edep89.t = 3.13813;
particleInput0_edep89.e = 0.852992;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep89));
supera::EDep particleInput0_edep90;
particleInput0_edep90.x = 0.2;
particleInput0_edep90.y = 1.4;
particleInput0_edep90.z = -52.2;
particleInput0_edep90.t = 3.155;
particleInput0_edep90.e = 1.63652;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep90));
supera::EDep particleInput0_edep91;
particleInput0_edep91.x = 0.2;
particleInput0_edep91.y = 1.4;
particleInput0_edep91.z = -52.2;
particleInput0_edep91.t = 3.17188;
particleInput0_edep91.e = 0.155056;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep91));
supera::EDep particleInput0_edep92;
particleInput0_edep92.x = 0.2;
particleInput0_edep92.y = 1.4;
particleInput0_edep92.z = -52.6;
particleInput0_edep92.t = 3.18149;
particleInput0_edep92.e = 2.01283;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep92));
supera::EDep particleInput0_edep93;
particleInput0_edep93.x = 0.2;
particleInput0_edep93.y = 1.4;
particleInput0_edep93.z = -53;
particleInput0_edep93.t = 3.19111;
particleInput0_edep93.e = 0.150277;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep93));
supera::EDep particleInput0_edep94;
particleInput0_edep94.x = 0.2;
particleInput0_edep94.y = 1.4;
particleInput0_edep94.z = -53;
particleInput0_edep94.t = 3.20073;
particleInput0_edep94.e = 1.79246;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep94));
supera::EDep particleInput0_edep95;
particleInput0_edep95.x = 0.2;
particleInput0_edep95.y = 1.4;
particleInput0_edep95.z = -53.4;
particleInput0_edep95.t = 3.21299;
particleInput0_edep95.e = 0.00952561;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep95));
supera::EDep particleInput0_edep96;
particleInput0_edep96.x = 0.2;
particleInput0_edep96.y = 1.4;
particleInput0_edep96.z = -53.4;
particleInput0_edep96.t = 3.22524;
particleInput0_edep96.e = 2.3022;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep96));
supera::EDep particleInput0_edep97;
particleInput0_edep97.x = 0.2;
particleInput0_edep97.y = 1.4;
particleInput0_edep97.z = -53.4;
particleInput0_edep97.t = 3.24654;
particleInput0_edep97.e = 0.568874;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep97));
supera::EDep particleInput0_edep98;
particleInput0_edep98.x = 0.2;
particleInput0_edep98.y = 1.4;
particleInput0_edep98.z = -53.8;
particleInput0_edep98.t = 3.25523;
particleInput0_edep98.e = 0.933877;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep98));
supera::EDep particleInput0_edep99;
particleInput0_edep99.x = 0.2;
particleInput0_edep99.y = 1.4;
particleInput0_edep99.z = -53.8;
particleInput0_edep99.t = 3.26393;
particleInput0_edep99.e = 1.54622;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep99));
supera::EDep particleInput0_edep100;
particleInput0_edep100.x = 0.2;
particleInput0_edep100.y = 1.4;
particleInput0_edep100.z = -53.8;
particleInput0_edep100.t = 3.2789;
particleInput0_edep100.e = 0.471395;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep100));
supera::EDep particleInput0_edep101;
particleInput0_edep101.x = 0.2;
particleInput0_edep101.y = 1.4;
particleInput0_edep101.z = -54.2;
particleInput0_edep101.t = 3.28521;
particleInput0_edep101.e = 0.655071;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep101));
supera::EDep particleInput0_edep102;
particleInput0_edep102.x = 0.2;
particleInput0_edep102.y = 1.4;
particleInput0_edep102.z = -54.2;
particleInput0_edep102.t = 3.29152;
particleInput0_edep102.e = 1.17253;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep102));
supera::EDep particleInput0_edep103;
particleInput0_edep103.x = 0.2;
particleInput0_edep103.y = 1.4;
particleInput0_edep103.z = -54.2;
particleInput0_edep103.t = 3.30252;
particleInput0_edep103.e = 0.923561;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep103));
supera::EDep particleInput0_edep104;
particleInput0_edep104.x = 0.2;
particleInput0_edep104.y = 1.4;
particleInput0_edep104.z = -54.2;
particleInput0_edep104.t = 3.31192;
particleInput0_edep104.e = 1.02956;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep104));
supera::EDep particleInput0_edep105;
particleInput0_edep105.x = 0.2;
particleInput0_edep105.y = 1.4;
particleInput0_edep105.z = -54.2;
particleInput0_edep105.t = 3.32013;
particleInput0_edep105.e = 0.29446;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep105));
supera::EDep particleInput0_edep106;
particleInput0_edep106.x = 0.2;
particleInput0_edep106.y = 1.4;
particleInput0_edep106.z = -54.6;
particleInput0_edep106.t = 3.32363;
particleInput0_edep106.e = 0.710265;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep106));
supera::EDep particleInput0_edep107;
particleInput0_edep107.x = 0.2;
particleInput0_edep107.y = 1.4;
particleInput0_edep107.z = -54.6;
particleInput0_edep107.t = 3.32712;
particleInput0_edep107.e = 0.660482;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep107));
supera::EDep particleInput0_edep108;
particleInput0_edep108.x = 0.2;
particleInput0_edep108.y = 1.4;
particleInput0_edep108.z = -54.6;
particleInput0_edep108.t = 3.33302;
particleInput0_edep108.e = 0.748718;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep108));
supera::EDep particleInput0_edep109;
particleInput0_edep109.x = 0.2;
particleInput0_edep109.y = 1.4;
particleInput0_edep109.z = -54.6;
particleInput0_edep109.t = 3.33827;
particleInput0_edep109.e = 0.747764;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep109));
supera::EDep particleInput0_edep110;
particleInput0_edep110.x = 0.2;
particleInput0_edep110.y = 1.4;
particleInput0_edep110.z = -54.6;
particleInput0_edep110.t = 3.34288;
particleInput0_edep110.e = 1.52848;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep110));
supera::EDep particleInput0_edep111;
particleInput0_edep111.x = 0.2;
particleInput0_edep111.y = 1.4;
particleInput0_edep111.z = -54.6;
particleInput0_edep111.t = 3.35049;
particleInput0_edep111.e = 0.79379;

particleInput0.pcloud.emplace_back(std::move(particleInput0_edep111));
particleInput0.valid = 1;
particleInput0.type = static_cast<supera::ProcessType>(0);
supera::Particle particleInput0_particle;
particleInput0_particle.id = 1;
particleInput0_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput0_particle.trackid = 0;
particleInput0_particle.pdg = -13;
particleInput0_particle.px = 4.10617;
particleInput0_particle.py = 3.74819;
particleInput0_particle.pz = -176.354;
particleInput0_particle.vtx = {0, 0, 0, 1};
particleInput0_particle.end_pt = {0.14755, 1.49633, -54.5587, 10.8502};
particleInput0_particle.first_step = {0, 0, 0, 1};
particleInput0_particle.last_step = {0.14755, 1.49633, -54.5587, 10.8502};
particleInput0_particle.dist_travel = 54.6454;
particleInput0_particle.energy_init = 205.658;
particleInput0_particle.energy_deposit = 94.0345;
particleInput0_particle.process = "primary";
particleInput0_particle.parent_trackid = kINVALID_TRACKID;
particleInput0_particle.parent_pdg = 0;
particleInput0_particle.parent_vtx = {0, 0, 0, 0};
particleInput0_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput0_particle.ancestor_pdg = 0;
particleInput0_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput0_particle.ancestor_process = "";
particleInput0_particle.parent_process = "";
particleInput0_particle.parent_id = supera::kINVALID_INSTANCEID;
particleInput0_particle.children_id = {  };
particleInput0_particle.group_id = kINVALID_INSTANCEID;
particleInput0_particle.interaction_id = 0;
particleInput0.part = std::move(particleInput0_particle);

evInput.push_back(std::move(particleInput0));
supera::ParticleInput particleInput1;
particleInput1.pcloud.reserve(1);
supera::EDep particleInput1_edep0;
particleInput1_edep0.x = 0.6;
particleInput1_edep0.y = 0.6;
particleInput1_edep0.z = -25.8;
particleInput1_edep0.t = 2.01024;
particleInput1_edep0.e = 0.439832;

particleInput1.pcloud.emplace_back(std::move(particleInput1_edep0));
particleInput1.valid = 1;
particleInput1.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput1_particle;
particleInput1_particle.id = 2;
particleInput1_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput1_particle.trackid = 1;
particleInput1_particle.pdg = 11;
particleInput1_particle.px = 0.402134;
particleInput1_particle.py = -0.449941;
particleInput1_particle.pz = -0.528015;
particleInput1_particle.vtx = {0.603066, 0.602993, -25.9686, 2.01024};
particleInput1_particle.end_pt = {0.624791, 0.578685, -25.9971, 2.01195};
particleInput1_particle.first_step = {0.603066, 0.602993, -25.9686, 2.01024};
particleInput1_particle.last_step = {0.624791, 0.578685, -25.9971, 2.01195};
particleInput1_particle.dist_travel = 0.0433199;
particleInput1_particle.energy_init = 0.95083;
particleInput1_particle.energy_deposit = 0.439832;
particleInput1_particle.process = "muIoni";
particleInput1_particle.parent_trackid = 0;
particleInput1_particle.parent_pdg = -13;
particleInput1_particle.parent_vtx = {0, 0, 0, 0};
particleInput1_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput1_particle.ancestor_pdg = 0;
particleInput1_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput1_particle.ancestor_process = "";
particleInput1_particle.parent_process = "";
particleInput1_particle.parent_id = 1;
particleInput1_particle.children_id = {  };
particleInput1_particle.group_id = kINVALID_INSTANCEID;
particleInput1_particle.interaction_id = 0;
particleInput1.part = std::move(particleInput1_particle);

evInput.push_back(std::move(particleInput1));
supera::ParticleInput particleInput2;
particleInput2.pcloud.reserve(2);
supera::EDep particleInput2_edep0;
particleInput2_edep0.x = 0.6;
particleInput2_edep0.y = 0.6;
particleInput2_edep0.z = -26.2;
particleInput2_edep0.t = 2.0169;
particleInput2_edep0.e = 0.351733;

particleInput2.pcloud.emplace_back(std::move(particleInput2_edep0));
supera::EDep particleInput2_edep1;
particleInput2_edep1.x = 0.6;
particleInput2_edep1.y = 0.6;
particleInput2_edep1.z = -26.2;
particleInput2_edep1.t = 2.02032;
particleInput2_edep1.e = 0.423428;

particleInput2.pcloud.emplace_back(std::move(particleInput2_edep1));
particleInput2.valid = 1;
particleInput2.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput2_particle;
particleInput2_particle.id = 3;
particleInput2_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput2_particle.trackid = 2;
particleInput2_particle.pdg = 11;
particleInput2_particle.px = 0.11796;
particleInput2_particle.py = -0.715356;
particleInput2_particle.pz = -0.931364;
particleInput2_particle.vtx = {0.60061, 0.606036, -26.139, 2.0169};
particleInput2_particle.end_pt = {0.627759, 0.559361, -26.2481, 2.02201};
particleInput2_particle.first_step = {0.60061, 0.606036, -26.139, 2.0169};
particleInput2_particle.last_step = {0.627759, 0.559361, -26.2481, 2.02201};
particleInput2_particle.dist_travel = 0.121774;
particleInput2_particle.energy_init = 1.28616;
particleInput2_particle.energy_deposit = 0.775161;
particleInput2_particle.process = "muIoni";
particleInput2_particle.parent_trackid = 0;
particleInput2_particle.parent_pdg = -13;
particleInput2_particle.parent_vtx = {0, 0, 0, 0};
particleInput2_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput2_particle.ancestor_pdg = 0;
particleInput2_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput2_particle.ancestor_process = "";
particleInput2_particle.parent_process = "";
particleInput2_particle.parent_id = 1;
particleInput2_particle.children_id = {  };
particleInput2_particle.group_id = kINVALID_INSTANCEID;
particleInput2_particle.interaction_id = 0;
particleInput2.part = std::move(particleInput2_particle);

evInput.push_back(std::move(particleInput2));
supera::ParticleInput particleInput3;
particleInput3.pcloud.reserve(4);
supera::EDep particleInput3_edep0;
particleInput3_edep0.x = 0.2;
particleInput3_edep0.y = 0.6;
particleInput3_edep0.z = -32.6;
particleInput3_edep0.t = 2.26384;
particleInput3_edep0.e = 0.379489;

particleInput3.pcloud.emplace_back(std::move(particleInput3_edep0));
supera::EDep particleInput3_edep1;
particleInput3_edep1.x = 0.2;
particleInput3_edep1.y = 0.6;
particleInput3_edep1.z = -32.6;
particleInput3_edep1.t = 2.26804;
particleInput3_edep1.e = 0.185989;

particleInput3.pcloud.emplace_back(std::move(particleInput3_edep1));
supera::EDep particleInput3_edep2;
particleInput3_edep2.x = 0.2;
particleInput3_edep2.y = 1;
particleInput3_edep2.z = -32.6;
particleInput3_edep2.t = 2.26919;
particleInput3_edep2.e = 0.277691;

particleInput3.pcloud.emplace_back(std::move(particleInput3_edep2));
supera::EDep particleInput3_edep3;
particleInput3_edep3.x = 0.2;
particleInput3_edep3.y = 1;
particleInput3_edep3.z = -32.6;
particleInput3_edep3.t = 2.27034;
particleInput3_edep3.e = 0.0918982;

particleInput3.pcloud.emplace_back(std::move(particleInput3_edep3));
particleInput3.valid = 1;
particleInput3.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput3_particle;
particleInput3_particle.id = 4;
particleInput3_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput3_particle.trackid = 3;
particleInput3_particle.pdg = 11;
particleInput3_particle.px = -0.654473;
particleInput3_particle.py = 0.231597;
particleInput3_particle.pz = -1.16104;
particleInput3_particle.vtx = {0.378264, 0.762396, -32.4307, 2.26384};
particleInput3_particle.end_pt = {0.341846, 0.824651, -32.5702, 2.27055};
particleInput3_particle.first_step = {0.378264, 0.762396, -32.4307, 2.26384};
particleInput3_particle.last_step = {0.341846, 0.824651, -32.5702, 2.27055};
particleInput3_particle.dist_travel = 0.157116;
particleInput3_particle.energy_init = 1.44607;
particleInput3_particle.energy_deposit = 0.935067;
particleInput3_particle.process = "muIoni";
particleInput3_particle.parent_trackid = 0;
particleInput3_particle.parent_pdg = -13;
particleInput3_particle.parent_vtx = {0, 0, 0, 0};
particleInput3_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput3_particle.ancestor_pdg = 0;
particleInput3_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput3_particle.ancestor_process = "";
particleInput3_particle.parent_process = "";
particleInput3_particle.parent_id = 1;
particleInput3_particle.children_id = {  };
particleInput3_particle.group_id = kINVALID_INSTANCEID;
particleInput3_particle.interaction_id = 0;
particleInput3.part = std::move(particleInput3_particle);

evInput.push_back(std::move(particleInput3));
supera::ParticleInput particleInput4;
particleInput4.pcloud.reserve(1);
supera::EDep particleInput4_edep0;
particleInput4_edep0.x = 0.6;
particleInput4_edep0.y = 1;
particleInput4_edep0.z = -35.4;
particleInput4_edep0.t = 2.37607;
particleInput4_edep0.e = 0.413441;

particleInput4.pcloud.emplace_back(std::move(particleInput4_edep0));
particleInput4.valid = 1;
particleInput4.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput4_particle;
particleInput4_particle.id = 5;
particleInput4_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput4_particle.trackid = 4;
particleInput4_particle.pdg = 11;
particleInput4_particle.px = 0.143948;
particleInput4_particle.py = 0.58689;
particleInput4_particle.pz = -0.477817;
particleInput4_particle.vtx = {0.523069, 0.827538, -35.2015, 2.37607};
particleInput4_particle.end_pt = {0.529273, 0.86147, -35.2268, 2.3779};
particleInput4_particle.first_step = {0.523069, 0.827538, -35.2015, 2.37607};
particleInput4_particle.last_step = {0.529273, 0.86147, -35.2268, 2.3779};
particleInput4_particle.dist_travel = 0.0427534;
particleInput4_particle.energy_init = 0.92444;
particleInput4_particle.energy_deposit = 0.413441;
particleInput4_particle.process = "muIoni";
particleInput4_particle.parent_trackid = 0;
particleInput4_particle.parent_pdg = -13;
particleInput4_particle.parent_vtx = {0, 0, 0, 0};
particleInput4_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput4_particle.ancestor_pdg = 0;
particleInput4_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput4_particle.ancestor_process = "";
particleInput4_particle.parent_process = "";
particleInput4_particle.parent_id = 1;
particleInput4_particle.children_id = {  };
particleInput4_particle.group_id = kINVALID_INSTANCEID;
particleInput4_particle.interaction_id = 0;
particleInput4.part = std::move(particleInput4_particle);

evInput.push_back(std::move(particleInput4));
supera::ParticleInput particleInput5;
particleInput5.pcloud.reserve(6);
supera::EDep particleInput5_edep0;
particleInput5_edep0.x = 0.6;
particleInput5_edep0.y = 1;
particleInput5_edep0.z = -35.4;
particleInput5_edep0.t = 2.38999;
particleInput5_edep0.e = 0.228324;

particleInput5.pcloud.emplace_back(std::move(particleInput5_edep0));
supera::EDep particleInput5_edep1;
particleInput5_edep1.x = 0.6;
particleInput5_edep1.y = 1;
particleInput5_edep1.z = -35.8;
particleInput5_edep1.t = 2.39201;
particleInput5_edep1.e = 0.131414;

particleInput5.pcloud.emplace_back(std::move(particleInput5_edep1));
supera::EDep particleInput5_edep2;
particleInput5_edep2.x = 0.6;
particleInput5_edep2.y = 1;
particleInput5_edep2.z = -35.8;
particleInput5_edep2.t = 2.39403;
particleInput5_edep2.e = 0.0669339;

particleInput5.pcloud.emplace_back(std::move(particleInput5_edep2));
supera::EDep particleInput5_edep3;
particleInput5_edep3.x = 0.6;
particleInput5_edep3.y = 0.6;
particleInput5_edep3.z = -35.8;
particleInput5_edep3.t = 2.39477;
particleInput5_edep3.e = 0.157791;

particleInput5.pcloud.emplace_back(std::move(particleInput5_edep3));
supera::EDep particleInput5_edep4;
particleInput5_edep4.x = 0.6;
particleInput5_edep4.y = 0.6;
particleInput5_edep4.z = -35.4;
particleInput5_edep4.t = 2.39552;
particleInput5_edep4.e = 0.0753152;

particleInput5.pcloud.emplace_back(std::move(particleInput5_edep4));
supera::EDep particleInput5_edep5;
particleInput5_edep5.x = 0.6;
particleInput5_edep5.y = 0.6;
particleInput5_edep5.z = -35.4;
particleInput5_edep5.t = 2.39626;
particleInput5_edep5.e = 0.243123;

particleInput5.pcloud.emplace_back(std::move(particleInput5_edep5));
particleInput5.valid = 1;
particleInput5.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput5_particle;
particleInput5_particle.id = 6;
particleInput5_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput5_particle.trackid = 5;
particleInput5_particle.pdg = 11;
particleInput5_particle.px = 0.556158;
particleInput5_particle.py = -0.40698;
particleInput5_particle.pz = -1.12386;
particleInput5_particle.vtx = {0.535816, 0.842702, -35.5389, 2.38998};
particleInput5_particle.end_pt = {0.595339, 0.771291, -35.5786, 2.39705};
particleInput5_particle.first_step = {0.535816, 0.842702, -35.5389, 2.38998};
particleInput5_particle.last_step = {0.595339, 0.771291, -35.5786, 2.39705};
particleInput5_particle.dist_travel = 0.101102;
particleInput5_particle.energy_init = 1.4139;
particleInput5_particle.energy_deposit = 0.902902;
particleInput5_particle.process = "muIoni";
particleInput5_particle.parent_trackid = 0;
particleInput5_particle.parent_pdg = -13;
particleInput5_particle.parent_vtx = {0, 0, 0, 0};
particleInput5_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput5_particle.ancestor_pdg = 0;
particleInput5_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput5_particle.ancestor_process = "";
particleInput5_particle.parent_process = "";
particleInput5_particle.parent_id = 1;
particleInput5_particle.children_id = {  };
particleInput5_particle.group_id = kINVALID_INSTANCEID;
particleInput5_particle.interaction_id = 0;
particleInput5.part = std::move(particleInput5_particle);

evInput.push_back(std::move(particleInput5));
supera::ParticleInput particleInput6;
particleInput6.pcloud.reserve(1);
supera::EDep particleInput6_edep0;
particleInput6_edep0.x = 0.6;
particleInput6_edep0.y = 1;
particleInput6_edep0.z = -37.4;
particleInput6_edep0.t = 2.47178;
particleInput6_edep0.e = 0.382866;

particleInput6.pcloud.emplace_back(std::move(particleInput6_edep0));
particleInput6.valid = 1;
particleInput6.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput6_particle;
particleInput6_particle.id = 7;
particleInput6_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput6_particle.trackid = 6;
particleInput6_particle.pdg = 11;
particleInput6_particle.px = -0.170983;
particleInput6_particle.py = 0.509824;
particleInput6_particle.pz = -0.498717;
particleInput6_particle.vtx = {0.589031, 0.940689, -37.5124, 2.47178};
particleInput6_particle.end_pt = {0.582896, 0.964219, -37.5374, 2.47338};
particleInput6_particle.first_step = {0.589031, 0.940689, -37.5124, 2.47178};
particleInput6_particle.last_step = {0.582896, 0.964219, -37.5374, 2.47338};
particleInput6_particle.dist_travel = 0.0348241;
particleInput6_particle.energy_init = 0.893865;
particleInput6_particle.energy_deposit = 0.382866;
particleInput6_particle.process = "muIoni";
particleInput6_particle.parent_trackid = 0;
particleInput6_particle.parent_pdg = -13;
particleInput6_particle.parent_vtx = {0, 0, 0, 0};
particleInput6_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput6_particle.ancestor_pdg = 0;
particleInput6_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput6_particle.ancestor_process = "";
particleInput6_particle.parent_process = "";
particleInput6_particle.parent_id = 1;
particleInput6_particle.children_id = {  };
particleInput6_particle.group_id = kINVALID_INSTANCEID;
particleInput6_particle.interaction_id = 0;
particleInput6.part = std::move(particleInput6_particle);

evInput.push_back(std::move(particleInput6));
supera::ParticleInput particleInput7;
particleInput7.pcloud.reserve(2);
supera::EDep particleInput7_edep0;
particleInput7_edep0.x = 0.6;
particleInput7_edep0.y = 1;
particleInput7_edep0.z = -39;
particleInput7_edep0.t = 2.53473;
particleInput7_edep0.e = 0.452837;

particleInput7.pcloud.emplace_back(std::move(particleInput7_edep0));
supera::EDep particleInput7_edep1;
particleInput7_edep1.x = 0.6;
particleInput7_edep1.y = 1;
particleInput7_edep1.z = -39;
particleInput7_edep1.t = 2.53776;
particleInput7_edep1.e = 0.245143;

particleInput7.pcloud.emplace_back(std::move(particleInput7_edep1));
particleInput7.valid = 1;
particleInput7.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput7_particle;
particleInput7_particle.id = 8;
particleInput7_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput7_particle.trackid = 7;
particleInput7_particle.pdg = 11;
particleInput7_particle.px = 0.138722;
particleInput7_particle.py = -0.622817;
particleInput7_particle.pz = -0.890712;
particleInput7_particle.vtx = {0.568644, 0.883198, -39.0087, 2.53473};
particleInput7_particle.end_pt = {0.587566, 0.822412, -39.0821, 2.53856};
particleInput7_particle.first_step = {0.568644, 0.883198, -39.0087, 2.53473};
particleInput7_particle.last_step = {0.587566, 0.822412, -39.0821, 2.53856};
particleInput7_particle.dist_travel = 0.0971839;
particleInput7_particle.energy_init = 1.20898;
particleInput7_particle.energy_deposit = 0.697981;
particleInput7_particle.process = "muIoni";
particleInput7_particle.parent_trackid = 0;
particleInput7_particle.parent_pdg = -13;
particleInput7_particle.parent_vtx = {0, 0, 0, 0};
particleInput7_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput7_particle.ancestor_pdg = 0;
particleInput7_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput7_particle.ancestor_process = "";
particleInput7_particle.parent_process = "";
particleInput7_particle.parent_id = 1;
particleInput7_particle.children_id = {  };
particleInput7_particle.group_id = kINVALID_INSTANCEID;
particleInput7_particle.interaction_id = 0;
particleInput7.part = std::move(particleInput7_particle);

evInput.push_back(std::move(particleInput7));
supera::ParticleInput particleInput8;
particleInput8.pcloud.reserve(1);
supera::EDep particleInput8_edep0;
particleInput8_edep0.x = 0.6;
particleInput8_edep0.y = 1;
particleInput8_edep0.z = -43.4;
particleInput8_edep0.t = 2.7189;
particleInput8_edep0.e = 0.506628;

particleInput8.pcloud.emplace_back(std::move(particleInput8_edep0));
particleInput8.valid = 1;
particleInput8.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput8_particle;
particleInput8_particle.id = 9;
particleInput8_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput8_particle.trackid = 8;
particleInput8_particle.pdg = 11;
particleInput8_particle.px = -0.078967;
particleInput8_particle.py = -0.537535;
particleInput8_particle.pz = -0.69229;
particleInput8_particle.vtx = {0.423152, 0.855246, -43.2781, 2.7189};
particleInput8_particle.end_pt = {0.416474, 0.827817, -43.3183, 2.72129};
particleInput8_particle.first_step = {0.423152, 0.855246, -43.2781, 2.7189};
particleInput8_particle.last_step = {0.416474, 0.827817, -43.3183, 2.72129};
particleInput8_particle.dist_travel = 0.049106;
particleInput8_particle.energy_init = 1.01763;
particleInput8_particle.energy_deposit = 0.506628;
particleInput8_particle.process = "muIoni";
particleInput8_particle.parent_trackid = 0;
particleInput8_particle.parent_pdg = -13;
particleInput8_particle.parent_vtx = {0, 0, 0, 0};
particleInput8_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput8_particle.ancestor_pdg = 0;
particleInput8_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput8_particle.ancestor_process = "";
particleInput8_particle.parent_process = "";
particleInput8_particle.parent_id = 1;
particleInput8_particle.children_id = {  };
particleInput8_particle.group_id = kINVALID_INSTANCEID;
particleInput8_particle.interaction_id = 0;
particleInput8.part = std::move(particleInput8_particle);

evInput.push_back(std::move(particleInput8));
supera::ParticleInput particleInput9;
particleInput9.pcloud.reserve(1);
supera::EDep particleInput9_edep0;
particleInput9_edep0.x = 0.2;
particleInput9_edep0.y = 1;
particleInput9_edep0.z = -47.4;
particleInput9_edep0.t = 2.90967;
particleInput9_edep0.e = 0.423917;

particleInput9.pcloud.emplace_back(std::move(particleInput9_edep0));
particleInput9.valid = 1;
particleInput9.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput9_particle;
particleInput9_particle.id = 10;
particleInput9_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput9_particle.trackid = 9;
particleInput9_particle.pdg = 11;
particleInput9_particle.px = 0.242912;
particleInput9_particle.py = 0.43325;
particleInput9_particle.pz = -0.605175;
particleInput9_particle.vtx = {0.327773, 0.863539, -47.4392, 2.90967};
particleInput9_particle.end_pt = {0.341062, 0.883623, -47.4636, 2.91171};
particleInput9_particle.first_step = {0.327773, 0.863539, -47.4392, 2.90967};
particleInput9_particle.last_step = {0.341062, 0.883623, -47.4636, 2.91171};
particleInput9_particle.dist_travel = 0.0342698;
particleInput9_particle.energy_init = 0.934916;
particleInput9_particle.energy_deposit = 0.423917;
particleInput9_particle.process = "muIoni";
particleInput9_particle.parent_trackid = 0;
particleInput9_particle.parent_pdg = -13;
particleInput9_particle.parent_vtx = {0, 0, 0, 0};
particleInput9_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput9_particle.ancestor_pdg = 0;
particleInput9_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput9_particle.ancestor_process = "";
particleInput9_particle.parent_process = "";
particleInput9_particle.parent_id = 1;
particleInput9_particle.children_id = {  };
particleInput9_particle.group_id = kINVALID_INSTANCEID;
particleInput9_particle.interaction_id = 0;
particleInput9.part = std::move(particleInput9_particle);

evInput.push_back(std::move(particleInput9));
supera::ParticleInput particleInput10;
particleInput10.pcloud.reserve(1);
supera::EDep particleInput10_edep0;
particleInput10_edep0.x = 0.2;
particleInput10_edep0.y = 1.4;
particleInput10_edep0.z = -51.8;
particleInput10_edep0.t = 3.13813;
particleInput10_edep0.e = 0.431848;

particleInput10.pcloud.emplace_back(std::move(particleInput10_edep0));
particleInput10.valid = 1;
particleInput10.type = static_cast<supera::ProcessType>(6);
supera::Particle particleInput10_particle;
particleInput10_particle.id = 11;
particleInput10_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput10_particle.trackid = 10;
particleInput10_particle.pdg = 11;
particleInput10_particle.px = -0.119298;
particleInput10_particle.py = -0.0551964;
particleInput10_particle.pz = -0.781385;
particleInput10_particle.vtx = {0.164775, 1.35218, -51.8076, 3.13813};
particleInput10_particle.end_pt = {0.157001, 1.35323, -51.8492, 3.14006};
particleInput10_particle.first_step = {0.164775, 1.35218, -51.8076, 3.13813};
particleInput10_particle.last_step = {0.157001, 1.35323, -51.8492, 3.14006};
particleInput10_particle.dist_travel = 0.0423668;
particleInput10_particle.energy_init = 0.942847;
particleInput10_particle.energy_deposit = 0.431848;
particleInput10_particle.process = "muIoni";
particleInput10_particle.parent_trackid = 0;
particleInput10_particle.parent_pdg = -13;
particleInput10_particle.parent_vtx = {0, 0, 0, 0};
particleInput10_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput10_particle.ancestor_pdg = 0;
particleInput10_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput10_particle.ancestor_process = "";
particleInput10_particle.parent_process = "";
particleInput10_particle.parent_id = 1;
particleInput10_particle.children_id = {  };
particleInput10_particle.group_id = kINVALID_INSTANCEID;
particleInput10_particle.interaction_id = 0;
particleInput10.part = std::move(particleInput10_particle);

evInput.push_back(std::move(particleInput10));
supera::ParticleInput particleInput11;
particleInput11.pcloud.reserve(141);
supera::EDep particleInput11_edep0;
particleInput11_edep0.x = 0.2;
particleInput11_edep0.y = 1.4;
particleInput11_edep0.z = -54.6;
particleInput11_edep0.t = 10.8502;
particleInput11_edep0.e = 1.03118;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep0));
supera::EDep particleInput11_edep1;
particleInput11_edep1.x = 0.6;
particleInput11_edep1.y = 1.4;
particleInput11_edep1.z = -54.6;
particleInput11_edep1.t = 10.8569;
particleInput11_edep1.e = 0.195663;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep1));
supera::EDep particleInput11_edep2;
particleInput11_edep2.x = 0.6;
particleInput11_edep2.y = 1.4;
particleInput11_edep2.z = -54.6;
particleInput11_edep2.t = 10.8636;
particleInput11_edep2.e = 0.031152;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep2));
supera::EDep particleInput11_edep3;
particleInput11_edep3.x = 0.6;
particleInput11_edep3.y = 1.4;
particleInput11_edep3.z = -54.2;
particleInput11_edep3.t = 10.8682;
particleInput11_edep3.e = 0.281882;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep3));
supera::EDep particleInput11_edep4;
particleInput11_edep4.x = 0.6;
particleInput11_edep4.y = 1;
particleInput11_edep4.z = -54.2;
particleInput11_edep4.t = 10.8728;
particleInput11_edep4.e = 0.637591;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep4));
supera::EDep particleInput11_edep5;
particleInput11_edep5.x = 1;
particleInput11_edep5.y = 1;
particleInput11_edep5.z = -54.2;
particleInput11_edep5.t = 10.8775;
particleInput11_edep5.e = 0.148015;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep5));
supera::EDep particleInput11_edep6;
particleInput11_edep6.x = 1;
particleInput11_edep6.y = 1;
particleInput11_edep6.z = -54.2;
particleInput11_edep6.t = 10.8821;
particleInput11_edep6.e = 0.56887;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep6));
supera::EDep particleInput11_edep7;
particleInput11_edep7.x = 1;
particleInput11_edep7.y = 1;
particleInput11_edep7.z = -53.8;
particleInput11_edep7.t = 10.8898;
particleInput11_edep7.e = 0.206413;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep7));
supera::EDep particleInput11_edep8;
particleInput11_edep8.x = 1.4;
particleInput11_edep8.y = 1;
particleInput11_edep8.z = -53.8;
particleInput11_edep8.t = 10.8974;
particleInput11_edep8.e = 0.0382903;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep8));
supera::EDep particleInput11_edep9;
particleInput11_edep9.x = 1.4;
particleInput11_edep9.y = 0.6;
particleInput11_edep9.z = -53.8;
particleInput11_edep9.t = 10.9051;
particleInput11_edep9.e = 0.860127;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep9));
supera::EDep particleInput11_edep10;
particleInput11_edep10.x = 1.8;
particleInput11_edep10.y = 0.6;
particleInput11_edep10.z = -53.8;
particleInput11_edep10.t = 10.9128;
particleInput11_edep10.e = 0.435149;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep10));
supera::EDep particleInput11_edep11;
particleInput11_edep11.x = 1.8;
particleInput11_edep11.y = 0.6;
particleInput11_edep11.z = -53.4;
particleInput11_edep11.t = 10.9204;
particleInput11_edep11.e = 0.294646;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep11));
supera::EDep particleInput11_edep12;
particleInput11_edep12.x = 1.8;
particleInput11_edep12.y = 0.2;
particleInput11_edep12.z = -53.4;
particleInput11_edep12.t = 10.9281;
particleInput11_edep12.e = 0.168623;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep12));
supera::EDep particleInput11_edep13;
particleInput11_edep13.x = 2.2;
particleInput11_edep13.y = 0.2;
particleInput11_edep13.z = -53.4;
particleInput11_edep13.t = 10.9358;
particleInput11_edep13.e = 0.898417;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep13));
supera::EDep particleInput11_edep14;
particleInput11_edep14.x = 2.6;
particleInput11_edep14.y = 0.2;
particleInput11_edep14.z = -53.4;
particleInput11_edep14.t = 10.9435;
particleInput11_edep14.e = 0.138522;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep14));
supera::EDep particleInput11_edep15;
particleInput11_edep15.x = 2.6;
particleInput11_edep15.y = 0.2;
particleInput11_edep15.z = -53.4;
particleInput11_edep15.t = 10.9511;
particleInput11_edep15.e = 0.0352126;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep15));
supera::EDep particleInput11_edep16;
particleInput11_edep16.x = 2.6;
particleInput11_edep16.y = 0.2;
particleInput11_edep16.z = -53;
particleInput11_edep16.t = 10.9557;
particleInput11_edep16.e = 0.359229;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep16));
supera::EDep particleInput11_edep17;
particleInput11_edep17.x = 2.6;
particleInput11_edep17.y = -0.2;
particleInput11_edep17.z = -53;
particleInput11_edep17.t = 10.9603;
particleInput11_edep17.e = 0.481244;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep17));
supera::EDep particleInput11_edep18;
particleInput11_edep18.x = 3;
particleInput11_edep18.y = -0.2;
particleInput11_edep18.z = -53;
particleInput11_edep18.t = 10.9649;
particleInput11_edep18.e = 0.137659;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep18));
supera::EDep particleInput11_edep19;
particleInput11_edep19.x = 3;
particleInput11_edep19.y = -0.2;
particleInput11_edep19.z = -53;
particleInput11_edep19.t = 10.9695;
particleInput11_edep19.e = 0.358218;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep19));
supera::EDep particleInput11_edep20;
particleInput11_edep20.x = 3;
particleInput11_edep20.y = -0.2;
particleInput11_edep20.z = -52.6;
particleInput11_edep20.t = 10.9755;
particleInput11_edep20.e = 0.465983;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep20));
supera::EDep particleInput11_edep21;
particleInput11_edep21.x = 3.4;
particleInput11_edep21.y = -0.6;
particleInput11_edep21.z = -52.6;
particleInput11_edep21.t = 10.9815;
particleInput11_edep21.e = 0.0613435;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep21));
supera::EDep particleInput11_edep22;
particleInput11_edep22.x = 3.4;
particleInput11_edep22.y = -0.6;
particleInput11_edep22.z = -52.6;
particleInput11_edep22.t = 10.9874;
particleInput11_edep22.e = 0.779836;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep22));
supera::EDep particleInput11_edep23;
particleInput11_edep23.x = 3.4;
particleInput11_edep23.y = -0.6;
particleInput11_edep23.z = -52.2;
particleInput11_edep23.t = 10.9946;
particleInput11_edep23.e = 0.232786;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep23));
supera::EDep particleInput11_edep24;
particleInput11_edep24.x = 3.8;
particleInput11_edep24.y = -0.6;
particleInput11_edep24.z = -52.2;
particleInput11_edep24.t = 11.0018;
particleInput11_edep24.e = 0.525544;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep24));
supera::EDep particleInput11_edep25;
particleInput11_edep25.x = 3.8;
particleInput11_edep25.y = -1;
particleInput11_edep25.z = -52.2;
particleInput11_edep25.t = 11.009;
particleInput11_edep25.e = 0.556794;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep25));
supera::EDep particleInput11_edep26;
particleInput11_edep26.x = 4.2;
particleInput11_edep26.y = -1;
particleInput11_edep26.z = -52.2;
particleInput11_edep26.t = 11.0162;
particleInput11_edep26.e = 0.0230653;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep26));
supera::EDep particleInput11_edep27;
particleInput11_edep27.x = 4.2;
particleInput11_edep27.y = -1;
particleInput11_edep27.z = -51.8;
particleInput11_edep27.t = 11.0234;
particleInput11_edep27.e = 1.03475;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep27));
supera::EDep particleInput11_edep28;
particleInput11_edep28.x = 4.2;
particleInput11_edep28.y = -1.4;
particleInput11_edep28.z = -51.8;
particleInput11_edep28.t = 11.0306;
particleInput11_edep28.e = 0.0245181;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep28));
supera::EDep particleInput11_edep29;
particleInput11_edep29.x = 4.6;
particleInput11_edep29.y = -1.4;
particleInput11_edep29.z = -51.8;
particleInput11_edep29.t = 11.0378;
particleInput11_edep29.e = 0.0414425;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep29));
supera::EDep particleInput11_edep30;
particleInput11_edep30.x = 4.6;
particleInput11_edep30.y = -1.4;
particleInput11_edep30.z = -51.8;
particleInput11_edep30.t = 11.045;
particleInput11_edep30.e = 0.250303;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep30));
supera::EDep particleInput11_edep31;
particleInput11_edep31.x = 4.6;
particleInput11_edep31.y = -1.4;
particleInput11_edep31.z = -51.4;
particleInput11_edep31.t = 11.0531;
particleInput11_edep31.e = 0.656857;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep31));
supera::EDep particleInput11_edep32;
particleInput11_edep32.x = 5;
particleInput11_edep32.y = -1.4;
particleInput11_edep32.z = -51.4;
particleInput11_edep32.t = 11.0612;
particleInput11_edep32.e = 0.753624;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep32));
supera::EDep particleInput11_edep33;
particleInput11_edep33.x = 5;
particleInput11_edep33.y = -1.4;
particleInput11_edep33.z = -51;
particleInput11_edep33.t = 11.0693;
particleInput11_edep33.e = 0.189654;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep33));
supera::EDep particleInput11_edep34;
particleInput11_edep34.x = 5.4;
particleInput11_edep34.y = -1.4;
particleInput11_edep34.z = -51;
particleInput11_edep34.t = 11.0774;
particleInput11_edep34.e = 0.671832;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep34));
supera::EDep particleInput11_edep35;
particleInput11_edep35.x = 5.4;
particleInput11_edep35.y = -1.8;
particleInput11_edep35.z = -51;
particleInput11_edep35.t = 11.0855;
particleInput11_edep35.e = 0.271446;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep35));
supera::EDep particleInput11_edep36;
particleInput11_edep36.x = 5.8;
particleInput11_edep36.y = -1.8;
particleInput11_edep36.z = -51;
particleInput11_edep36.t = 11.0936;
particleInput11_edep36.e = 0.27755;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep36));
supera::EDep particleInput11_edep37;
particleInput11_edep37.x = 5.8;
particleInput11_edep37.y = -1.8;
particleInput11_edep37.z = -50.6;
particleInput11_edep37.t = 11.1017;
particleInput11_edep37.e = 0.665728;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep37));
supera::EDep particleInput11_edep38;
particleInput11_edep38.x = 6.2;
particleInput11_edep38.y = -1.8;
particleInput11_edep38.z = -50.6;
particleInput11_edep38.t = 11.1098;
particleInput11_edep38.e = 0.744753;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep38));
supera::EDep particleInput11_edep39;
particleInput11_edep39.x = 6.2;
particleInput11_edep39.y = -1.8;
particleInput11_edep39.z = -50.2;
particleInput11_edep39.t = 11.1179;
particleInput11_edep39.e = 0.198525;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep39));
supera::EDep particleInput11_edep40;
particleInput11_edep40.x = 6.6;
particleInput11_edep40.y = -1.8;
particleInput11_edep40.z = -50.2;
particleInput11_edep40.t = 11.126;
particleInput11_edep40.e = 0.471697;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep40));
supera::EDep particleInput11_edep41;
particleInput11_edep41.x = 6.6;
particleInput11_edep41.y = -2.2;
particleInput11_edep41.z = -50.2;
particleInput11_edep41.t = 11.1341;
particleInput11_edep41.e = 0.471581;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep41));
supera::EDep particleInput11_edep42;
particleInput11_edep42.x = 7;
particleInput11_edep42.y = -2.2;
particleInput11_edep42.z = -50.2;
particleInput11_edep42.t = 11.1422;
particleInput11_edep42.e = 0.268679;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep42));
supera::EDep particleInput11_edep43;
particleInput11_edep43.x = 7;
particleInput11_edep43.y = -2.2;
particleInput11_edep43.z = -49.8;
particleInput11_edep43.t = 11.1503;
particleInput11_edep43.e = 0.674599;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep43));
supera::EDep particleInput11_edep44;
particleInput11_edep44.x = 7.4;
particleInput11_edep44.y = -2.2;
particleInput11_edep44.z = -49.8;
particleInput11_edep44.t = 11.1584;
particleInput11_edep44.e = 0.274029;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep44));
supera::EDep particleInput11_edep45;
particleInput11_edep45.x = 7.4;
particleInput11_edep45.y = -2.2;
particleInput11_edep45.z = -49.8;
particleInput11_edep45.t = 11.1665;
particleInput11_edep45.e = 0.385824;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep45));
supera::EDep particleInput11_edep46;
particleInput11_edep46.x = 7.4;
particleInput11_edep46.y = -2.2;
particleInput11_edep46.z = -49.4;
particleInput11_edep46.t = 11.1723;
particleInput11_edep46.e = 0.351533;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep46));
supera::EDep particleInput11_edep47;
particleInput11_edep47.x = 7.8;
particleInput11_edep47.y = -2.2;
particleInput11_edep47.z = -49.4;
particleInput11_edep47.t = 11.1782;
particleInput11_edep47.e = 0.673454;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep47));
supera::EDep particleInput11_edep48;
particleInput11_edep48.x = 7.8;
particleInput11_edep48.y = -2.6;
particleInput11_edep48.z = -49.4;
particleInput11_edep48.t = 11.184;
particleInput11_edep48.e = 0.153304;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep48));
supera::EDep particleInput11_edep49;
particleInput11_edep49.x = 7.8;
particleInput11_edep49.y = -2.6;
particleInput11_edep49.z = -49;
particleInput11_edep49.t = 11.1898;
particleInput11_edep49.e = 0.109913;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep49));
supera::EDep particleInput11_edep50;
particleInput11_edep50.x = 7.8;
particleInput11_edep50.y = -2.6;
particleInput11_edep50.z = -49;
particleInput11_edep50.t = 11.1957;
particleInput11_edep50.e = 0.09209;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep50));
supera::EDep particleInput11_edep51;
particleInput11_edep51.x = 8.2;
particleInput11_edep51.y = -2.6;
particleInput11_edep51.z = -49;
particleInput11_edep51.t = 11.2012;
particleInput11_edep51.e = 0.932787;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep51));
supera::EDep particleInput11_edep52;
particleInput11_edep52.x = 8.6;
particleInput11_edep52.y = -2.6;
particleInput11_edep52.z = -49;
particleInput11_edep52.t = 11.2068;
particleInput11_edep52.e = 0.0380627;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep52));
supera::EDep particleInput11_edep53;
particleInput11_edep53.x = 8.6;
particleInput11_edep53.y = -2.6;
particleInput11_edep53.z = -48.6;
particleInput11_edep53.t = 11.2123;
particleInput11_edep53.e = 0.139865;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep53));
supera::EDep particleInput11_edep54;
particleInput11_edep54.x = 8.6;
particleInput11_edep54.y = -2.6;
particleInput11_edep54.z = -48.6;
particleInput11_edep54.t = 11.2179;
particleInput11_edep54.e = 0.746039;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep54));
supera::EDep particleInput11_edep55;
particleInput11_edep55.x = 9;
particleInput11_edep55.y = -2.6;
particleInput11_edep55.z = -48.6;
particleInput11_edep55.t = 11.2257;
particleInput11_edep55.e = 0.09976;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep55));
supera::EDep particleInput11_edep56;
particleInput11_edep56.x = 9;
particleInput11_edep56.y = -2.6;
particleInput11_edep56.z = -48.6;
particleInput11_edep56.t = 11.2335;
particleInput11_edep56.e = 0.344618;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep56));
supera::EDep particleInput11_edep57;
particleInput11_edep57.x = 9;
particleInput11_edep57.y = -2.6;
particleInput11_edep57.z = -48.2;
particleInput11_edep57.t = 11.2407;
particleInput11_edep57.e = 0.305673;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep57));
supera::EDep particleInput11_edep58;
particleInput11_edep58.x = 9;
particleInput11_edep58.y = -3;
particleInput11_edep58.z = -48.2;
particleInput11_edep58.t = 11.2478;
particleInput11_edep58.e = 0.393622;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep58));
supera::EDep particleInput11_edep59;
particleInput11_edep59.x = 9.4;
particleInput11_edep59.y = -3;
particleInput11_edep59.z = -48.2;
particleInput11_edep59.t = 11.255;
particleInput11_edep59.e = 0.69349;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep59));
supera::EDep particleInput11_edep60;
particleInput11_edep60.x = 9.4;
particleInput11_edep60.y = -3;
particleInput11_edep60.z = -48.2;
particleInput11_edep60.t = 11.2621;
particleInput11_edep60.e = 0.392185;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep60));
supera::EDep particleInput11_edep61;
particleInput11_edep61.x = 9.8;
particleInput11_edep61.y = -3;
particleInput11_edep61.z = -48.2;
particleInput11_edep61.t = 11.2678;
particleInput11_edep61.e = 0.387826;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep61));
supera::EDep particleInput11_edep62;
particleInput11_edep62.x = 9.8;
particleInput11_edep62.y = -3.4;
particleInput11_edep62.z = -48.2;
particleInput11_edep62.t = 11.2735;
particleInput11_edep62.e = 0.0505778;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep62));
supera::EDep particleInput11_edep63;
particleInput11_edep63.x = 9.8;
particleInput11_edep63.y = -3.4;
particleInput11_edep63.z = -47.8;
particleInput11_edep63.t = 11.2792;
particleInput11_edep63.e = 0.523848;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep63));
supera::EDep particleInput11_edep64;
particleInput11_edep64.x = 10.2;
particleInput11_edep64.y = -3.4;
particleInput11_edep64.z = -47.8;
particleInput11_edep64.t = 11.2849;
particleInput11_edep64.e = 0.533987;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep64));
supera::EDep particleInput11_edep65;
particleInput11_edep65.x = 10.2;
particleInput11_edep65.y = -3.4;
particleInput11_edep65.z = -47.8;
particleInput11_edep65.t = 11.2907;
particleInput11_edep65.e = 0.356801;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep65));
supera::EDep particleInput11_edep66;
particleInput11_edep66.x = 10.6;
particleInput11_edep66.y = -3.4;
particleInput11_edep66.z = -47.8;
particleInput11_edep66.t = 11.2979;
particleInput11_edep66.e = 0.413808;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep66));
supera::EDep particleInput11_edep67;
particleInput11_edep67.x = 10.6;
particleInput11_edep67.y = -3.8;
particleInput11_edep67.z = -47.8;
particleInput11_edep67.t = 11.3051;
particleInput11_edep67.e = 0.352425;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep67));
supera::EDep particleInput11_edep68;
particleInput11_edep68.x = 10.6;
particleInput11_edep68.y = -3.8;
particleInput11_edep68.z = -47.8;
particleInput11_edep68.t = 11.3124;
particleInput11_edep68.e = 0.0538345;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep68));
supera::EDep particleInput11_edep69;
particleInput11_edep69.x = 11;
particleInput11_edep69.y = -3.8;
particleInput11_edep69.z = -47.8;
particleInput11_edep69.t = 11.3182;
particleInput11_edep69.e = 0.729074;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep69));
supera::EDep particleInput11_edep70;
particleInput11_edep70.x = 11;
particleInput11_edep70.y = -3.8;
particleInput11_edep70.z = -47.4;
particleInput11_edep70.t = 11.3241;
particleInput11_edep70.e = 0.217936;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep70));
supera::EDep particleInput11_edep71;
particleInput11_edep71.x = 11;
particleInput11_edep71.y = -3.8;
particleInput11_edep71.z = -47.4;
particleInput11_edep71.t = 11.33;
particleInput11_edep71.e = 0.203823;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep71));
supera::EDep particleInput11_edep72;
particleInput11_edep72.x = 11.4;
particleInput11_edep72.y = -3.8;
particleInput11_edep72.z = -47.4;
particleInput11_edep72.t = 11.3369;
particleInput11_edep72.e = 0.72283;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep72));
supera::EDep particleInput11_edep73;
particleInput11_edep73.x = 11.4;
particleInput11_edep73.y = -3.8;
particleInput11_edep73.z = -47.4;
particleInput11_edep73.t = 11.3438;
particleInput11_edep73.e = 0.219946;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep73));
supera::EDep particleInput11_edep74;
particleInput11_edep74.x = 11.8;
particleInput11_edep74.y = -3.8;
particleInput11_edep74.z = -47.4;
particleInput11_edep74.t = 11.3473;
particleInput11_edep74.e = 0.192491;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep74));
supera::EDep particleInput11_edep75;
particleInput11_edep75.x = 11.8;
particleInput11_edep75.y = -4.2;
particleInput11_edep75.z = -47.4;
particleInput11_edep75.t = 11.3507;
particleInput11_edep75.e = 0.250804;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep75));
supera::EDep particleInput11_edep76;
particleInput11_edep76.x = 11.8;
particleInput11_edep76.y = -4.2;
particleInput11_edep76.z = -47.4;
particleInput11_edep76.t = 11.3541;
particleInput11_edep76.e = 0.392249;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep76));
supera::EDep particleInput11_edep77;
particleInput11_edep77.x = 11.8;
particleInput11_edep77.y = -3.8;
particleInput11_edep77.z = -47.4;
particleInput11_edep77.t = 11.3567;
particleInput11_edep77.e = 0.088949;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep77));
supera::EDep particleInput11_edep78;
particleInput11_edep78.x = 11.8;
particleInput11_edep78.y = -3.8;
particleInput11_edep78.z = -47.8;
particleInput11_edep78.t = 11.3593;
particleInput11_edep78.e = 0.31239;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep78));
supera::EDep particleInput11_edep79;
particleInput11_edep79.x = 11.8;
particleInput11_edep79.y = -3.8;
particleInput11_edep79.z = -47.8;
particleInput11_edep79.t = 11.3618;
particleInput11_edep79.e = 0.365169;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep79));
supera::EDep particleInput11_edep80;
particleInput11_edep80.x = 11.8;
particleInput11_edep80.y = -3.8;
particleInput11_edep80.z = -47.8;
particleInput11_edep80.t = 11.366;
particleInput11_edep80.e = 0.361148;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep80));
supera::EDep particleInput11_edep81;
particleInput11_edep81.x = 11.8;
particleInput11_edep81.y = -3.8;
particleInput11_edep81.z = -47.8;
particleInput11_edep81.t = 11.3683;
particleInput11_edep81.e = 0.162493;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep81));
supera::EDep particleInput11_edep82;
particleInput11_edep82.x = -0.6;
particleInput11_edep82.y = -32.2;
particleInput11_edep82.z = -54.6;
particleInput11_edep82.t = 12.6498;
particleInput11_edep82.e = 0.0031776;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep82));
supera::EDep particleInput11_edep83;
particleInput11_edep83.x = -0.6;
particleInput11_edep83.y = -32.2;
particleInput11_edep83.z = -54.6;
particleInput11_edep83.t = 12.6674;
particleInput11_edep83.e = 0.101437;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep83));
supera::EDep particleInput11_edep84;
particleInput11_edep84.x = -1;
particleInput11_edep84.y = -32.2;
particleInput11_edep84.z = -54.2;
particleInput11_edep84.t = 12.6498;
particleInput11_edep84.e = 0.0714654;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep84));
supera::EDep particleInput11_edep85;
particleInput11_edep85.x = 1.4;
particleInput11_edep85.y = -31;
particleInput11_edep85.z = -54.6;
particleInput11_edep85.t = 12.5599;
particleInput11_edep85.e = 0.0192375;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep85));
supera::EDep particleInput11_edep86;
particleInput11_edep86.x = 6.6;
particleInput11_edep86.y = -28.6;
particleInput11_edep86.z = -50.2;
particleInput11_edep86.t = 12.3063;
particleInput11_edep86.e = 0.0633362;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep86));
supera::EDep particleInput11_edep87;
particleInput11_edep87.x = 8.2;
particleInput11_edep87.y = -28.2;
particleInput11_edep87.z = -51;
particleInput11_edep87.t = 12.2492;
particleInput11_edep87.e = 0.138336;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep87));
supera::EDep particleInput11_edep88;
particleInput11_edep88.x = 14.6;
particleInput11_edep88.y = -10.2;
particleInput11_edep88.z = -48.2;
particleInput11_edep88.t = 11.5959;
particleInput11_edep88.e = 0.114009;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep88));
supera::EDep particleInput11_edep89;
particleInput11_edep89.x = -9.4;
particleInput11_edep89.y = -8.6;
particleInput11_edep89.z = -28.6;
particleInput11_edep89.t = 12.619;
particleInput11_edep89.e = 0.00207092;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep89));
supera::EDep particleInput11_edep90;
particleInput11_edep90.x = -9.8;
particleInput11_edep90.y = -8.6;
particleInput11_edep90.z = -28.6;
particleInput11_edep90.t = 12.6324;
particleInput11_edep90.e = 0.00110668;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep90));
supera::EDep particleInput11_edep91;
particleInput11_edep91.x = -9.8;
particleInput11_edep91.y = -8.6;
particleInput11_edep91.z = -28.6;
particleInput11_edep91.t = 12.6457;
particleInput11_edep91.e = 0.0707433;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep91));
supera::EDep particleInput11_edep92;
particleInput11_edep92.x = -9;
particleInput11_edep92.y = -8.6;
particleInput11_edep92.z = -28.6;
particleInput11_edep92.t = 12.619;
particleInput11_edep92.e = 0.0101893;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep92));
supera::EDep particleInput11_edep93;
particleInput11_edep93.x = -8.6;
particleInput11_edep93.y = -8.6;
particleInput11_edep93.z = -28.2;
particleInput11_edep93.t = 12.5986;
particleInput11_edep93.e = 0.014702;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep93));
supera::EDep particleInput11_edep94;
particleInput11_edep94.x = -7.4;
particleInput11_edep94.y = -5.8;
particleInput11_edep94.z = -29;
particleInput11_edep94.t = 12.4956;
particleInput11_edep94.e = 0.0480985;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep94));
supera::EDep particleInput11_edep95;
particleInput11_edep95.x = -8.2;
particleInput11_edep95.y = -6.2;
particleInput11_edep95.z = -29.4;
particleInput11_edep95.t = 12.4547;
particleInput11_edep95.e = 0.104007;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep95));
supera::EDep particleInput11_edep96;
particleInput11_edep96.x = -1;
particleInput11_edep96.y = -7.8;
particleInput11_edep96.z = -36.2;
particleInput11_edep96.t = 12.1133;
particleInput11_edep96.e = 0.0521691;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep96));
supera::EDep particleInput11_edep97;
particleInput11_edep97.x = 5;
particleInput11_edep97.y = -0.6;
particleInput11_edep97.z = -41.8;
particleInput11_edep97.t = 11.7555;
particleInput11_edep97.e = 0.0447433;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep97));
supera::EDep particleInput11_edep98;
particleInput11_edep98.x = 6.6;
particleInput11_edep98.y = -0.6;
particleInput11_edep98.z = -43.4;
particleInput11_edep98.t = 11.6809;
particleInput11_edep98.e = 0.00119876;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep98));
supera::EDep particleInput11_edep99;
particleInput11_edep99.x = 9.8;
particleInput11_edep99.y = -0.2;
particleInput11_edep99.z = -47;
particleInput11_edep99.t = 11.5111;
particleInput11_edep99.e = 0.0138404;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep99));
supera::EDep particleInput11_edep100;
particleInput11_edep100.x = 10.2;
particleInput11_edep100.y = -0.6;
particleInput11_edep100.z = -47.4;
particleInput11_edep100.t = 11.4971;
particleInput11_edep100.e = 0.14813;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep100));
supera::EDep particleInput11_edep101;
particleInput11_edep101.x = 63.8;
particleInput11_edep101.y = -39.4;
particleInput11_edep101.z = -75;
particleInput11_edep101.t = 14.2161;
particleInput11_edep101.e = 0.0031776;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep101));
supera::EDep particleInput11_edep102;
particleInput11_edep102.x = 63.8;
particleInput11_edep102.y = -39.4;
particleInput11_edep102.z = -75;
particleInput11_edep102.t = 14.2826;
particleInput11_edep102.e = 0.0636857;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep102));
supera::EDep particleInput11_edep103;
particleInput11_edep103.x = 64.2;
particleInput11_edep103.y = -41;
particleInput11_edep103.z = -74.6;
particleInput11_edep103.t = 14.2161;
particleInput11_edep103.e = 0.003659;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep103));
supera::EDep particleInput11_edep104;
particleInput11_edep104.x = 64.2;
particleInput11_edep104.y = -41;
particleInput11_edep104.z = -74.6;
particleInput11_edep104.t = 14.2147;
particleInput11_edep104.e = 0.0251838;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep104));
supera::EDep particleInput11_edep105;
particleInput11_edep105.x = 64.2;
particleInput11_edep105.y = -40.6;
particleInput11_edep105.z = -74.6;
particleInput11_edep105.t = 14.2001;
particleInput11_edep105.e = 0.0403871;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep105));
supera::EDep particleInput11_edep106;
particleInput11_edep106.x = 63.8;
particleInput11_edep106.y = -45;
particleInput11_edep106.z = -79.8;
particleInput11_edep106.t = 13.9747;
particleInput11_edep106.e = 0.0356491;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep106));
supera::EDep particleInput11_edep107;
particleInput11_edep107.x = 59.8;
particleInput11_edep107.y = -43;
particleInput11_edep107.z = -82.6;
particleInput11_edep107.t = 13.8022;
particleInput11_edep107.e = 0.161841;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep107));
supera::EDep particleInput11_edep108;
particleInput11_edep108.x = 60.2;
particleInput11_edep108.y = -43;
particleInput11_edep108.z = -82.6;
particleInput11_edep108.t = 13.8024;
particleInput11_edep108.e = 0.00031343;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep108));
supera::EDep particleInput11_edep109;
particleInput11_edep109.x = 60.2;
particleInput11_edep109.y = -43;
particleInput11_edep109.z = -82.6;
particleInput11_edep109.t = 13.8077;
particleInput11_edep109.e = 0.0277466;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep109));
supera::EDep particleInput11_edep110;
particleInput11_edep110.x = 59.8;
particleInput11_edep110.y = -43.8;
particleInput11_edep110.z = -80.2;
particleInput11_edep110.t = 13.7255;
particleInput11_edep110.e = 0.288414;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep110));
supera::EDep particleInput11_edep111;
particleInput11_edep111.x = 14.6;
particleInput11_edep111.y = -1.8;
particleInput11_edep111.z = -29;
particleInput11_edep111.t = 12.3629;
particleInput11_edep111.e = 0.0171902;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep111));
supera::EDep particleInput11_edep112;
particleInput11_edep112.x = 19.8;
particleInput11_edep112.y = -3.4;
particleInput11_edep112.z = -32.6;
particleInput11_edep112.t = 12.1454;
particleInput11_edep112.e = 0.00737417;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep112));
supera::EDep particleInput11_edep113;
particleInput11_edep113.x = 23.4;
particleInput11_edep113.y = -4.6;
particleInput11_edep113.z = -38.6;
particleInput11_edep113.t = 11.9102;
particleInput11_edep113.e = 0.00625953;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep113));
supera::EDep particleInput11_edep114;
particleInput11_edep114.x = 23.4;
particleInput11_edep114.y = -4.6;
particleInput11_edep114.z = -39.4;
particleInput11_edep114.t = 11.8784;
particleInput11_edep114.e = 0.236041;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep114));
supera::EDep particleInput11_edep115;
particleInput11_edep115.x = 19.8;
particleInput11_edep115.y = -2.6;
particleInput11_edep115.z = -38.6;
particleInput11_edep115.t = 11.7295;
particleInput11_edep115.e = 0.159273;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep115));
supera::EDep particleInput11_edep116;
particleInput11_edep116.x = 15.4;
particleInput11_edep116.y = -2.6;
particleInput11_edep116.z = -41;
particleInput11_edep116.t = 11.5702;
particleInput11_edep116.e = 0.0906444;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep116));
supera::EDep particleInput11_edep117;
particleInput11_edep117.x = 11;
particleInput11_edep117.y = -3;
particleInput11_edep117.z = -47;
particleInput11_edep117.t = 11.3221;
particleInput11_edep117.e = 0.0516445;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep117));
supera::EDep particleInput11_edep118;
particleInput11_edep118.x = 4.6;
particleInput11_edep118.y = -1.4;
particleInput11_edep118.z = -51.8;
particleInput11_edep118.t = 11.045;
particleInput11_edep118.e = 0.0129185;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep118));
supera::EDep particleInput11_edep119;
particleInput11_edep119.x = 13;
particleInput11_edep119.y = -3;
particleInput11_edep119.z = -58.6;
particleInput11_edep119.t = 11.554;
particleInput11_edep119.e = 0.0031776;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep119));
supera::EDep particleInput11_edep120;
particleInput11_edep120.x = 13;
particleInput11_edep120.y = -3;
particleInput11_edep120.z = -58.6;
particleInput11_edep120.t = 11.6253;
particleInput11_edep120.e = 0.0532757;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep120));
supera::EDep particleInput11_edep121;
particleInput11_edep121.x = 14.2;
particleInput11_edep121.y = -1;
particleInput11_edep121.z = -58.2;
particleInput11_edep121.t = 11.554;
particleInput11_edep121.e = 0.01029;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep121));
supera::EDep particleInput11_edep122;
particleInput11_edep122.x = 8.2;
particleInput11_edep122.y = -1.8;
particleInput11_edep122.z = -53;
particleInput11_edep122.t = 11.2709;
particleInput11_edep122.e = 0.00253165;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep122));
supera::EDep particleInput11_edep123;
particleInput11_edep123.x = 7.8;
particleInput11_edep123.y = -3.4;
particleInput11_edep123.z = -51.8;
particleInput11_edep123.t = 11.2077;
particleInput11_edep123.e = 0.00910878;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep123));
supera::EDep particleInput11_edep124;
particleInput11_edep124.x = 7.4;
particleInput11_edep124.y = -3.4;
particleInput11_edep124.z = -51.4;
particleInput11_edep124.t = 11.193;
particleInput11_edep124.e = 0.000196779;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep124));
supera::EDep particleInput11_edep125;
particleInput11_edep125.x = 6.2;
particleInput11_edep125.y = -2.2;
particleInput11_edep125.z = -50.6;
particleInput11_edep125.t = 11.1217;
particleInput11_edep125.e = 0.0102075;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep125));
supera::EDep particleInput11_edep126;
particleInput11_edep126.x = 3;
particleInput11_edep126.y = -0.2;
particleInput11_edep126.z = -53;
particleInput11_edep126.t = 10.9511;
particleInput11_edep126.e = 0.0031776;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep126));
supera::EDep particleInput11_edep127;
particleInput11_edep127.x = 3;
particleInput11_edep127.y = -0.2;
particleInput11_edep127.z = -53;
particleInput11_edep127.t = 10.9711;
particleInput11_edep127.e = 0.0358089;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep127));
supera::EDep particleInput11_edep128;
particleInput11_edep128.x = 1;
particleInput11_edep128.y = 1;
particleInput11_edep128.z = -54.2;
particleInput11_edep128.t = 10.8821;
particleInput11_edep128.e = 0.0213072;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep128));
supera::EDep particleInput11_edep129;
particleInput11_edep129.x = -2.6;
particleInput11_edep129.y = -10.2;
particleInput11_edep129.z = -69.8;
particleInput11_edep129.t = 12.8106;
particleInput11_edep129.e = 0.0031776;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep129));
supera::EDep particleInput11_edep130;
particleInput11_edep130.x = -2.6;
particleInput11_edep130.y = -10.2;
particleInput11_edep130.z = -69.8;
particleInput11_edep130.t = 12.8541;
particleInput11_edep130.e = 0.0695706;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep130));
supera::EDep particleInput11_edep131;
particleInput11_edep131.x = -2.6;
particleInput11_edep131.y = -9.4;
particleInput11_edep131.z = -69;
particleInput11_edep131.t = 12.8106;
particleInput11_edep131.e = 0.00674182;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep131));
supera::EDep particleInput11_edep132;
particleInput11_edep132.x = -1.8;
particleInput11_edep132.y = -9.4;
particleInput11_edep132.z = -68.2;
particleInput11_edep132.t = 12.7764;
particleInput11_edep132.e = 0.0199659;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep132));
supera::EDep particleInput11_edep133;
particleInput11_edep133.x = -1.8;
particleInput11_edep133.y = -6.2;
particleInput11_edep133.z = -68.6;
particleInput11_edep133.t = 12.6665;
particleInput11_edep133.e = 0.00237202;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep133));
supera::EDep particleInput11_edep134;
particleInput11_edep134.x = -0.6;
particleInput11_edep134.y = -2.6;
particleInput11_edep134.z = -71.4;
particleInput11_edep134.t = 12.5053;
particleInput11_edep134.e = 0.00389041;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep134));
supera::EDep particleInput11_edep135;
particleInput11_edep135.x = 3.4;
particleInput11_edep135.y = -0.2;
particleInput11_edep135.z = -78.6;
particleInput11_edep135.t = 12.2276;
particleInput11_edep135.e = 0.0620167;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep135));
supera::EDep particleInput11_edep136;
particleInput11_edep136.x = 5.4;
particleInput11_edep136.y = -13;
particleInput11_edep136.z = -63;
particleInput11_edep136.t = 11.5574;
particleInput11_edep136.e = 0.074798;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep136));
supera::EDep particleInput11_edep137;
particleInput11_edep137.x = 3.8;
particleInput11_edep137.y = -1;
particleInput11_edep137.z = -52.2;
particleInput11_edep137.t = 11.0245;
particleInput11_edep137.e = 0.225807;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep137));
supera::EDep particleInput11_edep138;
particleInput11_edep138.x = 39;
particleInput11_edep138.y = -83.8;
particleInput11_edep138.z = 80.6;
particleInput11_edep138.t = 29.1565;
particleInput11_edep138.e = 0.00650139;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep138));
supera::EDep particleInput11_edep139;
particleInput11_edep139.x = -3.8;
particleInput11_edep139.y = 81.8;
particleInput11_edep139.z = 82.2;
particleInput11_edep139.t = 20.7103;
particleInput11_edep139.e = 0.0156247;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep139));
supera::EDep particleInput11_edep140;
particleInput11_edep140.x = -3.4;
particleInput11_edep140.y = 81;
particleInput11_edep140.z = 82.2;
particleInput11_edep140.t = 20.6764;
particleInput11_edep140.e = 0.00635404;

particleInput11.pcloud.emplace_back(std::move(particleInput11_edep140));
particleInput11.valid = 1;
particleInput11.type = static_cast<supera::ProcessType>(0);
supera::Particle particleInput11_particle;
particleInput11_particle.id = 12;
particleInput11_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput11_particle.trackid = 11;
particleInput11_particle.pdg = -11;
particleInput11_particle.px = 34.9792;
particleInput11_particle.py = -25.4476;
particleInput11_particle.pz = 17.6283;
particleInput11_particle.vtx = {0.14755, 1.49633, -54.5587, 10.8502};
particleInput11_particle.end_pt = {11.8708, -3.89269, -47.7023, 11.3687};
particleInput11_particle.first_step = {0.14755, 1.49633, -54.5587, 10.8502};
particleInput11_particle.last_step = {11.8708, -3.89269, -47.7023, 11.3687};
particleInput11_particle.dist_travel = 15.3603;
particleInput11_particle.energy_init = 46.7135;
particleInput11_particle.energy_deposit = 33.5085;
particleInput11_particle.process = "Decay";
particleInput11_particle.parent_trackid = 0;
particleInput11_particle.parent_pdg = -13;
particleInput11_particle.parent_vtx = {0, 0, 0, 0};
particleInput11_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput11_particle.ancestor_pdg = 0;
particleInput11_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput11_particle.ancestor_process = "";
particleInput11_particle.parent_process = "";
particleInput11_particle.parent_id = 1;
particleInput11_particle.children_id = {  };
particleInput11_particle.group_id = kINVALID_INSTANCEID;
particleInput11_particle.interaction_id = 0;
particleInput11.part = std::move(particleInput11_particle);

evInput.push_back(std::move(particleInput11));
supera::ParticleInput particleInput12;
particleInput12.pcloud.reserve(2);
supera::EDep particleInput12_edep0;
particleInput12_edep0.x = 3;
particleInput12_edep0.y = -0.2;
particleInput12_edep0.z = -53;
particleInput12_edep0.t = 10.9695;
particleInput12_edep0.e = 0.318759;

particleInput12.pcloud.emplace_back(std::move(particleInput12_edep0));
supera::EDep particleInput12_edep1;
particleInput12_edep1.x = 3;
particleInput12_edep1.y = -0.2;
particleInput12_edep1.z = -53;
particleInput12_edep1.t = 10.9715;
particleInput12_edep1.e = 0.17598;

particleInput12.pcloud.emplace_back(std::move(particleInput12_edep1));
particleInput12.valid = 1;
particleInput12.type = static_cast<supera::ProcessType>(8);
supera::Particle particleInput12_particle;
particleInput12_particle.id = 13;
particleInput12_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput12_particle.trackid = 12;
particleInput12_particle.pdg = 11;
particleInput12_particle.px = -0.0882954;
particleInput12_particle.py = -0.236592;
particleInput12_particle.pz = 0.828624;
particleInput12_particle.vtx = {2.85318, -0.151726, -52.913, 10.9695};
particleInput12_particle.end_pt = {2.85425, -0.168625, -52.8561, 10.972};
particleInput12_particle.first_step = {2.85318, -0.151726, -52.913, 10.9695};
particleInput12_particle.last_step = {2.85425, -0.168625, -52.8561, 10.972};
particleInput12_particle.dist_travel = 0.0593731;
particleInput12_particle.energy_init = 1.00574;
particleInput12_particle.energy_deposit = 0.494739;
particleInput12_particle.process = "eIoni";
particleInput12_particle.parent_trackid = 11;
particleInput12_particle.parent_pdg = -11;
particleInput12_particle.parent_vtx = {0, 0, 0, 0};
particleInput12_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput12_particle.ancestor_pdg = 0;
particleInput12_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput12_particle.ancestor_process = "";
particleInput12_particle.parent_process = "";
particleInput12_particle.parent_id = 12;
particleInput12_particle.children_id = {  };
particleInput12_particle.group_id = kINVALID_INSTANCEID;
particleInput12_particle.interaction_id = 0;
particleInput12.part = std::move(particleInput12_particle);

evInput.push_back(std::move(particleInput12));
supera::ParticleInput particleInput13;
particleInput13.pcloud.reserve(0);
particleInput13.valid = 1;
particleInput13.type = static_cast<supera::ProcessType>(2);
supera::Particle particleInput13_particle;
particleInput13_particle.id = 14;
particleInput13_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput13_particle.trackid = 13;
particleInput13_particle.pdg = 22;
particleInput13_particle.px = 6.31304;
particleInput13_particle.py = -7.26467;
particleInput13_particle.pz = 5.28439;
particleInput13_particle.vtx = {9.43697, -3.0664, -48.0411, 11.2621};
particleInput13_particle.end_pt = {1.79769e+308, 1.79769e+308, 1.79769e+308, 1.79769e+308};
particleInput13_particle.first_step = {1.79769e+308, 1.79769e+308, 1.79769e+308, 1.79769e+308};
particleInput13_particle.last_step = {1.79769e+308, 1.79769e+308, 1.79769e+308, 1.79769e+308};
particleInput13_particle.dist_travel = 935.355;
particleInput13_particle.energy_init = 10.9797;
particleInput13_particle.energy_deposit = 0;
particleInput13_particle.process = "brem";
particleInput13_particle.parent_trackid = 11;
particleInput13_particle.parent_pdg = -11;
particleInput13_particle.parent_vtx = {0, 0, 0, 0};
particleInput13_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput13_particle.ancestor_pdg = 0;
particleInput13_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput13_particle.ancestor_process = "";
particleInput13_particle.parent_process = "";
particleInput13_particle.parent_id = 12;
particleInput13_particle.children_id = {  };
particleInput13_particle.group_id = kINVALID_INSTANCEID;
particleInput13_particle.interaction_id = 0;
particleInput13.part = std::move(particleInput13_particle);

evInput.push_back(std::move(particleInput13));
supera::ParticleInput particleInput14;
particleInput14.pcloud.reserve(2);
supera::EDep particleInput14_edep0;
particleInput14_edep0.x = 25;
particleInput14_edep0.y = -20.6;
particleInput14_edep0.z = -35;
particleInput14_edep0.t = 12.1551;
particleInput14_edep0.e = 0.296807;

particleInput14.pcloud.emplace_back(std::move(particleInput14_edep0));
supera::EDep particleInput14_edep1;
particleInput14_edep1.x = 25;
particleInput14_edep1.y = -20.6;
particleInput14_edep1.z = -35;
particleInput14_edep1.t = 12.1577;
particleInput14_edep1.e = 0.314272;

particleInput14.pcloud.emplace_back(std::move(particleInput14_edep1));
particleInput14.valid = 1;
particleInput14.type = static_cast<supera::ProcessType>(4);
supera::Particle particleInput14_particle;
particleInput14_particle.id = 15;
particleInput14_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput14_particle.trackid = 14;
particleInput14_particle.pdg = 11;
particleInput14_particle.px = 0.60628;
particleInput14_particle.py = 0.130031;
particleInput14_particle.pz = 0.783233;
particleInput14_particle.vtx = {24.8292, -20.7789, -35.1569, 12.1551};
particleInput14_particle.end_pt = {24.8921, -20.7547, -35.1034, 12.1588};
particleInput14_particle.first_step = {24.8292, -20.7789, -35.1569, 12.1551};
particleInput14_particle.last_step = {24.8921, -20.7547, -35.1034, 12.1588};
particleInput14_particle.dist_travel = 0.0860548;
particleInput14_particle.energy_init = 1.12208;
particleInput14_particle.energy_deposit = 0.611079;
particleInput14_particle.process = "compt";
particleInput14_particle.parent_trackid = 13;
particleInput14_particle.parent_pdg = 22;
particleInput14_particle.parent_vtx = {0, 0, 0, 0};
particleInput14_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput14_particle.ancestor_pdg = 0;
particleInput14_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput14_particle.ancestor_process = "";
particleInput14_particle.parent_process = "";
particleInput14_particle.parent_id = 14;
particleInput14_particle.children_id = {  };
particleInput14_particle.group_id = kINVALID_INSTANCEID;
particleInput14_particle.interaction_id = 0;
particleInput14.part = std::move(particleInput14_particle);

evInput.push_back(std::move(particleInput14));
supera::ParticleInput particleInput15;
particleInput15.pcloud.reserve(7);
supera::EDep particleInput15_edep0;
particleInput15_edep0.x = 11.4;
particleInput15_edep0.y = -3.4;
particleInput15_edep0.z = -46.6;
particleInput15_edep0.t = 11.3483;
particleInput15_edep0.e = 0.286416;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep0));
supera::EDep particleInput15_edep1;
particleInput15_edep1.x = 11.8;
particleInput15_edep1.y = -3.4;
particleInput15_edep1.z = -46.6;
particleInput15_edep1.t = 11.3513;
particleInput15_edep1.e = 0.0856212;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep1));
supera::EDep particleInput15_edep2;
particleInput15_edep2.x = 11.8;
particleInput15_edep2.y = -3.4;
particleInput15_edep2.z = -46.2;
particleInput15_edep2.t = 11.3542;
particleInput15_edep2.e = 0.232888;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep2));
supera::EDep particleInput15_edep3;
particleInput15_edep3.x = 11.8;
particleInput15_edep3.y = -3.4;
particleInput15_edep3.z = -46.2;
particleInput15_edep3.t = 11.3572;
particleInput15_edep3.e = 0.183931;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep3));
supera::EDep particleInput15_edep4;
particleInput15_edep4.x = 11.8;
particleInput15_edep4.y = -3.8;
particleInput15_edep4.z = -46.2;
particleInput15_edep4.t = 11.3604;
particleInput15_edep4.e = 0.574434;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep4));
supera::EDep particleInput15_edep5;
particleInput15_edep5.x = 11.8;
particleInput15_edep5.y = -3.8;
particleInput15_edep5.z = -46.2;
particleInput15_edep5.t = 11.3637;
particleInput15_edep5.e = 0.335199;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep5));
supera::EDep particleInput15_edep6;
particleInput15_edep6.x = 11.8;
particleInput15_edep6.y = -3.8;
particleInput15_edep6.z = -46.2;
particleInput15_edep6.t = 11.3666;
particleInput15_edep6.e = 0.35511;

particleInput15.pcloud.emplace_back(std::move(particleInput15_edep6));
particleInput15.valid = 1;
particleInput15.type = static_cast<supera::ProcessType>(4);
supera::Particle particleInput15_particle;
particleInput15_particle.id = 16;
particleInput15_particle.shape = static_cast<supera::SemanticType_t>(6);
particleInput15_particle.trackid = 15;
particleInput15_particle.pdg = 11;
particleInput15_particle.px = 1.66336;
particleInput15_particle.py = -0.340332;
particleInput15_particle.pz = 1.85296;
particleInput15_particle.vtx = {11.518, -3.53797, -46.5186, 11.3483};
particleInput15_particle.end_pt = {11.6891, -3.68936, -46.2493, 11.3679};
particleInput15_particle.first_step = {11.518, -3.53797, -46.5186, 11.3483};
particleInput15_particle.last_step = {11.6891, -3.68936, -46.2493, 11.3679};
particleInput15_particle.dist_travel = 0.353138;
particleInput15_particle.energy_init = 2.5646;
particleInput15_particle.energy_deposit = 2.0536;
particleInput15_particle.process = "compt";
particleInput15_particle.parent_trackid = 11;
particleInput15_particle.parent_pdg = -11;
particleInput15_particle.parent_vtx = {0, 0, 0, 0};
particleInput15_particle.ancestor_trackid = kINVALID_TRACKID;
particleInput15_particle.ancestor_pdg = 0;
particleInput15_particle.ancestor_vtx = {0, 0, 0, 0};
particleInput15_particle.ancestor_process = "";
particleInput15_particle.parent_process = "";
particleInput15_particle.parent_id = 12;
particleInput15_particle.children_id = {  };
particleInput15_particle.group_id = kINVALID_INSTANCEID;
particleInput15_particle.interaction_id = 0;
particleInput15.part = std::move(particleInput15_particle);

evInput.push_back(std::move(particleInput15));

supera::EventOutput evtOutput;
evtOutput.Particles().reserve(4);
supera::ParticleLabel evtOutput_part0_label;
supera::Particle evtOutput_part0_label_part;
evtOutput_part0_label_part.id = 0;
evtOutput_part0_label_part.shape = static_cast<supera::SemanticType_t>(1);
evtOutput_part0_label_part.trackid = 0;
evtOutput_part0_label_part.pdg = -13;
evtOutput_part0_label_part.px = 4.10617;
evtOutput_part0_label_part.py = 3.74819;
evtOutput_part0_label_part.pz = -176.354;
evtOutput_part0_label_part.vtx = {0, 0, 0, 1};
evtOutput_part0_label_part.end_pt = {0.14755, 1.49633, -54.5587, 10.8502};
evtOutput_part0_label_part.first_step = {0, 0, 0, 1};
evtOutput_part0_label_part.last_step = {0.14755, 1.49633, -54.5587, 10.8502};
evtOutput_part0_label_part.dist_travel = 54.6454;
evtOutput_part0_label_part.energy_init = 205.658;
evtOutput_part0_label_part.energy_deposit = 99.9441;
evtOutput_part0_label_part.process = "primary";
evtOutput_part0_label_part.parent_trackid = kINVALID_TRACKID;
evtOutput_part0_label_part.parent_pdg = 0;
evtOutput_part0_label_part.parent_vtx = {0, 0, 0, 0};
evtOutput_part0_label_part.ancestor_trackid = kINVALID_TRACKID;
evtOutput_part0_label_part.ancestor_pdg = 0;
evtOutput_part0_label_part.ancestor_vtx = {0, 0, 0, 0};
evtOutput_part0_label_part.ancestor_process = "";
evtOutput_part0_label_part.parent_process = "";
evtOutput_part0_label_part.parent_id = 0;
evtOutput_part0_label_part.children_id = { 1 };
evtOutput_part0_label_part.group_id = 0;
evtOutput_part0_label_part.interaction_id = 0;
evtOutput_part0_label.part = std::move(evtOutput_part0_label_part);
evtOutput_part0_label.valid = 1;
evtOutput_part0_label.add_to_parent = 0;
evtOutput_part0_label.type = static_cast<supera::ProcessType>(0);
evtOutput_part0_label.trackid_v = {  };
supera::VoxelSet evtOutput_part0_label_energyVoxSet;
evtOutput_part0_label_energyVoxSet.id(18446744073709551615ul);
evtOutput_part0_label_energyVoxSet.reserve(83);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox0;
evtOutput_part0_label_energyVoxSet_vox0.set(static_cast<supera::VoxelID_t>(6959383750), 5.1895);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox0), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox1;
evtOutput_part0_label_energyVoxSet_vox1.set(static_cast<supera::VoxelID_t>(6965633750), 4.07518);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox1), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox2;
evtOutput_part0_label_energyVoxSet_vox2.set(static_cast<supera::VoxelID_t>(6971883750), 2.95149);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox2), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox3;
evtOutput_part0_label_energyVoxSet_vox3.set(static_cast<supera::VoxelID_t>(6978133750), 2.8806);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox3), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox4;
evtOutput_part0_label_energyVoxSet_vox4.set(static_cast<supera::VoxelID_t>(6984383750), 1.94274);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox4), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox5;
evtOutput_part0_label_energyVoxSet_vox5.set(static_cast<supera::VoxelID_t>(6990633750), 2.01283);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox5), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox6;
evtOutput_part0_label_energyVoxSet_vox6.set(static_cast<supera::VoxelID_t>(6996883750), 1.79158);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox6), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox7;
evtOutput_part0_label_energyVoxSet_vox7.set(static_cast<supera::VoxelID_t>(7003133750), 2.26157);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox7), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox8;
evtOutput_part0_label_energyVoxSet_vox8.set(static_cast<supera::VoxelID_t>(7009383750), 1.79643);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox8), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox9;
evtOutput_part0_label_energyVoxSet_vox9.set(static_cast<supera::VoxelID_t>(7015633750), 1.69942);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox9), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox10;
evtOutput_part0_label_energyVoxSet_vox10.set(static_cast<supera::VoxelID_t>(7021881250), 1.40495);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox10), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox11;
evtOutput_part0_label_energyVoxSet_vox11.set(static_cast<supera::VoxelID_t>(7021883750), 0.198421);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox11), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox12;
evtOutput_part0_label_energyVoxSet_vox12.set(static_cast<supera::VoxelID_t>(7028131250), 1.57963);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox12), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox13;
evtOutput_part0_label_energyVoxSet_vox13.set(static_cast<supera::VoxelID_t>(7034381250), 1.44692);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox13), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox14;
evtOutput_part0_label_energyVoxSet_vox14.set(static_cast<supera::VoxelID_t>(7040631250), 1.44692);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox14), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox15;
evtOutput_part0_label_energyVoxSet_vox15.set(static_cast<supera::VoxelID_t>(7046881250), 1.40586);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox15), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox16;
evtOutput_part0_label_energyVoxSet_vox16.set(static_cast<supera::VoxelID_t>(7053131250), 1.26383);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox16), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox17;
evtOutput_part0_label_energyVoxSet_vox17.set(static_cast<supera::VoxelID_t>(7059381250), 1.26383);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox17), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox18;
evtOutput_part0_label_energyVoxSet_vox18.set(static_cast<supera::VoxelID_t>(7065631250), 1.26383);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox18), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox19;
evtOutput_part0_label_energyVoxSet_vox19.set(static_cast<supera::VoxelID_t>(7071881250), 1.74209);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox19), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox20;
evtOutput_part0_label_energyVoxSet_vox20.set(static_cast<supera::VoxelID_t>(7078131250), 1.35471);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox20), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox21;
evtOutput_part0_label_energyVoxSet_vox21.set(static_cast<supera::VoxelID_t>(7084381250), 1.35471);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox21), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox22;
evtOutput_part0_label_energyVoxSet_vox22.set(static_cast<supera::VoxelID_t>(7090631250), 1.35471);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox22), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox23;
evtOutput_part0_label_energyVoxSet_vox23.set(static_cast<supera::VoxelID_t>(7096881250), 1.34205);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox23), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox24;
evtOutput_part0_label_energyVoxSet_vox24.set(static_cast<supera::VoxelID_t>(7103131250), 0.854753);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox24), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox25;
evtOutput_part0_label_energyVoxSet_vox25.set(static_cast<supera::VoxelID_t>(7103131251), 0.2676);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox25), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox26;
evtOutput_part0_label_energyVoxSet_vox26.set(static_cast<supera::VoxelID_t>(7109381251), 1.12235);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox26), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox27;
evtOutput_part0_label_energyVoxSet_vox27.set(static_cast<supera::VoxelID_t>(7115631251), 1.12235);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox27), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox28;
evtOutput_part0_label_energyVoxSet_vox28.set(static_cast<supera::VoxelID_t>(7121881251), 1.12235);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox28), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox29;
evtOutput_part0_label_energyVoxSet_vox29.set(static_cast<supera::VoxelID_t>(7128131251), 1.12235);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox29), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox30;
evtOutput_part0_label_energyVoxSet_vox30.set(static_cast<supera::VoxelID_t>(7134381251), 1.62874);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox30), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox31;
evtOutput_part0_label_energyVoxSet_vox31.set(static_cast<supera::VoxelID_t>(7140631251), 1.1211);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox31), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox32;
evtOutput_part0_label_energyVoxSet_vox32.set(static_cast<supera::VoxelID_t>(7146881251), 1.1211);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox32), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox33;
evtOutput_part0_label_energyVoxSet_vox33.set(static_cast<supera::VoxelID_t>(7153131251), 1.07164);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox33), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox34;
evtOutput_part0_label_energyVoxSet_vox34.set(static_cast<supera::VoxelID_t>(7159381251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox34), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox35;
evtOutput_part0_label_energyVoxSet_vox35.set(static_cast<supera::VoxelID_t>(7165631251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox35), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox36;
evtOutput_part0_label_energyVoxSet_vox36.set(static_cast<supera::VoxelID_t>(7171881251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox36), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox37;
evtOutput_part0_label_energyVoxSet_vox37.set(static_cast<supera::VoxelID_t>(7178131251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox37), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox38;
evtOutput_part0_label_energyVoxSet_vox38.set(static_cast<supera::VoxelID_t>(7184381251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox38), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox39;
evtOutput_part0_label_energyVoxSet_vox39.set(static_cast<supera::VoxelID_t>(7190631251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox39), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox40;
evtOutput_part0_label_energyVoxSet_vox40.set(static_cast<supera::VoxelID_t>(7196881251), 1.04203);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox40), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox41;
evtOutput_part0_label_energyVoxSet_vox41.set(static_cast<supera::VoxelID_t>(7203131251), 1.74108);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox41), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox42;
evtOutput_part0_label_energyVoxSet_vox42.set(static_cast<supera::VoxelID_t>(7209381251), 1.04407);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox42), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox43;
evtOutput_part0_label_energyVoxSet_vox43.set(static_cast<supera::VoxelID_t>(7215631251), 1.04407);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox43), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox44;
evtOutput_part0_label_energyVoxSet_vox44.set(static_cast<supera::VoxelID_t>(7221881251), 1.04407);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox44), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox45;
evtOutput_part0_label_energyVoxSet_vox45.set(static_cast<supera::VoxelID_t>(7228131251), 1.32907);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox45), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox46;
evtOutput_part0_label_energyVoxSet_vox46.set(static_cast<supera::VoxelID_t>(7234381251), 0.918773);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox46), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox47;
evtOutput_part0_label_energyVoxSet_vox47.set(static_cast<supera::VoxelID_t>(7240631251), 0.918773);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox47), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox48;
evtOutput_part0_label_energyVoxSet_vox48.set(static_cast<supera::VoxelID_t>(7246881251), 0.918773);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox48), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox49;
evtOutput_part0_label_energyVoxSet_vox49.set(static_cast<supera::VoxelID_t>(7253128751), 0.157791);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox49), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox50;
evtOutput_part0_label_energyVoxSet_vox50.set(static_cast<supera::VoxelID_t>(7253131251), 1.11712);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox50), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox51;
evtOutput_part0_label_energyVoxSet_vox51.set(static_cast<supera::VoxelID_t>(7259378751), 0.318438);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox51), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox52;
evtOutput_part0_label_energyVoxSet_vox52.set(static_cast<supera::VoxelID_t>(7259381251), 1.65497);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox52), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox53;
evtOutput_part0_label_energyVoxSet_vox53.set(static_cast<supera::VoxelID_t>(7265631251), 0.893369);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox53), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox54;
evtOutput_part0_label_energyVoxSet_vox54.set(static_cast<supera::VoxelID_t>(7271881251), 0.893369);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox54), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox55;
evtOutput_part0_label_energyVoxSet_vox55.set(static_cast<supera::VoxelID_t>(7278128751), 0.0672888);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox55), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox56;
evtOutput_part0_label_energyVoxSet_vox56.set(static_cast<supera::VoxelID_t>(7278131251), 0.82608);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox56), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox57;
evtOutput_part0_label_energyVoxSet_vox57.set(static_cast<supera::VoxelID_t>(7284378751), 0.893369);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox57), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox58;
evtOutput_part0_label_energyVoxSet_vox58.set(static_cast<supera::VoxelID_t>(7290628751), 0.893369);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox58), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox59;
evtOutput_part0_label_energyVoxSet_vox59.set(static_cast<supera::VoxelID_t>(7296878750), 0.103992);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox59), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox60;
evtOutput_part0_label_energyVoxSet_vox60.set(static_cast<supera::VoxelID_t>(7296878751), 0.789377);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox60), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox61;
evtOutput_part0_label_energyVoxSet_vox61.set(static_cast<supera::VoxelID_t>(7303128750), 1.45704);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox61), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox62;
evtOutput_part0_label_energyVoxSet_vox62.set(static_cast<supera::VoxelID_t>(7303131250), 0.369589);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox62), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox63;
evtOutput_part0_label_energyVoxSet_vox63.set(static_cast<supera::VoxelID_t>(7309378750), 0.873483);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox63), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox64;
evtOutput_part0_label_energyVoxSet_vox64.set(static_cast<supera::VoxelID_t>(7315628750), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox64), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox65;
evtOutput_part0_label_energyVoxSet_vox65.set(static_cast<supera::VoxelID_t>(7321878750), 0.593697);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox65), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox66;
evtOutput_part0_label_energyVoxSet_vox66.set(static_cast<supera::VoxelID_t>(7321878751), 0.283858);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox66), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox67;
evtOutput_part0_label_energyVoxSet_vox67.set(static_cast<supera::VoxelID_t>(7328128751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox67), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox68;
evtOutput_part0_label_energyVoxSet_vox68.set(static_cast<supera::VoxelID_t>(7334378751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox68), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox69;
evtOutput_part0_label_energyVoxSet_vox69.set(static_cast<supera::VoxelID_t>(7340628751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox69), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox70;
evtOutput_part0_label_energyVoxSet_vox70.set(static_cast<supera::VoxelID_t>(7346878751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox70), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox71;
evtOutput_part0_label_energyVoxSet_vox71.set(static_cast<supera::VoxelID_t>(7353128751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox71), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox72;
evtOutput_part0_label_energyVoxSet_vox72.set(static_cast<supera::VoxelID_t>(7359378751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox72), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox73;
evtOutput_part0_label_energyVoxSet_vox73.set(static_cast<supera::VoxelID_t>(7365628751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox73), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox74;
evtOutput_part0_label_energyVoxSet_vox74.set(static_cast<supera::VoxelID_t>(7371878751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox74), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox75;
evtOutput_part0_label_energyVoxSet_vox75.set(static_cast<supera::VoxelID_t>(7378128751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox75), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox76;
evtOutput_part0_label_energyVoxSet_vox76.set(static_cast<supera::VoxelID_t>(7384378751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox76), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox77;
evtOutput_part0_label_energyVoxSet_vox77.set(static_cast<supera::VoxelID_t>(7390628751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox77), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox78;
evtOutput_part0_label_energyVoxSet_vox78.set(static_cast<supera::VoxelID_t>(7396878751), 0.877555);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox78), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox79;
evtOutput_part0_label_energyVoxSet_vox79.set(static_cast<supera::VoxelID_t>(7403128751), 1.60019);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox79), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox80;
evtOutput_part0_label_energyVoxSet_vox80.set(static_cast<supera::VoxelID_t>(7409378751), 1.41842);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox80), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox81;
evtOutput_part0_label_energyVoxSet_vox81.set(static_cast<supera::VoxelID_t>(7415628751), 1.00007);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox81), false);
supera::Voxel evtOutput_part0_label_energyVoxSet_vox82;
evtOutput_part0_label_energyVoxSet_vox82.set(static_cast<supera::VoxelID_t>(7421878751), 0.500036);
evtOutput_part0_label_energyVoxSet.emplace(std::move(evtOutput_part0_label_energyVoxSet_vox82), false);
evtOutput_part0_label.energy = std::move(evtOutput_part0_label_energyVoxSet);





























































































































































































































































supera::EDep evtOutput_part0_label_firstEdep;
evtOutput_part0_label_firstEdep.x = 1.79769e+308;
evtOutput_part0_label_firstEdep.y = 1.79769e+308;
evtOutput_part0_label_firstEdep.z = 1.79769e+308;
evtOutput_part0_label_firstEdep.t = 1.79769e+308;
evtOutput_part0_label_firstEdep.e = 1.79769e+308;

evtOutput_part0_label.first_pt = std::move(evtOutput_part0_label_firstEdep);
supera::EDep evtOutput_part0_label_lastEdep;
evtOutput_part0_label_lastEdep.x = 1.79769e+308;
evtOutput_part0_label_lastEdep.y = 1.79769e+308;
evtOutput_part0_label_lastEdep.z = 1.79769e+308;
evtOutput_part0_label_lastEdep.t = 1.79769e+308;
evtOutput_part0_label_lastEdep.e = 1.79769e+308;

evtOutput_part0_label.last_pt = std::move(evtOutput_part0_label_lastEdep);
evtOutput.Particles().push_back(std::move(evtOutput_part0_label));
supera::ParticleLabel evtOutput_part1_label;
supera::Particle evtOutput_part1_label_part;
evtOutput_part1_label_part.id = 1;
evtOutput_part1_label_part.shape = static_cast<supera::SemanticType_t>(4);
evtOutput_part1_label_part.trackid = 11;
evtOutput_part1_label_part.pdg = -11;
evtOutput_part1_label_part.px = 34.9792;
evtOutput_part1_label_part.py = -25.4476;
evtOutput_part1_label_part.pz = 17.6283;
evtOutput_part1_label_part.vtx = {0.14755, 1.49633, -54.5587, 10.8502};
evtOutput_part1_label_part.end_pt = {11.8708, -3.89269, -47.7023, 11.3687};
evtOutput_part1_label_part.first_step = {0.14755, 1.49633, -54.5587, 10.8502};
evtOutput_part1_label_part.last_step = {11.8708, -3.89269, -47.7023, 11.3687};
evtOutput_part1_label_part.dist_travel = 15.3603;
evtOutput_part1_label_part.energy_init = 46.7135;
evtOutput_part1_label_part.energy_deposit = 34.0317;
evtOutput_part1_label_part.process = "Decay";
evtOutput_part1_label_part.parent_trackid = 0;
evtOutput_part1_label_part.parent_pdg = -13;
evtOutput_part1_label_part.parent_vtx = {0, 0, 0, 0};
evtOutput_part1_label_part.ancestor_trackid = kINVALID_TRACKID;
evtOutput_part1_label_part.ancestor_pdg = 0;
evtOutput_part1_label_part.ancestor_vtx = {0, 0, 0, 0};
evtOutput_part1_label_part.ancestor_process = "";
evtOutput_part1_label_part.parent_process = "";
evtOutput_part1_label_part.parent_id = 0;
evtOutput_part1_label_part.children_id = { 3 };
evtOutput_part1_label_part.group_id = 0;
evtOutput_part1_label_part.interaction_id = 0;
evtOutput_part1_label.part = std::move(evtOutput_part1_label_part);
evtOutput_part1_label.valid = 1;
evtOutput_part1_label.add_to_parent = 0;
evtOutput_part1_label.type = static_cast<supera::ProcessType>(0);
evtOutput_part1_label.trackid_v = {  };
supera::VoxelSet evtOutput_part1_label_energyVoxSet;
evtOutput_part1_label_energyVoxSet.id(18446744073709551615ul);
evtOutput_part1_label_energyVoxSet.reserve(109);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox0;
evtOutput_part1_label_energyVoxSet_vox0.set(static_cast<supera::VoxelID_t>(6521606399), 0.161841);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox0), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox1;
evtOutput_part1_label_energyVoxSet_vox1.set(static_cast<supera::VoxelID_t>(6521606400), 0.02806);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox1), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox2;
evtOutput_part1_label_energyVoxSet_vox2.set(static_cast<supera::VoxelID_t>(6559101399), 0.288414);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox2), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox3;
evtOutput_part1_label_energyVoxSet_vox3.set(static_cast<supera::VoxelID_t>(6565343909), 0.0356491);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox3), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox4;
evtOutput_part1_label_energyVoxSet_vox4.set(static_cast<supera::VoxelID_t>(6584373758), 0.0620167);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox4), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox5;
evtOutput_part1_label_energyVoxSet_vox5.set(static_cast<supera::VoxelID_t>(6640378909), 0.0668633);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox5), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox6;
evtOutput_part1_label_energyVoxSet_vox6.set(static_cast<supera::VoxelID_t>(6646618910), 0.0288428);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox6), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox7;
evtOutput_part1_label_energyVoxSet_vox7.set(static_cast<supera::VoxelID_t>(6646621410), 0.0403871);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox7), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox8;
evtOutput_part1_label_energyVoxSet_vox8.set(static_cast<supera::VoxelID_t>(6696858748), 0.00389041);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox8), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox9;
evtOutput_part1_label_energyVoxSet_vox9.set(static_cast<supera::VoxelID_t>(6721811243), 0.0727482);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox9), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox10;
evtOutput_part1_label_energyVoxSet_vox10.set(static_cast<supera::VoxelID_t>(6734316243), 0.00674182);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox10), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox11;
evtOutput_part1_label_energyVoxSet_vox11.set(static_cast<supera::VoxelID_t>(6740586245), 0.00237202);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox11), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox12;
evtOutput_part1_label_energyVoxSet_vox12.set(static_cast<supera::VoxelID_t>(6746816245), 0.0199659);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox12), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox13;
evtOutput_part1_label_energyVoxSet_vox13.set(static_cast<supera::VoxelID_t>(6828043763), 0.074798);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox13), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox14;
evtOutput_part1_label_energyVoxSet_vox14.set(static_cast<supera::VoxelID_t>(6896856282), 0.0564533);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox14), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox15;
evtOutput_part1_label_energyVoxSet_vox15.set(static_cast<supera::VoxelID_t>(6903118785), 0.01029);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox15), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox16;
evtOutput_part1_label_energyVoxSet_vox16.set(static_cast<supera::VoxelID_t>(6959173748), 0.104615);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox16), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox17;
evtOutput_part1_label_energyVoxSet_vox17.set(static_cast<supera::VoxelID_t>(6959181253), 0.0192375);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox17), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox18;
evtOutput_part1_label_energyVoxSet_vox18.set(static_cast<supera::VoxelID_t>(6959383750), 1.03118);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox18), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox19;
evtOutput_part1_label_energyVoxSet_vox19.set(static_cast<supera::VoxelID_t>(6959383751), 0.226815);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox19), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox20;
evtOutput_part1_label_energyVoxSet_vox20.set(static_cast<supera::VoxelID_t>(6965423747), 0.0714654);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox20), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox21;
evtOutput_part1_label_energyVoxSet_vox21.set(static_cast<supera::VoxelID_t>(6965631251), 0.637591);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox21), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox22;
evtOutput_part1_label_energyVoxSet_vox22.set(static_cast<supera::VoxelID_t>(6965631252), 0.738192);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox22), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox23;
evtOutput_part1_label_energyVoxSet_vox23.set(static_cast<supera::VoxelID_t>(6965633751), 0.281882);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox23), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox24;
evtOutput_part1_label_energyVoxSet_vox24.set(static_cast<supera::VoxelID_t>(6971878753), 0.860127);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox24), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox25;
evtOutput_part1_label_energyVoxSet_vox25.set(static_cast<supera::VoxelID_t>(6971878754), 0.435149);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox25), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox26;
evtOutput_part1_label_energyVoxSet_vox26.set(static_cast<supera::VoxelID_t>(6971881252), 0.206413);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox26), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox27;
evtOutput_part1_label_energyVoxSet_vox27.set(static_cast<supera::VoxelID_t>(6971881253), 0.0382903);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox27), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox28;
evtOutput_part1_label_energyVoxSet_vox28.set(static_cast<supera::VoxelID_t>(6978126254), 0.168623);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox28), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox29;
evtOutput_part1_label_energyVoxSet_vox29.set(static_cast<supera::VoxelID_t>(6978126255), 0.898417);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox29), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox30;
evtOutput_part1_label_energyVoxSet_vox30.set(static_cast<supera::VoxelID_t>(6978126256), 0.173735);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox30), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox31;
evtOutput_part1_label_energyVoxSet_vox31.set(static_cast<supera::VoxelID_t>(6978128754), 0.294646);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox31), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox32;
evtOutput_part1_label_energyVoxSet_vox32.set(static_cast<supera::VoxelID_t>(6984363770), 0.00253165);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox32), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox33;
evtOutput_part1_label_energyVoxSet_vox33.set(static_cast<supera::VoxelID_t>(6984373756), 0.481244);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox33), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox34;
evtOutput_part1_label_energyVoxSet_vox34.set(static_cast<supera::VoxelID_t>(6984373757), 1.0296);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox34), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox35;
evtOutput_part1_label_energyVoxSet_vox35.set(static_cast<supera::VoxelID_t>(6984376256), 0.359229);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox35), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox36;
evtOutput_part1_label_energyVoxSet_vox36.set(static_cast<supera::VoxelID_t>(6990621258), 0.841179);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox36), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox37;
evtOutput_part1_label_energyVoxSet_vox37.set(static_cast<supera::VoxelID_t>(6990623757), 0.465983);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox37), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox38;
evtOutput_part1_label_energyVoxSet_vox38.set(static_cast<supera::VoxelID_t>(6996868759), 0.782601);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox38), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox39;
evtOutput_part1_label_energyVoxSet_vox39.set(static_cast<supera::VoxelID_t>(6996868760), 0.0230653);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox39), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox40;
evtOutput_part1_label_energyVoxSet_vox40.set(static_cast<supera::VoxelID_t>(6996871258), 0.232786);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox40), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox41;
evtOutput_part1_label_energyVoxSet_vox41.set(static_cast<supera::VoxelID_t>(6996871259), 0.525544);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox41), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox42;
evtOutput_part1_label_energyVoxSet_vox42.set(static_cast<supera::VoxelID_t>(7003103769), 0.00910878);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox42), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox43;
evtOutput_part1_label_energyVoxSet_vox43.set(static_cast<supera::VoxelID_t>(7003116260), 0.0245181);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox43), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox44;
evtOutput_part1_label_energyVoxSet_vox44.set(static_cast<supera::VoxelID_t>(7003116261), 0.304664);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox44), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox45;
evtOutput_part1_label_energyVoxSet_vox45.set(static_cast<supera::VoxelID_t>(7003118760), 1.03475);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox45), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox46;
evtOutput_part1_label_energyVoxSet_vox46.set(static_cast<supera::VoxelID_t>(7009353768), 0.000196779);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox46), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox47;
evtOutput_part1_label_energyVoxSet_vox47.set(static_cast<supera::VoxelID_t>(7009366261), 0.656857);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox47), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox48;
evtOutput_part1_label_energyVoxSet_vox48.set(static_cast<supera::VoxelID_t>(7009366262), 0.753624);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox48), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox49;
evtOutput_part1_label_energyVoxSet_vox49.set(static_cast<supera::VoxelID_t>(7015448770), 0.138336);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox49), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox50;
evtOutput_part1_label_energyVoxSet_vox50.set(static_cast<supera::VoxelID_t>(7015613763), 0.271446);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox50), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox51;
evtOutput_part1_label_energyVoxSet_vox51.set(static_cast<supera::VoxelID_t>(7015613764), 0.27755);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox51), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox52;
evtOutput_part1_label_energyVoxSet_vox52.set(static_cast<supera::VoxelID_t>(7015616262), 0.189654);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox52), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox53;
evtOutput_part1_label_energyVoxSet_vox53.set(static_cast<supera::VoxelID_t>(7015616263), 0.671832);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox53), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox54;
evtOutput_part1_label_energyVoxSet_vox54.set(static_cast<supera::VoxelID_t>(7021861265), 0.0102075);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox54), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox55;
evtOutput_part1_label_energyVoxSet_vox55.set(static_cast<supera::VoxelID_t>(7021863764), 0.665728);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox55), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox56;
evtOutput_part1_label_energyVoxSet_vox56.set(static_cast<supera::VoxelID_t>(7021863765), 0.744753);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox56), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox57;
evtOutput_part1_label_energyVoxSet_vox57.set(static_cast<supera::VoxelID_t>(7027946266), 0.0633362);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox57), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox58;
evtOutput_part1_label_energyVoxSet_vox58.set(static_cast<supera::VoxelID_t>(7028111266), 0.471581);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox58), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox59;
evtOutput_part1_label_energyVoxSet_vox59.set(static_cast<supera::VoxelID_t>(7028111267), 0.268679);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox59), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox60;
evtOutput_part1_label_energyVoxSet_vox60.set(static_cast<supera::VoxelID_t>(7028113765), 0.198525);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox60), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox61;
evtOutput_part1_label_energyVoxSet_vox61.set(static_cast<supera::VoxelID_t>(7028113766), 0.471697);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox61), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox62;
evtOutput_part1_label_energyVoxSet_vox62.set(static_cast<supera::VoxelID_t>(7034361267), 0.674599);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox62), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox63;
evtOutput_part1_label_energyVoxSet_vox63.set(static_cast<supera::VoxelID_t>(7034361268), 0.659853);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox63), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox64;
evtOutput_part1_label_energyVoxSet_vox64.set(static_cast<supera::VoxelID_t>(7040608769), 0.153304);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox64), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox65;
evtOutput_part1_label_energyVoxSet_vox65.set(static_cast<supera::VoxelID_t>(7040611268), 0.351533);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox65), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox66;
evtOutput_part1_label_energyVoxSet_vox66.set(static_cast<supera::VoxelID_t>(7040611269), 0.673454);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox66), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox67;
evtOutput_part1_label_energyVoxSet_vox67.set(static_cast<supera::VoxelID_t>(7046858769), 0.202003);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox67), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox68;
evtOutput_part1_label_energyVoxSet_vox68.set(static_cast<supera::VoxelID_t>(7046858770), 0.932787);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox68), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox69;
evtOutput_part1_label_energyVoxSet_vox69.set(static_cast<supera::VoxelID_t>(7046858771), 0.0380627);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox69), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox70;
evtOutput_part1_label_energyVoxSet_vox70.set(static_cast<supera::VoxelID_t>(7053108771), 0.885904);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox70), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox71;
evtOutput_part1_label_energyVoxSet_vox71.set(static_cast<supera::VoxelID_t>(7053108772), 0.444378);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox71), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox72;
evtOutput_part1_label_energyVoxSet_vox72.set(static_cast<supera::VoxelID_t>(7059311286), 0.114009);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox72), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox73;
evtOutput_part1_label_energyVoxSet_vox73.set(static_cast<supera::VoxelID_t>(7059353774), 0.0505778);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox73), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox74;
evtOutput_part1_label_energyVoxSet_vox74.set(static_cast<supera::VoxelID_t>(7059356272), 0.393622);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox74), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox75;
evtOutput_part1_label_energyVoxSet_vox75.set(static_cast<supera::VoxelID_t>(7059356273), 1.08568);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox75), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox76;
evtOutput_part1_label_energyVoxSet_vox76.set(static_cast<supera::VoxelID_t>(7059356274), 0.387826);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox76), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox77;
evtOutput_part1_label_energyVoxSet_vox77.set(static_cast<supera::VoxelID_t>(7059358772), 0.305673);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox77), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox78;
evtOutput_part1_label_energyVoxSet_vox78.set(static_cast<supera::VoxelID_t>(7065601276), 0.40626);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox78), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox79;
evtOutput_part1_label_energyVoxSet_vox79.set(static_cast<supera::VoxelID_t>(7065601277), 0.729074);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox79), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox80;
evtOutput_part1_label_energyVoxSet_vox80.set(static_cast<supera::VoxelID_t>(7065601279), 1.2012);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox80), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox81;
evtOutput_part1_label_energyVoxSet_vox81.set(static_cast<supera::VoxelID_t>(7065603774), 0.523848);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox81), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox82;
evtOutput_part1_label_energyVoxSet_vox82.set(static_cast<supera::VoxelID_t>(7065603775), 0.890788);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox82), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox83;
evtOutput_part1_label_energyVoxSet_vox83.set(static_cast<supera::VoxelID_t>(7065603776), 0.413808);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox83), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox84;
evtOutput_part1_label_energyVoxSet_vox84.set(static_cast<supera::VoxelID_t>(7071848779), 0.643053);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox84), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox85;
evtOutput_part1_label_energyVoxSet_vox85.set(static_cast<supera::VoxelID_t>(7071851277), 0.421759);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox85), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox86;
evtOutput_part1_label_energyVoxSet_vox86.set(static_cast<supera::VoxelID_t>(7071851278), 0.942776);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox86), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox87;
evtOutput_part1_label_energyVoxSet_vox87.set(static_cast<supera::VoxelID_t>(7071851279), 0.28144);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox87), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox88;
evtOutput_part1_label_energyVoxSet_vox88.set(static_cast<supera::VoxelID_t>(7071871275), 0.14813);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox88), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox89;
evtOutput_part1_label_energyVoxSet_vox89.set(static_cast<supera::VoxelID_t>(7078106277), 0.0516445);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox89), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox90;
evtOutput_part1_label_energyVoxSet_vox90.set(static_cast<supera::VoxelID_t>(7078123774), 0.0138404);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox90), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox91;
evtOutput_part1_label_energyVoxSet_vox91.set(static_cast<supera::VoxelID_t>(7134371266), 0.00119876);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox91), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox92;
evtOutput_part1_label_energyVoxSet_vox92.set(static_cast<supera::VoxelID_t>(7159371262), 0.0447433);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox92), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox93;
evtOutput_part1_label_energyVoxSet_vox93.set(static_cast<supera::VoxelID_t>(7171858788), 0.0906444);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox93), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox94;
evtOutput_part1_label_energyVoxSet_vox94.set(static_cast<supera::VoxelID_t>(7196846308), 0.236041);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox94), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox95;
evtOutput_part1_label_energyVoxSet_vox95.set(static_cast<supera::VoxelID_t>(7209346308), 0.00625953);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox95), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox96;
evtOutput_part1_label_energyVoxSet_vox96.set(static_cast<supera::VoxelID_t>(7209358799), 0.159273);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox96), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox97;
evtOutput_part1_label_energyVoxSet_vox97.set(static_cast<supera::VoxelID_t>(7246826247), 0.0521691);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox97), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox98;
evtOutput_part1_label_energyVoxSet_vox98.set(static_cast<supera::VoxelID_t>(7303103799), 0.00737417);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox98), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox99;
evtOutput_part1_label_energyVoxSet_vox99.set(static_cast<supera::VoxelID_t>(7353086229), 0.104007);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox99), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox100;
evtOutput_part1_label_energyVoxSet_vox100.set(static_cast<supera::VoxelID_t>(7359338731), 0.0480985);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox100), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox101;
evtOutput_part1_label_energyVoxSet_vox101.set(static_cast<supera::VoxelID_t>(7359363786), 0.0171902);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox101), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox102;
evtOutput_part1_label_energyVoxSet_vox102.set(static_cast<supera::VoxelID_t>(7365571225), 0.07185);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox102), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox103;
evtOutput_part1_label_energyVoxSet_vox103.set(static_cast<supera::VoxelID_t>(7365571226), 0.00207092);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox103), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox104;
evtOutput_part1_label_energyVoxSet_vox104.set(static_cast<supera::VoxelID_t>(7365571227), 0.0101893);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox104), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox105;
evtOutput_part1_label_energyVoxSet_vox105.set(static_cast<supera::VoxelID_t>(7371821228), 0.014702);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox105), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox106;
evtOutput_part1_label_energyVoxSet_vox106.set(static_cast<supera::VoxelID_t>(9071351347), 0.00650139);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox106), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox107;
evtOutput_part1_label_energyVoxSet_vox107.set(static_cast<supera::VoxelID_t>(9097381241), 0.00635404);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox107), false);
supera::Voxel evtOutput_part1_label_energyVoxSet_vox108;
evtOutput_part1_label_energyVoxSet_vox108.set(static_cast<supera::VoxelID_t>(9097386240), 0.0156247);
evtOutput_part1_label_energyVoxSet.emplace(std::move(evtOutput_part1_label_energyVoxSet_vox108), false);
evtOutput_part1_label.energy = std::move(evtOutput_part1_label_energyVoxSet);











































































































































































































































































































































supera::EDep evtOutput_part1_label_firstEdep;
evtOutput_part1_label_firstEdep.x = 1.79769e+308;
evtOutput_part1_label_firstEdep.y = 1.79769e+308;
evtOutput_part1_label_firstEdep.z = 1.79769e+308;
evtOutput_part1_label_firstEdep.t = 1.79769e+308;
evtOutput_part1_label_firstEdep.e = 1.79769e+308;

evtOutput_part1_label.first_pt = std::move(evtOutput_part1_label_firstEdep);
supera::EDep evtOutput_part1_label_lastEdep;
evtOutput_part1_label_lastEdep.x = 1.79769e+308;
evtOutput_part1_label_lastEdep.y = 1.79769e+308;
evtOutput_part1_label_lastEdep.z = 1.79769e+308;
evtOutput_part1_label_lastEdep.t = 1.79769e+308;
evtOutput_part1_label_lastEdep.e = 1.79769e+308;

evtOutput_part1_label.last_pt = std::move(evtOutput_part1_label_lastEdep);
evtOutput.Particles().push_back(std::move(evtOutput_part1_label));
supera::ParticleLabel evtOutput_part2_label;
supera::Particle evtOutput_part2_label_part;
evtOutput_part2_label_part.id = 2;
evtOutput_part2_label_part.shape = static_cast<supera::SemanticType_t>(4);
evtOutput_part2_label_part.trackid = 14;
evtOutput_part2_label_part.pdg = 11;
evtOutput_part2_label_part.px = 0.60628;
evtOutput_part2_label_part.py = 0.130031;
evtOutput_part2_label_part.pz = 0.783233;
evtOutput_part2_label_part.vtx = {24.8292, -20.7789, -35.1569, 12.1551};
evtOutput_part2_label_part.end_pt = {24.8921, -20.7547, -35.1034, 12.1588};
evtOutput_part2_label_part.first_step = {24.8292, -20.7789, -35.1569, 12.1551};
evtOutput_part2_label_part.last_step = {24.8921, -20.7547, -35.1034, 12.1588};
evtOutput_part2_label_part.dist_travel = 0.0860548;
evtOutput_part2_label_part.energy_init = 1.12208;
evtOutput_part2_label_part.energy_deposit = 0.611079;
evtOutput_part2_label_part.process = "compt";
evtOutput_part2_label_part.parent_trackid = 13;
evtOutput_part2_label_part.parent_pdg = 22;
evtOutput_part2_label_part.parent_vtx = {0, 0, 0, 0};
evtOutput_part2_label_part.ancestor_trackid = kINVALID_TRACKID;
evtOutput_part2_label_part.ancestor_pdg = 0;
evtOutput_part2_label_part.ancestor_vtx = {0, 0, 0, 0};
evtOutput_part2_label_part.ancestor_process = "";
evtOutput_part2_label_part.parent_process = "";
evtOutput_part2_label_part.parent_id = 2;
evtOutput_part2_label_part.children_id = {  };
evtOutput_part2_label_part.group_id = 2;
evtOutput_part2_label_part.interaction_id = 0;
evtOutput_part2_label.part = std::move(evtOutput_part2_label_part);
evtOutput_part2_label.valid = 1;
evtOutput_part2_label.add_to_parent = 0;
evtOutput_part2_label.type = static_cast<supera::ProcessType>(4);
evtOutput_part2_label.trackid_v = {  };
supera::VoxelSet evtOutput_part2_label_energyVoxSet;
evtOutput_part2_label_energyVoxSet.id(18446744073709551615ul);
evtOutput_part2_label_energyVoxSet.reserve(1);
supera::Voxel evtOutput_part2_label_energyVoxSet_vox0;
evtOutput_part2_label_energyVoxSet_vox0.set(static_cast<supera::VoxelID_t>(7265496312), 0.611079);
evtOutput_part2_label_energyVoxSet.emplace(std::move(evtOutput_part2_label_energyVoxSet_vox0), false);
evtOutput_part2_label.energy = std::move(evtOutput_part2_label_energyVoxSet);







supera::EDep evtOutput_part2_label_firstEdep;
evtOutput_part2_label_firstEdep.x = 1.79769e+308;
evtOutput_part2_label_firstEdep.y = 1.79769e+308;
evtOutput_part2_label_firstEdep.z = 1.79769e+308;
evtOutput_part2_label_firstEdep.t = 1.79769e+308;
evtOutput_part2_label_firstEdep.e = 1.79769e+308;

evtOutput_part2_label.first_pt = std::move(evtOutput_part2_label_firstEdep);
supera::EDep evtOutput_part2_label_lastEdep;
evtOutput_part2_label_lastEdep.x = 1.79769e+308;
evtOutput_part2_label_lastEdep.y = 1.79769e+308;
evtOutput_part2_label_lastEdep.z = 1.79769e+308;
evtOutput_part2_label_lastEdep.t = 1.79769e+308;
evtOutput_part2_label_lastEdep.e = 1.79769e+308;

evtOutput_part2_label.last_pt = std::move(evtOutput_part2_label_lastEdep);
evtOutput.Particles().push_back(std::move(evtOutput_part2_label));
supera::ParticleLabel evtOutput_part3_label;
supera::Particle evtOutput_part3_label_part;
evtOutput_part3_label_part.id = 3;
evtOutput_part3_label_part.shape = static_cast<supera::SemanticType_t>(4);
evtOutput_part3_label_part.trackid = 15;
evtOutput_part3_label_part.pdg = 11;
evtOutput_part3_label_part.px = 1.66336;
evtOutput_part3_label_part.py = -0.340332;
evtOutput_part3_label_part.pz = 1.85296;
evtOutput_part3_label_part.vtx = {11.518, -3.53797, -46.5186, 11.3483};
evtOutput_part3_label_part.end_pt = {11.6891, -3.68936, -46.2493, 11.3679};
evtOutput_part3_label_part.first_step = {11.518, -3.53797, -46.5186, 11.3483};
evtOutput_part3_label_part.last_step = {11.6891, -3.68936, -46.2493, 11.3679};
evtOutput_part3_label_part.dist_travel = 0.353138;
evtOutput_part3_label_part.energy_init = 2.5646;
evtOutput_part3_label_part.energy_deposit = 2.0536;
evtOutput_part3_label_part.process = "compt";
evtOutput_part3_label_part.parent_trackid = 11;
evtOutput_part3_label_part.parent_pdg = -11;
evtOutput_part3_label_part.parent_vtx = {0, 0, 0, 0};
evtOutput_part3_label_part.ancestor_trackid = kINVALID_TRACKID;
evtOutput_part3_label_part.ancestor_pdg = 0;
evtOutput_part3_label_part.ancestor_vtx = {0, 0, 0, 0};
evtOutput_part3_label_part.ancestor_process = "";
evtOutput_part3_label_part.parent_process = "";
evtOutput_part3_label_part.parent_id = 1;
evtOutput_part3_label_part.children_id = {  };
evtOutput_part3_label_part.group_id = 0;
evtOutput_part3_label_part.interaction_id = 0;
evtOutput_part3_label.part = std::move(evtOutput_part3_label_part);
evtOutput_part3_label.valid = 1;
evtOutput_part3_label.add_to_parent = 0;
evtOutput_part3_label.type = static_cast<supera::ProcessType>(4);
evtOutput_part3_label.trackid_v = {  };
supera::VoxelSet evtOutput_part3_label_energyVoxSet;
evtOutput_part3_label_energyVoxSet.id(18446744073709551615ul);
evtOutput_part3_label_energyVoxSet.reserve(4);
supera::Voxel evtOutput_part3_label_energyVoxSet_vox0;
evtOutput_part3_label_energyVoxSet_vox0.set(static_cast<supera::VoxelID_t>(7084353778), 0.286416);
evtOutput_part3_label_energyVoxSet.emplace(std::move(evtOutput_part3_label_energyVoxSet_vox0), false);
supera::Voxel evtOutput_part3_label_energyVoxSet_vox1;
evtOutput_part3_label_energyVoxSet_vox1.set(static_cast<supera::VoxelID_t>(7084353779), 0.0856212);
evtOutput_part3_label_energyVoxSet.emplace(std::move(evtOutput_part3_label_energyVoxSet_vox1), false);
supera::Voxel evtOutput_part3_label_energyVoxSet_vox2;
evtOutput_part3_label_energyVoxSet_vox2.set(static_cast<supera::VoxelID_t>(7090601279), 1.26474);
evtOutput_part3_label_energyVoxSet.emplace(std::move(evtOutput_part3_label_energyVoxSet_vox2), false);
supera::Voxel evtOutput_part3_label_energyVoxSet_vox3;
evtOutput_part3_label_energyVoxSet_vox3.set(static_cast<supera::VoxelID_t>(7090603779), 0.416819);
evtOutput_part3_label_energyVoxSet.emplace(std::move(evtOutput_part3_label_energyVoxSet_vox3), false);
evtOutput_part3_label.energy = std::move(evtOutput_part3_label_energyVoxSet);
















supera::EDep evtOutput_part3_label_firstEdep;
evtOutput_part3_label_firstEdep.x = 1.79769e+308;
evtOutput_part3_label_firstEdep.y = 1.79769e+308;
evtOutput_part3_label_firstEdep.z = 1.79769e+308;
evtOutput_part3_label_firstEdep.t = 1.79769e+308;
evtOutput_part3_label_firstEdep.e = 1.79769e+308;

evtOutput_part3_label.first_pt = std::move(evtOutput_part3_label_firstEdep);
supera::EDep evtOutput_part3_label_lastEdep;
evtOutput_part3_label_lastEdep.x = 1.79769e+308;
evtOutput_part3_label_lastEdep.y = 1.79769e+308;
evtOutput_part3_label_lastEdep.z = 1.79769e+308;
evtOutput_part3_label_lastEdep.t = 1.79769e+308;
evtOutput_part3_label_lastEdep.e = 1.79769e+308;

evtOutput_part3_label.last_pt = std::move(evtOutput_part3_label_lastEdep);
evtOutput.Particles().push_back(std::move(evtOutput_part3_label));