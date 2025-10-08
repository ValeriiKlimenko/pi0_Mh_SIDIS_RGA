double dppC(float Px, float Py, float Pz, int sec, int ivec, int corON) {
    
    // 'Px'/'Py'/'Pz'   ==> Corresponds to the Cartesian Components of the particle momentum being corrected
    // 'sec'            ==> Corresponds to the Forward Detector Sectors where the given particle is detected (6 total)
    // 'ivec'           ==> Corresponds to the particle being corrected (See below)    
        // (*) ivec = 0 --> Electron Corrections
        // (*) ivec = 1 --> Pi+ Corrections
    // 'corON' ==> Controls which version of the particle correction is used
        // Includes:
            // (*) Correction On/Off
            // (*) Pass Version
            // (*) Data Set (Experimental or Monte Carlo)
        // corON == 0 --> DOES NOT apply the momentum corrections (i.e., turns the corrections 'off')
        // corON == 1 --> Fall  2018 - Pass 1 (Experimental Data)
        // corON == 2 --> Applies the (Pass 1) momentum corrections for the Monte Carlo (simulated) data
        // corON == 3 --> Fall  2018 - Pass 2 (Experimental Data)
        // corON == 4 --> Applies the (Pass 2) momentum corrections for the Monte Carlo (simulated) data
            // Not Available as of 8/30/2024

    // Momentum Magnitude
    double pp = sqrt(Px*Px + Py*Py + Pz*Pz);

    // Initializing the correction factor
    double dp = 0;

    // Defining Phi Angle
    double Phi = (180/3.1415926)*atan2(Py, Px);

    // Central Detector Corrections Not Included (Yet)

    // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
    if(((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)){
        Phi += 360;
    }

    // Getting Local Phi Angle
    double PhiLocal = Phi - (sec - 1)*60;

    // Applying Shift Functions to Phi Angles (local shifted phi = phi)
    double phi = PhiLocal;

    // For Electron Shift
    if(ivec == 0){
        phi = PhiLocal - 30/pp;
    }


    if(corON == 3 &&ivec == 0){ // Fall 2018 - Pass 2 Corrections
                if(sec == 1){
                    dp            =                ((-9.82416e-06)*phi*phi +            (-2.29956e-05)*phi +  (0.00029664199999999996))*pp*pp +           ((0.0001113414)*phi*phi +  (-2.041300000000001e-05)*phi +            (-0.00862226))*pp +            ((-0.000281738)*phi*phi +             (0.00058712)*phi +              (0.0652737));
                    if(pp < 7){dp = dp +            ((-3.4001e-06)*phi*phi +             (-2.2885e-05)*phi +              (9.9705e-04))*pp*pp +             ((2.1840e-05)*phi*phi +              (2.4238e-04)*phi +             (-0.0091904))*pp +             ((-2.9180e-05)*phi*phi +            (-6.4496e-04)*phi +               (0.022505));}
                    else{      dp = dp +            ((-6.3656e-05)*phi*phi +              (1.7266e-04)*phi +              (-0.0017909))*pp*pp +                ((0.00104)*phi*phi +              (-0.0028401)*phi +                (0.02981))*pp +              ((-0.0041995)*phi*phi +               (0.011537)*phi +                (-0.1196));}
                    dp            = dp + ((3.2780000000000006e-07)*phi*phi +              (6.7084e-07)*phi +  (-4.390000000000004e-05))*pp*pp + ((-7.230999999999999e-06)*phi*phi +            (-2.37482e-05)*phi +  (0.0004909000000000007))*pp +   ((3.285299999999999e-05)*phi*phi +            (9.63723e-05)*phi +               (-0.00115));
                }
                if(sec == 2){
                    dp            =               ((-7.741952e-06)*phi*phi + (-2.2402167000000004e-05)*phi + (-0.00042652900000000004))*pp*pp +            ((7.54079e-05)*phi*phi + (-1.3333999999999984e-05)*phi +  (0.0002420100000000004))*pp +            ((-0.000147876)*phi*phi +             (0.00057905)*phi +              (0.0253551));
                    if(pp < 7){dp = dp +             ((5.3611e-06)*phi*phi +              (8.1979e-06)*phi +              (5.9789e-04))*pp*pp +            ((-4.8185e-05)*phi*phi +             (-1.5188e-04)*phi +             (-0.0084675))*pp +              ((9.2324e-05)*phi*phi +             (6.4420e-04)*phi +               (0.026792));}
                    else{      dp = dp +            ((-6.1139e-05)*phi*phi +              (5.4087e-06)*phi +              (-0.0021284))*pp*pp +              ((0.0010007)*phi*phi +              (9.3492e-05)*phi +               (0.039813))*pp +              ((-0.0040434)*phi*phi +             (-0.0010953)*phi +               (-0.18112));}
                    dp            = dp +           ((6.221217e-07)*phi*phi +  (1.9596000000000003e-06)*phi +              (-9.826e-05))*pp*pp +           ((-1.28576e-05)*phi*phi +            (-4.36589e-05)*phi +             (0.00130342))*pp +             ((5.80399e-05)*phi*phi +            (0.000215388)*phi + (-0.0040414000000000005));
                }
                if(sec == 3){
                    dp            =      ((-5.115364000000001e-06)*phi*phi + (-1.1983000000000004e-05)*phi +  (-0.0006832899999999999))*pp*pp +            ((4.52287e-05)*phi*phi +  (0.00020855000000000003)*phi +  (0.0034986999999999996))*pp +  ((-9.044610000000001e-05)*phi*phi +            (-0.00106657)*phi +   (0.017954199999999997));
                    if(pp < 7){dp = dp +             ((9.9281e-07)*phi*phi +              (3.4879e-06)*phi +               (0.0011673))*pp*pp +            ((-2.0071e-05)*phi*phi +             (-3.1362e-05)*phi +              (-0.012329))*pp +              ((6.9463e-05)*phi*phi +             (3.5102e-05)*phi +               (0.037505));}
                    else{      dp = dp +            ((-3.2178e-06)*phi*phi +              (4.0630e-05)*phi +               (-0.005209))*pp*pp +             ((2.0884e-05)*phi*phi +             (-6.8800e-04)*phi +               (0.086513))*pp +              ((3.9530e-05)*phi*phi +              (0.0029306)*phi +                (-0.3507));}
                    dp            = dp + ((-4.045999999999999e-07)*phi*phi + (-1.3115999999999994e-06)*phi +  (3.9510000000000006e-05))*pp*pp +              ((5.521e-06)*phi*phi +  (2.4436999999999997e-05)*phi +             (-0.0016887))*pp + ((-1.0962999999999997e-05)*phi*phi +           (-0.000151944)*phi +   (0.009313599999999998));
                }
                if(sec == 4){
                    dp            =     ((-3.9278116999999996e-06)*phi*phi +  (2.2289300000000004e-05)*phi +  (0.00012665000000000002))*pp*pp + ((4.8649299999999995e-05)*phi*phi +             (-0.00012554)*phi +  (-0.005955500000000001))*pp + ((-0.00014617199999999997)*phi*phi +            (-0.00028571)*phi +              (0.0606998));
                    if(pp < 7){dp = dp +            ((-4.8455e-06)*phi*phi +             (-1.2074e-05)*phi +               (0.0013221))*pp*pp +             ((3.2207e-05)*phi*phi +              (1.3144e-04)*phi +              (-0.010451))*pp +             ((-3.7365e-05)*phi*phi +            (-4.2344e-04)*phi +               (0.019952));}
                    else{      dp = dp +            ((-3.9554e-05)*phi*phi +              (5.5496e-06)*phi +              (-0.0058293))*pp*pp +             ((6.5077e-04)*phi*phi +              (2.6735e-05)*phi +               (0.095025))*pp +              ((-0.0026457)*phi*phi +            (-6.1394e-04)*phi +                (-0.3793));}
                    dp            = dp +          ((-4.593089e-07)*phi*phi +             (1.40673e-05)*phi +                (6.69e-05))*pp*pp +             ((4.0239e-06)*phi*phi +            (-0.000180863)*phi + (-0.0008272199999999999))*pp + ((-5.1310000000000005e-06)*phi*phi +             (0.00049748)*phi +             (0.00255231));
                }
                if(sec == 5){
                    dp            =       ((8.036599999999999e-07)*phi*phi +             (2.58072e-05)*phi +             (0.000360217))*pp*pp + ((-9.932400000000002e-06)*phi*phi +           (-0.0005168531)*phi +              (-0.010904))*pp +  ((1.8516299999999998e-05)*phi*phi +  (0.0015570900000000001)*phi +               (0.066493));
                    if(pp < 7){dp = dp +             ((7.7156e-07)*phi*phi +             (-3.9566e-05)*phi +             (-2.3589e-04))*pp*pp +            ((-9.8309e-06)*phi*phi +              (3.7353e-04)*phi +              (0.0020382))*pp +              ((2.9506e-05)*phi*phi +            (-8.0409e-04)*phi +             (-0.0045615));}
                    else{      dp = dp +            ((-3.2410e-05)*phi*phi +             (-4.3301e-05)*phi +              (-0.0028742))*pp*pp +             ((5.3787e-04)*phi*phi +              (6.8921e-04)*phi +               (0.049578))*pp +              ((-0.0021955)*phi*phi +             (-0.0027698)*phi +               (-0.21142));}
                    dp            = dp +            ((-1.2151e-06)*phi*phi +             (-8.5087e-06)*phi +               (4.968e-05))*pp*pp +            ((1.46998e-05)*phi*phi +             (0.000115047)*phi +            (-0.00039269))*pp + ((-4.0368600000000005e-05)*phi*phi +            (-0.00037078)*phi +             (0.00073998));
                }
                if(sec == 6){
                    dp            =     ((-1.9552099999999998e-06)*phi*phi +   (8.042199999999997e-06)*phi + (-2.1324000000000028e-05))*pp*pp + ((1.6969399999999997e-05)*phi*phi +  (-6.306600000000001e-05)*phi +            (-0.00485568))*pp +             ((-2.7723e-05)*phi*phi + (-6.828400000000003e-05)*phi +              (0.0447535));
                    if(pp < 7){dp = dp +            ((-8.2535e-07)*phi*phi +              (9.1433e-06)*phi +              (3.5395e-04))*pp*pp +            ((-3.4272e-06)*phi*phi +             (-1.3012e-04)*phi +             (-0.0030724))*pp +              ((4.9211e-05)*phi*phi +             (4.5807e-04)*phi +              (0.0058932));}
                    else{      dp = dp +            ((-4.9760e-05)*phi*phi +             (-7.2903e-05)*phi +              (-0.0020453))*pp*pp +             ((8.0918e-04)*phi*phi +               (0.0011688)*phi +               (0.037042))*pp +              ((-0.0032504)*phi*phi +             (-0.0046169)*phi +               (-0.16331));}
                    dp            = dp + ((-7.153000000000002e-07)*phi*phi +             (1.62859e-05)*phi +               (8.129e-05))*pp*pp + ((7.2249999999999994e-06)*phi*phi +            (-0.000178946)*phi + (-0.0009485399999999999))*pp + ((-1.3018000000000003e-05)*phi*phi + (0.00046643000000000005)*phi +             (0.00266508));
                }
            }
    return dp/pp;
        }