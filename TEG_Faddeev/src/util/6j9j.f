! 6 and 9 J symbol, useless
      SUBROUTINE SIXJ(XJ1,XJ2,XJ3,XL1,XL2,XL3,C6J)                      TRO01770
C      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FLOG(301)                                               TRO01780
C                                                                       TRO01790
C     MISE EN DATA DES LOG(FACTORIELLE) ET DE EPS                       TRO01800
C                                                                       TRO01810
      DATA(FLOG(I),I=2,31)/0.D0,.69314718D0,1.7917595D0,3.1780538D0,
     A4.7874917D0,6.
     15792511D0,8.5251613D0,10.604603D0,12.801827D0,15.104413D0,
     B17.502307D0,19.98721
     24D0,22.552163D0,25.191221D0,27.899271D0,30.671860D0,33.505072D0,
     C36.395445D0,39.3
     339884D0,42.335616D0,45.380139D0,48.471180D0,51.606674D0,54.
     D784729D0,58.003604D0,
     461.261702D0,64.557537D0,67.889743D0,71.257038D0,74.658235D0/   
      DATA(FLOG(I),I=32,61)/78.092223D0,81.557959D0,85.054466D0,88.
     A580827D0,92.1
     136175D0,95.719694D0,99.330612D0,102.96820D0,106.63176D0,110.
     B32064D0,114.03421D0,
     2117.77188D0,121.53308D0,125.31727D0,129.12393D0,132.95257D0,
     C136.80272D0,140.67
     3392D0,144.56574D0,148.47776D0,152.40959D0,156.36083D0,160.3311
     D2D0,164.32011D0,16
     48.32744D0,172.35279D0,176.39584D0,180.45629D0,184.53383D0,
     E188.62817D0/        
      DATA(FLOG(I),I=62,91)/192.73904D0,196.86618D0,201.00931D0,
     A205.16820D0,209.
     134258D0,213.53224D0,217.73693D0,221.95644D0,226.19054D0,
     B230.43904D0,234.70172D0,
     2238.97839D0,243.26885D0,247.57291D0,251.89040D0,256.22113D0,
     C260.56494D0,264.92
     3164D0,269.29110D0,273.67312D0,278.06757D0,282.47429D0,286.
     D89313D0,291.32394D0,29
     45.76659D0,300.22094D0,304.68685D0,309.16419D0,313.65283D0,
     E318.15264D0/        
      DATA(FLOG(I),I=92,121)/322.66349D0,327.18529D0,331.71788D0,
     A336.26118D0,340
     1.81505D0,345.37940D0,349.95411D0,354.53908D0,359.13420D0,363.
     B73937D0,368.35449
     2D0,372.97946D0,377.61419D0,382.25859D0,386.91255D0,391.57598D0,
     C396.24881D0,400.9
     33094D0,405.62230D0,410.32277D0,415.03230D0,419.75080D0,424.47819
     DD0,429.21439D0,4
     433.95932D0,438.71291D0,443.47508D0,448.24576D0,453.02489D0,
     E457.81238D0/       
      DATA(FLOG(I),I=122,151)/462.60817D0,467.41220D0,472.22438D0,477.
     A04466D0,48
     11.87298D0,486.70926D0,491.55345D0,496.40547D0,501.26529D0,506.
     B13282D0,511.0080
     22D0,515.89082D0,520.78117D0,525.67901D0,530.58428D0,535.49694D0,
     C540.41692D0,545.
     334417D0,550.27865D0,555.22029D0,560.16905D0,565.12488D0,570.08772
     DD0,575.05753D0,
     4580.03427D0,585.01787D0,590.00830D0,595.00552D0,600.00946D0,
     E605.02010D0/      
      DATA(FLOG(I),I=152,181)/610.03738D0,615.06126D0,620.09170D0,
     A625.12866D0,63
     10.17208D0,635.22193D0,640.27818D0,645.34077D0,650.40968D0,655.
     B48486D0,660.5662
     26D0,665.65385D0,670.74760D0,675.84747D0,680.95341D0,686.06541D0,
     C691.18340D0,696.
     330735D0,701.43726D0,706.57306D0,711.71472D0,716.86221D0,722.
     D01551D0,727.17456D0,
     4732.33934D0,737.50983D0,742.68598D0,747.86776D0,753.05516D0,
     E758.24811D0/      
      DATA(FLOG(I),I=182,211)/763.44661D0,768.65061D0,773.86010D0,
     A779.07503D0,78
     14.29539D0,789.52114D0,794.75224D0,799.98869D0,805.23044D0,810.
     B47747D0,815.7297
     23D0,820.98722D0,826.24991D0,831.51778D0,836.79078D0,842.06890D0,
     C847.35209D0,852.
     364036D0,857.93366D0,863.23199D0,868.53529D0,873.84356D0,879.
     D15676D0,884.47488D0,
     4889.79789D0,895.12577D0,900.45848D0,905.79603D0,911.13836D0,
     E916.48547D0/      
      DATA(FLOG(I),I=212,241)/921.83732D0,927.19391D0,932.55521D0,
     A937.92118D0,94
     13.29181D0,948.66710D0,954.04699D0,959.43148D0,964.82056D0,970.
     B21419D0,975.6123
     25D0,981.01503D0,986.42220D0,991.83385D0,997.24995D0,1002.6705D0,
     C1008.0954D0,1013
     3.5248D0,1018.9585D0,1024.3966D0,1029.8389D0,1035.2857D0,1040.
     D7367D0,1046.1920D0,
     41051.6516D0,1057.1155D0,1062.5836D0,1068.0558D0,1073.5323D0,
     E1079.0129D0/      
      DATA(FLOG(I),I=242,271)/1084.4977D0,1089.9866D0,1095.4797D0,1100.
     A9768D0,11
     106.4781D0,1111.9834D0,1117.4928D0,1123.0063D0,1128.5237D0,1134.
     B0452D0,1139.570
     26D0,1145.1001D0,1150.6335D0,1156.1708D0,1161.7120D0,1167.2573D0,
     C1172.8063D0,1178
     3.3593D0,1183.9161D0,1189.4768D0,1195.0413D0,1200.6097D0,1206.
     D1818D0,1211.7577D0,
     41217.3375D0,1222.9209D0,1228.5082D0,1234.0992D0,1239.6939D0,
     E1245.2924D0/      
      DATA(FLOG(I),I=272,301)/1250.8944D0,1256.5003D0,1262.1097D0,
     A1267.7228D0,12
     173.3396D0,1278.9600D0,1284.5840D0,1290.2117D0,1295.8429D0,
     B1301.4777D0,1307.116
     20D0,1312.7580D0,1318.4034D0,1324.0524D0,1329.7048D0,1335.3609D0,
     C1341.0203D0,1346
     3.6833D0,1352.3497D0,1358.0196D0,1363.6929D0,1369.3697D0,
     D1375.0499D0,1380.7334D0,
     41386.4204D0,1392.1107D0,1397.8045D0,1403.5016D0,1409.2020D0,
     E1414.9058D0/      
      DATA EPS1,EPS2/.1D0,-.2D0/                                            TRO02320
C                                                                       TRO02330
C     CALCUL DES COMBINAISONS J,L                                       TRO02340
C                                                                       TRO02350
      XN=-XJ1+XJ2+XJ3+EPS1                                              TRO02360
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02370
      N1=XN                                                             TRO02380
      XN=-XL1+XL2+XJ3+EPS1                                              TRO02390
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02400
      N2=XN                                                             TRO02410
      XN=-XL1+XJ2+XL3+EPS1                                              TRO02420
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02430
      N3=XN                                                             TRO02440
      XN=-XJ1+XL2+XL3+EPS1                                              TRO02450
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02460
      N4=XN                                                             TRO02470
      XN=XJ1-XJ2+XJ3+EPS1                                               TRO02480
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02490
      N5=XN                                                             TRO02500
      XN=XL1-XL2+XJ3+EPS1                                               TRO02510
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02520
      N6=XN                                                             TRO02530
      XN=XL1-XJ2+XL3+EPS1                                               TRO02540
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02550
      N7=XN                                                             TRO02560
      XN=XJ1-XL2+XL3+EPS1                                               TRO02570
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02580
      N8=XN                                                             TRO02590
      XN=XJ1+XJ2-XJ3+EPS1                                               TRO02600
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02610
      N9=XN                                                             TRO02620
      XN=XL1+XL2-XJ3+EPS1                                               TRO02630
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02640
      N10=XN                                                            TRO02650
      XN=XL1+XJ2-XL3+EPS1                                               TRO02660
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02670
      N11=XN                                                            TRO02680
      XN=XJ1+XL2-XL3+EPS1                                               TRO02690
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02700
      N12=XN                                                            TRO02710
      XN=-XJ1-XL1+XJ3+XL3+EPS1                                          TRO02720
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02730
      N13=XN                                                            TRO02740
      XN=-XJ2-XL2+XJ3+XL3+EPS1                                          TRO02750
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02760
      N14=XN                                                            TRO02770
      XN=XJ1+XL1+XJ2+XL2+EPS1                                           TRO02780
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02790
      N15=XN                                                            TRO02800
      N15=N15+1                                                         TRO02810
      XN=XJ1+XJ2+XJ3+EPS1                                               TRO02820
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02830
      N16=XN                                                            TRO02840
      N16=N16+1                                                         TRO02850
      XN=XL1+XL2+XJ3+EPS1                                               TRO02860
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02870
      N17=XN                                                            TRO02880
      N17=N17+1                                                         TRO02890
      XN=XL1+XJ2+XL3+EPS1                                               TRO02900
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02910
      N18=XN                                                            TRO02920
      N18=N18+1                                                         TRO02930
      XN=XJ1+XL2+XL3+EPS1                                               TRO02940
      IF (XN.LT.0.D0) XN=XN+EPS2                                          TRO02950
      N19=XN                                                            TRO02960
      N19=N19+1                                                         TRO02970
C                                                                       TRO02980
C     TEST SUR LES J ET L                                               TRO02990
C                                                                       TRO03000
      IF(N9.LT.0) GO TO 50                                              TRO03010
      IF(N5.LT.0) GO TO 50                                              TRO03020
      IF(N1.LT.0) GO TO 50                                              TRO03030
      IF(N10.LT.0) GO TO 50                                             TRO03040
      IF(N6.LT.0) GO TO 50                                              TRO03050
      IF(N2.LT.0) GO TO 50                                              TRO03060
      IF(N11.LT.0) GO TO 50                                             TRO03070
      IF(N7.LT.0) GO TO 50                                              TRO03080
      IF(N3.LT.0) GO TO 50                                              TRO03090
      IF(N12.LT.0) GO TO 50                                             TRO03100
      IF(N8.LT.0) GO TO 50                                              TRO03110
      IF(N4.LT.0) GO TO 50                                              TRO03120
      K=N17+N18+N19-3                                                   TRO03130
      IF(K.GT.0) GO TO 54                                               TRO03140
      C6J=1.D0                                                            TRO03150
      RETURN                                                            TRO03160
   50 C6J=0.D0                                                            TRO03170
      RETURN                                                            TRO03180
C                                                                       TRO03190
C      CALCUL DE LA SOMME ALTERNEE                                      TRO03200
C                                                                       TRO03210
   54 K=0                                                               TRO03220
      L=-N13                                                            TRO03230
      IF(L.GT.K) K=L                                                    TRO03240
      L=-N14                                                            TRO03250
      IF(L.GT.K) K=L                                                    TRO03260
      L=N9                                                              TRO03270
      IF(N10.LT.L) L=N10                                                TRO03280
      IF(N11.LT.L) L=N11                                                TRO03290
      IF(N12.LT.L) L=N12                                                TRO03300
      IF(N15.LT.L) L=N15                                                TRO03310
      F=1.D0                                                              TRO03320
      S=1.D0                                                              TRO03330
      I=K+1                                                             TRO03340
   62 IF(I.GT.L) GO TO 80                                               TRO03350
      IM1=I-1                                                           TRO03360
      NN=(N9-IM1)*(N10-IM1)*(N11-IM1)*(N12-IM1)                         TRO03370
      ND=I*(N13+I)*(N14+I)*(N15-IM1)                                    TRO03380
      F=-F*DFLOAT(NN)/DFLOAT(ND)                                          TRO03390
      S=S+F                                                             TRO03400
      I=I+1                                                             TRO03410
      GO TO 62                                                          TRO03420
C                                                                       TRO03430
C     CALCUL DE LA RACINE                                               TRO03440
C                                                                       TRO03450
   80 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N3+1)+FLOG(N4+1)+FLOG(N5+1)+FLOG(N6TRO03460
     1+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N9+1)+FLOG(N10+1)+FLOG(N11+1)+FLOG(TRO03470
     2N12+1)                                                            TRO03480
      C2N=.5D0*C2N                                                        TRO03490
      C2D=FLOG(N16+1)+FLOG(N17+1)+FLOG(N18+1)+FLOG(N19+1)               TRO03500
      C2D=.5D0*C2D                                                        TRO03510
      KM1=K-1                                                           TRO03520
      KP1=K+1                                                           TRO03530
      C2N=C2N+FLOG(N15-KM1)                                             TRO03540
      C2D=C2D+FLOG(KP1)+FLOG(N13+KP1)+FLOG(N14+KP1)+FLOG(N9-KM1)+FLOG(N1TRO03550
     10-KM1)+FLOG(N11-KM1)+FLOG(N12-KM1)                                TRO03560
C                                                                       TRO03570
C     CALCUL DU C6J SANS PHASE                                          TRO03580
C                                                                       TRO03590
      F=C2D-C2N                                                         TRO03600
      IF(F.GT.80.D0) GO TO 98                                             TRO03610
      F=C2N/C2D                                                         TRO03620
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 98                          TRO03630
      C6J=S*DEXP(C2N-C2D)                                                TRO03640
      GO TO 106                                                         TRO03650
   98 IF(S) 100,50,102                                                  TRO03660
  100 S=DLOG(-S)                                                        TRO03670
      C6J=-DEXP(S+C2N-C2D)                                               TRO03680
      GO TO 106                                                         TRO03690
  102 S=DLOG(S)                                                         TRO03700
      C6J=DEXP(S+C2N-C2D)                                                TRO03710
C                                                                       TRO03720
C     CALCUL DE LA PHASE                                                TRO03730
C                                                                       TRO03740
  106 L=N15+KM1                                                         TRO03750
      K=L/2                                                             TRO03760
      K=2*K                                                             TRO03770
      IF(L.NE.K) C6J=-C6J                                               TRO03780
      RETURN                                                            TRO03790
      END                                                               TRO03800
      SUBROUTINE NEUFJ(XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJ33,C9J)TRO03820
C      IMPLICIT REAL*4 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FLOG(301)                                               TRO03830
C                                                                       TRO03840
C     MISE EN DATA DES LOG(FACTORIELLE) ET DE EPS                       TRO03850
C                                                                       TRO03860
      DATA(FLOG(I),I=2,31)/0.D0,.69314718D0,1.7917595D0,3.1780538D0,
     A4.7874917D0,6.
     15792511D0,8.5251613D0,10.604603D0,12.801827D0,15.104413D0,
     B17.502307D0,19.98721
     24D0,22.552163D0,25.191221D0,27.899271D0,30.671860D0,33.505072D0,
     C36.395445D0,39.3
     339884D0,42.335616D0,45.380139D0,48.471180D0,51.606674D0,54.
     D784729D0,58.003604D0,
     461.261702D0,64.557537D0,67.889743D0,71.257038D0,74.658235D0/   
      DATA(FLOG(I),I=32,61)/78.092223D0,81.557959D0,85.054466D0,88.
     A580827D0,92.1
     136175D0,95.719694D0,99.330612D0,102.96820D0,106.63176D0,110.
     B32064D0,114.03421D0,
     2117.77188D0,121.53308D0,125.31727D0,129.12393D0,132.95257D0,
     C136.80272D0,140.67
     3392D0,144.56574D0,148.47776D0,152.40959D0,156.36083D0,160.3311
     D2D0,164.32011D0,16
     48.32744D0,172.35279D0,176.39584D0,180.45629D0,184.53383D0,
     E188.62817D0/        
      DATA(FLOG(I),I=62,91)/192.73904D0,196.86618D0,201.00931D0,
     A205.16820D0,209.
     134258D0,213.53224D0,217.73693D0,221.95644D0,226.19054D0,
     B230.43904D0,234.70172D0,
     2238.97839D0,243.26885D0,247.57291D0,251.89040D0,256.22113D0,
     C260.56494D0,264.92
     3164D0,269.29110D0,273.67312D0,278.06757D0,282.47429D0,286.
     D89313D0,291.32394D0,29
     45.76659D0,300.22094D0,304.68685D0,309.16419D0,313.65283D0,
     E318.15264D0/        
      DATA(FLOG(I),I=92,121)/322.66349D0,327.18529D0,331.71788D0,
     A336.26118D0,340
     1.81505D0,345.37940D0,349.95411D0,354.53908D0,359.13420D0,363.
     B73937D0,368.35449
     2D0,372.97946D0,377.61419D0,382.25859D0,386.91255D0,391.57598D0,
     C396.24881D0,400.9
     33094D0,405.62230D0,410.32277D0,415.03230D0,419.75080D0,424.47819
     DD0,429.21439D0,4
     433.95932D0,438.71291D0,443.47508D0,448.24576D0,453.02489D0,
     E457.81238D0/       
      DATA(FLOG(I),I=122,151)/462.60817D0,467.41220D0,472.22438D0,477.
     A04466D0,48
     11.87298D0,486.70926D0,491.55345D0,496.40547D0,501.26529D0,506.
     B13282D0,511.0080
     22D0,515.89082D0,520.78117D0,525.67901D0,530.58428D0,535.49694D0,
     C540.41692D0,545.
     334417D0,550.27865D0,555.22029D0,560.16905D0,565.12488D0,570.08772
     DD0,575.05753D0,
     4580.03427D0,585.01787D0,590.00830D0,595.00552D0,600.00946D0,
     E605.02010D0/      
      DATA(FLOG(I),I=152,181)/610.03738D0,615.06126D0,620.09170D0,
     A625.12866D0,63
     10.17208D0,635.22193D0,640.27818D0,645.34077D0,650.40968D0,655.
     B48486D0,660.5662
     26D0,665.65385D0,670.74760D0,675.84747D0,680.95341D0,686.06541D0,
     C691.18340D0,696.
     330735D0,701.43726D0,706.57306D0,711.71472D0,716.86221D0,722.
     D01551D0,727.17456D0,
     4732.33934D0,737.50983D0,742.68598D0,747.86776D0,753.05516D0,
     E758.24811D0/      
      DATA(FLOG(I),I=182,211)/763.44661D0,768.65061D0,773.86010D0,
     A779.07503D0,78
     14.29539D0,789.52114D0,794.75224D0,799.98869D0,805.23044D0,810.
     B47747D0,815.7297
     23D0,820.98722D0,826.24991D0,831.51778D0,836.79078D0,842.06890D0,
     C847.35209D0,852.
     364036D0,857.93366D0,863.23199D0,868.53529D0,873.84356D0,879.
     D15676D0,884.47488D0,
     4889.79789D0,895.12577D0,900.45848D0,905.79603D0,911.13836D0,
     E916.48547D0/      
      DATA(FLOG(I),I=212,241)/921.83732D0,927.19391D0,932.55521D0,
     A937.92118D0,94
     13.29181D0,948.66710D0,954.04699D0,959.43148D0,964.82056D0,970.
     B21419D0,975.6123
     25D0,981.01503D0,986.42220D0,991.83385D0,997.24995D0,1002.6705D0,
     C1008.0954D0,1013
     3.5248D0,1018.9585D0,1024.3966D0,1029.8389D0,1035.2857D0,1040.
     D7367D0,1046.1920D0,
     41051.6516D0,1057.1155D0,1062.5836D0,1068.0558D0,1073.5323D0,
     E1079.0129D0/      
      DATA(FLOG(I),I=242,271)/1084.4977D0,1089.9866D0,1095.4797D0,1100.
     A9768D0,11
     106.4781D0,1111.9834D0,1117.4928D0,1123.0063D0,1128.5237D0,1134.
     B0452D0,1139.570
     26D0,1145.1001D0,1150.6335D0,1156.1708D0,1161.7120D0,1167.2573D0,
     C1172.8063D0,1178
     3.3593D0,1183.9161D0,1189.4768D0,1195.0413D0,1200.6097D0,1206.
     D1818D0,1211.7577D0,
     41217.3375D0,1222.9209D0,1228.5082D0,1234.0992D0,1239.6939D0,
     E1245.2924D0/      
      DATA(FLOG(I),I=272,301)/1250.8944D0,1256.5003D0,1262.1097D0,
     A1267.7228D0,12
     173.3396D0,1278.9600D0,1284.5840D0,1290.2117D0,1295.8429D0,
     B1301.4777D0,1307.116
     20D0,1312.7580D0,1318.4034D0,1324.0524D0,1329.7048D0,1335.3609D0,
     C1341.0203D0,1346
     3.6833D0,1352.3497D0,1358.0196D0,1363.6929D0,1369.3697D0,
     D1375.0499D0,1380.7334D0,
     41386.4204D0,1392.1107D0,1397.8045D0,1403.5016D0,1409.2020D0,
     E1414.9058D0/      
      DATA EPS1,EPS2/.1D0,-.2D0/                                            TRO04370
C                                                                       TRO04380
C     CALCUL DES COMBINAISONS XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJTRO04390
C                                                                       TRO04400
      XN=-XJ11+XJ21+XJ31+EPS1                                           TRO04410
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04420
      N1=XN                                                             TRO04430
      XN=-XJ32+XJ33+XJ31+EPS1                                           TRO04440
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04450
      N2=XN                                                             TRO04460
      XN=XJ11-XJ21+XJ31+EPS1                                            TRO04470
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04480
      N5=XN                                                             TRO04490
      XN=XJ32-XJ33+XJ31+EPS1                                            TRO04500
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04510
      N6=XN                                                             TRO04520
      XN=XJ11+XJ21-XJ31+EPS1                                            TRO04530
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04540
      N9=XN                                                             TRO04550
      XN=XJ32+XJ33-XJ31+EPS1                                            TRO04560
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04570
      N10=XN                                                            TRO04580
      XN=XJ11+XJ32+XJ21+XJ33+EPS1                                       TRO04590
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04600
      N15=XN                                                            TRO04610
      N15=N15+1                                                         TRO04620
      XN=XJ11+XJ21+XJ31+EPS1                                            TRO04630
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04640
      N16=XN                                                            TRO04650
      N16=N16+1                                                         TRO04660
      XN=XJ32+XJ33+XJ31+EPS1                                            TRO04670
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04680
      N17=XN                                                            TRO04690
      N17=N17+1                                                         TRO04700
      XN=-XJ12+XJ22+XJ32+EPS1                                           TRO04710
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04720
      N21=XN                                                            TRO04730
      XN=-XJ21+XJ22+XJ23+EPS1                                           TRO04740
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04750
      N23=XN                                                            TRO04760
      XN=XJ12-XJ22+XJ32+EPS1                                            TRO04770
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04780
      N25=XN                                                            TRO04790
      XN=XJ21-XJ22+XJ23+EPS1                                            TRO04800
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04810
      N27=XN                                                            TRO04820
      XN=XJ12+XJ22-XJ32+EPS1                                            TRO04830
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04840
      N29=XN                                                            TRO04850
      XN=XJ21+XJ22-XJ23+EPS1                                            TRO04860
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04870
      N31=XN                                                            TRO04880
      XN=-XJ12-XJ21+XJ32+XJ23+EPS1                                      TRO04890
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04900
      N33=XN                                                            TRO04910
      XN=XJ12+XJ22+XJ32+EPS1                                            TRO04920
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04930
      N36=XN                                                            TRO04940
      N36=N36+1                                                         TRO04950
      XN=XJ21+XJ22+XJ23+EPS1                                            TRO04960
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO04970
      N38=XN                                                            TRO04980
      N38=N38+1                                                         TRO04990
      XN=-XJ13+XJ23+XJ33+EPS1                                           TRO05000
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05010
      N41=XN                                                            TRO05020
      XN=-XJ13+XJ11+XJ12+EPS1                                           TRO05030
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05040
      N44=XN                                                            TRO05050
      XN=XJ13-XJ23+XJ33+EPS1                                            TRO05060
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05070
      N45=XN                                                            TRO05080
      XN=XJ13-XJ11+XJ12+EPS1                                            TRO05090
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05100
      N48=XN                                                            TRO05110
      XN=XJ13+XJ23-XJ33+EPS1                                            TRO05120
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05130
      N49=XN                                                            TRO05140
      XN=XJ13+XJ11-XJ12+EPS1                                            TRO05150
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05160
      N52=XN                                                            TRO05170
      XN=-XJ23-XJ11+XJ33+XJ12+EPS1                                      TRO05180
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05190
      N54=XN                                                            TRO05200
      XN=XJ13+XJ23+XJ33+EPS1                                            TRO05210
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05220
      N56=XN                                                            TRO05230
      N56=N56+1                                                         TRO05240
      XN=XJ13+XJ11+XJ12+EPS1                                            TRO05250
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05260
      N59=XN                                                            TRO05270
      N59=N59+1                                                         TRO05280
C                                                                       TRO05290
C     TEST SUR LES XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJ33         TRO05300
C                                                                       TRO05310
      IF(N9.LT.0) GO TO 50                                              TRO05320
      IF(N5.LT.0) GO TO 50                                              TRO05330
      IF(N1.LT.0) GO TO 50                                              TRO05340
      IF(N10.LT.0) GO TO 50                                             TRO05350
      IF(N6.LT.0) GO TO 50                                              TRO05360
      IF(N2.LT.0) GO TO 50                                              TRO05370
      IF(N29.LT.0) GO TO 50                                             TRO05380
      IF(N25.LT.0) GO TO 50                                             TRO05390
      IF(N21.LT.0) GO TO 50                                             TRO05400
      IF(N31.LT.0) GO TO 50                                             TRO05410
      IF(N27.LT.0) GO TO 50                                             TRO05420
      IF(N23.LT.0) GO TO 50                                             TRO05430
      IF(N49.LT.0) GO TO 50                                             TRO05440
      IF(N45.LT.0) GO TO 50                                             TRO05450
      IF(N41.LT.0) GO TO 50                                             TRO05460
      IF(N52.LT.0) GO TO 50                                             TRO05470
      IF(N48.LT.0) GO TO 50                                             TRO05480
      IF(N44.LT.0) GO TO 50                                             TRO05490
      K=N1+N2+N5+N6+N9+N10+N21+N23+N25+N27+N29+N31+N41+N44+N45+N48+N49+NTRO05500
     152                                                                TRO05510
      IF(K.GT.0) GO TO 54                                               TRO05520
      C9J=1.D0                                                            TRO05530
      RETURN                                                            TRO05540
  50  C9J=0.D0                                                            TRO05550
      RETURN                                                            TRO05560
C                                                                       TRO05570
C     CALCUL DE LA SOMME SUR XJ                                         TRO05580
C                                                                       TRO05590
   54 XN=2.D0*(XJ21-XJ32)+EPS1                                            TRO05600
      IF(XN.LT.0.D0) XN=-(XN+EPS2)                                        TRO05610
      JMIN=XN                                                           TRO05620
      XN=2.D0*(XJ11-XJ33)+EPS1                                            TRO05630
      IF(XN.LT.0.D0) XN=-(XN+EPS2)                                        TRO05640
      N=XN                                                              TRO05650
      IF(N.GT.JMIN) JMIN=N                                              TRO05660
      XN=2.D0*(XJ12-XJ23)+EPS1                                            TRO05670
      IF(XN.LT.0.D0) XN=-(XN+EPS2)                                        TRO05680
      N=XN                                                              TRO05690
      IF(N.GT.JMIN) JMIN=N                                              TRO05700
      XN=2.D0*(XJ21+XJ32)+EPS1                                            TRO05710
      JMAX=XN                                                           TRO05720
      XN=2.D0*(XJ11+XJ33)+EPS1                                            TRO05730
      N=XN                                                              TRO05740
      IF(N.LT.JMAX) JMAX=N                                              TRO05750
      XN=2.D0*(XJ12+XJ23)+EPS1                                            TRO05760
      N=XN                                                              TRO05770
      IF(N.LT.JMAX) JMAX=N                                              TRO05780
      XJMIN=DFLOAT(JMIN)/2.D0                                              TRO05790
      XJMAX=DFLOAT(JMAX)/2.D0                                              TRO05800
      S=0.D0                                                              TRO05810
      XJ=XJMIN                                                          TRO05820
      XN=-XJ32+XJ21+XJ+EPS1                                             TRO05830
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05840
      N3=XN                                                             TRO05850
      XN=-XJ11+XJ33+XJ+EPS1                                             TRO05860
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05870
      N4=XN                                                             TRO05880
      XN=XJ32-XJ21+XJ+EPS1                                              TRO05890
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05900
      N7=XN                                                             TRO05910
      XN=XJ11-XJ33+XJ+EPS1                                              TRO05920
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05930
      N8=XN                                                             TRO05940
      XN=XJ32+XJ21-XJ+EPS1                                              TRO05950
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05960
      N11=XN                                                            TRO05970
      XN=XJ11+XJ33-XJ+EPS1                                              TRO05980
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO05990
      N12=XN                                                            TRO06000
      XN=-XJ11-XJ32+XJ31+XJ+EPS1                                        TRO06010
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06020
      N13=XN                                                            TRO06030
      XN=-XJ21-XJ33+XJ31+XJ+EPS1                                        TRO06040
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06050
      N14=XN                                                            TRO06060
      XN=XJ32+XJ21+XJ+EPS1                                              TRO06070
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06080
      N18=XN                                                            TRO06090
      N18=N18+1                                                         TRO06100
      XN=XJ11+XJ33+XJ+EPS1                                              TRO06110
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06120
      N19=XN                                                            TRO06130
      N19=N19+1                                                         TRO06140
      XN=-XJ21+XJ+XJ32+EPS1                                             TRO06150
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06160
      N22=XN                                                            TRO06170
      XN=-XJ12+XJ23+XJ+EPS1                                             TRO06180
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06190
      N24=XN                                                            TRO06200
      XN=XJ21-XJ+XJ32+EPS1                                              TRO06210
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06220
      N26=XN                                                            TRO06230
      XN=XJ12-XJ+XJ23+EPS1                                              TRO06240
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06250
      N28=XN                                                            TRO06260
      XN=XJ21+XJ-XJ32+EPS1                                              TRO06270
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06280
      N30=XN                                                            TRO06290
      XN=XJ12+XJ-XJ23+EPS1                                              TRO06300
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06310
      N32=XN                                                            TRO06320
      XN=-XJ22-XJ+XJ32+XJ23+EPS1                                        TRO06330
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06340
      N34=XN                                                            TRO06350
      XN=XJ12+XJ21+XJ22+XJ+EPS1                                         TRO06360
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06370
      N35=XN                                                            TRO06380
      N35=N35+1                                                         TRO06390
      XN=XJ21+XJ+XJ32+EPS1                                              TRO06400
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06410
      N37=XN                                                            TRO06420
      N37=N37+1                                                         TRO06430
      XN=XJ12+XJ+XJ23+EPS1                                              TRO06440
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06450
      N39=XN                                                            TRO06460
      N39=N39+1                                                         TRO06470
      XN=-XJ+XJ11+XJ33+EPS1                                             TRO06480
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06490
      N42=XN                                                            TRO06500
      XN=-XJ+XJ23+XJ12+EPS1                                             TRO06510
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06520
      N43=XN                                                            TRO06530
      XN=XJ-XJ11+XJ33+EPS1                                              TRO06540
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06550
      N46=XN                                                            TRO06560
      XN=XJ-XJ23+XJ12+EPS1                                              TRO06570
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06580
      N47=XN                                                            TRO06590
      XN=XJ+XJ11-XJ33+EPS1                                              TRO06600
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06610
      N50=XN                                                            TRO06620
      XN=XJ+XJ23-XJ12+EPS1                                              TRO06630
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06640
      N51=XN                                                            TRO06650
      XN=-XJ13-XJ+XJ33+XJ12+EPS1                                        TRO06660
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06670
      N53=XN                                                            TRO06680
      XN=XJ13+XJ+XJ23+XJ11+EPS1                                         TRO06690
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06700
      N55=XN                                                            TRO06710
      N55=N55+1                                                         TRO06720
      XN=XJ+XJ11+XJ33+EPS1                                              TRO06730
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06740
      N57=XN                                                            TRO06750
      N57=N57+1                                                         TRO06760
      XN=XJ+XJ23+XJ12+EPS1                                              TRO06770
      IF(XN.LT.0.D0) XN=XN+EPS2                                           TRO06780
      N58=XN                                                            TRO06790
      N58=N58+1                                                         TRO06800
      GO TO 10                                                          TRO06810
   52 IF(XJ.GT.XJMAX) GO TO 120                                         TRO06820
      N3=N3+1                                                           TRO06830
      N4=N4+1                                                           TRO06840
      N7=N7+1                                                           TRO06850
      N8=N8+1                                                           TRO06860
      N11=N11-1                                                         TRO06870
      N12=N12-1                                                         TRO06880
      N13=N13+1                                                         TRO06890
      N14=N14+1                                                         TRO06900
      N18=N18+1                                                         TRO06910
      N19=N19+1                                                         TRO06920
      N22=N22+1                                                         TRO06930
      N24=N24+1                                                         TRO06940
      N26=N26-1                                                         TRO06950
      N28=N28-1                                                         TRO06960
      N30=N30+1                                                         TRO06970
      N32=N32+1                                                         TRO06980
      N34=N34-1                                                         TRO06990
      N35=N35+1                                                         TRO07000
      N37=N37+1                                                         TRO07010
      N39=N39+1                                                         TRO07020
      N42=N42-1                                                         TRO07030
      N43=N43-1                                                         TRO07040
      N46=N46+1                                                         TRO07050
      N47=N47+1                                                         TRO07060
      N50=N50+1                                                         TRO07070
      N51=N51+1                                                         TRO07080
      N53=N53-1                                                         TRO07090
      N55=N55+1                                                         TRO07100
      N57=N57+1                                                         TRO07110
      N58=N58+1                                                         TRO07120
C                                                                       TRO07130
C     CALCUL DES SOMMES ALTERNEES S1(K1),S2(K2),S3(K3)                  TRO07140
C                                                                       TRO07150
   10 K1=0.D0                                                             TRO07160
      L1=-N13                                                           TRO07170
      IF(L1.GT.K1) K1=L1                                                TRO07180
      L1=-N14                                                           TRO07190
      IF(L1.GT.K1) K1=L1                                                TRO07200
      L1=N9                                                             TRO07210
      IF(N10.LT.L1) L1=N10                                              TRO07220
      IF(N11.LT.L1) L1=N11                                              TRO07230
      IF(N12.LT.L1) L1=N12                                              TRO07240
      IF(N15.LT.L1) L1=N15                                              TRO07250
      F1=1.D0                                                             TRO07260
      S1=1.D0                                                             TRO07270
      I1=K1+1                                                           TRO07280
  62  IF(I1.GT.L1) GO TO 64                                             TRO07290
      I1M1=I1-1                                                         TRO07300
      NN1=(N9-I1M1)*(N10-I1M1)*(N11-I1M1)*(N12-I1M1)                    TRO07310
      ND1=I1*(N13+I1)*(N14+I1)*(N15-I1M1)                               TRO07320
      F1=-F1*DFLOAT(NN1)/DFLOAT(ND1)                                      TRO07330
      S1=S1+F1                                                          TRO07340
      I1=I1+1                                                           TRO07350
      GO TO 62                                                          TRO07360
   64 K2=0                                                              TRO07370
      L2=-N33                                                           TRO07380
      IF(L2.GT.K2) K2=L2                                                TRO07390
      L2=-N34                                                           TRO07400
      IF(L2.GT.K2) K2=L2                                                TRO07410
      L2=N29                                                            TRO07420
      IF(N30.LT.L2) L2=N30                                              TRO07430
      IF(N31.LT.L2) L2=N31                                              TRO07440
      IF(N32.LT.L2) L2=N32                                              TRO07450
      IF(N35.LT.L2) L2=N35                                              TRO07460
      F2=1.D0                                                             TRO07470
      S2=1.D0                                                             TRO07480
      I2=K2+1                                                           TRO07490
  70  IF(I2.GT.L2) GO TO 80                                             TRO07500
      I2M2=I2-1                                                         TRO07510
      NN2=(N29-I2M2)*(N30-I2M2)*(N31-I2M2)*(N32-I2M2)                   TRO07520
      ND2=I2*(N33+I2)*(N34+I2)*(N35-I2M2)                               TRO07530
      F2=-F2*DFLOAT(NN2)/DFLOAT(ND2)                                      TRO07540
      S2=S2+F2                                                          TRO07550
      I2=I2+1                                                           TRO07560
      GO TO 70                                                          TRO07570
   80 K3=0                                                              TRO07580
      L3=-N53                                                           TRO07590
      IF(L3.GT.K3) K3=L3                                                TRO07600
      L3=-N54                                                           TRO07610
      IF(L3.GT.K3) K3=L3                                                TRO07620
      L3=N49                                                            TRO07630
      IF(N50.LT.L3) L3=N50                                              TRO07640
      IF(N51.LT.L3) L3=N51                                              TRO07650
      IF(N52.LT.L3) L3=N52                                              TRO07660
      IF(N55.LT.L3) L3=N55                                              TRO07670
      F3=1.D0                                                             TRO07680
      S3=1.D0                                                             TRO07690
      I3=K3+1                                                           TRO07700
   84 IF(I3.GT.L3) GO TO 90                                             TRO07710
      I3M3=I3-1                                                         TRO07720
      NN3=(N49-I3M3)*(N50-I3M3)*(N51-I3M3)*(N52-I3M3)                   TRO07730
      ND3=I3*(N53+I3)*(N54+I3)*(N55-I3M3)                               TRO07740
      F3=-F3*DFLOAT(NN3)/DFLOAT(ND3)                                      TRO07750
      S3=S3+F3                                                          TRO07760
      I3=I3+1                                                           TRO07770
      GO TO 84                                                          TRO07780
C                                                                       TRO07790
C     CALCUL DE LA RACINE D'UN TERME DE LA SOMME SUR J                  TRO07800
C                                                                       TRO07810
   90 S2N=FLOG(N3+1)+FLOG(N4+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N11+1)+FLOG(NTRO07820
     112+1)+FLOG(N22+1)+FLOG(N24+1)+FLOG(N26+1)+FLOG(N28+1)+FLOG(N30+1)+TRO07830
     2FLOG(N32+1)+FLOG(N42+1)+FLOG(N43+1)+FLOG(N46+1)+FLOG(N47+1)+FLOG(NTRO07840
     350+1)+FLOG(N51+1)                                                 TRO07850
      S2N=.5D0*S2N                                                        TRO07860
      S2D=FLOG(N18+1)+FLOG(N19+1)+FLOG(N37+1)+FLOG(N39+1)+FLOG(N57+1)+FLTRO07870
     1OG(N58+1)                                                         TRO07880
      S2D=.5D0*S2D                                                        TRO07890
      KM1=K1-1                                                          TRO07900
      KP1=K1+1                                                          TRO07910
      KM2=K2-1                                                          TRO07920
      KP2=K2+1                                                          TRO07930
      KM3=K3-1                                                          TRO07940
      KP3=K3+1                                                          TRO07950
      S2N=S2N+FLOG(N15-KM1)+FLOG(N35-KM2)+FLOG(N55-KM3)                 TRO07960
      S2D=S2D+FLOG(KP1)+FLOG(KP2)+FLOG(KP3)+FLOG(N9-KM1)+FLOG(N10-KM1)+FTRO07970
     1LOG(N11-KM1)+FLOG(N12-KM1)+FLOG(N13+KP1)+FLOG(N14+KP1)+FLOG(N29-KMTRO07980
     22)+FLOG(N30-KM2)+FLOG(N31-KM2)+FLOG(N32-KM2)+FLOG(N33+KP2)+FLOG(N3TRO07990
     34+KP2)+FLOG(N49-KM3)+FLOG(N50-KM3)+FLOG(N51-KM3)                  TRO08000
     4+FLOG(N52-KM3)+FLOG(N53+KP3)+FLOG(N54+KP3)                        TRO08010
      F=S2D-S2N                                                         TRO08020
      IF (F.GT.80.D0) GO TO 100                                            TRO08030
      F=S2N/S2D                                                         TRO08040
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 100                         TRO08050
      F=S1*S2*S3*DEXP(S2N-S2D)                                           TRO08060
      GO TO 110                                                         TRO08070
  100 F=S1*S2*S3                                                        TRO08080
      IF(F) 102,112,104                                                 TRO08090
  102 F=DLOG(-F)                                                        TRO08100
      F=-DEXP(F+S2N-S2D)                                                 TRO08110
      GO TO 110                                                         TRO08120
  104 F=DLOG(F)                                                         TRO08130
      F=DEXP(F+S2N-S2D)                                                  TRO08140
  110 F=F*(2.D0*XJ+1.D0)                                                    TRO08150
C                                                                       TRO08160
C     CALCUL DE LA PHASE D'UN TERME DE LA SOMME SUR J                   TRO08170
C                                                                       TRO08180
      L1=K1+K2+K3                                                       TRO08190
      K1=L1/2                                                           TRO08200
      K1=2*K1                                                           TRO08210
      IF(L1.NE.K1) F=-F                                                 TRO08220
      S=S+F                                                             TRO08230
  112 XJ=XJ+1.D0                                                          TRO08240
      GO TO 52                                                          TRO08250
C                                                                       TRO08260
C     CALCUL DU C9J SANS PHASE                                          TRO08270
C                                                                       TRO08280
  120 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N5+1)+FLOG(N6+1)+FLOG(N9+1)+FLOG(N1TRO08290
     10+1)+FLOG(N21+1)+FLOG(N23+1)+FLOG(N25+1)+FLOG(N27+1)+FLOG(N29+1)+FTRO08300
     2LOG(N31+1)+FLOG(N41+1)+FLOG(N44+1)+FLOG(N45+1)+FLOG(N48+1)+FLOG(N4TRO08310
     39+1)+FLOG(N52+1)                                                  TRO08320
      C2N=.5D0*C2N                                                        TRO08330
      C2D=FLOG(N16+1)+FLOG(N17+1)+FLOG(N36+1)+FLOG(N38+1)+FLOG(N56+1)+FLTRO08340
     1OG(N59+1)                                                         TRO08350
      C2D=.5D0*C2D                                                        TRO08360
      F=C2D-C2N                                                         TRO08370
      IF(F.GT.80.D0) GO TO 122                                            TRO08380
      F=C2N/C2D                                                         TRO08390
      IF((F.LT.1.01D0).AND.(F.GT.0.98D0)) GO TO 122                         TRO08400
      C9J=S*DEXP(C2N-C2D)                                                TRO08410
      GO TO 130                                                         TRO08420
  122 IF(S) 124,50,126                                                  TRO08430
  124 S=DLOG(-S)                                                        TRO08440
      C9J=-DEXP(S+C2N-C2D)                                               TRO08450
      GO TO 130                                                         TRO08460
  126 S=DLOG(S)                                                         TRO08470
      C9J=DEXP(S+C2N-C2D)                                                TRO08480
C                                                                       TRO08490
C     CALCUL DE LA PHASE                                                TRO08500
C                                                                       TRO08510
  130 K=N9+N16+N36+N56-1                                                TRO08520
      L=K/2                                                             TRO08530
      L=2*L                                                             TRO08540
      IF (L.NE.K) C9J=-C9J                                              TRO08550
      RETURN                                                            TRO08560
      END                                                               TRO08570

