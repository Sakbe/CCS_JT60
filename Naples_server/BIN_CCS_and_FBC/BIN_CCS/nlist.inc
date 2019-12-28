&Pccs
 rj_dist=-0.20, zj_dist=0.00, ella=0.40, ellb=1.20, delta=0.00,
 iauto=1, startelipt=0.5
/
&Pdevice
 device = 'JT-60SA',
 icoilc = 0, ifluxl = 0,
 rrange  =  0.0, 6.0, zrange  = -4.5, 4.5,
 rvrange =  1.5, 4.5, zvrange = -3.0, 3.0
/
&FluxLoop
! 34 flux loops expected in JT-60SA
  nflux =   1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
           11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
           21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
           31, 32, 33, 34
/
&Probe
! 45 magnetic probes expected in JT-60SA
  nprob =   1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
           11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
           21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
           31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
           41, 42, 43, 44, 45
/
&Based
 ibase = 1
/
&Snoise
! Noise level of flux loop [Wb] (0.0001666 Wb expected in JT-60SA)
 fnoise = 0.00
/
&Snoise2
! Noise level of magnetic probe [G] (0.3333 G expected in JT-60SA)
 pnoise = 0.0
/
&Snoise3
! Noise level of PF coil current [%] (0.15 % expected in JT-60SA)
 cnoise = 0.0
/
&Times
! If you use MECS, you don't need to change here
 iwave = 1, starttime = 20.6, endtime = 24.6, interval = 200, ieddy=2
/
&Delay
 delaymeas = 0.0
/
