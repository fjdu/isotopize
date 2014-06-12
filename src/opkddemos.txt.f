
==========Source and Output for LSODE Demonstration Program=====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODE package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The package is used to solve two simple problems,
c one with a full Jacobian, the other with a banded Jacobian,
c with all 8 of the appropriate values of mf in each case.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.  All output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1, f2, jac2
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, er, erm, ero, hu, rtol, rwork, t,
     1   tout, tout1, y
      dimension y(25), rwork(697), iwork(45)
      data lout/6/, tout1/1.39283880203d0/, dtout/2.214773875d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-6
      lrw = 697
      liw = 45
      iopt = 0
c
c First problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(/' Demonstration program for DLSODE package'///
     1  ' Problem 1:  Van der Pol oscillator:'/
     2  '  xdotdot - 3*(1 - x**2)*xdot + x = 0, ',
     3  '   x(0) = 2, xdot(0) = 0'/
     4  ' neq =',i2/
     5  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
c
      do 195 meth = 1,2
      do 190 miter = 0,3
      mf = 10*meth + miter
      write (lout,120) mf
 120  format(///' Solution with mf =',i3//
     1     5x,'t               x               xdot       nq      h'//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      tout = tout1
      ero = 0.0d0
      do 170 iout = 1,nout
        call dlsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac1,mf)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,140) t,y(1),y(2),nqu,hu
 140    format(d15.5,d16.5,d14.3,i5,d14.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 170
        er = abs(y(1))/atol
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,150)
 150      format(//' Warning: error exceeds 1000 * tolerance'//)
          nerr = nerr + 1
        endif
 170    tout = tout + dtout
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - neq*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' error overrun =',d10.2)
 190  continue
 195  continue
c
c Second problem
c
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      nout = 5
      write (lout,210) neq,ml,mu,itol,rtol,atol
 210  format(///70('-')///
     1  ' Problem 2: ydot = A * y , where',
     2  '  A is a banded lower triangular matrix'/
     3     12x, 'derived from 2-D advection PDE'/
     4  ' neq =',i3,'   ml =',i2,'   mu =',i2/
     5  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
      do 295 meth = 1,2
      do 290 miter = 0,5
      if (miter .eq. 1 .or. miter .eq. 2) go to 290
      mf = 10*meth + miter
      write (lout,220) mf
 220  format(///' Solution with mf =',i3//
     1       5x,'t             max.err.     nq      h'//)
      t = 0.0d0
      do 230 i = 2,neq
 230    y(i) = 0.0d0
      y(1) = 1.0d0
      itask = 1
      istate = 1
      tout = 0.01d0
      ero = 0.0d0
      do 270 iout = 1,nout
        call dlsode(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,mf)
        call edit2(y,t,erm)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,240) t,erm,nqu,hu
 240    format(d15.5,d14.3,i5,d14.3)
        if (istate .lt. 0) go to 275
        er = erm/atol
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,150)
          nerr = nerr + 1
        endif
 270    tout = tout*10.0d0
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 5) nfea = nfe - mband*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 290  continue
 295  continue
      write (lout,300) nerr
 300  format(////' Number of errors encountered =',i3)
      stop
      end

      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      ydot(1) = y(2)
      ydot(2) = 3.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end

      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(neq), pd(nrowpd,neq)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
      return
      end

      subroutine f2 (neq, t, y, ydot)
      integer neq, i, j, k, ng
      double precision t, y, ydot, alph1, alph2, d
      dimension y(neq), ydot(neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end

      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(neq), pd(nrowpd,neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end

      subroutine edit2 (y, t, erm)
      integer i, j, k, ng
      double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
      dimension y(25)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      erm = 0.0d0
      if (t .eq. 0.0d0) return
      ex = 0.0d0
      if (t .le. 30.0d0) ex = exp(-2.0d0*t)
      a2 = 1.0d0
      do 60 j = 1,ng
        a1 = 1.0d0
        do 50 i = 1,ng
          k = i + (j - 1)*ng
          yt = t**(i+j-2)*ex*a1*a2
          er = abs(y(k)-yt)
          erm = max(erm,er)
          a1 = a1*alph1/i
 50       continue
        a2 = a2*alph2/j
 60     continue
      return
      end

................................................................................


 Demonstration program for DLSODE package


 Problem 1:  Van der Pol oscillator:
  xdotdot - 3*(1 - x**2)*xdot + x = 0,    x(0) = 2, xdot(0) = 0
 neq = 2
 itol =  1   rtol =   0.0E+00   atol =   0.1E-05





 Solution with mf = 10

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    3     0.123E+00
    0.36076E+01    -0.77986E-04    -0.317E+01    5     0.217E-01
    0.58224E+01    -0.16801E+01     0.291E+00    3     0.475E-01
    0.80372E+01     0.11669E-03     0.317E+01    5     0.234E-01


 Final statistics for this run:
 rwork size =  52   iwork size =  20
 number of steps =  297
 number of f-s   =  352
 (excluding J-s) =  352
 number of J-s   =    0
 error overrun =  0.12E+03



 Solution with mf = 11

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.121E+00
    0.36076E+01    -0.17732E-04    -0.317E+01    5     0.187E-01
    0.58224E+01    -0.16801E+01     0.291E+00    6     0.963E-01
    0.80372E+01     0.25894E-04     0.317E+01    5     0.190E-01


 Final statistics for this run:
 rwork size =  58   iwork size =  22
 number of steps =  203
 number of f-s   =  281
 (excluding J-s) =  281
 number of J-s   =   29
 error overrun =  0.26E+02



 Solution with mf = 12

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.121E+00
    0.36076E+01    -0.17732E-04    -0.317E+01    5     0.187E-01
    0.58224E+01    -0.16801E+01     0.291E+00    6     0.963E-01
    0.80372E+01     0.25894E-04     0.317E+01    5     0.190E-01


 Final statistics for this run:
 rwork size =  58   iwork size =  22
 number of steps =  203
 number of f-s   =  339
 (excluding J-s) =  281
 number of J-s   =   29
 error overrun =  0.26E+02



 Solution with mf = 13

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.739E-01
    0.36076E+01     0.34401E-04    -0.317E+01    6     0.260E-01
    0.58224E+01    -0.16801E+01     0.291E+00    4     0.133E+00
    0.80372E+01    -0.59053E-04     0.317E+01    5     0.205E-01


 Final statistics for this run:
 rwork size =  56   iwork size =  20
 number of steps =  198
 number of f-s   =  315
 (excluding J-s) =  289
 number of J-s   =   26
 error overrun =  0.59E+02



 Solution with mf = 20

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.549E-01
    0.36076E+01    -0.56579E-04    -0.317E+01    5     0.143E-01
    0.58224E+01    -0.16801E+01     0.291E+00    4     0.583E-01
    0.80372E+01     0.10387E-03     0.317E+01    5     0.149E-01


 Final statistics for this run:
 rwork size =  38   iwork size =  20
 number of steps =  289
 number of f-s   =  321
 (excluding J-s) =  321
 number of J-s   =    0
 error overrun =  0.10E+03



 Solution with mf = 21

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.676E-01
    0.36076E+01    -0.48977E-04    -0.317E+01    5     0.141E-01
    0.58224E+01    -0.16801E+01     0.291E+00    5     0.126E+00
    0.80372E+01     0.96867E-04     0.317E+01    5     0.142E-01


 Final statistics for this run:
 rwork size =  44   iwork size =  22
 number of steps =  262
 number of f-s   =  345
 (excluding J-s) =  345
 number of J-s   =   30
 error overrun =  0.97E+02



 Solution with mf = 22

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.676E-01
    0.36076E+01    -0.48977E-04    -0.317E+01    5     0.141E-01
    0.58224E+01    -0.16801E+01     0.291E+00    5     0.126E+00
    0.80372E+01     0.96867E-04     0.317E+01    5     0.142E-01


 Final statistics for this run:
 rwork size =  44   iwork size =  22
 number of steps =  262
 number of f-s   =  405
 (excluding J-s) =  345
 number of J-s   =   30
 error overrun =  0.97E+02



 Solution with mf = 23

     t               x               xdot       nq      h


    0.13928E+01     0.16801E+01    -0.291E+00    5     0.709E-01
    0.36076E+01    -0.46705E-04    -0.317E+01    5     0.139E-01
    0.58224E+01    -0.16801E+01     0.291E+00    3     0.719E-01
    0.80372E+01     0.54700E-04     0.317E+01    5     0.154E-01


 Final statistics for this run:
 rwork size =  42   iwork size =  20
 number of steps =  271
 number of f-s   =  414
 (excluding J-s) =  383
 number of J-s   =   31
 error overrun =  0.55E+02



----------------------------------------------------------------------


 Problem 2: ydot = A * y , where  A is a banded lower triangular matrix
            derived from 2-D advection PDE
 neq = 25   ml = 5   mu = 0
 itol =  1   rtol =   0.0E+00   atol =   0.1E-05





 Solution with mf = 10

     t             max.err.     nq      h


    0.10000E-01     0.556E-06    2     0.766E-02
    0.10000E+00     0.655E-05    3     0.249E-01
    0.10000E+01     0.274E-05    4     0.520E-01
    0.10000E+02     0.114E-05    3     0.117E+00
    0.10000E+03     0.221E-05    2     0.262E+00


 Final statistics for this run:
 rwork size = 420   iwork size =  20
 number of steps =  524
 number of f-s   =  552
 (excluding J-s) =  552
 number of J-s   =    0
 error overrun =  0.65E+01



 Solution with mf = 13

     t             max.err.     nq      h


    0.10000E-01     0.839E-06    2     0.949E-02
    0.10000E+00     0.208E-05    3     0.250E-01
    0.10000E+01     0.127E-03    3     0.168E-01
    0.10000E+02     0.113E-04    3     0.385E+00
    0.10000E+03     0.145E-05    2     0.149E+02


 Final statistics for this run:
 rwork size = 447   iwork size =  20
 number of steps =  129
 number of f-s   =  235
 (excluding J-s) =  201
 number of J-s   =   34
 error overrun =  0.13E+03



 Solution with mf = 14

     t             max.err.     nq      h


    0.10000E-01     0.877E-06    2     0.965E-02
    0.10000E+00     0.206E-05    3     0.250E-01
    0.10000E+01     0.126E-05    5     0.935E-01
    0.10000E+02     0.311E-06    6     0.442E+00
    0.10000E+03     0.159E-07    2     0.291E+02


 Final statistics for this run:
 rwork size = 697   iwork size =  45
 number of steps =   92
 number of f-s   =  113
 (excluding J-s) =  113
 number of J-s   =   18
 error overrun =  0.21E+01



 Solution with mf = 15

     t             max.err.     nq      h


    0.10000E-01     0.877E-06    2     0.965E-02
    0.10000E+00     0.206E-05    3     0.250E-01
    0.10000E+01     0.126E-05    5     0.935E-01
    0.10000E+02     0.311E-06    6     0.442E+00
    0.10000E+03     0.160E-07    2     0.291E+02


 Final statistics for this run:
 rwork size = 697   iwork size =  45
 number of steps =   92
 number of f-s   =  221
 (excluding J-s) =  113
 number of J-s   =   18
 error overrun =  0.21E+01



 Solution with mf = 20

     t             max.err.     nq      h


    0.10000E-01     0.465E-06    2     0.483E-02
    0.10000E+00     0.131E-05    3     0.148E-01
    0.10000E+01     0.427E-05    5     0.635E-01
    0.10000E+02     0.192E-05    4     0.351E+00
    0.10000E+03     0.929E-07    1     0.455E+00


 Final statistics for this run:
 rwork size = 245   iwork size =  20
 number of steps =  330
 number of f-s   =  530
 (excluding J-s) =  530
 number of J-s   =    0
 error overrun =  0.43E+01



 Solution with mf = 23

     t             max.err.     nq      h


    0.10000E-01     0.101E-05    2     0.598E-02
    0.10000E+00     0.446E-06    3     0.146E-01
    0.10000E+01     0.153E-05    5     0.738E-01
    0.10000E+02     0.578E-06    4     0.324E+00
    0.10000E+03     0.908E-08    1     0.992E+02


 Final statistics for this run:
 rwork size = 272   iwork size =  20
 number of steps =  180
 number of f-s   =  325
 (excluding J-s) =  274
 number of J-s   =   51
 error overrun =  0.15E+01



 Solution with mf = 24

     t             max.err.     nq      h


    0.10000E-01     0.104E-05    2     0.608E-02
    0.10000E+00     0.463E-06    3     0.146E-01
    0.10000E+01     0.247E-05    5     0.666E-01
    0.10000E+02     0.828E-06    5     0.391E+00
    0.10000E+03     0.384E-09    1     0.108E+03


 Final statistics for this run:
 rwork size = 522   iwork size =  45
 number of steps =  118
 number of f-s   =  136
 (excluding J-s) =  136
 number of J-s   =   18
 error overrun =  0.25E+01



 Solution with mf = 25

     t             max.err.     nq      h


    0.10000E-01     0.104E-05    2     0.608E-02
    0.10000E+00     0.463E-06    3     0.146E-01
    0.10000E+01     0.247E-05    5     0.666E-01
    0.10000E+02     0.828E-06    5     0.391E+00
    0.10000E+03     0.384E-09    1     0.108E+03


 Final statistics for this run:
 rwork size = 522   iwork size =  45
 number of steps =  118
 number of f-s   =  244
 (excluding J-s) =  136
 number of J-s   =   18
 error overrun =  0.25E+01




 Number of errors encountered =  0

==========Source and Output for LSODES Demonstration Program====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODES package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The package is used for each of the relevant values of mf to solve
c the problem ydot = A * y, where A is the 9 by 9 sparse matrix
c
c               -4  1     1
c                1 -4  1     1
c                   1 -4        1
c                        -4  1     1
c       A =               1 -4  1     1
c                            1 -4        1
c                                 -4  1
c                                  1 -4  1
c                                     1 -4
c
c The initial conditions are  y(0) = (1, 2, 3, ..., 9).
c Output is printed at t = 1, 2, and 3.
c Each case is solved first with nominal (large) values of lrw and liw,
c and then with values given by lenrw and leniw (optional outputs)
c on the first run, as a check on these computed work array lengths.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.
c All output is on unit lout, which is data-loaded to 6 below.
c-----------------------------------------------------------------------
      external fdem, jdem
      integer i, ia, igrid, iopt, iout, irun, istate, itask, itol,
     1  iwork, j, ja, k, l, leniw, lenrw, liw, lout, lrw,
     2  m, meth, mf, miter, moss, neq, nerr, nfe, nfea,
     3  ngp, nje, nlu, nnz, nout, nqu, nst, nzl, nzu
      double precision atol, erm, ero, hu, rtol, rwork, t, tout, y
      dimension y(9), ia(10), ja(50), iwork(90), rwork(1000)
      equivalence (ia(1),iwork(31)), (ja(1),iwork(41))
      data lout/6/
c
c Write heading and set fixed parameters.
      write(lout,10)
 10   format(/'Demonstration problem for the DLSODES package'//)
      nerr = 0
      igrid = 3
      neq = igrid**2
      t = 0.0d0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-5
      itask = 1
      iopt = 0
      do 20 i = 1,neq
 20     y(i) = i
      ia(1) = 1
      k = 1
      do 60 m = 1,igrid
        do 50 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .gt. 1) then
            ja(k) = j - igrid
            k = k + 1
          endif
 30       if (l .gt. 1) then
            ja(k) = j - 1
            k = k + 1
          endif
 35       ja(k) = j
          k = k + 1
          if (l .lt. igrid) then
            ja(k) = j + 1
            k = k + 1
          endif
 40       ia(j+1) = k
 50       continue
 60     continue
      write (lout,80)neq,t,rtol,atol,(y(i),i=1,neq)
 80   format(' neq =',i4,5x,'t0 =',f4.1,5x,'rtol =',d12.3,5x,
     1   'atol =',d12.3//' Initial y vector =  ',9f5.1)
c
c Loop over all relevant values of mf.
      do 193 moss = 0,2
      do 192 meth = 1,2
      do 191 miter = 0,3
      if ( (miter.eq.0 .or. miter.eq.3) .and. moss.ne.0) go to 191
      mf = 100*moss + 10*meth + miter
      write (lout,100)
 100  format(//80('*'))
c First run: nominal work array lengths, 3 output points.
      irun = 1
      lrw = 1000
      liw = 90
      nout = 3
 110  continue
      write (lout,120)mf,lrw,liw
 120  format(//'Run with mf =',i4,'.',5x,
     1       'Input work lengths lrw, liw =',2i6/)
      do 125 i = 1,neq
 125    y(i) = i
      t = 0.0d0
      tout = 1.0d0
      istate = 1
      ero = 0.0d0
c Loop over output points.  Do output and accuracy check at each.
      do 170 iout = 1,nout
        call dlsodes (fdem, neq, y, t, tout, itol, rtol, atol, itask,
     1                istate, iopt, rwork, lrw, iwork, liw, jdem, mf)
        nst = iwork(11)
        hu = rwork(11)
        nqu = iwork(14)
        call edit (y, iout, erm)
        write(lout,140)t,nst,hu,nqu,erm,(y(i),i=1,neq)
 140    format('At t =',f5.1,3x,'nst =',i4,3x,'hu =',d12.3,3x,
     1    'nqu =',i3,3x,' max. err. =',d11.3/
     2    '  y array =    ',4d15.6/5d15.6)
        if (istate .lt. 0) go to 175
        erm = erm/atol
        ero = max(ero,erm)
        if (erm .gt. 100.0d0) then
          write (lout,150)
 150      format(//' Warning: error exceeds 100 * tolerance'//)
          nerr = nerr + 1
        endif
        tout = tout + 1.0d0
 170    continue
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      if (irun .eq. 2) go to 191
c Print final statistics (first run only)
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nnz = iwork(19)
      ngp = iwork(20)
      nlu = iwork(21)
      nzl = iwork(25)
      nzu = iwork(26)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - ngp*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(/'Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' error overrun =',d10.2)
      if (miter .eq. 1 .or. miter .eq. 2)
     1   write (lout,185)nnz,ngp,nlu,nzl,nzu
 185  format(' number of nonzeros in J = ',i5/
     1  ' number of J index groups =',i5/
     2  ' number of LU decomp-s    =',i5/
     3  ' nonzeros in strict lower factor =',i5/
     4  ' nonzeros in strict upper factor =',i5)
      if (istate .lt. 0) go to 191
      if (miter .eq. 1 .or. miter .eq. 2)
     1   call ssout (neq, rwork(21), iwork, lout)
c Return for second run: minimal work array lengths, 1 output point.
      irun = irun + 1
      lrw = lenrw
      liw = leniw
      nout = 1
      go to 110
 191  continue
 192  continue
 193  continue
c
      write (lout,100)
      write (lout,200) nerr
 200  format(//'Number of errors encountered =',i3)
      stop
      end

      subroutine fdem (neq, t, y, ydot)
      integer neq,  i, igrid, j, l, m
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      data igrid/3/
      do 5 i = 1,neq
 5      ydot(i) = 0.0d0
      do 20 m = 1,igrid
        do 10 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .ne. 1) ydot(j-igrid) = ydot(j-igrid) + y(j)
          if (l .ne. 1) ydot(j-1) = ydot(j-1) + y(j)
          ydot(j) = ydot(j) - 4.0d0*y(j)
          if (l .ne. igrid) ydot(j+1) = ydot(j+1) + y(j)
 10       continue
 20     continue
      return
      end

      subroutine jdem (neq, t, y, j, ia, ja, pdj)
      integer neq, j, ia, ja,  igrid, l, m
      double precision t, y, pdj
      dimension y(neq), ia(*), ja(*), pdj(neq)
      data igrid/3/
      m = (j - 1)/igrid + 1
      l = j - (m - 1)*igrid
      pdj(j) = -4.0d0
      if (m .ne. 1) pdj(j-igrid) = 1.0d0
      if (l .ne. 1) pdj(j-1) = 1.0d0
      if (l .ne. igrid) pdj(j+1) = 1.0d0
      return
      end

      subroutine edit (y, iout, erm)
      integer iout,  i, neq
      double precision y, erm,   er, yex
      dimension y(*),yex(9,3)
      data neq /9/
      data yex /6.687279d-01, 9.901910d-01, 7.603061d-01,
     1   8.077979d-01, 1.170226e+00, 8.810605d-01, 5.013331d-01,
     2   7.201389d-01, 5.379644d-01, 1.340488d-01, 1.917157d-01,
     3   1.374034d-01, 1.007882d-01, 1.437868d-01, 1.028010d-01,
     4   3.844343d-02, 5.477593d-02, 3.911435d-02, 1.929166d-02,
     5   2.735444d-02, 1.939611d-02, 1.055981d-02, 1.496753d-02,
     6   1.060897d-02, 2.913689d-03, 4.128975d-03, 2.925977d-03/
      erm = 0.0d0
      do 10 i = 1,neq
        er = abs(y(i) - yex(i,iout))
 10     erm = max(erm,er)
      return
      end

      subroutine ssout (neq, iwk, iwork, lout)
      integer neq, iwk, iwork, lout
      integer i, i1, i2, ipian, ipjan, nnz
      dimension iwk(*), iwork(*)
      ipian = iwork(23)
      ipjan = iwork(24)
      nnz = iwork(19)
      i1 = ipian
      i2 = i1 + neq
      write (lout,10)(iwk(i),i=i1,i2)
 10   format(/' structure descriptor array ian ='/(20i4))
      i1 = ipjan
      i2 = i1 + nnz - 1
      write (lout,20)(iwk(i),i=i1,i2)
 20   format(/' structure descriptor array jan ='/(20i4))
      return
      end

................................................................................


Demonstration problem for the DLSODES package


 neq =   9     t0 = 0.0     rtol =   0.000E+00     atol =   0.100E-04

 Initial y vector =    1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0


********************************************************************************


Run with mf =  10.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  42   hu =   0.241E-01   nqu =  4    max. err. =  0.394E-05
  y array =       0.668730E+00   0.990188E+00   0.760308E+00   0.807799E+00
   0.117023E+01   0.881061E+00   0.501332E+00   0.720142E+00   0.537962E+00
At t =  2.0   nst =  71   hu =   0.680E-01   nqu =  3    max. err. =  0.273E-04
  y array =       0.134047E+00   0.191717E+00   0.137407E+00   0.100802E+00
   0.143808E+00   0.102820E+00   0.384623E-01   0.548033E-01   0.391361E-01
At t =  3.0   nst =  90   hu =   0.455E-01   nqu =  3    max. err. =  0.121E-04
  y array =       0.193008E-01   0.273568E-01   0.194059E-01   0.105663E-01
   0.149796E-01   0.106158E-01   0.291803E-02   0.413489E-02   0.293048E-02

Final statistics for this run:
 rwork size = 164   iwork size =  30
 number of steps =   90
 number of f-s   =   98
 (excluding J-s) =   98
 number of J-s   =    0
 error overrun =  0.27E+01


Run with mf =  10.     Input work lengths lrw, liw =   164    30

At t =  1.0   nst =  42   hu =   0.241E-01   nqu =  4    max. err. =  0.394E-05
  y array =       0.668730E+00   0.990188E+00   0.760308E+00   0.807799E+00
   0.117023E+01   0.881061E+00   0.501332E+00   0.720142E+00   0.537962E+00


********************************************************************************


Run with mf =  11.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00
At t =  2.0   nst =  37   hu =   0.115E+00   nqu =  5    max. err. =  0.593E-05
  y array =       0.134047E+00   0.191713E+00   0.137403E+00   0.100789E+00
   0.143787E+00   0.102803E+00   0.384477E-01   0.547819E-01   0.391199E-01
At t =  3.0   nst =  44   hu =   0.155E+00   nqu =  5    max. err. =  0.386E-05
  y array =       0.192920E-01   0.273552E-01   0.193970E-01   0.105624E-01
   0.149714E-01   0.106120E-01   0.291609E-02   0.413247E-02   0.292857E-02

Final statistics for this run:
 rwork size = 308   iwork size =  67
 number of steps =   44
 number of f-s   =   56
 (excluding J-s) =   56
 number of J-s   =    1
 error overrun =  0.99E+00
 number of nonzeros in J =    27
 number of J index groups =    0
 number of LU decomp-s    =    9
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf =  11.     Input work lengths lrw, liw =   308    67

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00


********************************************************************************


Run with mf =  12.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00
At t =  2.0   nst =  37   hu =   0.115E+00   nqu =  5    max. err. =  0.593E-05
  y array =       0.134047E+00   0.191713E+00   0.137403E+00   0.100789E+00
   0.143787E+00   0.102803E+00   0.384477E-01   0.547819E-01   0.391199E-01
At t =  3.0   nst =  44   hu =   0.155E+00   nqu =  5    max. err. =  0.386E-05
  y array =       0.192920E-01   0.273552E-01   0.193970E-01   0.105624E-01
   0.149714E-01   0.106120E-01   0.291609E-02   0.413247E-02   0.292857E-02

Final statistics for this run:
 rwork size = 315   iwork size =  67
 number of steps =   44
 number of f-s   =   60
 (excluding J-s) =   56
 number of J-s   =    1
 error overrun =  0.99E+00
 number of nonzeros in J =    27
 number of J index groups =    4
 number of LU decomp-s    =    9
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf =  12.     Input work lengths lrw, liw =   315    67

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00


********************************************************************************


Run with mf =  13.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  28   hu =   0.825E-01   nqu =  5    max. err. =  0.591E-05
  y array =       0.668729E+00   0.990192E+00   0.760304E+00   0.807797E+00
   0.117022E+01   0.881060E+00   0.501334E+00   0.720136E+00   0.537970E+00
At t =  2.0   nst =  41   hu =   0.825E-01   nqu =  5    max. err. =  0.468E-05
  y array =       0.134048E+00   0.191713E+00   0.137405E+00   0.100784E+00
   0.143782E+00   0.102798E+00   0.384467E-01   0.547752E-01   0.391129E-01
At t =  3.0   nst =  58   hu =   0.536E-01   nqu =  4    max. err. =  0.801E-04
  y array =       0.193021E-01   0.273573E-01   0.194762E-01   0.105617E-01
   0.149669E-01   0.106567E-01   0.291302E-02   0.412818E-02   0.292528E-02

Final statistics for this run:
 rwork size = 175   iwork size =  30
 number of steps =   58
 number of f-s   =   90
 (excluding J-s) =   80
 number of J-s   =   10
 error overrun =  0.80E+01


Run with mf =  13.     Input work lengths lrw, liw =   175    30

At t =  1.0   nst =  28   hu =   0.825E-01   nqu =  5    max. err. =  0.591E-05
  y array =       0.668729E+00   0.990192E+00   0.760304E+00   0.807797E+00
   0.117022E+01   0.881060E+00   0.501334E+00   0.720136E+00   0.537970E+00


********************************************************************************


Run with mf =  20.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  39   hu =   0.549E-01   nqu =  5    max. err. =  0.378E-04
  y array =       0.668726E+00   0.990193E+00   0.760309E+00   0.807791E+00
   0.117020E+01   0.881056E+00   0.501356E+00   0.720124E+00   0.538002E+00
At t =  2.0   nst =  53   hu =   0.677E-01   nqu =  5    max. err. =  0.113E-04
  y array =       0.134039E+00   0.191719E+00   0.137397E+00   0.100792E+00
   0.143779E+00   0.102808E+00   0.384485E-01   0.547872E-01   0.391222E-01
At t =  3.0   nst =  64   hu =   0.123E+00   nqu =  5    max. err. =  0.869E-05
  y array =       0.192944E-01   0.273518E-01   0.193999E-01   0.105634E-01
   0.149762E-01   0.106132E-01   0.291807E-02   0.413507E-02   0.293054E-02

Final statistics for this run:
 rwork size = 101   iwork size =  30
 number of steps =   64
 number of f-s   =   77
 (excluding J-s) =   77
 number of J-s   =    0
 error overrun =  0.38E+01


Run with mf =  20.     Input work lengths lrw, liw =   101    30

At t =  1.0   nst =  39   hu =   0.549E-01   nqu =  5    max. err. =  0.378E-04
  y array =       0.668726E+00   0.990193E+00   0.760309E+00   0.807791E+00
   0.117020E+01   0.881056E+00   0.501356E+00   0.720124E+00   0.538002E+00


********************************************************************************


Run with mf =  21.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00
At t =  2.0   nst =  52   hu =   0.105E+00   nqu =  5    max. err. =  0.132E-04
  y array =       0.134044E+00   0.191709E+00   0.137402E+00   0.100788E+00
   0.143786E+00   0.102805E+00   0.384531E-01   0.547890E-01   0.391276E-01
At t =  3.0   nst =  61   hu =   0.132E+00   nqu =  5    max. err. =  0.134E-04
  y array =       0.192907E-01   0.273543E-01   0.193977E-01   0.105672E-01
   0.149788E-01   0.106186E-01   0.292280E-02   0.414233E-02   0.293619E-02

Final statistics for this run:
 rwork size = 245   iwork size =  67
 number of steps =   61
 number of f-s   =   71
 (excluding J-s) =   71
 number of J-s   =    2
 error overrun =  0.25E+01
 number of nonzeros in J =    27
 number of J index groups =    0
 number of LU decomp-s    =    8
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf =  21.     Input work lengths lrw, liw =   245    67

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00


********************************************************************************


Run with mf =  22.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00
At t =  2.0   nst =  52   hu =   0.105E+00   nqu =  5    max. err. =  0.132E-04
  y array =       0.134044E+00   0.191709E+00   0.137402E+00   0.100788E+00
   0.143786E+00   0.102805E+00   0.384531E-01   0.547890E-01   0.391276E-01
At t =  3.0   nst =  61   hu =   0.132E+00   nqu =  5    max. err. =  0.134E-04
  y array =       0.192907E-01   0.273543E-01   0.193977E-01   0.105672E-01
   0.149788E-01   0.106186E-01   0.292280E-02   0.414233E-02   0.293619E-02

Final statistics for this run:
 rwork size = 252   iwork size =  67
 number of steps =   61
 number of f-s   =   79
 (excluding J-s) =   71
 number of J-s   =    2
 error overrun =  0.25E+01
 number of nonzeros in J =    27
 number of J index groups =    4
 number of LU decomp-s    =    8
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf =  22.     Input work lengths lrw, liw =   252    67

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00


********************************************************************************


Run with mf =  23.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  39   hu =   0.564E-01   nqu =  5    max. err. =  0.334E-04
  y array =       0.668727E+00   0.990191E+00   0.760307E+00   0.807793E+00
   0.117019E+01   0.881055E+00   0.501339E+00   0.720129E+00   0.537977E+00
At t =  2.0   nst =  53   hu =   0.720E-01   nqu =  5    max. err. =  0.268E-04
  y array =       0.134056E+00   0.191713E+00   0.137402E+00   0.100785E+00
   0.143781E+00   0.102798E+00   0.384437E-01   0.548027E-01   0.391161E-01
At t =  3.0   nst =  80   hu =   0.200E-01   nqu =  4    max. err. =  0.490E-04
  y array =       0.192918E-01   0.273431E-01   0.194010E-01   0.105628E-01
   0.149652E-01   0.106029E-01   0.290678E-02   0.412271E-02   0.297494E-02

Final statistics for this run:
 rwork size = 112   iwork size =  30
 number of steps =   80
 number of f-s   =  137
 (excluding J-s) =  122
 number of J-s   =   15
 error overrun =  0.49E+01


Run with mf =  23.     Input work lengths lrw, liw =   112    30

At t =  1.0   nst =  39   hu =   0.564E-01   nqu =  5    max. err. =  0.334E-04
  y array =       0.668727E+00   0.990191E+00   0.760307E+00   0.807793E+00
   0.117019E+01   0.881055E+00   0.501339E+00   0.720129E+00   0.537977E+00


********************************************************************************


Run with mf = 111.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00
At t =  2.0   nst =  37   hu =   0.115E+00   nqu =  5    max. err. =  0.593E-05
  y array =       0.134047E+00   0.191713E+00   0.137403E+00   0.100789E+00
   0.143787E+00   0.102803E+00   0.384477E-01   0.547819E-01   0.391199E-01
At t =  3.0   nst =  44   hu =   0.155E+00   nqu =  5    max. err. =  0.386E-05
  y array =       0.192920E-01   0.273552E-01   0.193970E-01   0.105624E-01
   0.149714E-01   0.106120E-01   0.291609E-02   0.413247E-02   0.292857E-02

Final statistics for this run:
 rwork size = 308   iwork size =  30
 number of steps =   44
 number of f-s   =   56
 (excluding J-s) =   56
 number of J-s   =    1
 error overrun =  0.99E+00
 number of nonzeros in J =    27
 number of J index groups =    0
 number of LU decomp-s    =    9
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 111.     Input work lengths lrw, liw =   308    30

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00


********************************************************************************


Run with mf = 112.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00
At t =  2.0   nst =  37   hu =   0.115E+00   nqu =  5    max. err. =  0.593E-05
  y array =       0.134047E+00   0.191713E+00   0.137403E+00   0.100789E+00
   0.143787E+00   0.102803E+00   0.384477E-01   0.547819E-01   0.391199E-01
At t =  3.0   nst =  44   hu =   0.155E+00   nqu =  5    max. err. =  0.386E-05
  y array =       0.192920E-01   0.273552E-01   0.193970E-01   0.105624E-01
   0.149714E-01   0.106120E-01   0.291609E-02   0.413247E-02   0.292857E-02

Final statistics for this run:
 rwork size = 315   iwork size =  30
 number of steps =   44
 number of f-s   =   60
 (excluding J-s) =   56
 number of J-s   =    1
 error overrun =  0.99E+00
 number of nonzeros in J =    27
 number of J index groups =    4
 number of LU decomp-s    =    9
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 112.     Input work lengths lrw, liw =   315    30

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00


********************************************************************************


Run with mf = 121.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00
At t =  2.0   nst =  52   hu =   0.105E+00   nqu =  5    max. err. =  0.132E-04
  y array =       0.134044E+00   0.191709E+00   0.137402E+00   0.100788E+00
   0.143786E+00   0.102805E+00   0.384531E-01   0.547890E-01   0.391276E-01
At t =  3.0   nst =  61   hu =   0.132E+00   nqu =  5    max. err. =  0.134E-04
  y array =       0.192907E-01   0.273543E-01   0.193977E-01   0.105672E-01
   0.149788E-01   0.106186E-01   0.292280E-02   0.414233E-02   0.293619E-02

Final statistics for this run:
 rwork size = 245   iwork size =  30
 number of steps =   61
 number of f-s   =   71
 (excluding J-s) =   71
 number of J-s   =    2
 error overrun =  0.25E+01
 number of nonzeros in J =    27
 number of J index groups =    0
 number of LU decomp-s    =    8
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 121.     Input work lengths lrw, liw =   245    30

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00


********************************************************************************


Run with mf = 122.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00
At t =  2.0   nst =  52   hu =   0.105E+00   nqu =  5    max. err. =  0.132E-04
  y array =       0.134044E+00   0.191709E+00   0.137402E+00   0.100788E+00
   0.143786E+00   0.102805E+00   0.384531E-01   0.547890E-01   0.391276E-01
At t =  3.0   nst =  61   hu =   0.132E+00   nqu =  5    max. err. =  0.134E-04
  y array =       0.192907E-01   0.273543E-01   0.193977E-01   0.105672E-01
   0.149788E-01   0.106186E-01   0.292280E-02   0.414233E-02   0.293619E-02

Final statistics for this run:
 rwork size = 252   iwork size =  30
 number of steps =   61
 number of f-s   =   79
 (excluding J-s) =   71
 number of J-s   =    2
 error overrun =  0.25E+01
 number of nonzeros in J =    27
 number of J index groups =    4
 number of LU decomp-s    =    8
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 122.     Input work lengths lrw, liw =   252    30

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00


********************************************************************************


Run with mf = 211.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00
At t =  2.0   nst =  37   hu =   0.115E+00   nqu =  5    max. err. =  0.593E-05
  y array =       0.134047E+00   0.191713E+00   0.137403E+00   0.100789E+00
   0.143787E+00   0.102803E+00   0.384477E-01   0.547819E-01   0.391199E-01
At t =  3.0   nst =  44   hu =   0.155E+00   nqu =  5    max. err. =  0.386E-05
  y array =       0.192920E-01   0.273552E-01   0.193970E-01   0.105624E-01
   0.149714E-01   0.106120E-01   0.291609E-02   0.413247E-02   0.292857E-02

Final statistics for this run:
 rwork size = 308   iwork size =  30
 number of steps =   44
 number of f-s   =   56
 (excluding J-s) =   56
 number of J-s   =    1
 error overrun =  0.99E+00
 number of nonzeros in J =    27
 number of J index groups =    0
 number of LU decomp-s    =    9
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 211.     Input work lengths lrw, liw =   308    30

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00


********************************************************************************


Run with mf = 212.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00
At t =  2.0   nst =  37   hu =   0.115E+00   nqu =  5    max. err. =  0.593E-05
  y array =       0.134047E+00   0.191713E+00   0.137403E+00   0.100789E+00
   0.143787E+00   0.102803E+00   0.384477E-01   0.547819E-01   0.391199E-01
At t =  3.0   nst =  44   hu =   0.155E+00   nqu =  5    max. err. =  0.386E-05
  y array =       0.192920E-01   0.273552E-01   0.193970E-01   0.105624E-01
   0.149714E-01   0.106120E-01   0.291609E-02   0.413247E-02   0.292857E-02

Final statistics for this run:
 rwork size = 315   iwork size =  30
 number of steps =   44
 number of f-s   =   60
 (excluding J-s) =   56
 number of J-s   =    1
 error overrun =  0.99E+00
 number of nonzeros in J =    27
 number of J index groups =    4
 number of LU decomp-s    =    9
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 212.     Input work lengths lrw, liw =   315    30

At t =  1.0   nst =  27   hu =   0.862E-01   nqu =  5    max. err. =  0.988E-05
  y array =       0.668727E+00   0.990190E+00   0.760307E+00   0.807796E+00
   0.117022E+01   0.881060E+00   0.501339E+00   0.720139E+00   0.537974E+00


********************************************************************************


Run with mf = 221.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00
At t =  2.0   nst =  52   hu =   0.105E+00   nqu =  5    max. err. =  0.132E-04
  y array =       0.134044E+00   0.191709E+00   0.137402E+00   0.100788E+00
   0.143786E+00   0.102805E+00   0.384531E-01   0.547890E-01   0.391276E-01
At t =  3.0   nst =  61   hu =   0.132E+00   nqu =  5    max. err. =  0.134E-04
  y array =       0.192907E-01   0.273543E-01   0.193977E-01   0.105672E-01
   0.149788E-01   0.106186E-01   0.292280E-02   0.414233E-02   0.293619E-02

Final statistics for this run:
 rwork size = 245   iwork size =  30
 number of steps =   61
 number of f-s   =   71
 (excluding J-s) =   71
 number of J-s   =    2
 error overrun =  0.25E+01
 number of nonzeros in J =    27
 number of J index groups =    0
 number of LU decomp-s    =    8
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 221.     Input work lengths lrw, liw =   245    30

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00


********************************************************************************


Run with mf = 222.     Input work lengths lrw, liw =  1000    90

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00
At t =  2.0   nst =  52   hu =   0.105E+00   nqu =  5    max. err. =  0.132E-04
  y array =       0.134044E+00   0.191709E+00   0.137402E+00   0.100788E+00
   0.143786E+00   0.102805E+00   0.384531E-01   0.547890E-01   0.391276E-01
At t =  3.0   nst =  61   hu =   0.132E+00   nqu =  5    max. err. =  0.134E-04
  y array =       0.192907E-01   0.273543E-01   0.193977E-01   0.105672E-01
   0.149788E-01   0.106186E-01   0.292280E-02   0.414233E-02   0.293619E-02

Final statistics for this run:
 rwork size = 252   iwork size =  30
 number of steps =   61
 number of f-s   =   79
 (excluding J-s) =   71
 number of J-s   =    2
 error overrun =  0.25E+01
 number of nonzeros in J =    27
 number of J index groups =    4
 number of LU decomp-s    =    8
 nonzeros in strict lower factor =    8
 nonzeros in strict upper factor =   14

 structure descriptor array ian =
   1   3   6   8  11  15  18  21  25  28

 structure descriptor array jan =
   1   2   3   1   2   3   2   1   5   4   2   6   5   4   3   6   5   7   4   8
   9   7   5   8   9   6   8


Run with mf = 222.     Input work lengths lrw, liw =   252    30

At t =  1.0   nst =  38   hu =   0.573E-01   nqu =  5    max. err. =  0.254E-04
  y array =       0.668726E+00   0.990191E+00   0.760308E+00   0.807793E+00
   0.117021E+01   0.881059E+00   0.501348E+00   0.720133E+00   0.537990E+00


********************************************************************************


Number of errors encountered =  0

==========Source and Output for LSODA Demonstration Program=====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODA package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The package is used to solve two simple problems,
c one with a full Jacobian, the other with a banded Jacobian,
c with the 2 appropriate values of jt in each case.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.  All output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1, f2, jac2
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   jt, leniw, lenrw, liw, lout, lrw, mband, mused,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, dtout0, dtout1, er, erm, ero, hu,
     1     rtol, rwork, t, tout, tout1, tsw, y
      dimension y(25), rwork(522), iwork(45)
      data lout/6/, tout1/16.921743d0/, dtout/17.341162d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-8
      lrw = 522
      liw = 45
      iopt = 0
c
c First problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(/'Demonstration program for DLSODA package'////
     1  ' Problem 1:   Van der Pol oscillator:'/
     2  '              xdotdot - 20*(1 - x**2)*xdot + x = 0, ',
     3  '   x(0) = 2, xdot(0) = 0'/' neq =',i2/
     4  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
c
      do 190 jt = 1,2
      write (lout,120) jt
 120  format(//' Solution with jt =',i3//
     1       '  t               x               xdot       meth',
     2       '   nq     h           tsw'//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      dtout0 = 0.5d0*tout1
      dtout1 = 0.5d0*dtout
      tout = dtout0
      ero = 0.0d0
      do 170 iout = 1,nout
        call dlsoda(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1              iopt,rwork,lrw,iwork,liw,jac1,jt)
        hu = rwork(11)
        tsw = rwork(15)
        nqu = iwork(14)
        mused = iwork(19)
        write (lout,140) t,y(1),y(2),mused,nqu,hu,tsw
 140    format(d12.5,d16.5,d14.3,2i6,2d13.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 160
        er = abs(y(1))
        ero = max(ero,er)
        if (er .gt. 1.0d-2) then
          write (lout,150)
 150      format(//' Warning: value at root exceeds 1.0d-2'//)
          nerr = nerr + 1
        endif
 160    if (iout .eq. 1) tout = tout + dtout0
        if (iout .gt. 1) tout = tout + dtout1
 170    continue
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' max. error at root =',d10.2)
 190  continue
c
c Second problem
c
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      atol = 1.0d-6
      nout = 5
      write (lout,210) neq,ml,mu,itol,rtol,atol
 210  format(///80('-')///
     1  ' Problem 2: ydot = A * y , where',
     2  '  A is a banded lower triangular matrix'/
     2  '            derived from 2-D advection PDE'/
     3  ' neq =',i3,'   ml =',i2,'   mu =',i2/
     4  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
      do 290 jt = 4,5
      write (lout,220) jt
 220  format(//' Solution with jt =',i3//
     1       '     t             max.err.     meth   ',
     2       'nq      h            tsw'//)
      t = 0.0d0
      do 230 i = 2,neq
 230    y(i) = 0.0d0
      y(1) = 1.0d0
      itask = 1
      istate = 1
      tout = 0.01d0
      ero = 0.0d0
      do 270 iout = 1,nout
        call dlsoda(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1              iopt,rwork,lrw,iwork,liw,jac2,jt)
        call edit2(y,t,erm)
        hu = rwork(11)
        tsw = rwork(15)
        nqu = iwork(14)
        mused = iwork(19)
        write (lout,240) t,erm,mused,nqu,hu,tsw
 240    format(d15.5,d14.3,2i6,2d14.3)
        if (istate .lt. 0) go to 275
        er = erm/atol
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,150)
          nerr = nerr + 1
        endif
 270    tout = tout*10.0d0
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 5) nfea = nfe - mband*nje
      write (lout,280) lenrw,leniw,nst,nfe,nfea,nje,ero
 280  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' error overrun =',d10.2)
 290  continue
      write (lout,300) nerr
 300  format(///' Number of errors encountered =',i3)
      stop
      end

      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      ydot(1) = y(2)
      ydot(2) = 20.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end

      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(neq), pd(nrowpd,neq)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -40.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 20.0d0*(1.0d0 - y(1)*y(1))
      return
      end

      subroutine f2 (neq, t, y, ydot)
      integer neq, i, j, k, ng
      double precision t, y, ydot, alph1, alph2, d
      dimension y(neq), ydot(neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end

      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(neq), pd(nrowpd,neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end

      subroutine edit2 (y, t, erm)
      integer i, j, k, ng
      double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
      dimension y(*)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      erm = 0.0d0
      if (t .eq. 0.0d0) return
      ex = 0.0d0
      if (t .le. 30.0d0) ex = exp(-2.0d0*t)
      a2 = 1.0d0
      do 60 j = 1,ng
        a1 = 1.0d0
        do 50 i = 1,ng
          k = i + (j - 1)*ng
          yt = t**(i+j-2)*ex*a1*a2
          er = abs(y(k)-yt)
          erm = max(erm,er)
          a1 = a1*alph1/i
 50       continue
        a2 = a2*alph2/j
 60     continue
      return
      end

................................................................................


Demonstration program for DLSODA package



 Problem 1:   Van der Pol oscillator:
              xdotdot - 20*(1 - x**2)*xdot + x = 0,    x(0) = 2, xdot(0) = 0
 neq = 2
 itol =  1   rtol =   0.0E+00   atol =   0.1E-07




 Solution with jt =  1

  t               x               xdot       meth   nq     h           tsw


 0.84609E+01     0.16731E+01    -0.464E-01     2     4    0.209E+00    0.311E+00
 0.16922E+02    -0.11574E-03    -0.141E+02     1     7    0.206E-02    0.158E+02
 0.25592E+02    -0.16828E+01     0.459E-01     2     4    0.240E+00    0.174E+02
 0.34263E+02     0.21448E-03     0.141E+02     1     8    0.293E-02    0.332E+02


 Final statistics for this run:
 rwork size =  52   iwork size =  22
 number of steps =  695
 number of f-s   = 1305
 (excluding J-s) = 1305
 number of J-s   =   30
 max. error at root =  0.21E-03


 Solution with jt =  2

  t               x               xdot       meth   nq     h           tsw


 0.84609E+01     0.16731E+01    -0.464E-01     2     4    0.209E+00    0.311E+00
 0.16922E+02    -0.11574E-03    -0.141E+02     1     7    0.206E-02    0.158E+02
 0.25592E+02    -0.16828E+01     0.459E-01     2     4    0.240E+00    0.174E+02
 0.34263E+02     0.21448E-03     0.141E+02     1     8    0.293E-02    0.332E+02


 Final statistics for this run:
 rwork size =  52   iwork size =  22
 number of steps =  695
 number of f-s   = 1365
 (excluding J-s) = 1305
 number of J-s   =   30
 max. error at root =  0.21E-03



--------------------------------------------------------------------------------


 Problem 2: ydot = A * y , where  A is a banded lower triangular matrix
            derived from 2-D advection PDE
 neq = 25   ml = 5   mu = 0
 itol =  1   rtol =   0.0E+00   atol =   0.1E-05




 Solution with jt =  4

     t             max.err.     meth   nq      h            tsw


    0.10000E-01     0.476E-06     1     2     0.714E-02     0.000E+00
    0.10000E+00     0.988E-06     1     4     0.343E-01     0.000E+00
    0.10000E+01     0.431E-06     1     5     0.724E-01     0.000E+00
    0.10000E+02     0.558E-07     1     3     0.323E+00     0.000E+00
    0.10000E+03     0.127E-11     2     1     0.239E+03     0.170E+02


 Final statistics for this run:
 rwork size = 522   iwork size =  45
 number of steps =  105
 number of f-s   =  207
 (excluding J-s) =  207
 number of J-s   =    3
 error overrun =  0.99E+00


 Solution with jt =  5

     t             max.err.     meth   nq      h            tsw


    0.10000E-01     0.476E-06     1     2     0.714E-02     0.000E+00
    0.10000E+00     0.988E-06     1     4     0.343E-01     0.000E+00
    0.10000E+01     0.431E-06     1     5     0.724E-01     0.000E+00
    0.10000E+02     0.558E-07     1     3     0.323E+00     0.000E+00
    0.10000E+03     0.127E-11     2     1     0.239E+03     0.170E+02


 Final statistics for this run:
 rwork size = 522   iwork size =  45
 number of steps =  105
 number of f-s   =  225
 (excluding J-s) =  207
 number of J-s   =    3
 error overrun =  0.99E+00



 Number of errors encountered =  0


==========Source and Output for LSODAR Demonstration Program====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODAR package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The DLSODAR package is used to solve two simple problems,
c one nonstiff and one intermittently stiff.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.  All output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, gr1, f2, jac2, gr2
      integer iopt, iout, istate, itask, itol, iwork, jroot, jt,
     1   kroot, leniw, lenrw, liw, lrw, lout, neq, nerr, ng,
     2   nfe, nfea, nge, nje, nst
      double precision atol, er, ero, errt, rtol, rwork,
     1   t, tout, tzero, y, yt
      dimension y(2), atol(2), rwork(57), iwork(22), jroot(2)
      data lout/6/
c
      nerr = 0
c-----------------------------------------------------------------------
c First problem.
c The initial value problem is:
c   dy/dt = ((2*log(y) + 8)/t - 5)*y,  y(1) = 1,  1 .le. t .le. 6
c The solution is  y(t) = exp(-t**2 + 5*t - 4)
c The two root functions are:
c   g1 = ((2*log(y)+8)/t - 5)*y (= dy/dt)  (with root at t = 2.5),
c   g2 = log(y) - 2.2491  (with roots at t = 2.47 and 2.53)
c-----------------------------------------------------------------------
c Set all input parameters and print heading.
      neq = 1
      y(1) = 1.0d0
      t = 1.0d0
      tout = 2.0d0
      itol = 1
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      itask = 1
      istate = 1
      iopt = 0
      lrw = 44
      liw = 21
      jt = 2
      ng = 2
      write (lout,110) itol,rtol,atol(1),jt
 110  format(/' Demonstration program for DLSODAR package'////
     1  ' First problem'///
     2  ' Problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1'//
     3  ' Solution is  y(t) = exp(-t**2 + 5*t - 4)'//
     4  ' Root functions are:'/
     5  10x,' g1 = dy/dt  (root at t = 2.5)'/
     6  10x,' g2 = log(y) - 2.2491  (roots at t = 2.47 and t = 2.53)'//
     7  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//
     8  ' jt =',i3///)
c
c Call DLSODAR in loop over tout values 2,3,4,5,6.
      ero = 0.0d0
      do 180 iout = 1,5
 120    continue
        call dlsodar(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jdum,jt,gr1,ng,jroot)
c
c Print y and error in y, and print warning if error too large.
        yt = exp(-t*t + 5.0d0*t - 4.0d0)
        er = y(1) - yt
        write (lout,130) t,y(1),er
 130    format(' At t =',d15.7,5x,'y =',d15.7,5x,'error =',d12.4)
        if (istate .lt. 0) go to 185
        er = abs(er)/(rtol*abs(y(1)) + atol(1))
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,140)
 140      format(//' Warning: error exceeds 1000 * tolerance'//)
          nerr = nerr + 1
          endif
        if (istate .ne. 3) go to 175
c
c If a root was found, write results and check root location.
c Then reset istate to 2 and return to DLSODAR call.
        write (lout,150) t,jroot(1),jroot(2)
 150    format(/' Root found at t =',d15.7,5x,'jroot =',2i5)
        if (jroot(1) .eq. 1) errt = t - 2.5d0
        if (jroot(2) .eq. 1 .and. t .le. 2.5d0) errt = t - 2.47d0
        if (jroot(2) .eq. 1 .and. t .gt. 2.5d0) errt = t - 2.53d0
        write (lout,160) errt
 160    format(' Error in t location of root is',d12.4/)
        if (abs(errt) .gt. 1.0d-3) then
          write (lout,170)
 170      format(//' Warning: root error exceeds 1.0d-3'//)
          nerr = nerr + 1
          endif
        istate = 2
        go to 120
c
c If no root found, increment tout and loop back.
 175    tout = tout + 1.0d0
 180    continue
c
c Problem complete.  Print final statistics.
 185  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,190) lenrw,leniw,nst,nfe,nfea,nje,nge,ero
 190  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding j-s) =',i5/
     5  ' number of j-s   =',i5/
     6  ' number of g-s   =',i5/
     7  ' error overrun =',d10.2)
c
c-----------------------------------------------------------------------
c Second problem (Van der Pol oscillator).
c The initial value problem is (after reduction of 2nd order ODE):
c   dy1/dt = y2,  dy2/dt = 100*(1 - y1**2)*y2 - y1,
c   y1(0) = 2,  y2(0) = 0,  0 .le. t .le. 200
c The root function is  g = y1.
c An analytic solution is not known, but the zeros of y1 are known
c to 15 figures for purposes of checking the accuracy.
c-----------------------------------------------------------------------
c Set tolerance parameters and print heading.
      itol = 2
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      atol(2) = 1.0d-4
      write (lout,200) itol,rtol,atol(1),atol(2)
 200  format(////80('*')//' Second problem (Van der Pol oscillator)'//
     1  ' Problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1'/
     2  '            y1(0) = 2,  y2(0) = 0'//
     3  ' Root function is  g = y1'//
     4  ' itol =',i3,'   rtol =',d10.1,'   atol =',2d10.1)
c
c Loop over jt = 1, 2.  Set remaining parameters and print jt.
      do 290 jt = 1,2
      neq = 2
      y(1) = 2.0d0
      y(2) = 0.0d0
      t = 0.0d0
      tout = 20.0d0
      itask = 1
      istate = 1
      iopt = 0
      lrw = 57
      liw = 22
      ng = 1
      write (lout,210) jt
 210  format(///' Solution with jt =',i2//)
c
c Call DLSODAR in loop over tout values 20,40,...,200.
      do 270 iout = 1,10
 220    continue
        call dlsodar(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,jt,gr2,ng,jroot)
c
c Print y1 and y2.
        write (lout,230) t,y(1),y(2)
 230    format(' At t =',d15.7,5x,'y1 =',d15.7,5x,'y2 =',d15.7)
        if (istate .lt. 0) go to 275
        if (istate .ne. 3) go to 265
c
c If a root was found, write results and check root location.
c Then reset istate to 2 and return to DLSODAR call.
        write (lout,240) t
 240    format(/' Root found at t =',d15.7)
        kroot = int(t/81.2d0 + 0.5d0)
        tzero = 81.17237787055d0 + (kroot-1)*81.41853556212d0
        errt = t - tzero
        write (lout,250) errt
 250    format(' Error in t location of root is',d12.4//)
        if (abs(errt) .gt. 1.0d-1) then
          write (lout,260)
 260      format(//' Warning: root error exceeds 1.0d-1'//)
          nerr = nerr + 1
          endif
        istate = 2
        go to 220
c
c If no root found, increment tout and loop back.
 265    tout = tout + 20.0d0
 270    continue
c
c Problem complete.  Print final statistics.
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,280) lenrw,leniw,nst,nfe,nfea,nje,nge
 280  format(//' Final statistics for this run:'/
     1  '  rwork size =',i4,'   iwork size =',i4/
     2  '  number of steps =',i5/
     3  '  number of f-s   =',i5/
     4  '  (excluding j-s) =',i5/
     5  '  number of j-s   =',i5/
     6  '  number of g-s   =',i5)
 290  continue
c
      write (lout,300) nerr
 300  format(///' Total number of errors encountered =',i3)
      stop
      end

      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(1), ydot(1)
      ydot(1) = ((2.0d0*log(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
      return
      end

      subroutine gr1 (neq, t, y, ng, groot)
      integer neq, ng
      double precision t, y, groot
      dimension y(1), groot(2)
      groot(1) = ((2.0d0*log(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
      groot(2) = log(y(1)) - 2.2491d0
      return
      end

      subroutine f2 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(2), ydot(2)
      ydot(1) = y(2)
      ydot(2) = 100.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end

      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(2), pd(nrowpd,2)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -200.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 100.0d0*(1.0d0 - y(1)*y(1))
      return
      end

      subroutine gr2 (neq, t, y, ng, groot)
      integer neq, ng
      double precision t, y, groot
      dimension y(2), groot(1)
      groot(1) = y(1)
      return
      end

................................................................................


 Demonstration program for DLSODAR package



 First problem


 Problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1

 Solution is  y(t) = exp(-t**2 + 5*t - 4)

 Root functions are:
           g1 = dy/dt  (root at t = 2.5)
           g2 = log(y) - 2.2491  (roots at t = 2.47 and t = 2.53)

 itol =  1   rtol =   0.1E-05   atol =   0.1E-05

 jt =  2



 At t =  0.2000000E+01     y =  0.7389071E+01     error =  0.1534E-04
 At t =  0.2469971E+01     y =  0.9479201E+01     error =  0.1631E-04

 Root found at t =  0.2469971E+01     jroot =    0    1
 Error in t location of root is -0.2865E-04

 At t =  0.2500001E+01     y =  0.9487752E+01     error =  0.1613E-04

 Root found at t =  0.2500001E+01     jroot =    1    0
 Error in t location of root is  0.6800E-06

 At t =  0.2530028E+01     y =  0.9479201E+01     error =  0.1601E-04

 Root found at t =  0.2530028E+01     jroot =    0    1
 Error in t location of root is  0.2813E-04

 At t =  0.3000000E+01     y =  0.7389083E+01     error =  0.2667E-04
 At t =  0.4000000E+01     y =  0.1000007E+01     error =  0.6703E-05
 At t =  0.5000000E+01     y =  0.1831582E-01     error =  0.1814E-06
 At t =  0.6000000E+01     y =  0.4536004E-04     error = -0.3989E-07


 Final statistics for this run:
 rwork size =  42   iwork size =  21
 number of steps =   71
 number of f-s   =  147
 (excluding j-s) =  147
 number of j-s   =    0
 number of g-s   =  114
 error overrun =  0.34E+01




********************************************************************************

 Second problem (Van der Pol oscillator)

 Problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1
            y1(0) = 2,  y2(0) = 0

 Root function is  g = y1

 itol =  2   rtol =   0.1E-05   atol =   0.1E-05   0.1E-03



 Solution with jt = 1


 At t =  0.2000000E+02     y1 =  0.1858228E+01     y2 = -0.7575094E-02
 At t =  0.4000000E+02     y1 =  0.1693230E+01     y2 = -0.9068584E-02
 At t =  0.6000000E+02     y1 =  0.1484608E+01     y2 = -0.1232742E-01
 At t =  0.8000000E+02     y1 =  0.1086291E+01     y2 = -0.5840716E-01
 At t =  0.8116520E+02     y1 = -0.1308482E-12     y2 = -0.6713980E+02

 Root found at t =  0.8116520E+02
 Error in t location of root is -0.7180E-02


 At t =  0.1000000E+03     y1 = -0.1868862E+01     y2 =  0.7497304E-02
 At t =  0.1200000E+03     y1 = -0.1705927E+01     y2 =  0.8930077E-02
 At t =  0.1400000E+03     y1 = -0.1501740E+01     y2 =  0.1196163E-01
 At t =  0.1600000E+03     y1 = -0.1148800E+01     y2 =  0.3568399E-01
 At t =  0.1625761E+03     y1 =  0.1153602E-11     y2 =  0.6713972E+02

 Root found at t =  0.1625761E+03
 Error in t location of root is -0.1485E-01


 At t =  0.1800000E+03     y1 =  0.1879384E+01     y2 = -0.7422067E-02
 At t =  0.2000000E+03     y1 =  0.1718431E+01     y2 = -0.8798201E-02


 Final statistics for this run:
  rwork size =  55   iwork size =  22
  number of steps =  478
  number of f-s   =  931
  (excluding j-s) =  931
  number of j-s   =   42
  number of g-s   =  513



 Solution with jt = 2


 At t =  0.2000000E+02     y1 =  0.1858228E+01     y2 = -0.7575094E-02
 At t =  0.4000000E+02     y1 =  0.1693230E+01     y2 = -0.9068584E-02
 At t =  0.6000000E+02     y1 =  0.1484608E+01     y2 = -0.1232742E-01
 At t =  0.8000000E+02     y1 =  0.1086291E+01     y2 = -0.5840716E-01
 At t =  0.8116520E+02     y1 = -0.8500767E-12     y2 = -0.6713980E+02

 Root found at t =  0.8116520E+02
 Error in t location of root is -0.7180E-02


 At t =  0.1000000E+03     y1 = -0.1868862E+01     y2 =  0.7497304E-02
 At t =  0.1200000E+03     y1 = -0.1705927E+01     y2 =  0.8930077E-02
 At t =  0.1400000E+03     y1 = -0.1501740E+01     y2 =  0.1196163E-01
 At t =  0.1600000E+03     y1 = -0.1148800E+01     y2 =  0.3568399E-01
 At t =  0.1625761E+03     y1 =  0.1057454E-11     y2 =  0.6713972E+02

 Root found at t =  0.1625761E+03
 Error in t location of root is -0.1485E-01


 At t =  0.1800000E+03     y1 =  0.1879384E+01     y2 = -0.7422067E-02
 At t =  0.2000000E+03     y1 =  0.1718431E+01     y2 = -0.8798201E-02


 Final statistics for this run:
  rwork size =  55   iwork size =  22
  number of steps =  478
  number of f-s   = 1015
  (excluding j-s) =  931
  number of j-s   =   42
  number of g-s   =  516



 Total number of errors encountered =  0

==========Source and Output for LSODPK Demonstration Program====================

c-----------------------------------------------------------------------
c Demonstration program for DLSODPK.
c ODE system from ns-species interaction pde in 2 dimensions.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c-----------------------------------------------------------------------
c This program solves a stiff ODE system that arises from a system
c of partial differential equations.  The PDE system is a food web
c population model, with predator-prey interaction and diffusion on
c the unit square in two dimensions.  The dependent variable vector is
c
c         1   2        ns
c   c = (c , c , ..., c  )
c
c and the PDEs are as follows:
c
c     i               i      i
c   dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
c                     xx     yy        i
c
c where
c                  i          ns         j
c   f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
c    i                       j=1
c
c The number of species is ns = 2*np, with the first np being prey and
c the last np being predators.  The coefficients a(i,j), b(i), d(i) are:
c
c   a(i,i) = -a  (all i)
c   a(i,j) = -g  (i .le. np, j .gt. np)
c   a(i,j) =  e  (i .gt. np, j .le. np)
c   b(i) =  b*(1 + alpha*x*y)  (i .le. np)
c   b(i) = -b*(1 + alpha*x*y)  (i .gt. np)
c   d(i) = dprey  (i .le. np)
c   d(i) = dpred  (i .gt. np)
c
c The various scalar parameters are set in subroutine setpar.
c
c The boundary conditions are: normal derivative = 0.
c A polynomial in x and y is used to set the initial conditions.
c
c The PDEs are discretized by central differencing on a mx by my mesh.
c
c The ODE system is solved by DLSODPK using method flag values
c mf = 10, 21, 22, 23, 24, 29.  The final time is tmax = 10, except
c that for mf = 10 it is tmax = 1.0d-3 because the problem is stiff,
c and for mf = 23 and 24 it is tmax = 2 because the lack of symmetry
c in the problem makes these methods more costly.
c
c Two preconditioner matrices are used.  One uses a fixed number of
c Gauss-Seidel iterations based on the diffusion terms only.
c The other preconditioner is a block-diagonal matrix based on
c the partial derivatives of the interaction terms f only, using
c block-grouping (computing only a subset of the ns by ns blocks).
c For mf = 21 and 22, these two preconditioners are applied on
c the left and right, respectively, and for mf = 23 and 24 the product
c of the two is used as the one preconditioner matrix.
c For mf = 29, the inverse of the product is applied.
c
c Two output files are written: one with the problem description and
c and performance statistics on unit 6, and one with solution profiles
c at selected output times (for mf = 22 only) on unit 8.
c-----------------------------------------------------------------------
c Note: In addition to the main program and 10 subroutines
c given below, this program requires the LINPACK subroutines
c DGEFA and DGESL, and the BLAS routine DAXPY.
c-----------------------------------------------------------------------
c Reference:
c     Peter N. Brown and Alan C. Hindmarsh,
c     Reduced Storage Matrix Methods in Stiff ODE Systems,
c     J. Appl. Math. & Comp., 31 (1989), pp. 40-91;
c     Also LLNL Report UCRL-95088, Rev. 1, June 1987.
c-----------------------------------------------------------------------
      external fweb, jacbg, solsbg
      integer ns, mx, my, mxns,
     1        mp, mq, mpsq, itmax,
     2        meshx,meshy,ngx,ngy,ngrp,mxmp,jgx,jgy,jigx,jigy,jxr,jyr
      integer i, imod3, iopt, iout, istate, itask, itol, iwork,
     1   jacflg, jpre, leniw, lenrw, liw, lrw, mf,
     2   ncfl, ncfn, neq, nfe, nfldif, nfndif, nli, nlidif, nni, nnidif,
     3   nout, npe, nps, nqu, nsdif, nst
      double precision aa, ee, gg, bb, dprey, dpred,
     1     ax, ay, acoef, bcoef, dx, dy, alph, diff, cox, coy,
     2     uround, srur
      double precision avdim, atol, cc, hu, rcfl, rcfn, dumach,
     1   rtol, rwork, t, tout
c
c The problem Common blocks below allow for up to 20 species,
c up to a 50x50 mesh, and up to a 20x20 group structure.
      common /pcom0/ aa, ee, gg, bb, dprey, dpred
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
      common /pcom3/ meshx, meshy, ngx, ngy, ngrp, mxmp,
     2   jgx(21), jgy(21), jigx(50), jigy(50), jxr(20), jyr(20)
c
c The dimension of cc below must be .ge. 2*neq, where neq = ns*mx*my.
c The dimension lrw of rwork must be .ge. 17*neq + ns*ns*ngrp + 61,
c and the dimension liw of iwork must be .ge. 35 + ns*ngrp.
      dimension cc(576), rwork(5213), iwork(67)
      data lrw/5213/, liw/67/
c
      open (unit=6, file='demout', status='new')
      open (unit=8, file='ccout', status='new')
c
      ax = 1.0d0
      ay = 1.0d0
c
c Call setpar to set problem parameters.
      call setpar
c
c Set remaining problem parameters.
      neq = ns*mx*my
      mxns = mx*ns
      dx = ax/(mx-1)
      dy = ay/(my-1)
      do 10 i = 1,ns
        cox(i) = diff(i)/dx**2
 10     coy(i) = diff(i)/dy**2
c
c Write heading.
      write(6,20)ns, mx,my,neq
 20   format(' Demonstration program for DLSODPK package'//
     1   ' Food web problem with ns species, ns =',i4/
     2   ' Predator-prey interaction and diffusion on a 2-d square'//
     3   ' Mesh dimensions (mx,my) =',2i4/
     4   ' Total system size is neq =',i7//)
      write(6,25) aa,ee,gg,bb,dprey,dpred,alph
 25   format(' Matrix parameters:  a =',d12.4,'   e =',d12.4,
     1   '   g =',d12.4/20x,' b =',d12.4//
     2   ' Diffusion coefficients: dprey =',d12.4,'   dpred =',d12.4/
     3   ' Rate parameter alpha =',d12.4//)
c
c Set remaining method parameters.
      jpre = 3
      jacflg = 1
      iwork(3) = jpre
      iwork(4) = jacflg
      iopt = 0
      mp = ns
      mq = mx*my
      mpsq = ns*ns
      uround = dumach()
      srur = sqrt(uround)
      meshx = mx
      meshy = my
      mxmp = meshx*mp
      ngx = 2
      ngy = 2
      ngrp = ngx*ngy
      call gset (meshx, ngx, jgx, jigx, jxr)
      call gset (meshy, ngy, jgy, jigy, jyr)
      iwork(1) = mpsq*ngrp
      iwork(2) = mp*ngrp
      itmax = 5
      itol = 1
      rtol = 1.0d-5
      atol = rtol
      itask = 1
      write(6,30)ngrp,ngx,ngy,itmax,rtol,atol
 30   format(' Preconditioning uses interaction-only block-diagonal',
     1   ' matrix'/' with block-grouping, and Gauss-Seidel iterations'//
     2   ' Number of diagonal block groups = ngrp =',i4,
     3   '   (ngx by ngy, ngx =',i2,'  ngy =',i2,' )'//
     4   ' G-S preconditioner uses itmax iterations, itmax =',i3//
     5   ' Tolerance parameters: rtol =',d10.2,'   atol =',d10.2)
c
c
c Loop over mf values 10, 21, ..., 24, 29.
c
      do 90 mf = 10,29
      if (mf .gt. 10 .and. mf .lt. 21) go to 90
      if (mf .gt. 24 .and. mf .lt. 29) go to 90
      write(6,40)mf
 40   format(//80('-')//' Solution with mf =',i3//
     1   '   t       nstep  nfe  nni  nli  npe  nq',
     2   4x,'h          avdim    ncf rate    lcf rate')
c
      t = 0.0d0
      tout = 1.0d-8
      nout = 18
      if (mf .eq. 10) nout = 6
      if (mf .eq. 23 .or. mf .eq. 24) nout = 10
      call cinit (cc)
      if (mf .eq. 22) call outweb (t, cc, ns, mx, my, 8)
      istate = 1
      nli = 0
      nni = 0
      ncfn = 0
      ncfl = 0
      nst = 0
c
c Loop over output times and call DLSODPK.
c
      do 70 iout = 1,nout
        call dlsodpk (fweb, neq, cc, t, tout, itol, rtol, atol, itask,
     1         istate, iopt, rwork, lrw, iwork, liw, jacbg, solsbg, mf)
        nsdif = iwork(11) - nst
        nst = iwork(11)
        nfe = iwork(12)
        npe = iwork(13)
        nnidif = iwork(19) - nni
        nni = iwork(19)
        nlidif = iwork(20) - nli
        nli = iwork(20)
        nfndif = iwork(22) - ncfn
        ncfn = iwork(22)
        nfldif = iwork(23) - ncfl
        ncfl = iwork(23)
        nqu = iwork(14)
        hu = rwork(11)
        avdim = 0.0d0
        rcfn = 0.0d0
        rcfl = 0.0d0
        if (nnidif .gt. 0) avdim = real(nlidif)/real(nnidif)
        if (nsdif .gt. 0) rcfn = real(nfndif)/real(nsdif)
        if (nnidif .gt. 0) rcfl = real(nfldif)/real(nnidif)
        write(6,50)t,nst,nfe,nni,nli,npe,nqu,hu,avdim,rcfn,rcfl
 50     format(d10.2,i5,i6,3i5,i4,2d11.2,d10.2,d12.2)
        imod3 = iout - 3*(iout/3)
        if (mf .eq. 22 .and. imod3 .eq. 0) call outweb (t,cc,ns,mx,my,8)
        if (istate .eq. 2) go to 65
        write(6,60)t
 60     format(//' final time reached =',d12.4//)
        go to 75
 65     continue
        if (tout .gt. 0.9d0) tout = tout + 1.0d0
        if (tout .lt. 0.9d0) tout = tout*10.0d0
 70     continue
c
 75   continue
      nst = iwork(11)
      nfe = iwork(12)
      npe = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nni = iwork(19)
      nli = iwork(20)
      nps = iwork(21)
      if (nni .gt. 0) avdim = real(nli)/real(nni)
      ncfn = iwork(22)
      ncfl = iwork(23)
      write (6,80) lenrw,leniw,nst,nfe,npe,nps,nni,nli,avdim,
     1               ncfn,ncfl
 80   format(//' Final statistics for this run:'/
     1   ' rwork size =',i8,'   iwork size =',i6/
     2   ' number of time steps            =',i5/
     3   ' number of f evaluations         =',i5/
     4   ' number of preconditioner evals. =',i5/
     4   ' number of preconditioner solves =',i5/
     5   ' number of nonlinear iterations  =',i5/
     5   ' number of linear iterations     =',i5/
     6   ' average subspace dimension  =',f8.4/
     7   i5,' nonlinear conv. failures,',i5,' linear conv. failures')
c
 90   continue
      stop
c------  end of main program for DLSODPK demonstration program ----------
      end

      subroutine setpar
c-----------------------------------------------------------------------
c This routine sets the problem parameters.
c It set ns, mx, my, and problem coefficients acoef, bcoef, diff, alph,
c using parameters np, aa, ee, gg, bb, dprey, dpred.
c-----------------------------------------------------------------------
      integer ns, mx, my, mxns
      integer i, j, np
      double precision aa, ee, gg, bb, dprey, dpred,
     1     ax, ay, acoef, bcoef, dx, dy, alph, diff, cox, coy
      common /pcom0/ aa, ee, gg, bb, dprey, dpred
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      np = 3
      mx = 6
      my = 6
      aa = 1.0d0
      ee = 1.0d4
      gg = 0.5d-6
      bb = 1.0d0
      dprey = 1.0d0
      dpred = 0.5d0
      alph = 1.0d0
      ns = 2*np
      do 70 j = 1,np
        do 60 i = 1,np
          acoef(np+i,j) = ee
          acoef(i,np+j) = -gg
 60       continue
        acoef(j,j) = -aa
        acoef(np+j,np+j) = -aa
        bcoef(j) = bb
        bcoef(np+j) = -bb
        diff(j) = dprey
        diff(np+j) = dpred
 70     continue
c
      return
c------------  end of subroutine setpar  -------------------------------
      end

      subroutine gset (m, ng, jg, jig, jr)
c-----------------------------------------------------------------------
c This routine sets arrays jg, jig, and jr describing
c a uniform partition of (1,2,...,m) into ng groups.
c-----------------------------------------------------------------------
      integer m, ng, jg, jig, jr
      dimension jg(*), jig(*), jr(*)
      integer ig, j, len1, mper, ngm1
c
      mper = m/ng
      do 10 ig = 1,ng
 10     jg(ig) = 1 + (ig - 1)*mper
      jg(ng+1) = m + 1
c
      ngm1 = ng - 1
      len1 = ngm1*mper
      do 20 j = 1,len1
 20     jig(j) = 1 + (j-1)/mper
      len1 = len1 + 1
      do 25 j = len1,m
 25     jig(j) = ng
c
      do 30 ig = 1,ngm1
 30     jr(ig) = 0.5d0 + (ig - 0.5d0)*mper
      jr(ng) = 0.5d0*(1 + ngm1*mper + m)
c
      return
c------------  end of subroutine gset  ---------------------------------
      end

      subroutine cinit (cc)
c-----------------------------------------------------------------------
c This routine computes and loads the vector of initial values.
c-----------------------------------------------------------------------
      double precision cc
      dimension cc(*)
      integer ns, mx, my, mxns
      integer i, ici, ioff, iyoff, jx, jy
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision argx, argy, x, y
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
        do 20 jy = 1,my
          y = (jy-1)*dy
          argy = 16.0d0*y*y*(ay-y)*(ay-y)
          iyoff = mxns*(jy-1)
          do 10 jx = 1,mx
            x = (jx-1)*dx
            argx = 16.0d0*x*x*(ax-x)*(ax-x)
            ioff = iyoff + ns*(jx-1)
            do 5 i = 1,ns
              ici = ioff + i
              cc(ici) = 10.0d0 + i*argx*argy
  5           continue
 10         continue
 20       continue
      return
c------------  end of subroutine cinit  --------------------------------
      end

      subroutine outweb (t, c, ns, mx, my, lun)
c-----------------------------------------------------------------------
c This routine prints the values of the individual species densities
c at the current time t.  The write statements use unit lun.
c-----------------------------------------------------------------------
      integer ns, mx, my, lun
      double precision t, c
      dimension c(ns,mx,my)
      integer i, jx, jy
c
      write(lun,10) t
 10   format(/80('-')/30x,'At time t = ',d16.8/80('-') )
c
      do 40 i = 1,ns
        write(lun,20) i
 20     format(' the species c(',i2,') values are:')
        do 30 jy = my,1,-1
          write(lun,25) (c(i,jx,jy),jx=1,mx)
 25       format(6(1x,g12.6))
 30       continue
        write(lun,35)
 35     format(80('-'),/)
 40     continue
c
      return
c------------  end of subroutine outweb  -------------------------------
      end

      subroutine fweb (neq, t, cc, cdot)
c-----------------------------------------------------------------------
c This routine computes the derivative of cc and returns it in cdot.
c The interaction rates are computed by calls to webr, and these are
c saved in cc(neq+1),...,cc(2*neq) for use in preconditioning.
c-----------------------------------------------------------------------
      integer neq
      double precision t, cc, cdot
      dimension cc(neq), cdot(neq)
      integer ns, mx, my, mxns
      integer i, ic, ici, idxl, idxu, idyl, idyu, iyoff, jx, jy
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision dcxli, dcxui, dcyli, dcyui, x, y
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      do 100 jy = 1,my
        y = (jy-1)*dy
        iyoff = mxns*(jy-1)
        idyu = mxns
        if (jy .eq. my) idyu = -mxns
        idyl = mxns
        if (jy .eq. 1) idyl = -mxns
        do 90 jx = 1,mx
          x = (jx-1)*dx
          ic = iyoff + ns*(jx-1) + 1
c Get interaction rates at one point (x,y).
          call webr (x, y, t, cc(ic), cc(neq+ic))
          idxu = ns
          if (jx .eq. mx) idxu = -ns
          idxl = ns
          if (jx .eq. 1) idxl = -ns
          do 80 i = 1,ns
            ici = ic + i - 1
c Do differencing in y.
            dcyli = cc(ici) - cc(ici-idyl)
            dcyui = cc(ici+idyu) - cc(ici)
c Do differencing in x.
            dcxli = cc(ici) - cc(ici-idxl)
            dcxui = cc(ici+idxu) - cc(ici)
c Collect terms and load cdot elements.
            cdot(ici) = coy(i)*(dcyui - dcyli) + cox(i)*(dcxui - dcxli)
     1                  + cc(neq+ici)
 80         continue
 90       continue
 100    continue
      return
c------------  end of subroutine fweb  ---------------------------------
      end

      subroutine webr (x, y, t, c, rate)
c-----------------------------------------------------------------------
c This routine computes the interaction rates for the species
c c(1),...,c(ns), at one spatial point and at time t.
c-----------------------------------------------------------------------
      double precision x, y, t, c, rate
      dimension c(*), rate(*)
      integer ns, mx, my, mxns
      integer i
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision fac
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      do 10 i = 1,ns
 10     rate(i) = 0.0d0
      do 15 i = 1,ns
        call daxpy (ns, c(i), acoef(1,i), 1, rate, 1)
 15     continue
      fac = 1.0d0 + alph*x*y
      do 20 i = 1,ns
 20     rate(i) = c(i)*(bcoef(i)*fac + rate(i))
      return
c------------  end of subroutine webr  ---------------------------------
      end

      subroutine jacbg (f, neq, t, cc, ccsv, rewt, f0, f1, hl0,
     1                  bd, ipbd, ier)
c-----------------------------------------------------------------------
c This routine generates part of the block-diagonal part of the
c Jacobian, multiplies by -hl0, adds the identity matrix,
c and calls DGEFA to do LU decomposition of each diagonal block.
c The computation of the diagonal blocks uses the block and grouping
c information in /pcom1/ and /pcom2/.  One block per group is computed.
c The Jacobian elements are generated by difference quotients
c using calls to the routine fbg.
c-----------------------------------------------------------------------
c The two Common blocks below are used for internal communication.
c The variables used are:
c   mp     = size of blocks in block-diagonal preconditioning matrix.
c   mq     = number of blocks in each direction (neq = mp*mq).
c   mpsq   = mp*mp.
c   uround = unit roundoff, generated by a call uround = dumach().
c   srur   = sqrt(uround).
c   meshx  = x mesh size
c   meshy  = y mesh size (mesh is meshx by meshy)
c   ngx    = no. groups in x direction in block-grouping scheme.
c   ngy    = no. groups in y direction in block-grouping scheme.
c   ngrp   = total number of groups = ngx*ngy.
c   mxmp   = meshx*mp.
c   jgx    = length ngx+1 array of group boundaries in x direction.
c            group igx has x indices jx = jgx(igx),...,jgx(igx+1)-1.
c   jigx   = length meshx array of x group indices vs x node index.
c            x node index jx is in x group jigx(jx).
c   jxr    = length ngx array of x indices representing the x groups.
c            the index for x group igx is jx = jxr(igx).
c   jgy, jigy, jyr = analogous arrays for grouping in y direction.
c-----------------------------------------------------------------------
      external f
      integer neq, ipbd, ier
      double precision t, cc, ccsv, rewt, f0, f1, hl0, bd
      dimension cc(neq), ccsv(neq), rewt(neq), f0(neq), f1(neq),
     1          bd(*), ipbd(*)
      integer mp, mq, mpsq, itmax,
     2        meshx,meshy,ngx,ngy,ngrp,mxmp,jgx,jgy,jigx,jigy,jxr,jyr
      integer i, ibd, idiag, if0, if00, ig, igx, igy, iip,
     1   j, jj, jx, jy, n
      double precision uround, srur
      double precision fac, r, r0, dvnorm
c
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
      common /pcom3/ meshx, meshy, ngx, ngy, ngrp, mxmp,
     1   jgx(21), jgy(21), jigx(50), jigy(50), jxr(20), jyr(20)
c
      n = neq
c
c-----------------------------------------------------------------------
c Make mp calls to fbg to approximate each diagonal block of Jacobian.
c Here cc(neq+1),...,cc(2*neq) contains the base fb value.
c r0 is a minimum increment factor for the difference quotient.
c-----------------------------------------------------------------------
 200  fac = dvnorm (n, f0, rewt)
      r0 = 1000.0d0*abs(hl0)*uround*n*fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      ibd = 0
      do 240 igy = 1,ngy
        jy = jyr(igy)
        if00 = (jy - 1)*mxmp
        do 230 igx = 1,ngx
          jx = jxr(igx)
          if0 = if00 + (jx - 1)*mp
          do 220 j = 1,mp
            jj = if0 + j
            r = max(srur*abs(cc(jj)),r0/rewt(jj))
            cc(jj) = cc(jj) + r
            fac = -hl0/r
            call fbg (neq, t, cc, jx, jy, f1)
            do 210 i = 1,mp
 210          bd(ibd+i) = (f1(i) - cc(neq+if0+i))*fac
            cc(jj) = ccsv(jj)
            ibd = ibd + mp
 220        continue
 230      continue
 240    continue
c
c Add identity matrix and do LU decompositions on blocks. --------------
 260  continue
      ibd = 1
      iip = 1
      do 280 ig = 1,ngrp
        idiag = ibd
        do 270 i = 1,mp
          bd(idiag) = bd(idiag) + 1.0d0
 270      idiag = idiag + (mp + 1)
        call dgefa (bd(ibd), mp, mp, ipbd(iip), ier)
        if (ier .ne. 0) go to 290
        ibd = ibd + mpsq
        iip = iip + mp
 280    continue
 290  return
c------------  end of subroutine jacbg  --------------------------------
      end

      subroutine fbg (neq, t, cc, jx, jy, cdot)
c-----------------------------------------------------------------------
c This routine computes one block of the interaction terms of the
c system, namely block (jx,jy), for use in preconditioning.
c-----------------------------------------------------------------------
      integer neq, jx, jy
      double precision t, cc, cdot
      dimension cc(neq), cdot(neq)
      integer ns, mx, my, mxns
      integer iblok, ic
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision x, y
c
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      iblok = jx + (jy-1)*mx
        y = (jy-1)*dy
          x = (jx-1)*dx
          ic = ns*(iblok-1) + 1
          call webr (x, y, t, cc(ic), cdot)
      return
c------------  end of subroutine fbg  ----------------------------------
      end

      subroutine solsbg (n, t, cc, f0, wk, hl0, bd, ipbd, v, lr, ier)
c-----------------------------------------------------------------------
c This routine applies one or two inverse preconditioner matrices
c to the array v, using the interaction-only block-diagonal Jacobian
c with block-grouping, and Gauss-Seidel applied to the diffusion terms.
c When lr = 1 or 3, it calls gs for a Gauss-Seidel approximation
c to ((I-hl0*Jd)-inverse)*v, and stores the result in v.
c When lr = 2 or 3, it computes ((I-hl0*dg/dc)-inverse)*v, using LU
c factors of the blocks in bd, and pivot information in ipbd.
c In both cases, the array v is overwritten with the solution.
c-----------------------------------------------------------------------
      integer n, ipbd, lr, ier
      double precision t, cc, f0, wk, hl0, bd, v
      dimension cc(n), f0(n), wk(n), bd(*), ipbd(*), v(n)
      integer mp, mq, mpsq, itmax,
     2        meshx,meshy,ngx,ngy,ngrp,mxmp,jgx,jgy,jigx,jigy,jxr,jyr
      integer ibd, ig0, igm1, igx, igy, iip, iv, jx, jy
      double precision uround, srur
c
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
      common /pcom3/ meshx, meshy, ngx, ngy, ngrp, mxmp,
     1   jgx(21), jgy(21), jigx(50), jigy(50), jxr(20), jyr(20)
c
      ier = 0
c
      if (lr.eq.0 .or. lr.eq.1 .or. lr.eq.3) call gs (n, hl0, v, wk)
      if (lr.eq.0 .or. lr.eq.2 .or. lr.eq.3) then
        iv = 1
        do 20 jy = 1,meshy
          igy = jigy(jy)
          ig0 = (igy - 1)*ngx
          do 10 jx = 1,meshx
            igx = jigx(jx)
            igm1 = igx - 1 + ig0
            ibd = 1 + igm1*mpsq
            iip = 1 + igm1*mp
            call dgesl (bd(ibd), mp, mp, ipbd(iip), v(iv), 0)
            iv = iv + mp
 10         continue
 20       continue
        endif
c
      return
c------------  end of subroutine solsbg  -------------------------------
      end

      subroutine gs (n, hl0, z, x)
c-----------------------------------------------------------------------
c This routine performs itmax Gauss-Seidel iterations
c to compute an approximation to P-inverse*z, where P = I - hl0*Jd, and
c Jd represents the diffusion contributions to the Jacobian.
c z contains the answer on return.
c The dimensions below assume ns .le. 20.
c-----------------------------------------------------------------------
      integer n
      double precision hl0, z, x
      dimension z(n), x(n)
      integer ns, mx, my, mxns,
     1        mp, mq, mpsq, itmax
      integer i, ic, ici, iter, iyoff, jx, jy
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy,
     2     uround, srur
      double precision beta,beta2,cof1,elamda,gamma,gamma2
      dimension beta(20), gamma(20), beta2(20), gamma2(20), cof1(20)
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
c
c-----------------------------------------------------------------------
c Write matrix as P = D - L - U.
c Load local arrays beta, beta2, gamma, gamma2, and cof1.
c-----------------------------------------------------------------------
      do 10 i = 1,ns
        elamda = 1.d0/(1.d0 + 2.d0*hl0*(cox(i) + coy(i)))
        beta(i) = hl0*cox(i)*elamda
        beta2(i) = 2.d0*beta(i)
        gamma(i) = hl0*coy(i)*elamda
        gamma2(i) = 2.d0*gamma(i)
        cof1(i) = elamda
 10     continue
c-----------------------------------------------------------------------
c Begin iteration loop.
c Load array x with (D-inverse)*z for first iteration.
c-----------------------------------------------------------------------
      iter = 1
c
      do 50 jy = 1,my
        iyoff = mxns*(jy-1)
        do 40 jx = 1,mx
          ic = iyoff + ns*(jx-1)
          do 30 i = 1,ns
            ici = ic + i
            x(ici) = cof1(i)*z(ici)
            z(ici) = 0.d0
 30         continue
 40       continue
 50     continue
      go to 160
c-----------------------------------------------------------------------
c Calculate (D-inverse)*U*x.
c-----------------------------------------------------------------------
 70   continue
      iter = iter + 1
      jy = 1
        jx = 1
        ic = ns*(jx-1)
        do 75 i = 1,ns
          ici = ic + i
 75       x(ici) = beta2(i)*x(ici+ns) + gamma2(i)*x(ici+mxns)
        do 85 jx = 2,mx-1
          ic = ns*(jx-1)
          do 80 i = 1,ns
            ici = ic + i
 80         x(ici) = beta(i)*x(ici+ns) + gamma2(i)*x(ici+mxns)
 85       continue
        jx = mx
        ic = ns*(jx-1)
        do 90 i = 1,ns
          ici = ic + i
 90       x(ici) = gamma2(i)*x(ici+mxns)
      do 115 jy = 2,my-1
        iyoff = mxns*(jy-1)
          jx = 1
          ic = iyoff
          do 95 i = 1,ns
            ici = ic + i
 95         x(ici) = beta2(i)*x(ici+ns) + gamma(i)*x(ici+mxns)
          do 105 jx = 2,mx-1
            ic = iyoff + ns*(jx-1)
            do 100 i = 1,ns
              ici = ic + i
 100          x(ici) = beta(i)*x(ici+ns) + gamma(i)*x(ici+mxns)
 105        continue
          jx = mx
          ic = iyoff + ns*(jx-1)
          do 110 i = 1,ns
            ici = ic + i
 110        x(ici) = gamma(i)*x(ici+mxns)
 115      continue
      jy = my
      iyoff = mxns*(jy-1)
        jx = 1
        ic = iyoff
        do 120 i = 1,ns
          ici = ic + i
 120      x(ici) = beta2(i)*x(ici+ns)
        do 130 jx = 2,mx-1
          ic = iyoff + ns*(jx-1)
          do 125 i = 1,ns
            ici = ic + i
 125      x(ici) = beta(i)*x(ici+ns)
 130      continue
        jx = mx
        ic = iyoff + ns*(jx-1)
        do 135 i = 1,ns
          ici = ic + i
 135      x(ici) = 0.0d0
c-----------------------------------------------------------------------
c Calculate (I - (D-inverse)*L)-inverse * x.
c-----------------------------------------------------------------------
 160  continue
      jy = 1
        do 175 jx = 2,mx-1
          ic = ns*(jx-1)
          do 170 i = 1,ns
            ici = ic + i
 170        x(ici) = x(ici) + beta(i)*x(ici-ns)
 175      continue
        jx = mx
        ic = ns*(jx-1)
        do 180 i = 1,ns
          ici = ic + i
 180      x(ici) = x(ici) + beta2(i)*x(ici-ns)
      do 210 jy = 2,my-1
        iyoff = mxns*(jy-1)
          jx = 1
          ic = iyoff
          do 185 i = 1,ns
            ici = ic + i
 185        x(ici) = x(ici) + gamma(i)*x(ici-mxns)
          do 200 jx = 2,mx-1
            ic = iyoff + ns*(jx-1)
            do 195 i = 1,ns
              ici = ic + i
              x(ici) = (x(ici) + beta(i)*x(ici-ns))
     1             + gamma(i)*x(ici-mxns)
 195          continue
 200        continue
            jx = mx
            ic = iyoff + ns*(jx-1)
            do 205 i = 1,ns
              ici = ic + i
              x(ici) = (x(ici) + beta2(i)*x(ici-ns))
     1             + gamma(i)*x(ici-mxns)
 205          continue
 210        continue
      jy = my
      iyoff = mxns*(jy-1)
        jx = 1
        ic = iyoff
        do 215 i = 1,ns
          ici = ic + i
 215      x(ici) = x(ici) + gamma2(i)*x(ici-mxns)
        do 225 jx = 2,mx-1
          ic = iyoff + ns*(jx-1)
          do 220 i = 1,ns
            ici = ic + i
            x(ici) = (x(ici) + beta(i)*x(ici-ns))
     1           + gamma2(i)*x(ici-mxns)
 220        continue
 225      continue
        jx = mx
        ic = iyoff + ns*(jx-1)
        do 230 i = 1,ns
          ici = ic + i
          x(ici) = (x(ici) + beta2(i)*x(ici-ns))
     1         + gamma2(i)*x(ici-mxns)
 230      continue
c-----------------------------------------------------------------------
c Add increment x to z.
c-----------------------------------------------------------------------
      do 300 i = 1,n
 300    z(i) = z(i) + x(i)
c
      if (iter .lt. itmax) go to 70
      return
c------------  end of subroutine gs  -----------------------------------
      end

................................................................................

 Demonstration program for DLSODPK package

 Food web problem with ns species, ns =   6
 Predator-prey interaction and diffusion on a 2-d square

 Mesh dimensions (mx,my) =   6   6
 Total system size is neq =    216


 Matrix parameters:  a =  0.1000E+01   e =  0.1000E+05   g =  0.5000E-06
                     b =  0.1000E+01

 Diffusion coefficients: dprey =  0.1000E+01   dpred =  0.5000E+00
 Rate parameter alpha =  0.1000E+01


 Preconditioning uses interaction-only block-diagonal matrix
 with block-grouping, and Gauss-Seidel iterations

 Number of diagonal block groups = ngrp =   4   (ngx by ngy, ngx = 2  ngy = 2 )

 G-S preconditioner uses itmax iterations, itmax =  5

 Tolerance parameters: rtol =  0.10E-04   atol =  0.10E-04


--------------------------------------------------------------------------------

 Solution with mf = 10

   t       nstep  nfe  nni  nli  npe  nq    h          avdim    ncf rate    lcf rate
  0.10E-07    3     4    0    0    0   2   0.10E-06   0.00E+00  0.00E+00    0.00E+00
  0.10E-06    3     4    0    0    0   2   0.10E-06   0.00E+00  0.00E+00    0.00E+00
  0.10E-05    7    13    0    0    0   3   0.35E-06   0.00E+00  0.00E+00    0.00E+00
  0.10E-04   22    41    0    0    0   6   0.10E-05   0.00E+00  0.00E+00    0.00E+00
  0.10E-03   81   155    0    0    0   2   0.24E-05   0.00E+00  0.34E-01    0.00E+00
  0.10E-02  397   899    0    0    0   2   0.43E-05   0.00E+00  0.22E+00    0.00E+00


 Final statistics for this run:
 rwork size =    3476   iwork size =    30
 number of time steps            =  397
 number of f evaluations         =  899
 number of preconditioner evals. =    0
 number of preconditioner solves =    0
 number of nonlinear iterations  =    0
 number of linear iterations     =    0
 average subspace dimension  =  0.0000
   73 nonlinear conv. failures,    0 linear conv. failures


--------------------------------------------------------------------------------

 Solution with mf = 21

   t       nstep  nfe  nni  nli  npe  nq    h          avdim    ncf rate    lcf rate
  0.10E-07    3     5    3    1    2   2   0.65E-07   0.33E+00  0.00E+00    0.00E+00
  0.10E-06    4     8    5    2    2   2   0.65E-07   0.50E+00  0.00E+00    0.00E+00
  0.10E-05   10    22   13    8    4   3   0.37E-06   0.75E+00  0.00E+00    0.00E+00
  0.10E-04   33    78   43   34    6   5   0.51E-06   0.87E+00  0.00E+00    0.00E+00
  0.10E-03  112   270  134  135   14   5   0.60E-05   0.11E+01  0.00E+00    0.00E+00
  0.10E-02  131   321  156  164   18   2   0.35E-03   0.13E+01  0.00E+00    0.00E+00
  0.10E-01  138   343  167  175   20   3   0.18E-02   0.10E+01  0.00E+00    0.00E+00
  0.10E+00  162   402  194  207   23   4   0.71E-02   0.12E+01  0.00E+00    0.00E+00
  0.10E+01  206   540  242  297   27   4   0.47E-01   0.19E+01  0.00E+00    0.00E+00
  0.20E+01  221   608  259  348   29   4   0.96E-01   0.30E+01  0.00E+00    0.00E+00
  0.30E+01  230   656  269  386   30   4   0.16E+00   0.38E+01  0.00E+00    0.00E+00
  0.40E+01  236   694  276  417   31   4   0.21E+00   0.44E+01  0.00E+00    0.29E+00
  0.50E+01  240   718  280  437   31   4   0.27E+00   0.50E+01  0.00E+00    0.75E+00
  0.60E+01  244   742  284  457   31   4   0.27E+00   0.50E+01  0.00E+00    0.75E+00
  0.70E+01  247   766  288  477   32   4   0.36E+00   0.50E+01  0.00E+00    0.10E+01
  0.80E+01  250   790  292  497   33   4   0.48E+00   0.50E+01  0.00E+00    0.10E+01
  0.90E+01  252   802  294  507   33   4   0.48E+00   0.50E+01  0.00E+00    0.10E+01
  0.10E+02  254   814  296  517   33   4   0.48E+00   0.50E+01  0.00E+00    0.10E+01


 Final statistics for this run:
 rwork size =    3861   iwork size =    59
 number of time steps            =  254
 number of f evaluations         =  814
 number of preconditioner evals. =   33
 number of preconditioner solves = 1554
 number of nonlinear iterations  =  296
 number of linear iterations     =  517
 average subspace dimension  =  1.7466
    0 nonlinear conv. failures,   20 linear conv. failures


--------------------------------------------------------------------------------

 Solution with mf = 22

   t       nstep  nfe  nni  nli  npe  nq    h          avdim    ncf rate    lcf rate
  0.10E-07    3     5    3    1    2   2   0.65E-07   0.33E+00  0.00E+00    0.00E+00
  0.10E-06    4     8    5    2    2   2   0.65E-07   0.50E+00  0.00E+00    0.00E+00
  0.10E-05   10    22   13    8    4   3   0.37E-06   0.75E+00  0.00E+00    0.00E+00
  0.10E-04   33    78   43   34    6   5   0.51E-06   0.87E+00  0.00E+00    0.00E+00
  0.10E-03  112   270  134  135   14   5   0.60E-05   0.11E+01  0.00E+00    0.00E+00
  0.10E-02  131   321  156  164   18   2   0.35E-03   0.13E+01  0.00E+00    0.00E+00
  0.10E-01  138   343  167  175   20   3   0.18E-02   0.10E+01  0.00E+00    0.00E+00
  0.10E+00  162   402  194  207   23   4   0.71E-02   0.12E+01  0.00E+00    0.00E+00
  0.10E+01  206   543  242  300   27   4   0.47E-01   0.19E+01  0.00E+00    0.00E+00
  0.20E+01  221   612  259  352   29   4   0.95E-01   0.31E+01  0.00E+00    0.00E+00
  0.30E+01  230   661  269  391   30   4   0.16E+00   0.39E+01  0.00E+00    0.10E+00
  0.40E+01  236   695  275  419   30   4   0.20E+00   0.47E+01  0.00E+00    0.67E+00
  0.50E+01  241   731  281  449   31   4   0.27E+00   0.50E+01  0.00E+00    0.10E+01
  0.60E+01  244   749  284  464   31   4   0.27E+00   0.50E+01  0.00E+00    0.10E+01
  0.70E+01  247   773  288  484   32   4   0.36E+00   0.50E+01  0.00E+00    0.10E+01
  0.80E+01  250   797  292  504   33   3   0.51E+00   0.50E+01  0.00E+00    0.10E+01
  0.90E+01  252   809  294  514   33   3   0.51E+00   0.50E+01  0.00E+00    0.10E+01
  0.10E+02  254   827  297  529   34   2   0.82E+00   0.50E+01  0.00E+00    0.10E+01


 Final statistics for this run:
 rwork size =    3877   iwork size =    54
 number of time steps            =  254
 number of f evaluations         =  827
 number of preconditioner evals. =   34
 number of preconditioner solves = 1580
 number of nonlinear iterations  =  297
 number of linear iterations     =  529
 average subspace dimension  =  1.7811
    0 nonlinear conv. failures,   27 linear conv. failures


--------------------------------------------------------------------------------

 Solution with mf = 23

   t       nstep  nfe  nni  nli  npe  nq    h          avdim    ncf rate    lcf rate
  0.10E-07    3     5    3    1    2   2   0.65E-07   0.33E+00  0.00E+00    0.00E+00
  0.10E-06    4     8    5    2    2   2   0.65E-07   0.50E+00  0.00E+00    0.00E+00
  0.10E-05   10    22   13    8    4   3   0.37E-06   0.75E+00  0.00E+00    0.00E+00
  0.10E-04   33    78   43   34    6   5   0.51E-06   0.87E+00  0.00E+00    0.00E+00
  0.10E-03  112   272  134  137   14   5   0.60E-05   0.11E+01  0.00E+00    0.00E+00
  0.10E-02  131   324  156  167   18   2   0.35E-03   0.14E+01  0.00E+00    0.00E+00
  0.10E-01  139   360  166  193   20   3   0.15E-02   0.26E+01  0.00E+00    0.40E+00
  0.10E+00  172   503  209  293   28   4   0.53E-02   0.23E+01  0.91E-01    0.26E+00
  0.10E+01  227   845  278  566   42   4   0.11E-01   0.40E+01  0.91E-01    0.52E+00
  0.20E+01  280  1287  360  926   70   2   0.17E-01   0.44E+01  0.25E+00    0.76E+00


 Final statistics for this run:
 rwork size =    3404   iwork size =    54
 number of time steps            =  280
 number of f evaluations         = 1287
 number of preconditioner evals. =   70
 number of preconditioner solves =  926
 number of nonlinear iterations  =  360
 number of linear iterations     =  926
 average subspace dimension  =  2.5722
   21 nonlinear conv. failures,  113 linear conv. failures


--------------------------------------------------------------------------------

 Solution with mf = 24

   t       nstep  nfe  nni  nli  npe  nq    h          avdim    ncf rate    lcf rate
  0.10E-07    3     5    3    1    2   2   0.65E-07   0.33E+00  0.00E+00    0.00E+00
  0.10E-06    4     8    5    2    2   2   0.65E-07   0.50E+00  0.00E+00    0.00E+00
  0.10E-05   10    22   13    8    4   3   0.37E-06   0.75E+00  0.00E+00    0.00E+00
  0.10E-04   33    78   43   34    6   5   0.51E-06   0.87E+00  0.00E+00    0.00E+00
  0.10E-03  112   270  134  135   14   5   0.60E-05   0.11E+01  0.00E+00    0.00E+00
  0.10E-02  131   321  156  164   18   2   0.35E-03   0.13E+01  0.00E+00    0.00E+00
  0.10E-01  139   356  167  188   20   3   0.16E-02   0.22E+01  0.00E+00    0.27E+00
  0.10E+00  162   412  193  218   23   4   0.71E-02   0.12E+01  0.00E+00    0.00E+00
  0.10E+01  223   780  274  505   38   4   0.23E-01   0.35E+01  0.82E-01    0.43E+00
  0.20E+01  263  1085  335  749   59   3   0.17E-01   0.40E+01  0.23E+00    0.56E+00


 Final statistics for this run:
 rwork size =    3404   iwork size =    54
 number of time steps            =  263
 number of f evaluations         = 1085
 number of preconditioner evals. =   59
 number of preconditioner solves =  749
 number of nonlinear iterations  =  335
 number of linear iterations     =  749
 average subspace dimension  =  2.2358
   14 nonlinear conv. failures,   72 linear conv. failures


--------------------------------------------------------------------------------

 Solution with mf = 29

   t       nstep  nfe  nni  nli  npe  nq    h          avdim    ncf rate    lcf rate
  0.10E-07    3     4    3    0    2   2   0.65E-07   0.00E+00  0.00E+00    0.00E+00
  0.10E-06    4     6    5    0    2   2   0.65E-07   0.00E+00  0.00E+00    0.00E+00
  0.10E-05   10    14   13    0    4   3   0.37E-06   0.00E+00  0.00E+00    0.00E+00
  0.10E-04   32    42   41    0    6   5   0.56E-06   0.00E+00  0.00E+00    0.00E+00
  0.10E-03  114   135  134    0   13   5   0.47E-05   0.00E+00  0.00E+00    0.00E+00
  0.10E-02  136   162  161    0   18   2   0.30E-03   0.00E+00  0.00E+00    0.00E+00
  0.10E-01  144   173  172    0   20   3   0.16E-02   0.00E+00  0.00E+00    0.00E+00
  0.10E+00  168   200  199    0   23   4   0.79E-02   0.00E+00  0.00E+00    0.00E+00
  0.10E+01  245   353  352    0   39   2   0.17E-01   0.00E+00  0.13E-01    0.00E+00
  0.20E+01  293   467  466    0   57   2   0.28E-01   0.00E+00  0.10E+00    0.00E+00
  0.30E+01  330   566  565    0   76   2   0.37E-01   0.00E+00  0.22E+00    0.00E+00
  0.40E+01  356   631  630    0   87   2   0.31E-01   0.00E+00  0.15E+00    0.00E+00
  0.50E+01  384   697  696    0   98   1   0.72E-01   0.00E+00  0.14E+00    0.00E+00
  0.60E+01  399   742  741    0  109   2   0.21E+00   0.00E+00  0.27E+00    0.00E+00
  0.70E+01  411   783  782    0  117   1   0.20E+00   0.00E+00  0.33E+00    0.00E+00
  0.80E+01  414   788  787    0  118   2   0.41E+00   0.00E+00  0.00E+00    0.00E+00
  0.90E+01  416   791  790    0  118   2   0.41E+00   0.00E+00  0.00E+00    0.00E+00
  0.10E+02  418   793  792    0  119   3   0.74E+00   0.00E+00  0.00E+00    0.00E+00


 Final statistics for this run:
 rwork size =    2756   iwork size =    54
 number of time steps            =  418
 number of f evaluations         =  793
 number of preconditioner evals. =  119
 number of preconditioner solves =  777
 number of nonlinear iterations  =  792
 number of linear iterations     =    0
 average subspace dimension  =  0.0000
   30 nonlinear conv. failures,    0 linear conv. failures

................................................................................


--------------------------------------------------------------------------------
                              At time t =   0.00000000E+00
--------------------------------------------------------------------------------
 the species c( 1) values are:
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
  10.0000      10.1678      10.3775      10.3775      10.1678      10.0000    
  10.0000      10.3775      10.8493      10.8493      10.3775      10.0000    
  10.0000      10.3775      10.8493      10.8493      10.3775      10.0000    
  10.0000      10.1678      10.3775      10.3775      10.1678      10.0000    
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
  10.0000      10.3355      10.7550      10.7550      10.3355      10.0000    
  10.0000      10.7550      11.6987      11.6987      10.7550      10.0000    
  10.0000      10.7550      11.6987      11.6987      10.7550      10.0000    
  10.0000      10.3355      10.7550      10.7550      10.3355      10.0000    
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
  10.0000      10.5033      11.1325      11.1325      10.5033      10.0000    
  10.0000      11.1325      12.5480      12.5480      11.1325      10.0000    
  10.0000      11.1325      12.5480      12.5480      11.1325      10.0000    
  10.0000      10.5033      11.1325      11.1325      10.5033      10.0000    
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
  10.0000      10.6711      11.5099      11.5099      10.6711      10.0000    
  10.0000      11.5099      13.3974      13.3974      11.5099      10.0000    
  10.0000      11.5099      13.3974      13.3974      11.5099      10.0000    
  10.0000      10.6711      11.5099      11.5099      10.6711      10.0000    
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
  10.0000      10.8389      11.8874      11.8874      10.8389      10.0000    
  10.0000      11.8874      14.2467      14.2467      11.8874      10.0000    
  10.0000      11.8874      14.2467      14.2467      11.8874      10.0000    
  10.0000      10.8389      11.8874      11.8874      10.8389      10.0000    
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
  10.0000      11.0066      12.2649      12.2649      11.0066      10.0000    
  10.0000      12.2649      15.0961      15.0961      12.2649      10.0000    
  10.0000      12.2649      15.0961      15.0961      12.2649      10.0000    
  10.0000      11.0066      12.2649      12.2649      11.0066      10.0000    
  10.0000      10.0000      10.0000      10.0000      10.0000      10.0000    
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
                              At time t =   0.10000000E-05
--------------------------------------------------------------------------------
 the species c( 1) values are:
  9.99991      9.99992      9.99993      9.99993      9.99993      9.99992    
  9.99992      10.1677      10.3774      10.3774      10.1677      9.99993    
  9.99993      10.3774      10.8492      10.8492      10.3774      9.99993    
  9.99993      10.3774      10.8492      10.8492      10.3774      9.99993    
  9.99992      10.1677      10.3774      10.3774      10.1677      9.99992    
  9.99991      9.99992      9.99993      9.99993      9.99992      9.99991    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  9.99991      9.99993      9.99995      9.99995      9.99993      9.99992    
  9.99993      10.3355      10.7549      10.7549      10.3355      9.99993    
  9.99995      10.7549      11.6985      11.6985      10.7549      9.99995    
  9.99995      10.7549      11.6985      11.6985      10.7549      9.99995    
  9.99993      10.3355      10.7549      10.7549      10.3355      9.99993    
  9.99991      9.99993      9.99995      9.99995      9.99993      9.99991    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  9.99991      9.99994      9.99997      9.99997      9.99994      9.99992    
  9.99994      10.5032      11.1323      11.1323      10.5032      9.99994    
  9.99997      11.1323      12.5478      12.5478      11.1323      9.99997    
  9.99997      11.1323      12.5478      12.5478      11.1323      9.99997    
  9.99994      10.5032      11.1323      11.1323      10.5032      9.99994    
  9.99991      9.99994      9.99997      9.99997      9.99994      9.99991    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  13.4987      13.4987      13.4987      13.4987      13.4987      13.4987    
  13.4987      14.5503      15.8929      15.8929      14.5503      13.4987    
  13.4987      15.8929      19.0303      19.0303      15.8929      13.4987    
  13.4987      15.8929      19.0303      19.0303      15.8929      13.4987    
  13.4987      14.5503      15.8929      15.8929      14.5503      13.4987    
  13.4987      13.4987      13.4987      13.4987      13.4987      13.4987    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  13.4987      13.4987      13.4987      13.4987      13.4987      13.4987    
  13.4987      14.7791      16.4141      16.4141      14.7791      13.4987    
  13.4987      16.4141      20.2367      20.2367      16.4141      13.4987    
  13.4987      16.4141      20.2367      20.2367      16.4141      13.4987    
  13.4987      14.7791      16.4141      16.4141      14.7791      13.4987    
  13.4987      13.4987      13.4987      13.4987      13.4987      13.4987    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  13.4987      13.4987      13.4988      13.4987      13.4987      13.4987    
  13.4987      15.0078      16.9353      16.9353      15.0078      13.4987    
  13.4988      16.9353      21.4431      21.4431      16.9353      13.4987    
  13.4988      16.9353      21.4431      21.4431      16.9353      13.4988    
  13.4987      15.0078      16.9353      16.9353      15.0078      13.4987    
  13.4987      13.4987      13.4988      13.4988      13.4987      13.4987    
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
                              At time t =   0.10000000E-02
--------------------------------------------------------------------------------
 the species c( 1) values are:
  9.90702      9.91664      9.92836      9.93033      9.92253      9.91674    
  9.91472      10.0746      10.2769      10.2785      10.0795      9.92253    
  9.92446      10.2748      10.7181      10.7194      10.2785      9.93033    
  9.92445      10.2744      10.7173      10.7181      10.2769      9.92836    
  9.91469      10.0734      10.2744      10.2748      10.0746      9.91664    
  9.90697      9.91469      9.92445      9.92446      9.91472      9.90702    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  9.90741      9.92474      9.94623      9.94820      9.93064      9.91713    
  9.92282      10.2412      10.6440      10.6457      10.2461      9.93064    
  9.94232      10.6419      11.5267      11.5281      10.6457      9.94820    
  9.94231      10.6415      11.5258      11.5267      10.6440      9.94623    
  9.92279      10.2400      10.6415      10.6419      10.2412      9.92474    
  9.90736      9.92279      9.94231      9.94232      9.92282      9.90741    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  9.90780      9.93284      9.96409      9.96606      9.93874      9.91752    
  9.93092      10.4078      11.0109      11.0127      10.4127      9.93874    
  9.96018      11.0088      12.3339      12.3354      11.0127      9.96606    
  9.96017      11.0083      12.3329      12.3339      11.0109      9.96409    
  9.93089      10.4065      11.0083      11.0088      10.4078      9.93284    
  9.90776      9.93089      9.96017      9.96018      9.93092      9.90780    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  297231.      297749.      298393.      298451.      297925.      297520.    
  297692.      307244.      319327.      319378.      307390.      297925.    
  298276.      319264.      345799.      345840.      319378.      298451.    
  298276.      319252.      345771.      345799.      319327.      298393.    
  297691.      307208.      319252.      319264.      307244.      297749.    
  297229.      297691.      298276.      298276.      297692.      297231.    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  297231.      297749.      298393.      298451.      297925.      297520.    
  297692.      307244.      319327.      319378.      307390.      297925.    
  298276.      319264.      345799.      345840.      319378.      298451.    
  298276.      319252.      345771.      345799.      319327.      298393.    
  297691.      307208.      319252.      319264.      307244.      297749.    
  297229.      297691.      298276.      298276.      297692.      297231.    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  297231.      297749.      298393.      298451.      297925.      297520.    
  297692.      307244.      319327.      319378.      307390.      297925.    
  298276.      319264.      345799.      345840.      319378.      298451.    
  298276.      319252.      345771.      345799.      319327.      298393.    
  297691.      307208.      319252.      319264.      307244.      297749.    
  297229.      297691.      298276.      298276.      297692.      297231.    
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
                              At time t =   0.10000000E+01
--------------------------------------------------------------------------------
 the species c( 1) values are:
  1.58846      1.59918      1.62146      1.64759      1.67030      1.68143    
  1.58527      1.59498      1.61542      1.63946      1.66027      1.67030    
  1.57751      1.58542      1.60234      1.62229      1.63946      1.64759    
  1.56815      1.57406      1.58700      1.60234      1.61542      1.62146    
  1.56043      1.56457      1.57406      1.58542      1.59498      1.59918    
  1.55727      1.56043      1.56815      1.57751      1.58527      1.58846    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  1.59061      1.60135      1.62365      1.64981      1.67255      1.68369    
  1.58742      1.59714      1.61761      1.64167      1.66251      1.67255    
  1.57965      1.58757      1.60451      1.62449      1.64167      1.64981    
  1.57028      1.57620      1.58916      1.60451      1.61761      1.62365    
  1.56255      1.56670      1.57620      1.58757      1.59714      1.60135    
  1.55939      1.56255      1.57028      1.57965      1.58742      1.59061    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  1.59265      1.60340      1.62572      1.65191      1.67468      1.68583    
  1.58946      1.59918      1.61967      1.64377      1.66462      1.67468    
  1.58168      1.58960      1.60656      1.62656      1.64377      1.65191    
  1.57230      1.57823      1.59119      1.60656      1.61967      1.62572    
  1.56456      1.56872      1.57823      1.58960      1.59918      1.60340    
  1.56140      1.56456      1.57230      1.58168      1.58946      1.59265    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  47716.6      48038.6      48707.3      49491.7      50173.7      50507.6    
  47621.0      47912.3      48526.1      49247.9      49872.6      50173.7    
  47388.0      47625.3      48133.2      48732.4      49247.9      49491.7    
  47106.8      47284.3      47672.8      48133.2      48526.1      48707.3    
  46874.9      46999.4      47284.3      47625.3      47912.3      48038.6    
  46780.1      46874.9      47106.8      47388.0      47621.0      47716.6    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  47716.6      48038.6      48707.3      49491.7      50173.7      50507.6    
  47621.0      47912.3      48526.1      49247.9      49872.6      50173.7    
  47388.0      47625.3      48133.2      48732.4      49247.9      49491.7    
  47106.8      47284.3      47672.8      48133.2      48526.1      48707.3    
  46874.9      46999.4      47284.3      47625.3      47912.3      48038.6    
  46780.1      46874.9      47106.8      47388.0      47621.0      47716.6    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  47716.6      48038.6      48707.3      49491.7      50173.7      50507.6    
  47621.0      47912.3      48526.1      49247.9      49872.6      50173.7    
  47388.0      47625.3      48133.2      48732.4      49247.9      49491.7    
  47106.8      47284.3      47672.8      48133.2      48526.1      48707.3    
  46874.9      46999.4      47284.3      47625.3      47912.3      48038.6    
  46780.1      46874.9      47106.8      47388.0      47621.0      47716.6    
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
                              At time t =   0.40000000E+01
--------------------------------------------------------------------------------
 the species c( 1) values are:
  1.19534      1.20367      1.22108      1.24156      1.25934      1.26799    
  1.19279      1.20034      1.21635      1.23521      1.25152      1.25934    
  1.18656      1.19272      1.20601      1.22172      1.23521      1.24156    
  1.17903      1.18367      1.19388      1.20601      1.21635      1.22108    
  1.17283      1.17611      1.18367      1.19272      1.20034      1.20367    
  1.17031      1.17283      1.17903      1.18656      1.19279      1.19534    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  1.19537      1.20370      1.22112      1.24159      1.25937      1.26802    
  1.19282      1.20037      1.21638      1.23525      1.25156      1.25937    
  1.18659      1.19276      1.20605      1.22175      1.23525      1.24159    
  1.17906      1.18370      1.19391      1.20605      1.21638      1.22112    
  1.17287      1.17615      1.18370      1.19276      1.20037      1.20370    
  1.17034      1.17287      1.17906      1.18659      1.19282      1.19537    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  1.19540      1.20373      1.22115      1.24163      1.25940      1.26806    
  1.19286      1.20040      1.21641      1.23528      1.25159      1.25940    
  1.18662      1.19279      1.20608      1.22179      1.23528      1.24163    
  1.17910      1.18373      1.19395      1.20608      1.21641      1.22115    
  1.17290      1.17618      1.18373      1.19279      1.20040      1.20373    
  1.17037      1.17290      1.17910      1.18662      1.19286      1.19540    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  35860.3      36109.8      36632.0      37246.1      37779.1      38038.3    
  35783.8      36010.0      36490.0      37055.9      37544.9      37779.1    
  35596.8      35781.6      36180.1      36651.2      37055.9      37246.1    
  35371.0      35510.0      35816.3      36180.1      36490.0      36632.0    
  35185.1      35283.4      35510.0      35781.6      36010.0      36109.8    
  35109.4      35185.1      35371.0      35596.8      35783.8      35860.3    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  35860.3      36109.8      36632.0      37246.1      37779.1      38038.3    
  35783.8      36010.0      36490.0      37055.9      37544.9      37779.1    
  35596.8      35781.6      36180.1      36651.2      37055.9      37246.1    
  35371.0      35510.0      35816.3      36180.1      36490.0      36632.0    
  35185.1      35283.4      35510.0      35781.6      36010.0      36109.8    
  35109.4      35185.1      35371.0      35596.8      35783.8      35860.3    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  35860.3      36109.8      36632.0      37246.1      37779.1      38038.3    
  35783.8      36010.0      36490.0      37055.9      37544.9      37779.1    
  35596.8      35781.6      36180.1      36651.2      37055.9      37246.1    
  35371.0      35510.0      35816.3      36180.1      36490.0      36632.0    
  35185.1      35283.4      35510.0      35781.6      36010.0      36109.8    
  35109.4      35185.1      35371.0      35596.8      35783.8      35860.3    
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
                              At time t =   0.70000000E+01
--------------------------------------------------------------------------------
 the species c( 1) values are:
  1.18854      1.19682      1.21415      1.23453      1.25221      1.26082    
  1.18600      1.19351      1.20944      1.22821      1.24444      1.25221    
  1.17980      1.18593      1.19916      1.21479      1.22821      1.23453    
  1.17231      1.17692      1.18708      1.19916      1.20944      1.21415    
  1.16614      1.16940      1.17692      1.18593      1.19351      1.19682    
  1.16363      1.16614      1.17231      1.17980      1.18600      1.18854    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  1.18854      1.19682      1.21415      1.23453      1.25221      1.26082    
  1.18600      1.19351      1.20944      1.22821      1.24444      1.25221    
  1.17980      1.18593      1.19916      1.21479      1.22821      1.23453    
  1.17231      1.17692      1.18708      1.19916      1.20944      1.21415    
  1.16614      1.16940      1.17692      1.18593      1.19351      1.19682    
  1.16363      1.16614      1.17231      1.17980      1.18600      1.18854    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  1.18854      1.19682      1.21415      1.23453      1.25221      1.26082    
  1.18600      1.19351      1.20944      1.22821      1.24444      1.25221    
  1.17980      1.18593      1.19916      1.21479      1.22821      1.23453    
  1.17231      1.17692      1.18709      1.19916      1.20944      1.21415    
  1.16614      1.16940      1.17692      1.18593      1.19351      1.19682    
  1.16363      1.16614      1.17231      1.17980      1.18600      1.18854    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  35655.3      35903.5      36423.1      37034.0      37564.3      37822.2    
  35579.2      35804.2      36281.9      36844.8      37331.4      37564.3    
  35393.1      35576.9      35973.5      36442.2      36844.8      37034.0    
  35168.3      35306.6      35611.4      35973.5      36281.9      36423.1    
  34983.2      35081.1      35306.6      35576.9      35804.2      35903.5    
  34907.9      34983.2      35168.3      35393.1      35579.2      35655.3    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  35655.3      35903.5      36423.1      37034.0      37564.3      37822.2    
  35579.2      35804.2      36281.9      36844.8      37331.4      37564.3    
  35393.1      35576.9      35973.5      36442.2      36844.8      37034.0    
  35168.3      35306.6      35611.4      35973.5      36281.9      36423.1    
  34983.2      35081.1      35306.6      35576.9      35804.2      35903.5    
  34907.9      34983.2      35168.3      35393.1      35579.2      35655.3    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  35655.3      35903.5      36423.1      37034.0      37564.3      37822.2    
  35579.2      35804.2      36281.9      36844.8      37331.4      37564.3    
  35393.1      35576.9      35973.5      36442.2      36844.8      37034.0    
  35168.3      35306.6      35611.4      35973.5      36281.9      36423.1    
  34983.2      35081.1      35306.6      35576.9      35804.2      35903.5    
  34907.9      34983.2      35168.3      35393.1      35579.2      35655.3    
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
                              At time t =   0.10000000E+02
--------------------------------------------------------------------------------
 the species c( 1) values are:
  1.18838      1.19666      1.21399      1.23436      1.25205      1.26065    
  1.18584      1.19335      1.20928      1.22805      1.24428      1.25205    
  1.17964      1.18577      1.19899      1.21463      1.22805      1.23436    
  1.17215      1.17676      1.18692      1.19899      1.20928      1.21399    
  1.16598      1.16924      1.17676      1.18577      1.19335      1.19666    
  1.16347      1.16598      1.17215      1.17964      1.18584      1.18838    
--------------------------------------------------------------------------------

 the species c( 2) values are:
  1.18838      1.19666      1.21399      1.23436      1.25205      1.26065    
  1.18584      1.19335      1.20928      1.22805      1.24428      1.25205    
  1.17964      1.18577      1.19899      1.21462      1.22805      1.23436    
  1.17215      1.17676      1.18692      1.19899      1.20928      1.21399    
  1.16598      1.16924      1.17676      1.18577      1.19335      1.19666    
  1.16347      1.16598      1.17215      1.17964      1.18584      1.18838    
--------------------------------------------------------------------------------

 the species c( 3) values are:
  1.18838      1.19666      1.21399      1.23436      1.25205      1.26065    
  1.18584      1.19335      1.20928      1.22805      1.24428      1.25205    
  1.17964      1.18577      1.19899      1.21462      1.22805      1.23436    
  1.17215      1.17676      1.18692      1.19899      1.20928      1.21399    
  1.16598      1.16924      1.17676      1.18577      1.19335      1.19666    
  1.16347      1.16598      1.17215      1.17964      1.18584      1.18838    
--------------------------------------------------------------------------------

 the species c( 4) values are:
  35650.5      35898.7      36418.2      37029.1      37559.4      37817.2    
  35574.4      35799.4      36277.0      36839.9      37326.5      37559.4    
  35388.3      35572.1      35968.6      36437.3      36839.9      37029.1    
  35163.6      35301.8      35606.6      35968.6      36277.0      36418.2    
  34978.5      35076.4      35301.8      35572.1      35799.4      35898.7    
  34903.1      34978.5      35163.6      35388.3      35574.4      35650.5    
--------------------------------------------------------------------------------

 the species c( 5) values are:
  35650.5      35898.7      36418.2      37029.1      37559.4      37817.2    
  35574.4      35799.4      36277.0      36839.9      37326.5      37559.4    
  35388.3      35572.1      35968.6      36437.3      36839.9      37029.1    
  35163.6      35301.8      35606.6      35968.6      36277.0      36418.2    
  34978.5      35076.4      35301.8      35572.1      35799.4      35898.7    
  34903.1      34978.5      35163.6      35388.3      35574.4      35650.5    
--------------------------------------------------------------------------------

 the species c( 6) values are:
  35650.5      35898.7      36418.2      37029.1      37559.4      37817.2    
  35574.4      35799.4      36277.0      36839.9      37326.5      37559.4    
  35388.3      35572.1      35968.6      36437.3      36839.9      37029.1    
  35163.6      35301.8      35606.6      35968.6      36277.0      36418.2    
  34978.5      35076.4      35301.8      35572.1      35799.4      35898.7    
  34903.1      34978.5      35163.6      35388.3      35574.4      35650.5    
--------------------------------------------------------------------------------


==========Source and Output for LSODKR Demonstration Program====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODKR package.
c This is the version of 27 April 2005.
c
c This version is in double precision.
c
c An ODE system is generated from the following 2-species diurnal
c kinetics advection-diffusion PDE system in 2 space dimensions:
c
c dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
c                 + Ri(c1,c2,t)      for i = 1,2,   where
c   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
c   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
c   Kv(z) = Kv0*exp(z/5) ,
c Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
c vary diurnally.  The species are oxygen singlet and ozone.
c The problem is posed on the square
c   0 .le. x .le. 20,    30 .le. z .le. 50   (all in km),
c with homogeneous Neumann boundary conditions, and for time t in
c   0 .le. t .le. 86400 sec (1 day).
c
c The PDE system is treated by central differences on a uniform
c 10 x 10 mesh, with simple polynomial initial profiles.
c
c The problem is solved with DLSODKR, with the BDF/GMRES method and
c the block-diagonal part of the Jacobian as a left preconditioner.
c At intervals of 7200 sec (2 hrs), output includes sample values
c of c1, c2, and c2tot = total of all c2 values.
c
c Roots of the function g = d(c2tot)/dt are found, i.e. the points at
c which the total c2 (ozone) is stationary.
c
c Note: The preconditioner routines call LINPACK routines DGEFA, DGESL,
c and the BLAS routines DCOPY and DSCAL.
c-----------------------------------------------------------------------
      external fdem, jacbd, solbd, gdem
      integer mx, mz, mm, iout, istate, iwork, jacflg, jpre, jroot,
     1   jx, jz, leniw, lenrw, liw, lrw, mf, neq, nst, nsfi, nfe,
     2   nge, npe, njev, nps, nni, nli, ncfn, ncfl   
      double precision q1,q2,q3,q4, a3,a4, om, c3, dz, hdco, vdco, haco,
     1   dkh, vel, dkv0, halfda, pi, twohr, rtol, floor,
     2   dx, atol, t, tout, x, cx, z, cz, y, rwork, c2tot, avdim
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
      dimension y(2,10,10), rwork(4264), iwork(235)
      data dkh/4.0d-6/, vel/0.001d0/, dkv0/1.0d-8/, halfda/4.32d4/,
     1  pi/3.1415926535898d0/, twohr/7200.0d0/, rtol/1.0d-5/,
     2  floor/100.0d0/, lrw/4264/, liw/235/, mf/22/, jpre/1/, jacflg/1/
c
c Load Common block of problem parameters.
      mx = 10
      mz = 10
      mm = mx*mz
      q1 = 1.63d-16
      q2 = 4.66d-16
      a3 = 22.62d0
      a4 = 7.601d0
      om = pi/halfda
      c3 = 3.7d16
      dx = 20.0d0/(mx - 1.0d0)
      dz = 20.0d0/(mz - 1.0d0)
      hdco = dkh/dx**2
      haco = vel/(2.0d0*dx)
      vdco = (1.0d0/dz**2)*dkv0
c Set other input arguments.
      atol = rtol*floor
      neq = 2*mx*mz
      iwork(1) = 8*mx*mz
      iwork(2) = neq
      iwork(3) = jpre
      iwork(4) = jacflg
      t = 0.0d0
      tout = 0.0d0
      istate = 1
c Initialize values for tout = 0 output.
      iwork(11) = 0
      iwork(14) = 0
      rwork(11) = 0.0d0
c Set initial profiles.
      do 20 jz = 1,mz
        z = 30.0d0 + (jz - 1.0d0)*dz
        cz = (0.1d0*(z - 40.0d0))**2
        cz = 1.0d0 - cz + 0.5d0*cz**2
        do 10 jx = 1,mx
          x = (jx - 1.0d0)*dx
          cx = (0.1d0*(x - 10.0d0))**2
          cx = 1.0d0 - cx + 0.5d0*cx**2
          y(1,jx,jz) = 1.0d6*cx*cz
          y(2,jx,jz) = 1.0d12*cx*cz
 10       continue
 20     continue
c
c Write heading, problem parameters, solution parameters.
      write(6,30) mx, mz, mf, rtol, atol
 30   format('Demonstration program for DLSODKR package'//
     1       '2D diurnal kinetics-transport PDE system with 2 species'/
     2       'Spatial mesh is',i3,' by',i3/'Method flag is mf =',i3,
     3       '   Tolerances are rtol =',d8.1,'   atol =',d8.1/
     4       'Left preconditioner uses block-diagonal part of Jacobian'/
     5       'Root function finds stationary points of total ozone,'/
     6       '  i.e. roots of (d/dt)(sum of c2 over all mesh points)'/)
c
c Loop over output points, call DLSODKR, print sample solution values.
      do 70 iout = 1,13
 40     call dlsodkr (fdem, neq, y, t, tout, 1, rtol, atol, 1, istate,
     1      0, rwork, lrw, iwork, liw, jacbd, solbd, mf, gdem, 1, jroot)
        write(6,50) t,iwork(11),iwork(14),rwork(11)
 50     format(/' t =',d10.3,4x,'no. steps =',i5,
     1          '   order =',i2,'   stepsize =',d10.3)
        call c2sum (y, mx, mz, c2tot)
        write(6,60) y(1,1,1), y(1,5,5), y(1,10,10),
     1              y(2,1,1), y(2,5,5), y(2,10,10)
 60     format('   c1 (bot.left/middle/top rt.) =',3d12.3/
     1         '   c2 (bot.left/middle/top rt.) =',3d12.3)
        write(6,62)c2tot,jroot
 62     format('   total c2 =',d15.6,
     1         '   jroot =',i2' (1 = root found, 0 = no root)')
        if (istate .lt. 0) then
          write(6,65)istate
 65       format('DLSODKR returned istate = ',i3)
          go to 75
        endif
        if (istate .eq. 3) then
          istate = 2
          go to 40
        endif
        tout = tout + twohr
 70     continue
c
c Print final statistics.
 75   lenrw = iwork(17)
      leniw = iwork(18)
      nst = iwork(11)
      nsfi = iwork(24)
      nfe = iwork(12)
      nge = iwork(10)
      npe = iwork(13)
      njev = iwork(25)
      nps = iwork(21)
      nni = iwork(19)
      nli = iwork(20)
      avdim = real(nli)/real(nni)
      ncfn = iwork(22)
      ncfl = iwork(23)
      write (6,80) lenrw,leniw,nst,nsfi,nfe,nge,npe,njev,nps,nni,nli,
     1             avdim,ncfn,ncfl
 80   format(//' Final statistics:'/
     1 ' rwork size =',i5,5x,' iwork size =',i4/
     2 ' number of steps        =',i5,5x,'no. fnal. iter. steps  =',i5/
     3 ' number of f evals.     =',i5,5x,'number of g evals.     =',i5/
     4 ' number of prec. evals. =',i5,5x,'number of Jac. evals.  =',i5/
     5 ' number of prec. solves =',i5/
     6 ' number of nonl. iters. =',i5,5x,'number of lin. iters.  =',i5/
     7 ' average Krylov subspace dimension (nli/nni)  =',f8.4/
     8 ' number of conv. failures:  nonlinear =',i3,'  linear =',i3)
      stop
      end

      subroutine fdem (neq, t, y, ydot)
      integer neq,  mx, mz, mm,
     1   iblok, iblok0, idn, iup, ileft, iright, jx, jz
      double precision t, y(2,*), ydot(2,*),
     1   q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco,
     2   c1, c1dn, c1up, c1lt, c1rt, c2, c2dn, c2up, c2lt, c2rt, 
     3   czdn, czup, horad1, hord1, horad2, hord2, qq1, qq2, qq3, qq4, 
     4   rkin1, rkin2, s, vertd1, vertd2, zdn, zup
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c
c Set diurnal rate coefficients.
      s = sin(om*t)
      if (s .gt. 0.0d0) then
        q3 = exp(-a3/s)
        q4 = exp(-a4/s)
      else
        q3 = 0.0d0
        q4 = 0.0d0
      endif
c Loop over all grid points.
      do 20 jz = 1,mz
        zdn = 30.0d0 + (jz - 1.5d0)*dz
        zup = zdn + dz
        czdn = vdco*exp(0.2d0*zdn)
        czup = vdco*exp(0.2d0*zup)
        iblok0 = (jz-1)*mx
        idn = -mx
        if (jz .eq. 1) idn = mx
        iup = mx
        if (jz .eq. mz) iup = -mx
        do 10 jx = 1,mx
          iblok = iblok0 + jx
          c1 = y(1,iblok)
          c2 = y(2,iblok)
c Set kinetic rate terms.
          qq1 = q1*c1*c3
          qq2 = q2*c1*c2
          qq3 = q3*c3
          qq4 = q4*c2
          rkin1 = -qq1 - qq2 + 2.0d0*qq3 + qq4
          rkin2 = qq1 - qq2 - qq4
c Set vertical diffusion terms.
          c1dn = y(1,iblok+idn)
          c2dn = y(2,iblok+idn)
          c1up = y(1,iblok+iup)
          c2up = y(2,iblok+iup)
          vertd1 = czup*(c1up - c1) - czdn*(c1 - c1dn)
          vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn)
c Set horizontal diffusion and advection terms.
          ileft = -1
          if (jx .eq. 1) ileft = 1
          iright = 1
          if (jx .eq. mx) iright = -1
          c1lt = y(1,iblok+ileft)
          c2lt = y(2,iblok+ileft)
          c1rt = y(1,iblok+iright)
          c2rt = y(2,iblok+iright)
          hord1 = hdco*(c1rt - 2.0d0*c1 + c1lt)
          hord2 = hdco*(c2rt - 2.0d0*c2 + c2lt)
          horad1 = haco*(c1rt - c1lt)
          horad2 = haco*(c2rt - c2lt)
c Load all terms into ydot.
          ydot(1,iblok) = vertd1 + hord1 + horad1 + rkin1
          ydot(2,iblok) = vertd2 + hord2 + horad2 + rkin2
 10       continue
 20     continue
      return
      end

      subroutine gdem (neq, t, y, ng, gout)
      integer neq, ng,  mx, mz, mm,
     1   iblok, iblok0, idn, iup, ileft, iright, jx, jz
      double precision t, y(2,*), gout,
     1   q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco,
     2   c1, c2, c2dn, c2up, c2lt, c2rt, c2dot, czdn, czup, horad2,
     3   hord2, qq1, qq2, qq4, rkin2, s, sum, vertd2, zdn, zup
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c
c This routine computes the rates for c2 and adds them.
c
c Set diurnal rate coefficient q4.
      s = sin(om*t)
      if (s .gt. 0.0d0) then
        q4 = exp(-a4/s)
      else
        q4 = 0.0d0
      endif
      sum = 0.0d0
c Loop over all grid points.
      do 20 jz = 1,mz
        zdn = 30.0d0 + (jz - 1.5d0)*dz
        zup = zdn + dz
        czdn = vdco*exp(0.2d0*zdn)
        czup = vdco*exp(0.2d0*zup)
        iblok0 = (jz-1)*mx
        idn = -mx
        if (jz .eq. 1) idn = mx
        iup = mx
        if (jz .eq. mz) iup = -mx
        do 10 jx = 1,mx
          iblok = iblok0 + jx
          c1 = y(1,iblok)
          c2 = y(2,iblok)
c Set kinetic rate term for c2.
          qq1 = q1*c1*c3
          qq2 = q2*c1*c2
          qq4 = q4*c2
          rkin2 = qq1 - qq2 - qq4
c Set vertical diffusion terms for c2.
          c2dn = y(2,iblok+idn)
          c2up = y(2,iblok+iup)
          vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn)
c Set horizontal diffusion and advection terms for c2.
          ileft = -1
          if (jx .eq. 1) ileft = 1
          iright = 1
          if (jx .eq. mx) iright = -1
          c2lt = y(2,iblok+ileft)
          c2rt = y(2,iblok+iright)
          hord2 = hdco*(c2rt - 2.0d0*c2 + c2lt)
          horad2 = haco*(c2rt - c2lt)
c Load all terms into c2dot and sum.
          c2dot = vertd2 + hord2 + horad2 + rkin2
          sum = sum + c2dot
 10       continue
 20     continue
      gout = sum
      return
      end

      subroutine jacbd (f, neq, t, y, ysv, rewt, f0, f1, hl0, jok,
     1    bd, ipbd, ier)
      external f
      integer neq, jok, ipbd(2,*), ier, mx, mz, mm,
     1    iblok, iblok0, jx, jz, lenbd
      double precision t, y(2,*), ysv(neq), rewt(neq), f0(neq), f1(neq),
     1   hl0, bd(2,2,*),
     2   q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco,
     3   c1, c2, czdn, czup, diag, temp, zdn, zup
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c
        lenbd = 4*mm
c If jok = 1, copy saved block-diagonal approximate Jacobian into bd.
      if (jok .eq. 1) then
        call dcopy (lenbd, bd(1,1,mm+1), 1, bd, 1)
        go to 30
        endif
c
c If jok = -1, compute and save diagonal Jacobian blocks
c  (using q3 and q4 values computed on last f call).
      do 20 jz = 1,mz
        zdn = 30.0d0 + (jz - 1.5d0)*dz
        zup = zdn + dz
        czdn = vdco*exp(0.2d0*zdn)
        czup = vdco*exp(0.2d0*zup)
        diag = -(czdn + czup + 2.0d0*hdco)
        iblok0 = (jz-1)*mx
        do 10 jx = 1,mx
          iblok = iblok0 + jx
          c1 = y(1,iblok)
          c2 = y(2,iblok)
          bd(1,1,iblok) = (-q1*c3 - q2*c2) + diag
          bd(1,2,iblok) = -q2*c1 + q4
          bd(2,1,iblok) = q1*c3 - q2*c2
          bd(2,2,iblok) = (-q2*c1 - q4) + diag
 10       continue
 20     continue
      call dcopy (lenbd, bd, 1, bd(1,1,mm+1), 1)
c Scale by -hl0, add identity matrix and LU-decompose blocks.
 30   temp = -hl0
      call dscal (lenbd, temp, bd, 1)
      do 40 iblok = 1,mm
        bd(1,1,iblok) = bd(1,1,iblok) + 1.0d0
        bd(2,2,iblok) = bd(2,2,iblok) + 1.0d0
        call dgefa (bd(1,1,iblok), 2, 2, ipbd(1,iblok), ier)
        if (ier .ne. 0) return
 40     continue
      return
      end

      subroutine solbd (neq, t, y, f0, wk, hl0, bd, ipbd, v, lr, ier)
      integer neq, ipbd(2,*), lr, ier,  mx, mz, mm,  i
      double precision t, y(neq), f0(neq), wk(neq), hl0, bd(2,2,*),
     2   v(2,*),  q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c Solve the block-diagonal system Px = v using LU factors stored in bd
c and pivot data in ipbd, and return the solution in v.
      ier = 0
      do 10 i = 1,mm
        call dgesl (bd(1,1,i), 2, 2, ipbd(1,i), v(1,i), 0)
 10     continue
      return
      end

      subroutine c2sum (y, mx, mz, c2tot)
      integer mx, mz, jx, jz
      double precision y(2,mx,mz), c2tot, sum
c Sum the c2 values.
      sum = 0.0d0
      do 20 jz = 1,mz
        do 20 jx = 1,mx
 20       sum = sum + y(2,jx,jz)
      c2tot = sum
      return
      end

................................................................................

Demonstration program for DLSODKR package

2D diurnal kinetics-transport PDE system with 2 species
Spatial mesh is 10 by 10
Method flag is mf = 22   Tolerances are rtol = 0.1E-04   atol = 0.1E-02
Left preconditioner uses block-diagonal part of Jacobian
Root function finds stationary points of total ozone,
  i.e. roots of (d/dt)(sum of c2 over all mesh points)


 t = 0.000E+00    no. steps =    0   order = 0   stepsize = 0.000E+00
   c1 (bot.left/middle/top rt.) =   0.250E+06   0.976E+06   0.250E+06
   c2 (bot.left/middle/top rt.) =   0.250E+12   0.976E+12   0.250E+12
   total c2 =   0.547546E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.184E+03    no. steps =  132   order = 2   stepsize = 0.138E+03
   c1 (bot.left/middle/top rt.) =   0.438E-07   0.171E-06   0.440E-07
   c2 (bot.left/middle/top rt.) =   0.250E+12   0.979E+12   0.251E+12
   total c2 =   0.547575E+14   jroot = 1 (1 = root found, 0 = no root)

 t = 0.720E+04    no. steps =  192   order = 5   stepsize = 0.158E+03
   c1 (bot.left/middle/top rt.) =   0.105E+05   0.296E+05   0.112E+05
   c2 (bot.left/middle/top rt.) =   0.253E+12   0.715E+12   0.270E+12
   total c2 =   0.503689E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.144E+05    no. steps =  222   order = 5   stepsize = 0.352E+03
   c1 (bot.left/middle/top rt.) =   0.666E+07   0.532E+07   0.730E+07
   c2 (bot.left/middle/top rt.) =   0.258E+12   0.206E+12   0.283E+12
   total c2 =   0.405024E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.183E+05    no. steps =  255   order = 5   stepsize = 0.270E+03
   c1 (bot.left/middle/top rt.) =   0.187E+08   0.551E+07   0.206E+08
   c2 (bot.left/middle/top rt.) =   0.269E+12   0.694E+11   0.297E+12
   total c2 =   0.379640E+14   jroot = 1 (1 = root found, 0 = no root)

 t = 0.216E+05    no. steps =  265   order = 5   stepsize = 0.409E+03
   c1 (bot.left/middle/top rt.) =   0.266E+08   0.104E+08   0.293E+08
   c2 (bot.left/middle/top rt.) =   0.299E+12   0.103E+12   0.331E+12
   total c2 =   0.397443E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.252E+05    no. steps =  274   order = 5   stepsize = 0.409E+03
   c1 (bot.left/middle/top rt.) =   0.218E+08   0.189E+08   0.241E+08
   c2 (bot.left/middle/top rt.) =   0.331E+12   0.284E+12   0.366E+12
   total c2 =   0.417179E+14   jroot = 1 (1 = root found, 0 = no root)

 t = 0.288E+05    no. steps =  295   order = 5   stepsize = 0.142E+03
   c1 (bot.left/middle/top rt.) =   0.870E+07   0.129E+08   0.965E+07
   c2 (bot.left/middle/top rt.) =   0.338E+12   0.503E+12   0.375E+12
   total c2 =   0.395883E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.360E+05    no. steps =  325   order = 5   stepsize = 0.119E+03
   c1 (bot.left/middle/top rt.) =   0.140E+05   0.203E+05   0.156E+05
   c2 (bot.left/middle/top rt.) =   0.339E+12   0.489E+12   0.377E+12
   total c2 =   0.303138E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.405E+05    no. steps =  371   order = 5   stepsize = 0.123E+03
   c1 (bot.left/middle/top rt.) =   0.116E-05   0.750E-06   0.130E-05
   c2 (bot.left/middle/top rt.) =   0.339E+12   0.218E+12   0.378E+12
   total c2 =   0.278867E+14   jroot = 1 (1 = root found, 0 = no root)

 t = 0.432E+05    no. steps =  379   order = 5   stepsize = 0.499E+03
   c1 (bot.left/middle/top rt.) =   0.181E-07   0.331E-06   0.211E-07
   c2 (bot.left/middle/top rt.) =   0.338E+12   0.136E+12   0.380E+12
   total c2 =   0.287486E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.504E+05    no. steps =  401   order = 5   stepsize = 0.483E+03
   c1 (bot.left/middle/top rt.) =   0.643E-08   0.423E-06   0.198E-07
   c2 (bot.left/middle/top rt.) =   0.336E+12   0.493E+12   0.386E+12
   total c2 =   0.361416E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.576E+05    no. steps =  414   order = 5   stepsize = 0.334E+03
   c1 (bot.left/middle/top rt.) =  -0.972E-11  -0.176E-09   0.200E-12
   c2 (bot.left/middle/top rt.) =   0.332E+12   0.965E+12   0.391E+12
   total c2 =   0.446354E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.648E+05    no. steps =  427   order = 5   stepsize = 0.789E+03
   c1 (bot.left/middle/top rt.) =  -0.987E-10   0.244E-08  -0.985E-10
   c2 (bot.left/middle/top rt.) =   0.331E+12   0.892E+12   0.396E+12
   total c2 =   0.492882E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.720E+05    no. steps =  436   order = 5   stepsize = 0.789E+03
   c1 (bot.left/middle/top rt.) =   0.187E-11   0.395E-10  -0.109E-12
   c2 (bot.left/middle/top rt.) =   0.333E+12   0.619E+12   0.404E+12
   total c2 =   0.529683E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.792E+05    no. steps =  445   order = 5   stepsize = 0.789E+03
   c1 (bot.left/middle/top rt.) =   0.794E-12  -0.147E-10   0.680E-12
   c2 (bot.left/middle/top rt.) =   0.333E+12   0.667E+12   0.412E+12
   total c2 =   0.605274E+14   jroot = 0 (1 = root found, 0 = no root)

 t = 0.824E+05    no. steps =  449   order = 5   stepsize = 0.789E+03
   c1 (bot.left/middle/top rt.) =  -0.178E-12   0.308E-11  -0.148E-12
   c2 (bot.left/middle/top rt.) =   0.334E+12   0.804E+12   0.415E+12
   total c2 =   0.619212E+14   jroot = 1 (1 = root found, 0 = no root)

 t = 0.864E+05    no. steps =  454   order = 5   stepsize = 0.789E+03
   c1 (bot.left/middle/top rt.) =   0.157E-14  -0.419E-12   0.104E-13
   c2 (bot.left/middle/top rt.) =   0.335E+12   0.911E+12   0.416E+12
   total c2 =   0.597970E+14   jroot = 0 (1 = root found, 0 = no root)


 Final statistics:
 rwork size = 4264      iwork size = 230
 number of steps        =  454     no. fnal. iter. steps  =  112
 number of f evals.     = 1547     number of g evals.     =  504
 number of prec. evals. =   52     number of Jac. evals.  =   11
 number of prec. solves = 1231
 number of nonl. iters. =  461     number of lin. iters.  =  847
 average Krylov subspace dimension (nli/nni)  =  1.8373
 number of conv. failures:  nonlinear =  0  linear =  0

==========Source and Output for LSODI Demonstration Program=====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODI package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
C this program solves a semi-discretized form of the Burgers equation,
c
c     u  = -(u*u/2)  + eta * u
c      t           x          xx
c
c for a = -1 .le. x .le. 1 = b, t .ge. 0.
c Here eta = 0.05.
c Boundary conditions: u(-1,t) = u(1,t) = 0.
c Initial profile: square wave
c     u(0,x) = 0    for 1/2 .lt. abs(x) .le. 1
c     u(0,x) = 1/2  for abs(x) = 1/2
c     u(0,x) = 1    for 0 .le. abs(x) .lt. 1/2
c
c An ODE system is generated by a simplified Galerkin treatment
c of the spatial variable x.
c
c Reference:
c R. C. Y. Chin, G. W. Hedstrom, and K. E. Karlsson,
c A Simplified Galerkin Method for Hyperbolic Equations,
c Math. Comp., vol. 33, no. 146 (April 1979), pp. 647-658.
c
c The problem is run with the DLSODI package with a 10-point mesh
c and a 100-point mesh.  In each case, it is run with two tolerances
c and for various appropriate values of the method flag mf.
c Output is on unit lout, set to 6 in a data statement below.
c-----------------------------------------------------------------------
      external res, addabd, addafl, jacbd, jacfl
      integer i, io, istate, itol, iwork, j,
     1   lout, liw, lrw, meth, miter, mf, ml, mu,
     2   n, nout, npts, nerr,
     3   nptsm1, n14, n34, n14m1, n14p1, n34m1, n34p1
      integer nm1
      double precision a, b, eta, delta,
     1   zero, fourth, half, one, hun,
     2   t, tout, tlast, tinit, errfac,
     3   atol, rtol, rwork, y, ydoti, elkup
      double precision eodsq, r4d
      dimension y(99), ydoti(99), tout(4), atol(2), rtol(2)
      dimension rwork(2002), iwork(125)
c Pass problem parameters in the Common block test1.
      common /test1/ r4d, eodsq, nm1
c
c Set problem parameters and run parameters
      data eta/0.05d0/, a/-1.0d0/, b/1.0d0/
      data zero/0.0d0/, fourth/0.25d0/, half/.5d0/, one/1.0d0/,
     1   hun/100.0d0/
      data tinit/0.0d0/, tlast/0.4d0/
      data tout/.10d0,.20d0,.30d0,.40d0/
      data ml/1/, mu/1/, lout/6/
      data nout/4/, lrw/2002/, liw/125/
      data itol/1/, rtol/1.0d-3, 1.0d-6/, atol/1.0d-3, 1.0d-6/
c
      iwork(1) = ml
      iwork(2) = mu
      nerr = 0
c
c Loop over two values of npts.
      do 300  npts = 10, 100, 90
c
c Compute the mesh width delta and other parameters.
      delta = (b - a)/npts
      r4d = fourth/delta
      eodsq = eta/delta**2
      nptsm1 = npts - 1
      n14 = npts/4
      n34 = 3 * n14
      n14m1 = n14 - 1
      n14p1 = n14m1 + 2
      n34m1 = n34 - 1
      n34p1 = n34m1 + 2
      n = nptsm1
      nm1 = n - 1
c
c Set the initial profile (for output purposes only).
c
      do 10 i = 1,n14m1
   10   y(i) = zero
      y(n14) = half
      do 20 i = n14p1,n34m1
   20   y(i) = one
      y(n34) = half
      do 30 i = n34p1,nptsm1
   30   y(i) = zero
c
      if (npts .gt. 10) write (lout,1010)
      write (lout,1000)
      write (lout,1100) eta,a,b,tinit,tlast,ml,mu,n
      write (lout,1200) zero, (y(i), i=1,n), zero
c
c The j loop is over error tolerances.
c
      do 200 j = 1,2
c
c Loop over method flag loop (for demonstration).
c
      do 100 meth = 1,2
       do 100 miter = 1,5
        if (miter .eq. 3)  go to 100
        if (miter .le. 2 .and. npts .gt. 10)  go to 100
        if (miter .eq. 5 .and. npts .lt. 100)  go to 100
        mf = 10*meth + miter
c
c Set the initial profile.
c
        do 40 i = 1,n14m1
   40     y(i) = zero
        y(n14) = half
        do 50 i = n14p1,n34m1
   50     y(i) = one
        y(n34) = half
        do 60 i = n34p1,nptsm1
   60     y(i) = zero
c
        t = tinit
        istate = 0
c
        write (lout,1500) rtol(j), atol(j), mf, npts
c
c  Output loop for each case
c
        do 80 io = 1,nout
c
c         call DLSODI
          if (miter .le. 2) call dlsodi (res, addafl, jacfl, n, y,
     1                      ydoti, t, tout(io), itol, rtol(j), atol(j),
     2                      1, istate, 0, rwork, lrw, iwork, liw, mf)
          if (miter .ge. 4) call dlsodi (res, addabd, jacbd, n, y,
     1                      ydoti, t, tout(io), itol, rtol(j), atol(j),
     2                      1, istate, 0, rwork, lrw, iwork, liw, mf)
          write (lout,2000) t, rwork(11), iwork(14),(y(i), i=1,n)
c
c If istate is not 2 on return, print message and loop.
          if (istate .ne. 2) then
            write (lout,4000) mf, t, istate
            nerr = nerr + 1
            go to 100
          endif
c
   80     continue
c
        write (lout,3000) mf, iwork(11), iwork(12), iwork(13),
     1                iwork(17), iwork(18)
c
c Estimate final error and print result.
        errfac = elkup( n, y, rwork(21), itol, rtol(j), atol(j) )
        if (errfac .gt. hun) then
          write (lout,5001)  errfac
          nerr = nerr + 1
        else
          write (lout,5000)  errfac
        endif
  100   continue
  200 continue
  300 continue
c
      write (lout,6000) nerr
      stop
c
 1000 format(20x,' Demonstration Problem for DLSODI')
 1010 format(///80('*')///)
 1100 format(/10x,' Simplified Galerkin Solution of Burgers Equation'//
     1       13x,'Diffusion coefficient is eta =',d10.2/
     2       13x,'Uniform mesh on interval',d12.3,' to ',d12.3/
     3       13x,'Zero boundary conditions'/
     4       13x,'Time limits: t0 = ',d12.5,'   tlast = ',d12.5/
     5       13x,'Half-bandwidths ml = ',i2,'   mu = ',i2/
     6       13x,'System size neq = ',i3/)
c
 1200 format('Initial profile:'/17(6d12.4/))
c
 1500 format(///80('-')///'Run with rtol =',d12.2,'  atol =',d12.2,
     1       '   mf =',i3,'   npts =',i4,':'//)
c
 2000 format('Output for time t = ',d12.5,'   current h =',
     1       d12.5,'   current order =',i2,':'/17(6d12.4/))
c
 3000 format(//'Final statistics for mf = ',i2,':'/
     1       i4,' steps,',i5,' res,',i4,' Jacobians,',
     2       '   rwork size =',i6,',   iwork size =',i6)
c
 4000 format(///80('*')//20x,'Final time reached for mf = ',i2,
     1       ' was t = ',d12.5/25x,'at which istate = ',i2////80('*'))
 5000 format('  Final output is correct to within ',d8.1,
     1       '  times local error tolerance')
 5001 format('  Final output is wrong by ',d8.1,
     1       '  times local error tolerance')
 6000 format(//80('*')//
     1       'Run completed.  Number of errors encountered =',i3)
c
c end of main program for the DLSODI demonstration problem.
      end

      subroutine res (n, t, y, v, r, ires)
c This subroutine computes the residual vector
c   r = g(t,y) - A(t,y)*v .
c It uses nm1 = n - 1 from Common.
c If ires = -1, only g(t,y) is returned in r, since A(t,y) does
c not depend on y.
c
      integer i, ires, n, nm1
      double precision t, y, v, r, r4d, eodsq, one, four, six,
     1   fact1, fact4
      dimension y(n), v(n), r(n)
      common /test1/ r4d, eodsq, nm1
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      call gfun (n, t, y, r)
      if (ires .eq. -1) return
c
      fact1 = one/six
      fact4 = four/six
      r(1) = r(1) - (fact4*v(1) + fact1*v(2))
      do 10 i = 2, nm1
  10   r(i) = r(i) - (fact1*v(i-1) + fact4*v(i) + fact1*v(i+1))
      r(n) = r(n) - (fact1*v(nm1) + fact4*v(n))
      return
c end of subroutine res for the DLSODI demonstration problem.
      end

      subroutine gfun (n, t, y, g)
c This subroutine computes the right-hand side function g(y,t).
c It uses r4d = 1/(4*delta), eodsq = eta/delta**2, and nm1 = n - 1
c from the Common block test1.
c
      integer i, n, nm1
      double precision t, y, g, r4d, eodsq, two
      dimension g(n), y(n)
      common /test1/ r4d, eodsq, nm1
      data two/2.0d0/
c
      g(1) = -r4d*y(2)**2 + eodsq*(y(2) - two*y(1))
c
      do 20 i = 2,nm1
        g(i) = r4d*(y(i-1)**2 - y(i+1)**2)
     1        + eodsq*(y(i+1) - two*y(i) + y(i-1))
   20   continue
c
      g(n) = r4d*y(nm1)**2 + eodsq*(y(nm1) - two*y(n))
c
      return
c end of subroutine gfun for the DLSODI demonstration problem.
      end

      subroutine addabd (n, t, y, ml, mu, pa, m0)
c This subroutine computes the matrix A in band form, adds it to pa,
c and returns the sum in pa.   The matrix A is tridiagonal, of order n,
c with nonzero elements (reading across) of  1/6, 4/6, 1/6.
c
      integer i, n, m0, ml, mu, mup1, mup2
      double precision t, y, pa, fact1, fact4, one, four, six
      dimension y(n), pa(m0,n)
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
c Set the pointers.
      mup1 = mu + 1
      mup2 = mu + 2
c Compute the elements of A.
      fact1 = one/six
      fact4 = four/six
c Add the matrix A to the matrix pa (banded).
      do 10 i = 1,n
        pa(mu,i) = pa(mu,i) + fact1
        pa(mup1,i) = pa(mup1,i) + fact4
        pa(mup2,i) = pa(mup2,i) + fact1
   10   continue
      return
c end of subroutine addabd for the DLSODI demonstration problem.
      end

      subroutine addafl (n, t, y, ml, mu, pa, m0)
c This subroutine computes the matrix A in full form, adds it to
c pa, and returns the sum in pa.
c It uses nm1 = n - 1 from Common.
c The matrix A is tridiagonal, of order n, with nonzero elements
c (reading across) of  1/6, 4/6, 1/6.
c
      integer i, n, m0, ml, mu, nm1
      double precision t, y, pa, r4d, eodsq, one, four, six,
     1   fact1, fact4
      dimension y(n), pa(m0,n)
      common /test1/ r4d, eodsq, nm1
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
c Compute the elements of A.
      fact1 = one/six
      fact4 = four/six
c
c Add the matrix A to the matrix pa (full).
c
      do 110  i = 2, nm1
         pa(i,i+1) = pa(i,i+1) + fact1
         pa(i,i) = pa(i,i) + fact4
         pa(i,i-1) = pa(i,i-1) + fact1
  110    continue
      pa(1,2) = pa(1,2) + fact1
      pa(1,1) = pa(1,1) + fact4
      pa(n,n) = pa(n,n) + fact4
      pa(n,nm1) = pa(n,nm1) + fact1
      return
c end of subroutine addafl for the DLSODI demonstration problem.
      end

      subroutine jacbd (n, t, y, s, ml, mu, pa, m0)
c This subroutine computes the Jacobian dg/dy = d(g-a*s)/dy
c and stores elements
c   i   j
c dg /dy   in  pa(i-j+mu+1,j)  in band matrix format.
c It uses r4d = 1/(4*delta), eodsq = eta/delta**2, and nm1 = n - 1
c from the Common block test1.
c
      integer i, n, m0, ml, mu, mup1, mup2, nm1
      double precision t, y, s, pa, diag, r4d, eodsq, two, r2d
      dimension y(n), s(n), pa(m0,n)
      common /test1/ r4d, eodsq, nm1
      data two/2.0d0/
c
      mup1 = mu + 1
      mup2 = mu + 2
      diag = -two*eodsq
      r2d = two*r4d
c                     1   1
c Compute and store dg /dy
      pa(mup1,1) = diag
c
c                     1   2
c Compute and store dg /dy
      pa(mu,2) = -r2d*y(2) + eodsq
c
      do 20 i = 2,nm1
c
c                     i   i-1
c Compute and store dg /dy
        pa(mup2,i-1) = r2d*y(i-1) + eodsq
c
c                     i   i
c Compute and store dg /dy
      pa(mup1,i) = diag
c
c                     i   i+1
c Compute and store dg /dy
        pa(mu,i+1) = -r2d*y(i+1) + eodsq
   20   continue
c
c                     n   n-1
c Compute and store dg /dy
      pa(mup2,nm1) = r2d*y(nm1) + eodsq
c
c                     n   n
c Compute and store dg /dy
      pa(mup1,n) = diag
c
      return
c end of subroutine jacbd for the DLSODI demonstration problem.
      end

      subroutine jacfl (n, t, y, s, ml, mu, pa, m0)
c This subroutine computes the Jacobian dg/dy = d(g-a*s)/dy
c and stores elements
c   i   j
c dg /dy   in  pa(i,j) in full matrix format.
c It uses r4d = 1/(4*delta), eodsq = eta/delta**2, and nm1 = n - 1
c from the Common block test1.
c
      integer i, n, m0, ml, mu, nm1
      double precision t, y, s, pa, diag, r4d, eodsq, two, r2d
      dimension y(n), s(n), pa(m0,n)
      common /test1/ r4d, eodsq, nm1
      data two/2.0d0/
c
      diag = -two*eodsq
      r2d = two*r4d
c
c                     1   1
c Compute and store dg /dy
      pa(1,1) = diag
c
c                     1   2
c Compute and store dg /dy
      pa(1,2) = -r2d*y(2) + eodsq
c
      do 120  i = 2,nm1
c
c                     i   i-1
c Compute and store dg /dy
        pa(i,i-1) = r2d*y(i-1) + eodsq
c
c                     i   i
c Compute and store dg /dy
      pa(i,i) = diag
c
c                     i   i+1
c Compute and store dg /dy
        pa(i,i+1) = -r2d*y(i+1) + eodsq
  120   continue
c
c                     n   n-1
c Compute and store dg /dy
      pa(n,nm1) = r2d*y(nm1) + eodsq
c
c                     n   n
c Compute and store dg /dy
      pa(n,n) = diag
c
      return
c end of subroutine jacfl for the DLSODI demonstration problem.
      end

      double precision function elkup (n, y, ewt, itol, rtol, atol)
c This routine looks up approximately correct values of y at t = 0.4,
c ytrue = y9 or y99 depending on whether n = 9 or 99.  These were
c obtained by running DLSODI with very tight tolerances.
c The returned value is
c elkup  =  norm of  ( y - ytrue ) / ( rtol*abs(ytrue) + atol ).
c
      integer n, itol, i
      double precision y, ewt, rtol, atol, y9, y99, y99a, y99b, y99c,
     1   y99d, y99e, y99f, y99g, dvnorm
      dimension y(n), ewt(n), y9(9), y99(99)
      dimension y99a(16), y99b(16), y99c(16), y99d(16), y99e(16),
     1          y99f(16), y99g(3)
      equivalence (y99a(1),y99(1)), (y99b(1),y99(17)),
     1      (y99c(1),y99(33)), (y99d(1),y99(49)), (y99e(1),y99(65)),
     1      (y99f(1),y99(81)), (y99g(1),y99(97))
      data y9 /
     1 1.07001457d-01, 2.77432492d-01, 5.02444616d-01, 7.21037157d-01,
     1 9.01670441d-01, 8.88832048d-01, 4.96572850d-01, 9.46924362d-02,
     1-6.90855199d-03 /
      data y99a /
     1 2.05114384d-03, 4.19527452d-03, 6.52533872d-03, 9.13412751d-03,
     1 1.21140191d-02, 1.55565301d-02, 1.95516488d-02, 2.41869487d-02,
     1 2.95465081d-02, 3.57096839d-02, 4.27498067d-02, 5.07328729d-02,
     1 5.97163151d-02, 6.97479236d-02, 8.08649804d-02, 9.30936515d-02 /
      data y99b /
     1 1.06448659d-01, 1.20933239d-01, 1.36539367d-01, 1.53248227d-01,
     1 1.71030869d-01, 1.89849031d-01, 2.09656044d-01, 2.30397804d-01,
     1 2.52013749d-01, 2.74437805d-01, 2.97599285d-01, 3.21423708d-01,
     1 3.45833531d-01, 3.70748792d-01, 3.96087655d-01, 4.21766871d-01 /
      data y99c /
     1 4.47702161d-01, 4.73808532d-01, 5.00000546d-01, 5.26192549d-01,
     1 5.52298887d-01, 5.78234121d-01, 6.03913258d-01, 6.29252015d-01,
     1 6.54167141d-01, 6.78576790d-01, 7.02400987d-01, 7.25562165d-01,
     1 7.47985803d-01, 7.69601151d-01, 7.90342031d-01, 8.10147715d-01 /
      data y99d /
     1 8.28963844d-01, 8.46743353d-01, 8.63447369d-01, 8.79046021d-01,
     1 8.93519106d-01, 9.06856541d-01, 9.19058529d-01, 9.30135374d-01,
     1 9.40106872d-01, 9.49001208d-01, 9.56853318d-01, 9.63702661d-01,
     1 9.69590361d-01, 9.74555682d-01, 9.78631814d-01, 9.81840924d-01 /
      data y99e /
     1 9.84188430d-01, 9.85656465d-01, 9.86196496d-01, 9.85721098d-01,
     1 9.84094964d-01, 9.81125395d-01, 9.76552747d-01, 9.70041743d-01,
     1 9.61175143d-01, 9.49452051d-01, 9.34294085d-01, 9.15063568d-01,
     1 8.91098383d-01, 8.61767660d-01, 8.26550038d-01, 7.85131249d-01 /
      data y99f /
     1 7.37510044d-01, 6.84092540d-01, 6.25748369d-01, 5.63802368d-01,
     1 4.99946558d-01, 4.36077986d-01, 3.74091566d-01, 3.15672765d-01,
     1 2.62134958d-01, 2.14330497d-01, 1.72640946d-01, 1.37031155d-01,
     1 1.07140815d-01, 8.23867920d-02, 6.20562432d-02, 4.53794321d-02 /
      data y99g / 3.15789227d-02, 1.98968820d-02, 9.60472135d-03 /
c
      if (n .eq. 99) go to 99
c
c Compute local error tolerance using correct y (n = 9).
c
      call dewset( n, itol, rtol, atol, y9, ewt )
c
c Invert ewt and replace y by the error, y - ytrue.
c
      do 20  i = 1, 9
        ewt(i) = 1.0d0/ewt(i)
 20     y(i) = y(i) - y9(i)
      go to 200
c
c Compute local error tolerance using correct y (n = 99).
c
 99   call dewset( n, itol, rtol, atol, y99, ewt )
c
c Invert ewt and replace y by the error, y - ytrue.
c
      do 120  i = 1, 99
        ewt(i) = 1.0d0/ewt(i)
 120    y(i) = y(i) - y99(i)
c
c Find weighted norm of the error and return.
c
 200  elkup = dvnorm (n, y, ewt)
      return
c end of function elkup for the DLSODI demonstration program.
      end

................................................................................

                     Demonstration Problem for DLSODI

           Simplified Galerkin Solution of Burgers Equation

             Diffusion coefficient is eta =  0.50E-01
             Uniform mesh on interval  -0.100E+01 to    0.100E+01
             Zero boundary conditions
             Time limits: t0 =  0.00000E+00   tlast =  0.40000E+00
             Half-bandwidths ml =  1   mu =  1
             System size neq =   9

Initial profile:
  0.0000E+00  0.0000E+00  0.5000E+00  0.1000E+01  0.1000E+01  0.1000E+01
  0.5000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 11   npts =  10:


Output for time t =  0.10000E+00   current h = 0.42496E-01   current order = 2:
  0.6081E-01  0.3900E+00  0.8101E+00  0.9897E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1388E-01  0.1294E-02
Output for time t =  0.20000E+00   current h = 0.66617E-01   current order = 2:
  0.8935E-01  0.3342E+00  0.6690E+00  0.9069E+00  0.1011E+01  0.7573E+00
  0.2383E+00 -0.1250E-02 -0.3190E-02
Output for time t =  0.30000E+00   current h = 0.66617E-01   current order = 2:
  0.1022E+00  0.3006E+00  0.5707E+00  0.8091E+00  0.9703E+00  0.8439E+00
  0.3682E+00  0.3492E-01 -0.7198E-02
Output for time t =  0.40000E+00   current h = 0.98382E-01   current order = 2:
  0.1073E+00  0.2768E+00  0.5009E+00  0.7212E+00  0.9043E+00  0.8896E+00
  0.4949E+00  0.9423E-01 -0.6444E-02


Final statistics for mf = 11:
   9 steps,   13 res,   4 Jacobians,   rwork size =   247,   iwork size =    29
  Final output is correct to within  0.8E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 12   npts =  10:


Output for time t =  0.10000E+00   current h = 0.42496E-01   current order = 2:
  0.6081E-01  0.3900E+00  0.8101E+00  0.9897E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1388E-01  0.1294E-02
Output for time t =  0.20000E+00   current h = 0.66617E-01   current order = 2:
  0.8935E-01  0.3342E+00  0.6690E+00  0.9069E+00  0.1011E+01  0.7573E+00
  0.2383E+00 -0.1250E-02 -0.3190E-02
Output for time t =  0.30000E+00   current h = 0.66617E-01   current order = 2:
  0.1022E+00  0.3006E+00  0.5707E+00  0.8091E+00  0.9703E+00  0.8439E+00
  0.3682E+00  0.3492E-01 -0.7198E-02
Output for time t =  0.40000E+00   current h = 0.98382E-01   current order = 2:
  0.1073E+00  0.2768E+00  0.5009E+00  0.7212E+00  0.9043E+00  0.8896E+00
  0.4949E+00  0.9423E-01 -0.6444E-02


Final statistics for mf = 12:
   9 steps,   53 res,   4 Jacobians,   rwork size =   247,   iwork size =    29
  Final output is correct to within  0.8E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 14   npts =  10:


Output for time t =  0.10000E+00   current h = 0.42496E-01   current order = 2:
  0.6081E-01  0.3900E+00  0.8101E+00  0.9897E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1388E-01  0.1294E-02
Output for time t =  0.20000E+00   current h = 0.66617E-01   current order = 2:
  0.8935E-01  0.3342E+00  0.6690E+00  0.9069E+00  0.1011E+01  0.7573E+00
  0.2383E+00 -0.1250E-02 -0.3190E-02
Output for time t =  0.30000E+00   current h = 0.66617E-01   current order = 2:
  0.1022E+00  0.3006E+00  0.5707E+00  0.8091E+00  0.9703E+00  0.8439E+00
  0.3682E+00  0.3492E-01 -0.7198E-02
Output for time t =  0.40000E+00   current h = 0.98382E-01   current order = 2:
  0.1073E+00  0.2768E+00  0.5009E+00  0.7212E+00  0.9043E+00  0.8896E+00
  0.4949E+00  0.9423E-01 -0.6444E-02


Final statistics for mf = 14:
   9 steps,   13 res,   4 Jacobians,   rwork size =   202,   iwork size =    29
  Final output is correct to within  0.8E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 21   npts =  10:


Output for time t =  0.10000E+00   current h = 0.50696E-01   current order = 3:
  0.6089E-01  0.3899E+00  0.8101E+00  0.9900E+00  0.1012E+01  0.6408E+00
  0.1147E+00 -0.1391E-01  0.1328E-02
Output for time t =  0.20000E+00   current h = 0.50696E-01   current order = 3:
  0.8864E-01  0.3358E+00  0.6687E+00  0.9046E+00  0.1013E+01  0.7569E+00
  0.2384E+00 -0.1103E-02 -0.3436E-02
Output for time t =  0.30000E+00   current h = 0.57639E-01   current order = 3:
  0.1014E+00  0.3024E+00  0.5717E+00  0.8058E+00  0.9704E+00  0.8448E+00
  0.3684E+00  0.3507E-01 -0.7538E-02
Output for time t =  0.40000E+00   current h = 0.57639E-01   current order = 3:
  0.1067E+00  0.2782E+00  0.5030E+00  0.7188E+00  0.9018E+00  0.8904E+00
  0.4961E+00  0.9453E-01 -0.6919E-02


Final statistics for mf = 21:
  11 steps,   14 res,   3 Jacobians,   rwork size =   184,   iwork size =    29
  Final output is correct to within  0.6E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 22   npts =  10:


Output for time t =  0.10000E+00   current h = 0.50696E-01   current order = 3:
  0.6089E-01  0.3899E+00  0.8101E+00  0.9900E+00  0.1012E+01  0.6408E+00
  0.1147E+00 -0.1391E-01  0.1328E-02
Output for time t =  0.20000E+00   current h = 0.50696E-01   current order = 3:
  0.8864E-01  0.3358E+00  0.6687E+00  0.9046E+00  0.1013E+01  0.7569E+00
  0.2384E+00 -0.1103E-02 -0.3436E-02
Output for time t =  0.30000E+00   current h = 0.57639E-01   current order = 3:
  0.1014E+00  0.3024E+00  0.5717E+00  0.8058E+00  0.9704E+00  0.8448E+00
  0.3684E+00  0.3507E-01 -0.7538E-02
Output for time t =  0.40000E+00   current h = 0.57639E-01   current order = 3:
  0.1067E+00  0.2782E+00  0.5030E+00  0.7188E+00  0.9018E+00  0.8904E+00
  0.4961E+00  0.9453E-01 -0.6919E-02


Final statistics for mf = 22:
  11 steps,   44 res,   3 Jacobians,   rwork size =   184,   iwork size =    29
  Final output is correct to within  0.6E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 24   npts =  10:


Output for time t =  0.10000E+00   current h = 0.50696E-01   current order = 3:
  0.6089E-01  0.3899E+00  0.8101E+00  0.9900E+00  0.1012E+01  0.6408E+00
  0.1147E+00 -0.1391E-01  0.1328E-02
Output for time t =  0.20000E+00   current h = 0.50696E-01   current order = 3:
  0.8864E-01  0.3358E+00  0.6687E+00  0.9046E+00  0.1013E+01  0.7569E+00
  0.2384E+00 -0.1103E-02 -0.3436E-02
Output for time t =  0.30000E+00   current h = 0.57639E-01   current order = 3:
  0.1014E+00  0.3024E+00  0.5717E+00  0.8058E+00  0.9704E+00  0.8448E+00
  0.3684E+00  0.3507E-01 -0.7538E-02
Output for time t =  0.40000E+00   current h = 0.57639E-01   current order = 3:
  0.1067E+00  0.2782E+00  0.5030E+00  0.7188E+00  0.9018E+00  0.8904E+00
  0.4961E+00  0.9453E-01 -0.6919E-02


Final statistics for mf = 24:
  11 steps,   14 res,   3 Jacobians,   rwork size =   139,   iwork size =    29
  Final output is correct to within  0.6E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 11   npts =  10:


Output for time t =  0.10000E+00   current h = 0.16073E-01   current order = 4:
  0.6054E-01  0.3907E+00  0.8099E+00  0.9886E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1382E-01  0.1185E-02
Output for time t =  0.20000E+00   current h = 0.24112E-01   current order = 5:
  0.8890E-01  0.3352E+00  0.6696E+00  0.9050E+00  0.1011E+01  0.7581E+00
  0.2384E+00 -0.1212E-02 -0.3381E-02
Output for time t =  0.30000E+00   current h = 0.33274E-01   current order = 5:
  0.1018E+00  0.3014E+00  0.5720E+00  0.8078E+00  0.9686E+00  0.8445E+00
  0.3690E+00  0.3491E-01 -0.7446E-02
Output for time t =  0.40000E+00   current h = 0.33274E-01   current order = 5:
  0.1070E+00  0.2774E+00  0.5024E+00  0.7210E+00  0.9017E+00  0.8888E+00
  0.4966E+00  0.9469E-01 -0.6909E-02


Final statistics for mf = 11:
  24 steps,   32 res,   6 Jacobians,   rwork size =   247,   iwork size =    29
  Final output is correct to within  0.4E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 12   npts =  10:


Output for time t =  0.10000E+00   current h = 0.16073E-01   current order = 4:
  0.6054E-01  0.3907E+00  0.8099E+00  0.9886E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1382E-01  0.1185E-02
Output for time t =  0.20000E+00   current h = 0.24112E-01   current order = 5:
  0.8890E-01  0.3352E+00  0.6696E+00  0.9050E+00  0.1011E+01  0.7581E+00
  0.2384E+00 -0.1212E-02 -0.3381E-02
Output for time t =  0.30000E+00   current h = 0.33274E-01   current order = 5:
  0.1018E+00  0.3014E+00  0.5720E+00  0.8078E+00  0.9686E+00  0.8445E+00
  0.3690E+00  0.3491E-01 -0.7446E-02
Output for time t =  0.40000E+00   current h = 0.33274E-01   current order = 5:
  0.1070E+00  0.2774E+00  0.5024E+00  0.7210E+00  0.9017E+00  0.8888E+00
  0.4966E+00  0.9469E-01 -0.6909E-02


Final statistics for mf = 12:
  24 steps,   92 res,   6 Jacobians,   rwork size =   247,   iwork size =    29
  Final output is correct to within  0.4E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 14   npts =  10:


Output for time t =  0.10000E+00   current h = 0.16073E-01   current order = 4:
  0.6054E-01  0.3907E+00  0.8099E+00  0.9886E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1382E-01  0.1185E-02
Output for time t =  0.20000E+00   current h = 0.24112E-01   current order = 5:
  0.8890E-01  0.3352E+00  0.6696E+00  0.9050E+00  0.1011E+01  0.7581E+00
  0.2384E+00 -0.1212E-02 -0.3381E-02
Output for time t =  0.30000E+00   current h = 0.33274E-01   current order = 5:
  0.1018E+00  0.3014E+00  0.5720E+00  0.8078E+00  0.9686E+00  0.8445E+00
  0.3690E+00  0.3491E-01 -0.7446E-02
Output for time t =  0.40000E+00   current h = 0.33274E-01   current order = 5:
  0.1070E+00  0.2774E+00  0.5024E+00  0.7210E+00  0.9017E+00  0.8888E+00
  0.4966E+00  0.9469E-01 -0.6909E-02


Final statistics for mf = 14:
  24 steps,   32 res,   6 Jacobians,   rwork size =   202,   iwork size =    29
  Final output is correct to within  0.4E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 21   npts =  10:


Output for time t =  0.10000E+00   current h = 0.10895E-01   current order = 4:
  0.6054E-01  0.3907E+00  0.8099E+00  0.9886E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1382E-01  0.1185E-02
Output for time t =  0.20000E+00   current h = 0.16165E-01   current order = 5:
  0.8890E-01  0.3352E+00  0.6696E+00  0.9050E+00  0.1011E+01  0.7581E+00
  0.2384E+00 -0.1212E-02 -0.3381E-02
Output for time t =  0.30000E+00   current h = 0.26124E-01   current order = 5:
  0.1018E+00  0.3014E+00  0.5720E+00  0.8078E+00  0.9686E+00  0.8445E+00
  0.3690E+00  0.3491E-01 -0.7446E-02
Output for time t =  0.40000E+00   current h = 0.26124E-01   current order = 5:
  0.1070E+00  0.2774E+00  0.5024E+00  0.7210E+00  0.9017E+00  0.8888E+00
  0.4966E+00  0.9469E-01 -0.6909E-02


Final statistics for mf = 21:
  35 steps,   44 res,   7 Jacobians,   rwork size =   184,   iwork size =    29
  Final output is correct to within  0.1E+01  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 22   npts =  10:


Output for time t =  0.10000E+00   current h = 0.10895E-01   current order = 4:
  0.6054E-01  0.3907E+00  0.8099E+00  0.9886E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1382E-01  0.1185E-02
Output for time t =  0.20000E+00   current h = 0.16165E-01   current order = 5:
  0.8890E-01  0.3352E+00  0.6696E+00  0.9050E+00  0.1011E+01  0.7581E+00
  0.2384E+00 -0.1212E-02 -0.3381E-02
Output for time t =  0.30000E+00   current h = 0.26124E-01   current order = 5:
  0.1018E+00  0.3014E+00  0.5720E+00  0.8078E+00  0.9686E+00  0.8445E+00
  0.3690E+00  0.3491E-01 -0.7446E-02
Output for time t =  0.40000E+00   current h = 0.26124E-01   current order = 5:
  0.1070E+00  0.2774E+00  0.5024E+00  0.7210E+00  0.9017E+00  0.8888E+00
  0.4966E+00  0.9469E-01 -0.6909E-02


Final statistics for mf = 22:
  35 steps,  114 res,   7 Jacobians,   rwork size =   184,   iwork size =    29
  Final output is correct to within  0.1E+01  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 24   npts =  10:


Output for time t =  0.10000E+00   current h = 0.10895E-01   current order = 4:
  0.6054E-01  0.3907E+00  0.8099E+00  0.9886E+00  0.1013E+01  0.6407E+00
  0.1147E+00 -0.1382E-01  0.1185E-02
Output for time t =  0.20000E+00   current h = 0.16165E-01   current order = 5:
  0.8890E-01  0.3352E+00  0.6696E+00  0.9050E+00  0.1011E+01  0.7581E+00
  0.2384E+00 -0.1212E-02 -0.3381E-02
Output for time t =  0.30000E+00   current h = 0.26124E-01   current order = 5:
  0.1018E+00  0.3014E+00  0.5720E+00  0.8078E+00  0.9686E+00  0.8445E+00
  0.3690E+00  0.3491E-01 -0.7446E-02
Output for time t =  0.40000E+00   current h = 0.26124E-01   current order = 5:
  0.1070E+00  0.2774E+00  0.5024E+00  0.7210E+00  0.9017E+00  0.8888E+00
  0.4966E+00  0.9469E-01 -0.6909E-02


Final statistics for mf = 24:
  35 steps,   44 res,   7 Jacobians,   rwork size =   139,   iwork size =    29
  Final output is correct to within  0.1E+01  times local error tolerance



********************************************************************************



                     Demonstration Problem for DLSODI

           Simplified Galerkin Solution of Burgers Equation

             Diffusion coefficient is eta =  0.50E-01
             Uniform mesh on interval  -0.100E+01 to    0.100E+01
             Zero boundary conditions
             Time limits: t0 =  0.00000E+00   tlast =  0.40000E+00
             Half-bandwidths ml =  1   mu =  1
             System size neq =  99

Initial profile:
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.5000E+00  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.5000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 14   npts = 100:


Output for time t =  0.10000E+00   current h = 0.32596E-01   current order = 2:
  0.8566E-06  0.2312E-05  0.5353E-05  0.1197E-04  0.2629E-04  0.5678E-04
  0.1203E-03  0.2492E-03  0.5031E-03  0.9860E-03  0.1869E-02  0.3419E-02
  0.6021E-02  0.1020E-01  0.1661E-01  0.2604E-01  0.3930E-01  0.5724E-01
  0.8055E-01  0.1098E+00  0.1452E+00  0.1869E+00  0.2344E+00  0.2872E+00
  0.3444E+00  0.4050E+00  0.4676E+00  0.5310E+00  0.5938E+00  0.6547E+00
  0.7123E+00  0.7656E+00  0.8136E+00  0.8556E+00  0.8913E+00  0.9206E+00
  0.9437E+00  0.9614E+00  0.9743E+00  0.9835E+00  0.9896E+00  0.9937E+00
  0.9963E+00  0.9978E+00  0.9988E+00  0.9993E+00  0.9996E+00  0.9998E+00
  0.9999E+00  0.9999E+00  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9791E+00  0.9647E+00  0.9427E+00  0.9107E+00
  0.8661E+00  0.8071E+00  0.7334E+00  0.6464E+00  0.5502E+00  0.4506E+00
  0.3541E+00  0.2668E+00  0.1928E+00  0.1336E+00  0.8887E-01  0.5683E-01
  0.3497E-01  0.2071E-01  0.1183E-01  0.6516E-02  0.3469E-02  0.1789E-02
  0.8957E-03  0.4369E-03  0.2082E-03  0.9727E-04  0.4467E-04  0.2021E-04
  0.8997E-05  0.3873E-05  0.1432E-05
Output for time t =  0.20000E+00   current h = 0.57739E-01   current order = 2:
  0.1580E-03  0.3506E-03  0.6171E-03  0.1007E-02  0.1584E-02  0.2433E-02
  0.3662E-02  0.5405E-02  0.7830E-02  0.1113E-01  0.1553E-01  0.2126E-01
  0.2858E-01  0.3774E-01  0.4899E-01  0.6253E-01  0.7855E-01  0.9717E-01
  0.1185E+00  0.1425E+00  0.1692E+00  0.1984E+00  0.2301E+00  0.2641E+00
  0.3000E+00  0.3377E+00  0.3767E+00  0.4169E+00  0.4578E+00  0.4991E+00
  0.5405E+00  0.5815E+00  0.6219E+00  0.6613E+00  0.6993E+00  0.7357E+00
  0.7700E+00  0.8022E+00  0.8319E+00  0.8589E+00  0.8831E+00  0.9045E+00
  0.9230E+00  0.9389E+00  0.9521E+00  0.9631E+00  0.9719E+00  0.9789E+00
  0.9844E+00  0.9886E+00  0.9917E+00  0.9941E+00  0.9958E+00  0.9971E+00
  0.9980E+00  0.9986E+00  0.9990E+00  0.9993E+00  0.9994E+00  0.9994E+00
  0.9993E+00  0.9990E+00  0.9985E+00  0.9976E+00  0.9961E+00  0.9939E+00
  0.9905E+00  0.9855E+00  0.9782E+00  0.9679E+00  0.9535E+00  0.9341E+00
  0.9082E+00  0.8748E+00  0.8328E+00  0.7816E+00  0.7213E+00  0.6528E+00
  0.5782E+00  0.5002E+00  0.4222E+00  0.3475E+00  0.2789E+00  0.2184E+00
  0.1671E+00  0.1249E+00  0.9146E-01  0.6562E-01  0.4617E-01  0.3190E-01
  0.2165E-01  0.1444E-01  0.9475E-02  0.6114E-02  0.3877E-02  0.2407E-02
  0.1446E-02  0.8090E-03  0.3610E-03
Output for time t =  0.30000E+00   current h = 0.57739E-01   current order = 2:
  0.9060E-03  0.1896E-02  0.3055E-02  0.4476E-02  0.6255E-02  0.8498E-02
  0.1132E-01  0.1483E-01  0.1916E-01  0.2444E-01  0.3078E-01  0.3832E-01
  0.4715E-01  0.5739E-01  0.6911E-01  0.8238E-01  0.9724E-01  0.1137E+00
  0.1318E+00  0.1515E+00  0.1728E+00  0.1955E+00  0.2197E+00  0.2452E+00
  0.2720E+00  0.2998E+00  0.3287E+00  0.3584E+00  0.3888E+00  0.4197E+00
  0.4511E+00  0.4828E+00  0.5146E+00  0.5464E+00  0.5779E+00  0.6092E+00
  0.6399E+00  0.6700E+00  0.6993E+00  0.7277E+00  0.7549E+00  0.7810E+00
  0.8056E+00  0.8288E+00  0.8505E+00  0.8704E+00  0.8887E+00  0.9052E+00
  0.9200E+00  0.9330E+00  0.9444E+00  0.9543E+00  0.9627E+00  0.9699E+00
  0.9758E+00  0.9807E+00  0.9847E+00  0.9879E+00  0.9904E+00  0.9923E+00
  0.9937E+00  0.9947E+00  0.9951E+00  0.9950E+00  0.9944E+00  0.9931E+00
  0.9910E+00  0.9878E+00  0.9831E+00  0.9767E+00  0.9679E+00  0.9561E+00
  0.9407E+00  0.9207E+00  0.8954E+00  0.8639E+00  0.8255E+00  0.7799E+00
  0.7270E+00  0.6676E+00  0.6028E+00  0.5346E+00  0.4652E+00  0.3970E+00
  0.3322E+00  0.2728E+00  0.2200E+00  0.1743E+00  0.1359E+00  0.1043E+00
  0.7895E-01  0.5891E-01  0.4336E-01  0.3144E-01  0.2240E-01  0.1557E-01
  0.1037E-01  0.6312E-02  0.2981E-02
Output for time t =  0.40000E+00   current h = 0.90345E-01   current order = 2:
  0.2030E-02  0.4154E-02  0.6468E-02  0.9065E-02  0.1204E-01  0.1548E-01
  0.1949E-01  0.2414E-01  0.2952E-01  0.3572E-01  0.4279E-01  0.5081E-01
  0.5982E-01  0.6989E-01  0.8103E-01  0.9327E-01  0.1066E+00  0.1211E+00
  0.1367E+00  0.1534E+00  0.1711E+00  0.1899E+00  0.2096E+00  0.2302E+00
  0.2518E+00  0.2741E+00  0.2971E+00  0.3208E+00  0.3451E+00  0.3700E+00
  0.3952E+00  0.4208E+00  0.4466E+00  0.4726E+00  0.4988E+00  0.5249E+00
  0.5511E+00  0.5770E+00  0.6028E+00  0.6283E+00  0.6533E+00  0.6780E+00
  0.7020E+00  0.7255E+00  0.7482E+00  0.7702E+00  0.7913E+00  0.8114E+00
  0.8305E+00  0.8486E+00  0.8655E+00  0.8812E+00  0.8958E+00  0.9091E+00
  0.9212E+00  0.9321E+00  0.9419E+00  0.9505E+00  0.9581E+00  0.9646E+00
  0.9702E+00  0.9749E+00  0.9787E+00  0.9816E+00  0.9838E+00  0.9851E+00
  0.9855E+00  0.9849E+00  0.9833E+00  0.9803E+00  0.9758E+00  0.9694E+00
  0.9606E+00  0.9490E+00  0.9340E+00  0.9149E+00  0.8910E+00  0.8618E+00
  0.8266E+00  0.7851E+00  0.7374E+00  0.6838E+00  0.6252E+00  0.5632E+00
  0.4992E+00  0.4353E+00  0.3734E+00  0.3150E+00  0.2615E+00  0.2138E+00
  0.1722E+00  0.1367E+00  0.1069E+00  0.8223E-01  0.6196E-01  0.4533E-01
  0.3156E-01  0.1989E-01  0.9605E-02


Final statistics for mf = 14:
  27 steps,   40 res,  11 Jacobians,   rwork size =  2002,   iwork size =   119
  Final output is correct to within  0.5E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 15   npts = 100:


Output for time t =  0.10000E+00   current h = 0.32596E-01   current order = 2:
  0.8566E-06  0.2312E-05  0.5353E-05  0.1197E-04  0.2629E-04  0.5678E-04
  0.1203E-03  0.2492E-03  0.5031E-03  0.9860E-03  0.1869E-02  0.3419E-02
  0.6021E-02  0.1020E-01  0.1661E-01  0.2604E-01  0.3930E-01  0.5724E-01
  0.8055E-01  0.1098E+00  0.1452E+00  0.1869E+00  0.2344E+00  0.2872E+00
  0.3444E+00  0.4050E+00  0.4676E+00  0.5310E+00  0.5938E+00  0.6547E+00
  0.7123E+00  0.7656E+00  0.8136E+00  0.8556E+00  0.8913E+00  0.9206E+00
  0.9437E+00  0.9614E+00  0.9743E+00  0.9835E+00  0.9896E+00  0.9937E+00
  0.9963E+00  0.9978E+00  0.9988E+00  0.9993E+00  0.9996E+00  0.9998E+00
  0.9999E+00  0.9999E+00  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9791E+00  0.9647E+00  0.9427E+00  0.9107E+00
  0.8661E+00  0.8071E+00  0.7334E+00  0.6464E+00  0.5502E+00  0.4506E+00
  0.3541E+00  0.2668E+00  0.1928E+00  0.1336E+00  0.8887E-01  0.5683E-01
  0.3497E-01  0.2071E-01  0.1183E-01  0.6516E-02  0.3469E-02  0.1789E-02
  0.8957E-03  0.4369E-03  0.2082E-03  0.9727E-04  0.4467E-04  0.2021E-04
  0.8997E-05  0.3873E-05  0.1432E-05
Output for time t =  0.20000E+00   current h = 0.57739E-01   current order = 2:
  0.1580E-03  0.3506E-03  0.6171E-03  0.1007E-02  0.1584E-02  0.2433E-02
  0.3662E-02  0.5405E-02  0.7830E-02  0.1113E-01  0.1553E-01  0.2126E-01
  0.2858E-01  0.3774E-01  0.4899E-01  0.6253E-01  0.7855E-01  0.9717E-01
  0.1185E+00  0.1425E+00  0.1692E+00  0.1984E+00  0.2301E+00  0.2641E+00
  0.3000E+00  0.3377E+00  0.3767E+00  0.4169E+00  0.4578E+00  0.4991E+00
  0.5405E+00  0.5815E+00  0.6219E+00  0.6613E+00  0.6993E+00  0.7357E+00
  0.7700E+00  0.8022E+00  0.8319E+00  0.8589E+00  0.8831E+00  0.9045E+00
  0.9230E+00  0.9389E+00  0.9521E+00  0.9631E+00  0.9719E+00  0.9789E+00
  0.9844E+00  0.9886E+00  0.9917E+00  0.9941E+00  0.9958E+00  0.9971E+00
  0.9980E+00  0.9986E+00  0.9990E+00  0.9993E+00  0.9994E+00  0.9994E+00
  0.9993E+00  0.9990E+00  0.9985E+00  0.9976E+00  0.9961E+00  0.9939E+00
  0.9905E+00  0.9855E+00  0.9782E+00  0.9679E+00  0.9535E+00  0.9341E+00
  0.9082E+00  0.8748E+00  0.8328E+00  0.7816E+00  0.7213E+00  0.6528E+00
  0.5782E+00  0.5002E+00  0.4222E+00  0.3475E+00  0.2789E+00  0.2184E+00
  0.1671E+00  0.1249E+00  0.9146E-01  0.6562E-01  0.4617E-01  0.3190E-01
  0.2165E-01  0.1444E-01  0.9475E-02  0.6114E-02  0.3877E-02  0.2407E-02
  0.1446E-02  0.8090E-03  0.3610E-03
Output for time t =  0.30000E+00   current h = 0.57739E-01   current order = 2:
  0.9060E-03  0.1896E-02  0.3055E-02  0.4476E-02  0.6255E-02  0.8498E-02
  0.1132E-01  0.1483E-01  0.1916E-01  0.2444E-01  0.3078E-01  0.3832E-01
  0.4715E-01  0.5739E-01  0.6911E-01  0.8238E-01  0.9724E-01  0.1137E+00
  0.1318E+00  0.1515E+00  0.1728E+00  0.1955E+00  0.2197E+00  0.2452E+00
  0.2720E+00  0.2998E+00  0.3287E+00  0.3584E+00  0.3888E+00  0.4197E+00
  0.4511E+00  0.4828E+00  0.5146E+00  0.5464E+00  0.5779E+00  0.6092E+00
  0.6399E+00  0.6700E+00  0.6993E+00  0.7277E+00  0.7549E+00  0.7810E+00
  0.8056E+00  0.8288E+00  0.8505E+00  0.8704E+00  0.8887E+00  0.9052E+00
  0.9200E+00  0.9330E+00  0.9444E+00  0.9543E+00  0.9627E+00  0.9699E+00
  0.9758E+00  0.9807E+00  0.9847E+00  0.9879E+00  0.9904E+00  0.9923E+00
  0.9937E+00  0.9947E+00  0.9951E+00  0.9950E+00  0.9944E+00  0.9931E+00
  0.9910E+00  0.9878E+00  0.9831E+00  0.9767E+00  0.9679E+00  0.9561E+00
  0.9407E+00  0.9207E+00  0.8954E+00  0.8639E+00  0.8255E+00  0.7799E+00
  0.7270E+00  0.6676E+00  0.6028E+00  0.5346E+00  0.4652E+00  0.3970E+00
  0.3322E+00  0.2728E+00  0.2200E+00  0.1743E+00  0.1359E+00  0.1043E+00
  0.7895E-01  0.5891E-01  0.4336E-01  0.3144E-01  0.2240E-01  0.1557E-01
  0.1037E-01  0.6312E-02  0.2981E-02
Output for time t =  0.40000E+00   current h = 0.90345E-01   current order = 2:
  0.2030E-02  0.4154E-02  0.6468E-02  0.9065E-02  0.1204E-01  0.1548E-01
  0.1949E-01  0.2414E-01  0.2952E-01  0.3572E-01  0.4279E-01  0.5081E-01
  0.5982E-01  0.6989E-01  0.8103E-01  0.9327E-01  0.1066E+00  0.1211E+00
  0.1367E+00  0.1534E+00  0.1711E+00  0.1899E+00  0.2096E+00  0.2302E+00
  0.2518E+00  0.2741E+00  0.2971E+00  0.3208E+00  0.3451E+00  0.3700E+00
  0.3952E+00  0.4208E+00  0.4466E+00  0.4726E+00  0.4988E+00  0.5249E+00
  0.5511E+00  0.5770E+00  0.6028E+00  0.6283E+00  0.6533E+00  0.6780E+00
  0.7020E+00  0.7255E+00  0.7482E+00  0.7702E+00  0.7913E+00  0.8114E+00
  0.8305E+00  0.8486E+00  0.8655E+00  0.8812E+00  0.8958E+00  0.9091E+00
  0.9212E+00  0.9321E+00  0.9419E+00  0.9505E+00  0.9581E+00  0.9646E+00
  0.9702E+00  0.9749E+00  0.9787E+00  0.9816E+00  0.9838E+00  0.9851E+00
  0.9855E+00  0.9849E+00  0.9833E+00  0.9803E+00  0.9758E+00  0.9694E+00
  0.9606E+00  0.9490E+00  0.9340E+00  0.9149E+00  0.8910E+00  0.8618E+00
  0.8266E+00  0.7851E+00  0.7374E+00  0.6838E+00  0.6252E+00  0.5632E+00
  0.4992E+00  0.4353E+00  0.3734E+00  0.3150E+00  0.2615E+00  0.2138E+00
  0.1722E+00  0.1367E+00  0.1069E+00  0.8223E-01  0.6196E-01  0.4533E-01
  0.3156E-01  0.1989E-01  0.9605E-02


Final statistics for mf = 15:
  27 steps,   84 res,  11 Jacobians,   rwork size =  2002,   iwork size =   119
  Final output is correct to within  0.5E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 24   npts = 100:


Output for time t =  0.10000E+00   current h = 0.20794E-01   current order = 3:
  0.8937E-06  0.2397E-05  0.5502E-05  0.1218E-04  0.2650E-04  0.5674E-04
  0.1194E-03  0.2461E-03  0.4957E-03  0.9720E-03  0.1849E-02  0.3399E-02
  0.6018E-02  0.1023E-01  0.1670E-01  0.2616E-01  0.3942E-01  0.5730E-01
  0.8054E-01  0.1097E+00  0.1451E+00  0.1867E+00  0.2342E+00  0.2871E+00
  0.3445E+00  0.4053E+00  0.4681E+00  0.5315E+00  0.5942E+00  0.6550E+00
  0.7124E+00  0.7654E+00  0.8131E+00  0.8550E+00  0.8907E+00  0.9202E+00
  0.9436E+00  0.9615E+00  0.9745E+00  0.9837E+00  0.9898E+00  0.9938E+00
  0.9963E+00  0.9978E+00  0.9988E+00  0.9993E+00  0.9996E+00  0.9998E+00
  0.9999E+00  0.9999E+00  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9791E+00  0.9648E+00  0.9428E+00  0.9109E+00
  0.8663E+00  0.8073E+00  0.7333E+00  0.6462E+00  0.5499E+00  0.4503E+00
  0.3540E+00  0.2668E+00  0.1928E+00  0.1337E+00  0.8901E-01  0.5693E-01
  0.3500E-01  0.2069E-01  0.1177E-01  0.6463E-02  0.3431E-02  0.1767E-02
  0.8858E-03  0.4335E-03  0.2078E-03  0.9780E-04  0.4530E-04  0.2068E-04
  0.9287E-05  0.4027E-05  0.1497E-05
Output for time t =  0.20000E+00   current h = 0.28924E-01   current order = 2:
  0.1546E-03  0.3436E-03  0.6066E-03  0.9934E-03  0.1569E-02  0.2420E-02
  0.3656E-02  0.5415E-02  0.7865E-02  0.1120E-01  0.1563E-01  0.2140E-01
  0.2874E-01  0.3789E-01  0.4909E-01  0.6256E-01  0.7848E-01  0.9700E-01
  0.1182E+00  0.1422E+00  0.1689E+00  0.1982E+00  0.2300E+00  0.2641E+00
  0.3003E+00  0.3381E+00  0.3774E+00  0.4178E+00  0.4588E+00  0.5002E+00
  0.5415E+00  0.5823E+00  0.6223E+00  0.6613E+00  0.6988E+00  0.7348E+00
  0.7688E+00  0.8008E+00  0.8305E+00  0.8576E+00  0.8821E+00  0.9038E+00
  0.9227E+00  0.9388E+00  0.9523E+00  0.9633E+00  0.9722E+00  0.9793E+00
  0.9847E+00  0.9889E+00  0.9920E+00  0.9943E+00  0.9960E+00  0.9972E+00
  0.9980E+00  0.9986E+00  0.9990E+00  0.9993E+00  0.9994E+00  0.9994E+00
  0.9993E+00  0.9990E+00  0.9985E+00  0.9976E+00  0.9961E+00  0.9939E+00
  0.9905E+00  0.9856E+00  0.9784E+00  0.9681E+00  0.9538E+00  0.9343E+00
  0.9084E+00  0.8749E+00  0.8328E+00  0.7814E+00  0.7210E+00  0.6524E+00
  0.5778E+00  0.4999E+00  0.4220E+00  0.3475E+00  0.2790E+00  0.2186E+00
  0.1673E+00  0.1251E+00  0.9161E-01  0.6570E-01  0.4620E-01  0.3187E-01
  0.2159E-01  0.1436E-01  0.9396E-02  0.6043E-02  0.3819E-02  0.2363E-02
  0.1416E-02  0.7900E-03  0.3520E-03
Output for time t =  0.30000E+00   current h = 0.49754E-01   current order = 3:
  0.9126E-03  0.1910E-02  0.3079E-02  0.4512E-02  0.6305E-02  0.8562E-02
  0.1139E-01  0.1492E-01  0.1925E-01  0.2452E-01  0.3085E-01  0.3836E-01
  0.4716E-01  0.5736E-01  0.6904E-01  0.8227E-01  0.9710E-01  0.1136E+00
  0.1316E+00  0.1514E+00  0.1726E+00  0.1955E+00  0.2197E+00  0.2454E+00
  0.2722E+00  0.3002E+00  0.3292E+00  0.3590E+00  0.3895E+00  0.4205E+00
  0.4520E+00  0.4836E+00  0.5154E+00  0.5471E+00  0.5786E+00  0.6097E+00
  0.6402E+00  0.6700E+00  0.6991E+00  0.7271E+00  0.7541E+00  0.7799E+00
  0.8044E+00  0.8276E+00  0.8492E+00  0.8692E+00  0.8876E+00  0.9043E+00
  0.9193E+00  0.9325E+00  0.9441E+00  0.9542E+00  0.9627E+00  0.9700E+00
  0.9760E+00  0.9809E+00  0.9850E+00  0.9882E+00  0.9907E+00  0.9926E+00
  0.9940E+00  0.9949E+00  0.9953E+00  0.9953E+00  0.9946E+00  0.9933E+00
  0.9912E+00  0.9879E+00  0.9833E+00  0.9769E+00  0.9681E+00  0.9563E+00
  0.9408E+00  0.9208E+00  0.8954E+00  0.8638E+00  0.8253E+00  0.7796E+00
  0.7268E+00  0.6674E+00  0.6027E+00  0.5346E+00  0.4653E+00  0.3972E+00
  0.3325E+00  0.2731E+00  0.2202E+00  0.1745E+00  0.1360E+00  0.1044E+00
  0.7900E-01  0.5892E-01  0.4333E-01  0.3139E-01  0.2233E-01  0.1550E-01
  0.1031E-01  0.6267E-02  0.2957E-02
Output for time t =  0.40000E+00   current h = 0.49754E-01   current order = 3:
  0.2055E-02  0.4206E-02  0.6545E-02  0.9167E-02  0.1216E-01  0.1562E-01
  0.1964E-01  0.2429E-01  0.2966E-01  0.3582E-01  0.4285E-01  0.5081E-01
  0.5976E-01  0.6975E-01  0.8082E-01  0.9300E-01  0.1063E+00  0.1208E+00
  0.1363E+00  0.1530E+00  0.1708E+00  0.1897E+00  0.2096E+00  0.2304E+00
  0.2521E+00  0.2747E+00  0.2979E+00  0.3219E+00  0.3463E+00  0.3713E+00
  0.3967E+00  0.4224E+00  0.4482E+00  0.4743E+00  0.5003E+00  0.5263E+00
  0.5522E+00  0.5779E+00  0.6034E+00  0.6284E+00  0.6530E+00  0.6772E+00
  0.7008E+00  0.7238E+00  0.7462E+00  0.7679E+00  0.7888E+00  0.8090E+00
  0.8282E+00  0.8464E+00  0.8636E+00  0.8796E+00  0.8945E+00  0.9082E+00
  0.9206E+00  0.9318E+00  0.9418E+00  0.9506E+00  0.9584E+00  0.9650E+00
  0.9707E+00  0.9754E+00  0.9792E+00  0.9822E+00  0.9844E+00  0.9857E+00
  0.9861E+00  0.9855E+00  0.9839E+00  0.9809E+00  0.9763E+00  0.9697E+00
  0.9608E+00  0.9491E+00  0.9339E+00  0.9146E+00  0.8906E+00  0.8613E+00
  0.8262E+00  0.7849E+00  0.7374E+00  0.6841E+00  0.6258E+00  0.5639E+00
  0.5000E+00  0.4361E+00  0.3741E+00  0.3156E+00  0.2620E+00  0.2141E+00
  0.1724E+00  0.1368E+00  0.1069E+00  0.8219E-01  0.6188E-01  0.4524E-01
  0.3147E-01  0.1982E-01  0.9566E-02


Final statistics for mf = 24:
  36 steps,   47 res,  13 Jacobians,   rwork size =  1309,   iwork size =   119
  Final output is correct to within  0.4E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-02  atol =    0.10E-02   mf = 25   npts = 100:


Output for time t =  0.10000E+00   current h = 0.20794E-01   current order = 3:
  0.8937E-06  0.2397E-05  0.5502E-05  0.1218E-04  0.2650E-04  0.5674E-04
  0.1194E-03  0.2461E-03  0.4957E-03  0.9720E-03  0.1849E-02  0.3399E-02
  0.6018E-02  0.1023E-01  0.1670E-01  0.2616E-01  0.3942E-01  0.5730E-01
  0.8054E-01  0.1097E+00  0.1451E+00  0.1867E+00  0.2342E+00  0.2871E+00
  0.3445E+00  0.4053E+00  0.4681E+00  0.5315E+00  0.5942E+00  0.6550E+00
  0.7124E+00  0.7654E+00  0.8131E+00  0.8550E+00  0.8907E+00  0.9202E+00
  0.9436E+00  0.9615E+00  0.9745E+00  0.9837E+00  0.9898E+00  0.9938E+00
  0.9963E+00  0.9978E+00  0.9988E+00  0.9993E+00  0.9996E+00  0.9998E+00
  0.9999E+00  0.9999E+00  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9791E+00  0.9648E+00  0.9428E+00  0.9109E+00
  0.8663E+00  0.8073E+00  0.7333E+00  0.6462E+00  0.5499E+00  0.4503E+00
  0.3540E+00  0.2668E+00  0.1928E+00  0.1337E+00  0.8901E-01  0.5693E-01
  0.3500E-01  0.2069E-01  0.1177E-01  0.6463E-02  0.3431E-02  0.1767E-02
  0.8858E-03  0.4335E-03  0.2078E-03  0.9780E-04  0.4530E-04  0.2068E-04
  0.9287E-05  0.4027E-05  0.1497E-05
Output for time t =  0.20000E+00   current h = 0.28924E-01   current order = 2:
  0.1546E-03  0.3436E-03  0.6066E-03  0.9934E-03  0.1569E-02  0.2420E-02
  0.3656E-02  0.5415E-02  0.7865E-02  0.1120E-01  0.1563E-01  0.2140E-01
  0.2874E-01  0.3789E-01  0.4909E-01  0.6256E-01  0.7848E-01  0.9700E-01
  0.1182E+00  0.1422E+00  0.1689E+00  0.1982E+00  0.2300E+00  0.2641E+00
  0.3003E+00  0.3381E+00  0.3774E+00  0.4178E+00  0.4588E+00  0.5002E+00
  0.5415E+00  0.5823E+00  0.6223E+00  0.6613E+00  0.6988E+00  0.7348E+00
  0.7688E+00  0.8008E+00  0.8305E+00  0.8576E+00  0.8821E+00  0.9038E+00
  0.9227E+00  0.9388E+00  0.9523E+00  0.9633E+00  0.9722E+00  0.9793E+00
  0.9847E+00  0.9889E+00  0.9920E+00  0.9943E+00  0.9960E+00  0.9972E+00
  0.9980E+00  0.9986E+00  0.9990E+00  0.9993E+00  0.9994E+00  0.9994E+00
  0.9993E+00  0.9990E+00  0.9985E+00  0.9976E+00  0.9961E+00  0.9939E+00
  0.9905E+00  0.9856E+00  0.9784E+00  0.9681E+00  0.9538E+00  0.9343E+00
  0.9084E+00  0.8749E+00  0.8328E+00  0.7814E+00  0.7210E+00  0.6524E+00
  0.5778E+00  0.4999E+00  0.4220E+00  0.3475E+00  0.2790E+00  0.2186E+00
  0.1673E+00  0.1251E+00  0.9161E-01  0.6570E-01  0.4620E-01  0.3187E-01
  0.2159E-01  0.1436E-01  0.9396E-02  0.6043E-02  0.3819E-02  0.2363E-02
  0.1416E-02  0.7900E-03  0.3520E-03
Output for time t =  0.30000E+00   current h = 0.49754E-01   current order = 3:
  0.9126E-03  0.1910E-02  0.3079E-02  0.4512E-02  0.6305E-02  0.8562E-02
  0.1139E-01  0.1492E-01  0.1925E-01  0.2452E-01  0.3085E-01  0.3836E-01
  0.4716E-01  0.5736E-01  0.6904E-01  0.8227E-01  0.9710E-01  0.1136E+00
  0.1316E+00  0.1514E+00  0.1726E+00  0.1955E+00  0.2197E+00  0.2454E+00
  0.2722E+00  0.3002E+00  0.3292E+00  0.3590E+00  0.3895E+00  0.4205E+00
  0.4520E+00  0.4836E+00  0.5154E+00  0.5471E+00  0.5786E+00  0.6097E+00
  0.6402E+00  0.6700E+00  0.6991E+00  0.7271E+00  0.7541E+00  0.7799E+00
  0.8044E+00  0.8276E+00  0.8492E+00  0.8692E+00  0.8876E+00  0.9043E+00
  0.9193E+00  0.9325E+00  0.9441E+00  0.9542E+00  0.9627E+00  0.9700E+00
  0.9760E+00  0.9809E+00  0.9850E+00  0.9882E+00  0.9907E+00  0.9926E+00
  0.9940E+00  0.9949E+00  0.9953E+00  0.9953E+00  0.9946E+00  0.9933E+00
  0.9912E+00  0.9879E+00  0.9833E+00  0.9769E+00  0.9681E+00  0.9563E+00
  0.9408E+00  0.9208E+00  0.8954E+00  0.8638E+00  0.8253E+00  0.7796E+00
  0.7268E+00  0.6674E+00  0.6027E+00  0.5346E+00  0.4653E+00  0.3972E+00
  0.3325E+00  0.2731E+00  0.2202E+00  0.1745E+00  0.1360E+00  0.1044E+00
  0.7900E-01  0.5892E-01  0.4333E-01  0.3139E-01  0.2233E-01  0.1550E-01
  0.1031E-01  0.6267E-02  0.2957E-02
Output for time t =  0.40000E+00   current h = 0.49754E-01   current order = 3:
  0.2055E-02  0.4206E-02  0.6545E-02  0.9167E-02  0.1216E-01  0.1562E-01
  0.1964E-01  0.2429E-01  0.2966E-01  0.3582E-01  0.4285E-01  0.5081E-01
  0.5976E-01  0.6975E-01  0.8082E-01  0.9300E-01  0.1063E+00  0.1208E+00
  0.1363E+00  0.1530E+00  0.1708E+00  0.1897E+00  0.2096E+00  0.2304E+00
  0.2521E+00  0.2747E+00  0.2979E+00  0.3219E+00  0.3463E+00  0.3713E+00
  0.3967E+00  0.4224E+00  0.4482E+00  0.4743E+00  0.5003E+00  0.5263E+00
  0.5522E+00  0.5779E+00  0.6034E+00  0.6284E+00  0.6530E+00  0.6772E+00
  0.7008E+00  0.7238E+00  0.7462E+00  0.7679E+00  0.7888E+00  0.8090E+00
  0.8282E+00  0.8464E+00  0.8636E+00  0.8796E+00  0.8945E+00  0.9082E+00
  0.9206E+00  0.9318E+00  0.9418E+00  0.9506E+00  0.9584E+00  0.9650E+00
  0.9707E+00  0.9754E+00  0.9792E+00  0.9822E+00  0.9844E+00  0.9857E+00
  0.9861E+00  0.9855E+00  0.9839E+00  0.9809E+00  0.9763E+00  0.9697E+00
  0.9608E+00  0.9491E+00  0.9339E+00  0.9146E+00  0.8906E+00  0.8613E+00
  0.8262E+00  0.7849E+00  0.7374E+00  0.6841E+00  0.6258E+00  0.5639E+00
  0.5000E+00  0.4361E+00  0.3741E+00  0.3156E+00  0.2620E+00  0.2141E+00
  0.1724E+00  0.1368E+00  0.1069E+00  0.8219E-01  0.6188E-01  0.4524E-01
  0.3147E-01  0.1982E-01  0.9566E-02


Final statistics for mf = 25:
  36 steps,   99 res,  13 Jacobians,   rwork size =  1309,   iwork size =   119
  Final output is correct to within  0.4E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 14   npts = 100:


Output for time t =  0.10000E+00   current h = 0.73068E-02   current order = 3:
  0.2343E-06  0.8142E-06  0.2455E-05  0.6875E-05  0.1806E-04  0.4468E-04
  0.1046E-03  0.2321E-03  0.4899E-03  0.9854E-03  0.1892E-02  0.3476E-02
  0.6116E-02  0.1033E-01  0.1675E-01  0.2615E-01  0.3935E-01  0.5718E-01
  0.8039E-01  0.1095E+00  0.1449E+00  0.1866E+00  0.2342E+00  0.2871E+00
  0.3446E+00  0.4054E+00  0.4682E+00  0.5318E+00  0.5946E+00  0.6554E+00
  0.7129E+00  0.7658E+00  0.8134E+00  0.8551E+00  0.8905E+00  0.9196E+00
  0.9428E+00  0.9606E+00  0.9738E+00  0.9832E+00  0.9897E+00  0.9939E+00
  0.9965E+00  0.9981E+00  0.9990E+00  0.9995E+00  0.9998E+00  0.9999E+00
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9792E+00  0.9648E+00  0.9429E+00  0.9109E+00
  0.8663E+00  0.8072E+00  0.7333E+00  0.6462E+00  0.5498E+00  0.4502E+00
  0.3538E+00  0.2667E+00  0.1928E+00  0.1337E+00  0.8909E-01  0.5707E-01
  0.3516E-01  0.2084E-01  0.1188E-01  0.6510E-02  0.3423E-02  0.1725E-02
  0.8322E-03  0.3833E-03  0.1682E-03  0.7012E-04  0.2769E-04  0.1032E-04
  0.3609E-05  0.1174E-05  0.3324E-06
Output for time t =  0.20000E+00   current h = 0.63398E-02   current order = 2:
  0.1565E-03  0.3492E-03  0.6187E-03  0.1016E-02  0.1605E-02  0.2469E-02
  0.3717E-02  0.5481E-02  0.7923E-02  0.1123E-01  0.1563E-01  0.2135E-01
  0.2864E-01  0.3776E-01  0.4895E-01  0.6243E-01  0.7838E-01  0.9695E-01
  0.1182E+00  0.1422E+00  0.1690E+00  0.1983E+00  0.2301E+00  0.2641E+00
  0.3002E+00  0.3380E+00  0.3772E+00  0.4175E+00  0.4586E+00  0.5000E+00
  0.5414E+00  0.5825E+00  0.6228E+00  0.6620E+00  0.6998E+00  0.7359E+00
  0.7699E+00  0.8017E+00  0.8310E+00  0.8578E+00  0.8818E+00  0.9031E+00
  0.9216E+00  0.9376E+00  0.9511E+00  0.9622E+00  0.9714E+00  0.9787E+00
  0.9844E+00  0.9888E+00  0.9921E+00  0.9945E+00  0.9963E+00  0.9975E+00
  0.9984E+00  0.9989E+00  0.9993E+00  0.9995E+00  0.9996E+00  0.9996E+00
  0.9994E+00  0.9991E+00  0.9986E+00  0.9976E+00  0.9962E+00  0.9940E+00
  0.9906E+00  0.9856E+00  0.9783E+00  0.9680E+00  0.9537E+00  0.9342E+00
  0.9083E+00  0.8748E+00  0.8327E+00  0.7814E+00  0.7210E+00  0.6525E+00
  0.5779E+00  0.5000E+00  0.4221E+00  0.3475E+00  0.2790E+00  0.2186E+00
  0.1673E+00  0.1252E+00  0.9168E-01  0.6579E-01  0.4630E-01  0.3197E-01
  0.2167E-01  0.1442E-01  0.9417E-02  0.6036E-02  0.3791E-02  0.2326E-02
  0.1378E-02  0.7604E-03  0.3358E-03
Output for time t =  0.30000E+00   current h = 0.17382E-01   current order = 4:
  0.9295E-03  0.1942E-02  0.3124E-02  0.4565E-02  0.6361E-02  0.8616E-02
  0.1144E-01  0.1495E-01  0.1926E-01  0.2451E-01  0.3082E-01  0.3831E-01
  0.4709E-01  0.5727E-01  0.6894E-01  0.8216E-01  0.9698E-01  0.1134E+00
  0.1316E+00  0.1513E+00  0.1726E+00  0.1954E+00  0.2197E+00  0.2454E+00
  0.2723E+00  0.3003E+00  0.3293E+00  0.3591E+00  0.3897E+00  0.4208E+00
  0.4523E+00  0.4841E+00  0.5159E+00  0.5477E+00  0.5792E+00  0.6103E+00
  0.6409E+00  0.6707E+00  0.6997E+00  0.7277E+00  0.7546E+00  0.7803E+00
  0.8046E+00  0.8274E+00  0.8487E+00  0.8684E+00  0.8866E+00  0.9030E+00
  0.9178E+00  0.9311E+00  0.9427E+00  0.9529E+00  0.9617E+00  0.9691E+00
  0.9754E+00  0.9806E+00  0.9849E+00  0.9883E+00  0.9910E+00  0.9930E+00
  0.9945E+00  0.9954E+00  0.9958E+00  0.9957E+00  0.9951E+00  0.9937E+00
  0.9915E+00  0.9882E+00  0.9835E+00  0.9770E+00  0.9682E+00  0.9564E+00
  0.9408E+00  0.9208E+00  0.8953E+00  0.8637E+00  0.8252E+00  0.7795E+00
  0.7267E+00  0.6673E+00  0.6027E+00  0.5346E+00  0.4654E+00  0.3973E+00
  0.3327E+00  0.2733E+00  0.2205E+00  0.1748E+00  0.1363E+00  0.1046E+00
  0.7914E-01  0.5902E-01  0.4339E-01  0.3142E-01  0.2235E-01  0.1550E-01
  0.1030E-01  0.6260E-02  0.2953E-02
Output for time t =  0.40000E+00   current h = 0.24310E-01   current order = 4:
  0.2051E-02  0.4195E-02  0.6525E-02  0.9134E-02  0.1211E-01  0.1556E-01
  0.1955E-01  0.2419E-01  0.2955E-01  0.3571E-01  0.4275E-01  0.5073E-01
  0.5972E-01  0.6975E-01  0.8087E-01  0.9309E-01  0.1064E+00  0.1209E+00
  0.1365E+00  0.1532E+00  0.1710E+00  0.1898E+00  0.2097E+00  0.2304E+00
  0.2520E+00  0.2744E+00  0.2976E+00  0.3214E+00  0.3458E+00  0.3707E+00
  0.3961E+00  0.4218E+00  0.4477E+00  0.4738E+00  0.5000E+00  0.5262E+00
  0.5523E+00  0.5782E+00  0.6039E+00  0.6293E+00  0.6542E+00  0.6786E+00
  0.7024E+00  0.7256E+00  0.7480E+00  0.7696E+00  0.7903E+00  0.8102E+00
  0.8290E+00  0.8467E+00  0.8634E+00  0.8790E+00  0.8935E+00  0.9069E+00
  0.9191E+00  0.9301E+00  0.9401E+00  0.9490E+00  0.9569E+00  0.9637E+00
  0.9696E+00  0.9746E+00  0.9786E+00  0.9818E+00  0.9842E+00  0.9857E+00
  0.9862E+00  0.9857E+00  0.9841E+00  0.9811E+00  0.9766E+00  0.9700E+00
  0.9612E+00  0.9495E+00  0.9343E+00  0.9151E+00  0.8911E+00  0.8618E+00
  0.8265E+00  0.7851E+00  0.7375E+00  0.6841E+00  0.6257E+00  0.5638E+00
  0.4999E+00  0.4361E+00  0.3741E+00  0.3157E+00  0.2621E+00  0.2143E+00
  0.1726E+00  0.1370E+00  0.1071E+00  0.8239E-01  0.6206E-01  0.4538E-01
  0.3158E-01  0.1990E-01  0.9605E-02


Final statistics for mf = 14:
  95 steps,  123 res,  20 Jacobians,   rwork size =  2002,   iwork size =   119
  Final output is correct to within  0.6E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 15   npts = 100:


Output for time t =  0.10000E+00   current h = 0.73068E-02   current order = 3:
  0.2343E-06  0.8142E-06  0.2455E-05  0.6875E-05  0.1806E-04  0.4468E-04
  0.1046E-03  0.2321E-03  0.4899E-03  0.9854E-03  0.1892E-02  0.3476E-02
  0.6116E-02  0.1033E-01  0.1675E-01  0.2615E-01  0.3935E-01  0.5718E-01
  0.8039E-01  0.1095E+00  0.1449E+00  0.1866E+00  0.2342E+00  0.2871E+00
  0.3446E+00  0.4054E+00  0.4682E+00  0.5318E+00  0.5946E+00  0.6554E+00
  0.7129E+00  0.7658E+00  0.8134E+00  0.8551E+00  0.8905E+00  0.9196E+00
  0.9428E+00  0.9606E+00  0.9738E+00  0.9832E+00  0.9897E+00  0.9939E+00
  0.9965E+00  0.9981E+00  0.9990E+00  0.9995E+00  0.9998E+00  0.9999E+00
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9792E+00  0.9648E+00  0.9429E+00  0.9109E+00
  0.8663E+00  0.8072E+00  0.7333E+00  0.6462E+00  0.5498E+00  0.4502E+00
  0.3538E+00  0.2667E+00  0.1928E+00  0.1337E+00  0.8909E-01  0.5707E-01
  0.3516E-01  0.2084E-01  0.1188E-01  0.6510E-02  0.3423E-02  0.1725E-02
  0.8322E-03  0.3833E-03  0.1682E-03  0.7012E-04  0.2769E-04  0.1032E-04
  0.3609E-05  0.1174E-05  0.3324E-06
Output for time t =  0.20000E+00   current h = 0.63398E-02   current order = 2:
  0.1565E-03  0.3492E-03  0.6187E-03  0.1016E-02  0.1605E-02  0.2469E-02
  0.3717E-02  0.5481E-02  0.7923E-02  0.1123E-01  0.1563E-01  0.2135E-01
  0.2864E-01  0.3776E-01  0.4895E-01  0.6243E-01  0.7838E-01  0.9695E-01
  0.1182E+00  0.1422E+00  0.1690E+00  0.1983E+00  0.2301E+00  0.2641E+00
  0.3002E+00  0.3380E+00  0.3772E+00  0.4175E+00  0.4586E+00  0.5000E+00
  0.5414E+00  0.5825E+00  0.6228E+00  0.6620E+00  0.6998E+00  0.7359E+00
  0.7699E+00  0.8017E+00  0.8310E+00  0.8578E+00  0.8818E+00  0.9031E+00
  0.9216E+00  0.9376E+00  0.9511E+00  0.9622E+00  0.9714E+00  0.9787E+00
  0.9844E+00  0.9888E+00  0.9921E+00  0.9945E+00  0.9963E+00  0.9975E+00
  0.9984E+00  0.9989E+00  0.9993E+00  0.9995E+00  0.9996E+00  0.9996E+00
  0.9994E+00  0.9991E+00  0.9986E+00  0.9976E+00  0.9962E+00  0.9940E+00
  0.9906E+00  0.9856E+00  0.9783E+00  0.9680E+00  0.9537E+00  0.9342E+00
  0.9083E+00  0.8748E+00  0.8327E+00  0.7814E+00  0.7210E+00  0.6525E+00
  0.5779E+00  0.5000E+00  0.4221E+00  0.3475E+00  0.2790E+00  0.2186E+00
  0.1673E+00  0.1252E+00  0.9168E-01  0.6579E-01  0.4630E-01  0.3197E-01
  0.2167E-01  0.1442E-01  0.9417E-02  0.6036E-02  0.3791E-02  0.2326E-02
  0.1378E-02  0.7604E-03  0.3358E-03
Output for time t =  0.30000E+00   current h = 0.17382E-01   current order = 4:
  0.9295E-03  0.1942E-02  0.3124E-02  0.4565E-02  0.6361E-02  0.8616E-02
  0.1144E-01  0.1495E-01  0.1926E-01  0.2451E-01  0.3082E-01  0.3831E-01
  0.4709E-01  0.5727E-01  0.6894E-01  0.8216E-01  0.9698E-01  0.1134E+00
  0.1316E+00  0.1513E+00  0.1726E+00  0.1954E+00  0.2197E+00  0.2454E+00
  0.2723E+00  0.3003E+00  0.3293E+00  0.3591E+00  0.3897E+00  0.4208E+00
  0.4523E+00  0.4841E+00  0.5159E+00  0.5477E+00  0.5792E+00  0.6103E+00
  0.6409E+00  0.6707E+00  0.6997E+00  0.7277E+00  0.7546E+00  0.7803E+00
  0.8046E+00  0.8274E+00  0.8487E+00  0.8684E+00  0.8866E+00  0.9030E+00
  0.9178E+00  0.9311E+00  0.9427E+00  0.9529E+00  0.9617E+00  0.9691E+00
  0.9754E+00  0.9806E+00  0.9849E+00  0.9883E+00  0.9910E+00  0.9930E+00
  0.9945E+00  0.9954E+00  0.9958E+00  0.9957E+00  0.9951E+00  0.9937E+00
  0.9915E+00  0.9882E+00  0.9835E+00  0.9770E+00  0.9682E+00  0.9564E+00
  0.9408E+00  0.9208E+00  0.8953E+00  0.8637E+00  0.8252E+00  0.7795E+00
  0.7267E+00  0.6673E+00  0.6027E+00  0.5346E+00  0.4654E+00  0.3973E+00
  0.3327E+00  0.2733E+00  0.2205E+00  0.1748E+00  0.1363E+00  0.1046E+00
  0.7914E-01  0.5902E-01  0.4339E-01  0.3142E-01  0.2235E-01  0.1550E-01
  0.1030E-01  0.6260E-02  0.2953E-02
Output for time t =  0.40000E+00   current h = 0.24310E-01   current order = 4:
  0.2051E-02  0.4195E-02  0.6525E-02  0.9134E-02  0.1211E-01  0.1556E-01
  0.1955E-01  0.2419E-01  0.2955E-01  0.3571E-01  0.4275E-01  0.5073E-01
  0.5972E-01  0.6975E-01  0.8087E-01  0.9309E-01  0.1064E+00  0.1209E+00
  0.1365E+00  0.1532E+00  0.1710E+00  0.1898E+00  0.2097E+00  0.2304E+00
  0.2520E+00  0.2744E+00  0.2976E+00  0.3214E+00  0.3458E+00  0.3707E+00
  0.3961E+00  0.4218E+00  0.4477E+00  0.4738E+00  0.5000E+00  0.5262E+00
  0.5523E+00  0.5782E+00  0.6039E+00  0.6293E+00  0.6542E+00  0.6786E+00
  0.7024E+00  0.7256E+00  0.7480E+00  0.7696E+00  0.7903E+00  0.8102E+00
  0.8290E+00  0.8467E+00  0.8634E+00  0.8790E+00  0.8935E+00  0.9069E+00
  0.9191E+00  0.9301E+00  0.9401E+00  0.9490E+00  0.9569E+00  0.9637E+00
  0.9696E+00  0.9746E+00  0.9786E+00  0.9818E+00  0.9842E+00  0.9857E+00
  0.9862E+00  0.9857E+00  0.9841E+00  0.9811E+00  0.9766E+00  0.9700E+00
  0.9612E+00  0.9495E+00  0.9343E+00  0.9151E+00  0.8911E+00  0.8618E+00
  0.8265E+00  0.7851E+00  0.7375E+00  0.6841E+00  0.6257E+00  0.5638E+00
  0.4999E+00  0.4361E+00  0.3741E+00  0.3157E+00  0.2621E+00  0.2143E+00
  0.1726E+00  0.1370E+00  0.1071E+00  0.8239E-01  0.6206E-01  0.4538E-01
  0.3158E-01  0.1990E-01  0.9605E-02


Final statistics for mf = 15:
  95 steps,  203 res,  20 Jacobians,   rwork size =  2002,   iwork size =   119
  Final output is correct to within  0.6E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 24   npts = 100:


Output for time t =  0.10000E+00   current h = 0.64687E-02   current order = 5:
  0.2331E-06  0.8110E-06  0.2449E-05  0.6866E-05  0.1805E-04  0.4468E-04
  0.1046E-03  0.2321E-03  0.4900E-03  0.9854E-03  0.1892E-02  0.3476E-02
  0.6116E-02  0.1033E-01  0.1675E-01  0.2615E-01  0.3935E-01  0.5718E-01
  0.8039E-01  0.1095E+00  0.1449E+00  0.1866E+00  0.2342E+00  0.2871E+00
  0.3446E+00  0.4054E+00  0.4682E+00  0.5318E+00  0.5946E+00  0.6554E+00
  0.7129E+00  0.7658E+00  0.8134E+00  0.8551E+00  0.8905E+00  0.9196E+00
  0.9428E+00  0.9606E+00  0.9738E+00  0.9832E+00  0.9897E+00  0.9939E+00
  0.9965E+00  0.9981E+00  0.9990E+00  0.9995E+00  0.9998E+00  0.9999E+00
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9792E+00  0.9648E+00  0.9429E+00  0.9109E+00
  0.8663E+00  0.8072E+00  0.7333E+00  0.6462E+00  0.5498E+00  0.4502E+00
  0.3538E+00  0.2667E+00  0.1928E+00  0.1337E+00  0.8909E-01  0.5707E-01
  0.3516E-01  0.2084E-01  0.1188E-01  0.6510E-02  0.3423E-02  0.1726E-02
  0.8323E-03  0.3833E-03  0.1682E-03  0.7010E-04  0.2767E-04  0.1030E-04
  0.3600E-05  0.1170E-05  0.3311E-06
Output for time t =  0.20000E+00   current h = 0.11780E-01   current order = 5:
  0.1566E-03  0.3492E-03  0.6188E-03  0.1016E-02  0.1605E-02  0.2470E-02
  0.3717E-02  0.5481E-02  0.7924E-02  0.1123E-01  0.1563E-01  0.2135E-01
  0.2864E-01  0.3776E-01  0.4895E-01  0.6243E-01  0.7838E-01  0.9695E-01
  0.1182E+00  0.1422E+00  0.1690E+00  0.1983E+00  0.2301E+00  0.2641E+00
  0.3002E+00  0.3380E+00  0.3772E+00  0.4175E+00  0.4586E+00  0.5000E+00
  0.5414E+00  0.5825E+00  0.6228E+00  0.6620E+00  0.6998E+00  0.7359E+00
  0.7699E+00  0.8017E+00  0.8310E+00  0.8578E+00  0.8818E+00  0.9030E+00
  0.9216E+00  0.9376E+00  0.9511E+00  0.9622E+00  0.9714E+00  0.9786E+00
  0.9844E+00  0.9888E+00  0.9921E+00  0.9945E+00  0.9963E+00  0.9975E+00
  0.9984E+00  0.9989E+00  0.9993E+00  0.9995E+00  0.9996E+00  0.9996E+00
  0.9994E+00  0.9991E+00  0.9986E+00  0.9976E+00  0.9962E+00  0.9940E+00
  0.9906E+00  0.9856E+00  0.9783E+00  0.9680E+00  0.9537E+00  0.9342E+00
  0.9083E+00  0.8748E+00  0.8327E+00  0.7814E+00  0.7210E+00  0.6525E+00
  0.5779E+00  0.5000E+00  0.4221E+00  0.3475E+00  0.2790E+00  0.2186E+00
  0.1673E+00  0.1252E+00  0.9168E-01  0.6580E-01  0.4630E-01  0.3197E-01
  0.2167E-01  0.1442E-01  0.9416E-02  0.6035E-02  0.3791E-02  0.2325E-02
  0.1378E-02  0.7601E-03  0.3357E-03
Output for time t =  0.30000E+00   current h = 0.15321E-01   current order = 5:
  0.9295E-03  0.1942E-02  0.3124E-02  0.4565E-02  0.6361E-02  0.8616E-02
  0.1144E-01  0.1495E-01  0.1926E-01  0.2451E-01  0.3082E-01  0.3831E-01
  0.4709E-01  0.5727E-01  0.6894E-01  0.8216E-01  0.9698E-01  0.1134E+00
  0.1316E+00  0.1513E+00  0.1726E+00  0.1954E+00  0.2197E+00  0.2454E+00
  0.2723E+00  0.3003E+00  0.3293E+00  0.3591E+00  0.3897E+00  0.4208E+00
  0.4523E+00  0.4841E+00  0.5159E+00  0.5477E+00  0.5792E+00  0.6103E+00
  0.6409E+00  0.6707E+00  0.6997E+00  0.7277E+00  0.7546E+00  0.7803E+00
  0.8046E+00  0.8274E+00  0.8487E+00  0.8685E+00  0.8866E+00  0.9030E+00
  0.9178E+00  0.9311E+00  0.9427E+00  0.9529E+00  0.9617E+00  0.9691E+00
  0.9754E+00  0.9806E+00  0.9849E+00  0.9883E+00  0.9910E+00  0.9930E+00
  0.9945E+00  0.9954E+00  0.9958E+00  0.9957E+00  0.9951E+00  0.9937E+00
  0.9915E+00  0.9882E+00  0.9835E+00  0.9770E+00  0.9682E+00  0.9564E+00
  0.9408E+00  0.9208E+00  0.8953E+00  0.8637E+00  0.8252E+00  0.7795E+00
  0.7267E+00  0.6673E+00  0.6027E+00  0.5346E+00  0.4654E+00  0.3973E+00
  0.3327E+00  0.2733E+00  0.2205E+00  0.1748E+00  0.1363E+00  0.1046E+00
  0.7914E-01  0.5902E-01  0.4339E-01  0.3142E-01  0.2235E-01  0.1550E-01
  0.1030E-01  0.6260E-02  0.2953E-02
Output for time t =  0.40000E+00   current h = 0.19409E-01   current order = 5:
  0.2051E-02  0.4195E-02  0.6525E-02  0.9134E-02  0.1211E-01  0.1556E-01
  0.1955E-01  0.2419E-01  0.2955E-01  0.3571E-01  0.4275E-01  0.5073E-01
  0.5972E-01  0.6975E-01  0.8087E-01  0.9309E-01  0.1064E+00  0.1209E+00
  0.1365E+00  0.1532E+00  0.1710E+00  0.1898E+00  0.2097E+00  0.2304E+00
  0.2520E+00  0.2744E+00  0.2976E+00  0.3214E+00  0.3458E+00  0.3707E+00
  0.3961E+00  0.4218E+00  0.4477E+00  0.4738E+00  0.5000E+00  0.5262E+00
  0.5523E+00  0.5782E+00  0.6039E+00  0.6292E+00  0.6542E+00  0.6786E+00
  0.7024E+00  0.7256E+00  0.7480E+00  0.7696E+00  0.7903E+00  0.8101E+00
  0.8290E+00  0.8467E+00  0.8635E+00  0.8790E+00  0.8935E+00  0.9069E+00
  0.9191E+00  0.9301E+00  0.9401E+00  0.9490E+00  0.9569E+00  0.9637E+00
  0.9696E+00  0.9746E+00  0.9786E+00  0.9818E+00  0.9842E+00  0.9857E+00
  0.9862E+00  0.9857E+00  0.9841E+00  0.9811E+00  0.9766E+00  0.9700E+00
  0.9612E+00  0.9495E+00  0.9343E+00  0.9151E+00  0.8911E+00  0.8618E+00
  0.8265E+00  0.7851E+00  0.7375E+00  0.6841E+00  0.6257E+00  0.5638E+00
  0.4999E+00  0.4361E+00  0.3741E+00  0.3157E+00  0.2621E+00  0.2143E+00
  0.1726E+00  0.1370E+00  0.1071E+00  0.8239E-01  0.6206E-01  0.4538E-01
  0.3158E-01  0.1990E-01  0.9605E-02


Final statistics for mf = 24:
 100 steps,  121 res,  18 Jacobians,   rwork size =  1309,   iwork size =   119
  Final output is correct to within  0.7E+00  times local error tolerance



--------------------------------------------------------------------------------


Run with rtol =    0.10E-05  atol =    0.10E-05   mf = 25   npts = 100:


Output for time t =  0.10000E+00   current h = 0.64687E-02   current order = 5:
  0.2331E-06  0.8110E-06  0.2449E-05  0.6866E-05  0.1805E-04  0.4468E-04
  0.1046E-03  0.2321E-03  0.4900E-03  0.9854E-03  0.1892E-02  0.3476E-02
  0.6116E-02  0.1033E-01  0.1675E-01  0.2615E-01  0.3935E-01  0.5718E-01
  0.8039E-01  0.1095E+00  0.1449E+00  0.1866E+00  0.2342E+00  0.2871E+00
  0.3446E+00  0.4054E+00  0.4682E+00  0.5318E+00  0.5946E+00  0.6554E+00
  0.7129E+00  0.7658E+00  0.8134E+00  0.8551E+00  0.8905E+00  0.9196E+00
  0.9428E+00  0.9606E+00  0.9738E+00  0.9832E+00  0.9897E+00  0.9939E+00
  0.9965E+00  0.9981E+00  0.9990E+00  0.9995E+00  0.9998E+00  0.9999E+00
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01  0.1000E+01
  0.9999E+00  0.9998E+00  0.9996E+00  0.9992E+00  0.9983E+00  0.9966E+00
  0.9935E+00  0.9881E+00  0.9792E+00  0.9648E+00  0.9429E+00  0.9109E+00
  0.8663E+00  0.8072E+00  0.7333E+00  0.6462E+00  0.5498E+00  0.4502E+00
  0.3538E+00  0.2667E+00  0.1928E+00  0.1337E+00  0.8909E-01  0.5707E-01
  0.3516E-01  0.2084E-01  0.1188E-01  0.6510E-02  0.3423E-02  0.1726E-02
  0.8323E-03  0.3833E-03  0.1682E-03  0.7010E-04  0.2767E-04  0.1030E-04
  0.3600E-05  0.1170E-05  0.3311E-06
Output for time t =  0.20000E+00   current h = 0.11780E-01   current order = 5:
  0.1566E-03  0.3492E-03  0.6188E-03  0.1016E-02  0.1605E-02  0.2470E-02
  0.3717E-02  0.5481E-02  0.7924E-02  0.1123E-01  0.1563E-01  0.2135E-01
  0.2864E-01  0.3776E-01  0.4895E-01  0.6243E-01  0.7838E-01  0.9695E-01
  0.1182E+00  0.1422E+00  0.1690E+00  0.1983E+00  0.2301E+00  0.2641E+00
  0.3002E+00  0.3380E+00  0.3772E+00  0.4175E+00  0.4586E+00  0.5000E+00
  0.5414E+00  0.5825E+00  0.6228E+00  0.6620E+00  0.6998E+00  0.7359E+00
  0.7699E+00  0.8017E+00  0.8310E+00  0.8578E+00  0.8818E+00  0.9030E+00
  0.9216E+00  0.9376E+00  0.9511E+00  0.9622E+00  0.9714E+00  0.9786E+00
  0.9844E+00  0.9888E+00  0.9921E+00  0.9945E+00  0.9963E+00  0.9975E+00
  0.9984E+00  0.9989E+00  0.9993E+00  0.9995E+00  0.9996E+00  0.9996E+00
  0.9994E+00  0.9991E+00  0.9986E+00  0.9976E+00  0.9962E+00  0.9940E+00
  0.9906E+00  0.9856E+00  0.9783E+00  0.9680E+00  0.9537E+00  0.9342E+00
  0.9083E+00  0.8748E+00  0.8327E+00  0.7814E+00  0.7210E+00  0.6525E+00
  0.5779E+00  0.5000E+00  0.4221E+00  0.3475E+00  0.2790E+00  0.2186E+00
  0.1673E+00  0.1252E+00  0.9168E-01  0.6580E-01  0.4630E-01  0.3197E-01
  0.2167E-01  0.1442E-01  0.9416E-02  0.6035E-02  0.3791E-02  0.2325E-02
  0.1378E-02  0.7601E-03  0.3357E-03
Output for time t =  0.30000E+00   current h = 0.15321E-01   current order = 5:
  0.9295E-03  0.1942E-02  0.3124E-02  0.4565E-02  0.6361E-02  0.8616E-02
  0.1144E-01  0.1495E-01  0.1926E-01  0.2451E-01  0.3082E-01  0.3831E-01
  0.4709E-01  0.5727E-01  0.6894E-01  0.8216E-01  0.9698E-01  0.1134E+00
  0.1316E+00  0.1513E+00  0.1726E+00  0.1954E+00  0.2197E+00  0.2454E+00
  0.2723E+00  0.3003E+00  0.3293E+00  0.3591E+00  0.3897E+00  0.4208E+00
  0.4523E+00  0.4841E+00  0.5159E+00  0.5477E+00  0.5792E+00  0.6103E+00
  0.6409E+00  0.6707E+00  0.6997E+00  0.7277E+00  0.7546E+00  0.7803E+00
  0.8046E+00  0.8274E+00  0.8487E+00  0.8685E+00  0.8866E+00  0.9030E+00
  0.9178E+00  0.9311E+00  0.9427E+00  0.9529E+00  0.9617E+00  0.9691E+00
  0.9754E+00  0.9806E+00  0.9849E+00  0.9883E+00  0.9910E+00  0.9930E+00
  0.9945E+00  0.9954E+00  0.9958E+00  0.9957E+00  0.9951E+00  0.9937E+00
  0.9915E+00  0.9882E+00  0.9835E+00  0.9770E+00  0.9682E+00  0.9564E+00
  0.9408E+00  0.9208E+00  0.8953E+00  0.8637E+00  0.8252E+00  0.7795E+00
  0.7267E+00  0.6673E+00  0.6027E+00  0.5346E+00  0.4654E+00  0.3973E+00
  0.3327E+00  0.2733E+00  0.2205E+00  0.1748E+00  0.1363E+00  0.1046E+00
  0.7914E-01  0.5902E-01  0.4339E-01  0.3142E-01  0.2235E-01  0.1550E-01
  0.1030E-01  0.6260E-02  0.2953E-02
Output for time t =  0.40000E+00   current h = 0.19409E-01   current order = 5:
  0.2051E-02  0.4195E-02  0.6525E-02  0.9134E-02  0.1211E-01  0.1556E-01
  0.1955E-01  0.2419E-01  0.2955E-01  0.3571E-01  0.4275E-01  0.5073E-01
  0.5972E-01  0.6975E-01  0.8087E-01  0.9309E-01  0.1064E+00  0.1209E+00
  0.1365E+00  0.1532E+00  0.1710E+00  0.1898E+00  0.2097E+00  0.2304E+00
  0.2520E+00  0.2744E+00  0.2976E+00  0.3214E+00  0.3458E+00  0.3707E+00
  0.3961E+00  0.4218E+00  0.4477E+00  0.4738E+00  0.5000E+00  0.5262E+00
  0.5523E+00  0.5782E+00  0.6039E+00  0.6292E+00  0.6542E+00  0.6786E+00
  0.7024E+00  0.7256E+00  0.7480E+00  0.7696E+00  0.7903E+00  0.8101E+00
  0.8290E+00  0.8467E+00  0.8635E+00  0.8790E+00  0.8935E+00  0.9069E+00
  0.9191E+00  0.9301E+00  0.9401E+00  0.9490E+00  0.9569E+00  0.9637E+00
  0.9696E+00  0.9746E+00  0.9786E+00  0.9818E+00  0.9842E+00  0.9857E+00
  0.9862E+00  0.9857E+00  0.9841E+00  0.9811E+00  0.9766E+00  0.9700E+00
  0.9612E+00  0.9495E+00  0.9343E+00  0.9151E+00  0.8911E+00  0.8618E+00
  0.8265E+00  0.7851E+00  0.7375E+00  0.6841E+00  0.6257E+00  0.5638E+00
  0.4999E+00  0.4361E+00  0.3741E+00  0.3157E+00  0.2621E+00  0.2143E+00
  0.1726E+00  0.1370E+00  0.1071E+00  0.8239E-01  0.6206E-01  0.4538E-01
  0.3158E-01  0.1990E-01  0.9605E-02


Final statistics for mf = 25:
 100 steps,  193 res,  18 Jacobians,   rwork size =  1309,   iwork size =   119
  Final output is correct to within  0.7E+00  times local error tolerance


********************************************************************************

Run completed.  Number of errors encountered =  0

==========Source and Output for LSOIBT Demonstration Program====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSOIBT package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c This program solves a semi-discretized form of the following system
c Of three PDEs (each similar to a Burgers equation):
c
c   u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)    (i=1,2,3),
c       t                           x                xx
c
c on the interval  -1 .le. x .le. 1, and with time t .ge. 0.
c The diffusion coefficients are eta(*) = .1, .02, .01.
c The boundary conditions are u(i) = 0 at x = -1 and x = 1 for all i.
c The initial profile for each u(i) is a square wave:
c     u(i) = 0         on 1/2 .lt. abs(x) .le. 1
c     u(i) = amp(i)/2  on abs(x) = 1/2
c     u(i) = amp(i)    on 0 .le. abs(x) .lt. 1/2
c where the amplitudes are amp(*) = .2, .3, .5.
c
c A simplified Galerkin treatment of the spatial variable x is used,
c with piecewise linear basis functions on a uniform mesh of 100
c intervals.  The result is a system of ODEs in the discrete values
c u(i,k) approximating u(i)  (i=1,2,3) at the interior points
c (k = 1,...,99).  The ODEs are:
c
c    .            .        .
c   (u(i,k-1) + 4 u(i,k) + u(i,k+1))/6  =
c
c     -(1/6dx) (c(k-1)dul(i) + 2c(k)(dul(i)+dur(i)) + c(k+1)dur(i))
c
c     + (eta(i)/dx**2) (dur(i) - dul(i))     (i=1,2,3,  k=1,...,99),
c
c where
c     c(j) = u(1,j)+u(2,j)+u(3,j),   dx = .02 = the interval size,
c     dul(i) = u(i,k) - u(i,k-1),   dur(i) = u(i,k+1) - u(i,k).
c Terms involving boundary values (subscripts 0 or 100) are dropped
c from the equations for k = 1 and k = 99 above.
c
c The problem is run for each of the 4 values of mf, and for two values
c of the tolerances.  Output is taken at t = .1, .2, .3, .4.
c Output is on unit lout, set to 6 in a data statement below.
c-----------------------------------------------------------------------
      external res, addabt, jacbt
      integer ncomp, nip, nm1
      integer i, io, istate, itol, iwork, jtol, lout, liw, lrw,
     1   meth, miter, mf, neq, nerr, nint, nout
      double precision eodsq, r6d
      double precision abermx, atol, dx, errfac, eta, hun, one,
     1   rtol, rwork, six, t, tinit, tlast, tout, tols, two, y, ydoti
      dimension eta(3), y(297), ydoti(297), tout(4), tols(2)
      dimension rwork(7447), iwork(317)
c Pass problem parameters in the common block par.
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
c
c Set problem parameters and run parameters
      data eta/0.1d0,0.02d0,0.01d0/, tinit/0.0d0/, tlast/0.4d0/
      data one/1.0d0/, two/2.0d0/, six/6.0d0/, hun/100.0d0/
      data tout/.10d0,.20d0,.30d0,.40d0/
      data lout/6/, nout/4/, lrw/7447/, liw/317/
      data itol/1/, tols/1.0d-3, 1.0d-6/
c
c Set mesh parameters nint, dxc etc.
      nint = 100
      ncomp = 3
      dx = two/nint
      r6d = one/(six*dx)
      do 10 i = 1,ncomp
 10     eodsq(i) = eta(i)/dx**2
      nip = nint - 1
      neq = ncomp*nip
      nm1 = nip - 1
      iwork(1) = ncomp
      iwork(2) = nip
c
      nerr = 0
c
c Set the initial conditions (for output purposes only).
      call setic (nint, ncomp, y)
c
      write (lout,1000)
      write (lout,1100) (eta(i),i=1,ncomp), tinit, tlast, nint,
     1      ncomp, nip, neq
      write (lout,1200)
      call edit (y, ncomp, nip, lout)
c
c The jtol loop is over error tolerances.
      do 200 jtol = 1,2
      rtol = tols(jtol)
      atol = rtol
c
c The meth/miter loops cover 4 values of method flag mf.
      do 100 meth = 1,2
       do 100 miter = 1,2
        mf = 10*meth + miter
c
c Set the initial conditions.
        call setic (nint, ncomp, y)
        t = tinit
        istate = 0
c
        write (lout,1500)  rtol, atol, mf
c
c Loop over output times for each case
        do 80 io = 1,nout
c
          call dlsoibt (res, addabt,jacbt, neq, y, ydoti, t, tout(io),
     1     itol,rtol,atol, 1, istate, 0, rwork,lrw,iwork,liw, mf)
c
          write (lout,2000) t, rwork(11), iwork(14), iwork(11)
          if (io .eq. nout) call edit (y, ncomp, nip, lout)
c
c If istate is not 2 on return, print message and go to next case.
          if (istate .ne. 2) then
            write (lout,4000) mf, t, istate
            nerr = nerr + 1
            go to 100
            endif
c
 80       continue
c
c Print final statistics.
        write (lout,3000) mf, iwork(11), iwork(12), iwork(13),
     1         iwork(17), iwork(18)
c
c Estimate final error and print result.
        call maxerr (y, ncomp, nip, abermx)
        errfac = abermx/tols(jtol)
        if (errfac .lt. hun) then
          write (lout,5000) errfac
        else  
          write (lout,5100) errfac
          nerr = nerr + 1
          endif
 100    continue
 200  continue
c
      write (lout,6000) nerr
      stop
c
 1000 format(/20x,' Demonstration Problem for DLSOIBT'//
     1   10x,'Galerkin method solution of system of 3 PDEs:'//
     2   10x,'  u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)',
     3   5x,'(i=1,2,3)',/16x,'t',27x,'x',16x,'xx'//
     4   10x,'x interval is -1 to 1,  zero boundary conditions'/
     5   10x,'x discretized using piecewise linear basis functions')
 1100 format(/10x,'Fixed parameters are as follows:'/
     1       13x,'Diffusion coefficients are eta =',3d10.2/
     2       13x,'t0 = ',d12.5/13x,'tlast = ',d12.5/
     3       13x,'Uniform mesh, number of intervals =',i4/
     4       13x,'Block size mb =',i2/13x,'Number of blocks nb =',i4/
     5       13x,'ODE system size neq =',i5//)
c
 1200 format(/'Initial profiles:'/)
c
 1500 format(////90('*')//'Run with rtol =',d9.1,'  atol =',d9.1,
     1       '   mf =',i3///)
c
 2000 format(' At time t =',d12.5,'  current h =',d12.5,
     1       '  current order =',i2,'  current nst =',i5/)
c
 3000 format(//'Final statistics for mf = ',i2,':',
     1       i5,' steps,',i6,' res,',i6,' jacobians,'/
     2       30x,       'rwork size =',i6,',  iwork size =',i6)
c
 4000 format(//20x,'Final time reached for mf = ',i2,
     1       ' was t = ',d12.5/25x,'at which istate = ',i2//)
 5000 format('Final output is correct to within ',d9.2,
     1       '  times local error tolerance. ')
 5100 format('Final output is wrong by ',d9.2,
     1       '  times local error tolerance.')
 6000 format(//90('*')//'Run completed: ',i3,' errors encountered')
c
c end of main program for the DLSOIBT demonstration problem
      end

      subroutine setic (nint, mb, y)
c This routine loads the y array with initial data based on a
c square wave profile for each of the mb PDE variables.
c
      integer nint, mb, i, k, nip, n14, n34
      double precision y,  amp, half, zero
      dimension y(mb,*), amp(3)
      data zero/0.0d0/, half/0.5d0/, amp/0.2d0,0.3d0,0.5d0/
c
      nip = nint - 1
      n14 = nint/4
      n34 = 3*n14
c
      do 15 k = 1,n14-1
        do 10 i = 1,mb
 10       y(i,k) = zero
 15     continue
c
      do 20 i = 1,mb
 20     y(i,n14) = half*amp(i)
c
      do 35 k = n14+1,n34-1
        do 30 i = 1,mb
 30       y(i,k) = amp(i)
 35     continue
c
      do 40 i = 1,mb
 40     y(i,n34) = half*amp(i)
c
      do 55 k = n34+1,nip
        do 50 i = 1,mb
 50       y(i,k) = zero
 55     continue
c
      return
c end of subroutine setic
      end

      subroutine res (n, t, y, v, r, ires)
c This subroutine computes the residual vector
c   r = g(t,y) - A(t,y)*v
c using routines gfun and subav.
c If ires = -1, only g(t,y) is returned in r, since A(t,y) does
c not depend on y.
c No changes need to be made to this routine if nip is changed.
c
      integer ires, n,  ncomp, nip, nm1
      double precision t, y, v, r, r6d, eodsq
      dimension y(n), v(n), r(n)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
c
      call gfun (t, y, r, ncomp)
      if (ires .eq. -1) return
c
      call subav (r, v, ncomp)
c
      return
c end of subroutine res
      end

      subroutine gfun (t, y, g, mb)
c This subroutine computes the right-hand side function g(y,t).
c It uses r6d = 1/(6*dx), eodsq(*) = eta(*)/dx**2, nip,
c and nm1 = nip - 1 from the Common block par.
c
      integer mb,  ncomp, nip, nm1,  i, k
      double precision t, y, g,  r6d, eodsq,  cc, cl, cr, dli, dri, two
      dimension g(mb,*), y(mb,*)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data two/2.0d0/
c
c left-most interior point (k = 1)
      cc = y(1,1) + y(2,1) + y(3,1)
      cr = y(1,2) + y(2,2) + y(3,2)
      do 10 i = 1,mb
        dri = y(i,2) - y(i,1)
        g(i,1) = -r6d*(two*cc*y(i,2) + cr*dri)
     1         + eodsq(i)*(dri - y(i,1))
 10     continue
c
c interior points k = 2 to nip-1
      do 20 k = 2,nm1
        cl = y(1,k-1) + y(2,k-1) + y(3,k-1)
        cc = y(1,k) + y(2,k) + y(3,k)
        cr = y(1,k+1) + y(2,k+1) + y(3,k+1)
        do 15 i = 1,mb
          dli = y(i,k) - y(i,k-1)
          dri = y(i,k+1) - y(i,k)
          g(i,k) = -r6d*(cl*dli + two*cc*(dli + dri) + cr*dri)
     1           + eodsq(i)*(dri - dli)
 15       continue
 20     continue
c
c right-most interior point (k = nip)
      cl = y(1,nm1) + y(2,nm1) + y(3,nm1)
      cc = y(1,nip) + y(2,nip) + y(3,nip)
      do 30 i = 1,mb
        dli = y(i,nip) - y(i,nm1)
        g(i,nip) = -r6d*(cl*dli - two*cc*y(i,nm1))
     1           - eodsq(i)*(y(i,nip) + dli)
 30     continue
c
      return
c end of subroutine gfun
      end

      subroutine subav (r, v, mb)
c This routine subtracts the matrix a time the vector v from r,
c in order to form the residual vector, stored in r.
c
      integer mb,  ncomp, nip, nm1,  i, k
      double precision r, v,  r6d, eodsq,  aa1, aa4, four, one, six
      dimension r(mb,*), v(mb,*)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      aa1 = one/six
      aa4 = four/six
c
      do 10 i = 1,mb
 10     r(i,1) = r(i,1) - (aa4*v(i,1) + aa1*v(i,2))
c
      do 20 k = 2,nm1
        do 15 i = 1,mb
 15       r(i,k) = r(i,k) - (aa1*v(i,k-1) + aa4*v(i,k) + aa1*v(i,k+1))
 20     continue
c
      do 30 i = 1,mb
 30     r(i,nip) = r(i,nip) - (aa1*v(i,nm1) + aa4*v(i,nip))
c
      return
c end of subroutine subav
      end

      subroutine addabt (n, t, y, mb, nb, pa, pb, pc)
c This subroutine computes the elements of the matrix A,
c and adds them to pa, pb, and pc in the appropriate manner.
c The matrix A is tridiagonal, of order n, with
c nonzero elements (reading across) of 1/6, 4/6, 1/6.
c
      integer n, mb, nb, i, k
      double precision pa, pb, pc, t, y,  aa1, aa4, four, one, six
      dimension y(mb,nb),pa(mb,mb,nb),pb(mb,mb,nb),pc(mb,mb,nb)
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
      aa1 = one/six
      aa4 = four/six
      do 50 k = 1,nb
        do 10 i = 1,mb
 10       pa(i,i,k) = pa(i,i,k) + aa4
        if (k .ne. nb) then
          do 20 i = 1,mb
 20         pb(i,i,k) = pb(i,i,k) + aa1
        endif
        if (k .ne. 1) then
          do 30 i = 1,mb
 30         pc(i,i,k) = pc(i,i,k) + aa1
        endif
 50     continue
c
      return
c end of subroutine addabt
      end

      subroutine jacbt (n, t, y, s, mb, nb, pa, pb, pc)
c This subroutine computes the Jacobian dg/dy = d(g-A*s)/dy
c which has block-tridiagonal structure.  The main, upper, and
c lower diagonals are stored in pa, pb, and pc, respectively.
c
      integer n, mb, nb,  ncomp, nip, nm1, i, j, k
      double precision t, y, s, pa, pb, pc,  r6d, eodsq,  cc, cl, cr,
     1   dlj, drj, paij, pbij, pcij, terma, termb, termc, two
      dimension y(mb,nb),s(n),pa(mb,mb,nb),pb(mb,mb,nb),pc(mb,mb,nb)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data two/2.0d0/
c
c left-most interior point (k = 1)
      cc = y(1,1) + y(2,1) + y(3,1)
      cr = y(1,2) + y(2,2) + y(3,2)
      terma = r6d*cr
      termb = -r6d*(two*cc + cr)
      do 20 j = 1,mb
        drj = y(j,2) - y(j,1)
        paij = -r6d*two*y(j,2)
        pbij = -r6d*drj
        do 10 i = 1,mb
          pa(i,j,1) = paij
 10       pb(i,j,1) = pbij
        pa(j,j,1) = pa(j,j,1) + terma - two*eodsq(j)
        pb(j,j,1) = pb(j,j,1) + termb + eodsq(j)
 20     continue
c
c interior points k = 2 to nip-1
      do 50 k = 2,nm1
        cl = y(1,k-1) + y(2,k-1) + y(3,k-1)
        cc = y(1,k) + y(2,k) + y(3,k)
        cr = y(1,k+1) + y(2,k+1) + y(3,k+1)
        terma = r6d*(cr - cl)
        termb = -r6d*(two*cc + cr)
        termc = r6d*(two*cc + cl)
        do 40 j = 1,mb
          dlj = y(j,k) - y(j,k-1)
          drj = y(j,k+1) - y(j,k)
          paij = -r6d*two*(dlj + drj)
          pbij = -r6d*drj
          pcij = -r6d*dlj
          do 30 i = 1,mb
            pa(i,j,k) = paij
            pb(i,j,k) = pbij
 30         pc(i,j,k) = pcij
          pa(j,j,k) = pa(j,j,k) + terma - two*eodsq(j)
          pb(j,j,k) = pb(j,j,k) + termb + eodsq(j)
          pc(j,j,k) = pc(j,j,k) + termc + eodsq(j)
 40       continue
 50     continue
c
c right-most interior point (k = nip)
      cl = y(1,nm1) + y(2,nm1) + y(3,nm1)
      cc = y(1,nip) + y(2,nip) + y(3,nip)
      terma = -r6d*cl
      termc = r6d*(two*cc + cl)
      do 70 j = 1,mb
        dlj = y(j,nip) - y(j,nm1)
        paij = r6d*two*y(j,nm1)
        pcij = -r6d*dlj
        do 60 i = 1,mb
          pa(i,j,nip) = paij
 60       pc(i,j,nip) = pcij
        pa(j,j,nip) = pa(j,j,nip) + terma - two*eodsq(j)
        pc(j,j,nip) = pc(j,j,nip) + termc + eodsq(j)
 70     continue
c
      return
c end of subroutine jacbt
      end

      subroutine edit (y, mb, nip, lout)
c This routine prints output.  For each of the mb PDE components, the
c values at the nip points are printed.  All output is on unit lout.
c
      integer mb, nip, lout,  i, k
      double precision y
      dimension y(mb,nip)
c
      do 10 i = 1,mb
 10      write (lout,20) i, (y(i,k),k=1,nip)
c
 20   format(' Values of PDE component i =',i3/15(7d12.4/))
c
      return
c end of subroutine edit
      end

      subroutine maxerr (y, mb, nb, abermx)
c This routine computes the maximum absolute error in the y array,
c as a computed solution at t = 0.4, using data-loaded values for
c accurate answers (from running with much smaller tolerances).
c
      integer mb, nb,  k
      double precision y, abermx,  ae1, ae2, ae3, u1, u2, u3, zero,
     1   u1a, u1b, u1c, u1d, u1e, u1f, u1g,
     2   u2a, u2b, u2c, u2d, u2e, u2f, u2g,
     3   u3a, u3b, u3c, u3d, u3e, u3f, u3g
      dimension y(mb,nb), u1(99), u2(99), u3(99)
      dimension u1a(16),u1b(16),u1c(16),u1d(16),u1e(16),u1f(16),u1g(3),
     2          u2a(16),u2b(16),u2c(16),u2d(16),u2e(16),u2f(16),u2g(3),
     3          u3a(16),u3b(16),u3c(16),u3d(16),u3e(16),u3f(16),u3g(3)
      equivalence (u1a(1),u1(1)), (u1b(1),u1(17)),
     1      (u1c(1),u1(33)), (u1d(1),u1(49)), (u1e(1),u1(65)),
     1      (u1f(1),u1(81)), (u1g(1),u1(97))
      equivalence (u2a(1),u2(1)), (u2b(1),u2(17)),
     1      (u2c(1),u2(33)), (u2d(1),u2(49)), (u2e(1),u2(65)),
     1      (u2f(1),u2(81)), (u2g(1),u2(97))
      equivalence (u3a(1),u3(1)), (u3b(1),u3(17)),
     1      (u3c(1),u3(33)), (u3d(1),u3(49)), (u3e(1),u3(65)),
     1      (u3f(1),u3(81)), (u3g(1),u3(97))
c
      data u1a /
     1  1.70956682d-03, 3.43398445d-03, 5.18783349d-03, 6.98515842d-03,
     1  8.83921016d-03, 1.07622016d-02, 1.27650806d-02, 1.48573251d-02,
     1  1.70467655d-02, 1.93394396d-02, 2.17394852d-02, 2.42490773d-02,
     1  2.68684152d-02, 2.95957660d-02, 3.24275691d-02, 3.53586054d-02/
      data u1b /
     1  3.83822285d-02, 4.14906520d-02, 4.46752791d-02, 4.79270545d-02,
     1  5.12368132d-02, 5.45956048d-02, 5.79949684d-02, 6.14271460d-02,
     1  6.48852271d-02, 6.83632267d-02, 7.18561029d-02, 7.53597274d-02,
     1  7.88708192d-02, 8.23868545d-02, 8.59059616d-02, 8.94268082d-02/
      data u1c /
     1  9.29484864d-02, 9.64703968d-02, 9.99921344d-02, 1.03513375d-01,
     1  1.07033760d-01, 1.10552783d-01, 1.14069668d-01, 1.17583246d-01,
     1  1.21091827d-01, 1.24593066d-01, 1.28083828d-01, 1.31560049d-01,
     1  1.35016617d-01, 1.38447256d-01, 1.41844451d-01, 1.45199401d-01/
      data u1d /
     1  1.48502033d-01, 1.51741065d-01, 1.54904135d-01, 1.57977973d-01,
     1  1.60948623d-01, 1.63801670d-01, 1.66522463d-01, 1.69096305d-01,
     1  1.71508595d-01, 1.73744902d-01, 1.75790974d-01, 1.77632682d-01,
     1  1.79255895d-01, 1.80646319d-01, 1.81789276d-01, 1.82669470d-01/
      data u1e /
     1  1.83270725d-01, 1.83575716d-01, 1.83565712d-01, 1.83220322d-01,
     1  1.82517279d-01, 1.81432251d-01, 1.79938706d-01, 1.78007835d-01,
     1  1.75608540d-01, 1.72707519d-01, 1.69269456d-01, 1.65257378d-01,
     1  1.60633244d-01, 1.55358941d-01, 1.49398029d-01, 1.42718981d-01/
      data u1f /
     1  1.35301474d-01, 1.27148627d-01, 1.18308730d-01, 1.08905085d-01,
     1  9.91559295d-02, 8.93515884d-02, 7.97824293d-02, 7.06663514d-02,
     1  6.21244732d-02, 5.41994827d-02, 4.68848207d-02, 4.01465202d-02,
     1  3.39357642d-02, 2.81954415d-02, 2.28635569d-02, 1.78750916d-02/
      data u1g / 1.31630892d-02, 8.65933391d-03, 4.29480447d-03/
      data u2a /
     1  7.17416019d-06, 1.70782645d-05, 3.31245126d-05, 6.01588363d-05,
     1  1.05339286d-04, 1.79174771d-04, 2.96719122d-04, 4.78862606d-04,
     1  7.53598916d-04, 1.15707860d-03, 1.73420412d-03, 2.53849668d-03,
     1  3.63099110d-03, 5.07800919d-03, 6.94782549d-03, 9.30645443d-03/
      data u2b /
     1  1.22130079d-02, 1.57152366d-02, 1.98459102d-02, 2.46205841d-02,
     1  3.00370492d-02, 3.60764461d-02, 4.27057301d-02, 4.98809820d-02,
     1  5.75510102d-02, 6.56607602d-02, 7.41541974d-02, 8.29764928d-02,
     1  9.20754824d-02, 1.01402468d-01, 1.10912474d-01, 1.20564094d-01/
      data u2c /
     1  1.30319039d-01, 1.40141489d-01, 1.49997326d-01, 1.59853293d-01,
     1  1.69676126d-01, 1.79431680d-01, 1.89084097d-01, 1.98595037d-01,
     1  2.07923034d-01, 2.17023055d-01, 2.25846345d-01, 2.34340694d-01,
     1  2.42451240d-01, 2.50121934d-01, 2.57297724d-01, 2.63927433d-01/
      data u2d /
     1  2.69967170d-01, 2.75383917d-01, 2.80158840d-01, 2.84289739d-01,
     1  2.87792167d-01, 2.90698875d-01, 2.93057586d-01, 2.94927384d-01,
     1  2.96374262d-01, 2.97466488d-01, 2.98270390d-01, 2.98847025d-01,
     1  2.99249945d-01, 2.99524080d-01, 2.99705593d-01, 2.99822450d-01/
      data u2e /
     1  2.99895431d-01, 2.99939301d-01, 2.99963931d-01, 2.99975129d-01,
     1  2.99974996d-01, 2.99961526d-01, 2.99927041d-01, 2.99854809d-01,
     1  2.99712769d-01, 2.99442742d-01, 2.98942676d-01, 2.98038511d-01,
     1  2.96441259d-01, 2.93684573d-01, 2.89040478d-01, 2.81421884d-01/
      data u2f /
     1  2.69315148d-01, 2.50874185d-01, 2.24457680d-01, 1.89885662d-01,
     1  1.49894358d-01, 1.09927672d-01, 7.54041273d-02, 4.90259517d-02,
     1  3.06080023d-02, 1.85165524d-02, 1.09104125d-02, 6.27726960d-03,
     1  3.53002680d-03, 1.94049735d-03, 1.04218859d-03, 5.45964314d-04/
      data u2g / 2.77379128d-04, 1.33343739d-04, 5.32660444d-05/
      data u3a /
     1  1.86765383d-10, 1.96772458d-09, 1.19111389d-08, 5.54964761d-08,
     1  2.18340713d-07, 7.55899524d-07, 2.35604385d-06, 6.70801745d-06,
     1  1.76224112d-05, 4.30351929d-05, 9.82592148d-05, 2.10736217d-04,
     1  4.26209304d-04, 8.15657041d-04, 1.48160943d-03, 2.56186555d-03/
      data u3b /
     1  4.22851247d-03, 6.68078970d-03, 1.01317466d-02, 1.47903961d-02,
     1  2.08424987d-02, 2.84336008d-02, 3.76573037d-02, 4.85502549d-02,
     1  6.10936693d-02, 7.52198901d-02, 9.08218891d-02, 1.07763660d-01,
     1  1.25889931d-01, 1.45034247d-01, 1.65025016d-01, 1.85689556d-01/
      data u3c /
     1  2.06856371d-01, 2.28356037d-01, 2.50021072d-01, 2.71685149d-01,
     1  2.93181998d-01, 3.14344301d-01, 3.35002907d-01, 3.54986687d-01,
     1  3.74123404d-01, 3.92241969d-01, 4.09176451d-01, 4.24772089d-01,
     1  4.38893320d-01, 4.51433444d-01, 4.62324969d-01, 4.71549073d-01/
      data u3d /
     1  4.79142163d-01, 4.85197409d-01, 4.89859810d-01, 4.93314543d-01,
     1  4.95770115d-01, 4.97439231d-01, 4.98520996d-01, 4.99187563d-01,
     1  4.99576941d-01, 4.99791928d-01, 4.99903753d-01, 4.99958343d-01,
     1  4.99983239d-01, 4.99993785d-01, 4.99997902d-01, 4.99999367d-01/
      data u3e /
     1  4.99999835d-01, 4.99999965d-01, 4.99999995d-01, 5.00000000d-01,
     1  5.00000000d-01, 4.99999997d-01, 4.99999976d-01, 4.99999863d-01,
     1  4.99999315d-01, 4.99996914d-01, 4.99987300d-01, 4.99951740d-01,
     1  4.99829328d-01, 4.99435130d-01, 4.98245007d-01, 4.94883400d-01/
      data u3f /
     1  4.86081966d-01, 4.65174923d-01, 4.21856650d-01, 3.47885738d-01,
     1  2.49649938d-01, 1.51648615d-01, 7.80173239d-02, 3.47983164d-02,
     1  1.38686441d-02, 5.05765688d-03, 1.71052539d-03, 5.38966324d-04,
     1  1.57923694d-04, 4.27352191d-05, 1.05512005d-05, 2.33068621d-06/
      data u3g / 4.45404604d-07, 6.88336884d-08, 7.23875975d-09/
      data zero/0.0d0/
c
      abermx = zero
      do 10 k = 1,99
        ae1 = abs(y(1,k) - u1(k))
        ae2 = abs(y(2,k) - u2(k))
        ae3 = abs(y(3,k) - u3(k))
        abermx = max(abermx, ae1, ae2, ae3)
 10     continue
c
      return
c end of subroutine maxerr
      end

................................................................................


                     Demonstration Problem for DLSOIBT

          Galerkin method solution of system of 3 PDEs:

            u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)     (i=1,2,3)
                t                           x                xx

          x interval is -1 to 1,  zero boundary conditions
          x discretized using piecewise linear basis functions

          Fixed parameters are as follows:
             Diffusion coefficients are eta =  0.10E+00  0.20E-01  0.10E-01
             t0 =  0.00000E+00
             tlast =  0.40000E+00
             Uniform mesh, number of intervals = 100
             Block size mb = 3
             Number of blocks nb =  99
             ODE system size neq =  297



Initial profiles:

 Values of PDE component i =  1
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.1000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00
  0.2000E+00  0.2000E+00  0.2000E+00  0.2000E+00  0.1000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00
 Values of PDE component i =  2
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.1500E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.1500E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00
 Values of PDE component i =  3
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.2500E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.2500E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
  0.0000E+00




******************************************************************************************

Run with rtol =  0.1E-02  atol =  0.1E-02   mf = 11



 At time t = 0.10000E+00  current h = 0.22309E-01  current order = 2  current nst =   17

 At time t = 0.20000E+00  current h = 0.41111E-01  current order = 2  current nst =   20

 At time t = 0.30000E+00  current h = 0.51073E-01  current order = 2  current nst =   22

 At time t = 0.40000E+00  current h = 0.58735E-01  current order = 2  current nst =   24

 Values of PDE component i =  1
  0.1714E-02  0.3444E-02  0.5203E-02  0.7006E-02  0.8866E-02  0.1079E-01  0.1280E-01
  0.1490E-01  0.1710E-01  0.1940E-01  0.2181E-01  0.2433E-01  0.2695E-01  0.2968E-01
  0.3252E-01  0.3546E-01  0.3848E-01  0.4161E-01  0.4476E-01  0.4807E-01  0.5130E-01
  0.5475E-01  0.5807E-01  0.6152E-01  0.6497E-01  0.6846E-01  0.7191E-01  0.7531E-01
  0.7899E-01  0.8240E-01  0.8592E-01  0.8937E-01  0.9290E-01  0.9639E-01  0.9989E-01
  0.1034E+00  0.1069E+00  0.1104E+00  0.1139E+00  0.1174E+00  0.1209E+00  0.1244E+00
  0.1279E+00  0.1313E+00  0.1348E+00  0.1383E+00  0.1417E+00  0.1451E+00  0.1484E+00
  0.1517E+00  0.1549E+00  0.1580E+00  0.1610E+00  0.1639E+00  0.1666E+00  0.1692E+00
  0.1717E+00  0.1739E+00  0.1760E+00  0.1778E+00  0.1794E+00  0.1808E+00  0.1819E+00
  0.1828E+00  0.1834E+00  0.1836E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1814E+00
  0.1799E+00  0.1779E+00  0.1755E+00  0.1725E+00  0.1691E+00  0.1650E+00  0.1603E+00
  0.1548E+00  0.1488E+00  0.1424E+00  0.1348E+00  0.1266E+00  0.1176E+00  0.1082E+00
  0.9850E-01  0.8871E-01  0.7919E-01  0.7011E-01  0.6161E-01  0.5373E-01  0.4646E-01
  0.3977E-01  0.3361E-01  0.2792E-01  0.2263E-01  0.1770E-01  0.1303E-01  0.8572E-02
  0.4251E-02
 Values of PDE component i =  2
  0.7744E-05  0.1814E-04  0.3450E-04  0.6154E-04  0.1063E-03  0.1791E-03  0.2948E-03
  0.4744E-03  0.7460E-03  0.1146E-02  0.1721E-02  0.2524E-02  0.3617E-02  0.5068E-02
  0.6946E-02  0.9315E-02  0.1223E-01  0.1575E-01  0.1989E-01  0.2467E-01  0.3009E-01
  0.3613E-01  0.4275E-01  0.4991E-01  0.5755E-01  0.6565E-01  0.7409E-01  0.8288E-01
  0.9193E-01  0.1012E+00  0.1106E+00  0.1202E+00  0.1299E+00  0.1396E+00  0.1494E+00
  0.1592E+00  0.1690E+00  0.1787E+00  0.1885E+00  0.1981E+00  0.2076E+00  0.2170E+00
  0.2260E+00  0.2348E+00  0.2431E+00  0.2510E+00  0.2582E+00  0.2649E+00  0.2709E+00
  0.2762E+00  0.2808E+00  0.2847E+00  0.2880E+00  0.2907E+00  0.2929E+00  0.2946E+00
  0.2960E+00  0.2970E+00  0.2978E+00  0.2984E+00  0.2989E+00  0.2992E+00  0.2994E+00
  0.2996E+00  0.2997E+00  0.2998E+00  0.2999E+00  0.2999E+00  0.2999E+00  0.2999E+00
  0.2999E+00  0.2998E+00  0.2997E+00  0.2995E+00  0.2990E+00  0.2983E+00  0.2967E+00
  0.2937E+00  0.2890E+00  0.2817E+00  0.2700E+00  0.2505E+00  0.2224E+00  0.1867E+00
  0.1467E+00  0.1076E+00  0.7410E-01  0.4856E-01  0.3063E-01  0.1876E-01  0.1121E-01
  0.6560E-02  0.3762E-02  0.2117E-02  0.1169E-02  0.6331E-03  0.3338E-03  0.1667E-03
  0.6863E-04
 Values of PDE component i =  3
  0.9946E-07  0.2111E-06  0.3530E-06  0.5634E-06  0.9413E-06  0.1748E-05  0.3651E-05
  0.8255E-05  0.1918E-04  0.4406E-04  0.9777E-04  0.2074E-03  0.4186E-03  0.8031E-03
  0.1465E-02  0.2544E-02  0.4215E-02  0.6676E-02  0.1014E-01  0.1480E-01  0.2085E-01
  0.2842E-01  0.3759E-01  0.4839E-01  0.6081E-01  0.7478E-01  0.9019E-01  0.1069E+00
  0.1248E+00  0.1437E+00  0.1634E+00  0.1839E+00  0.2049E+00  0.2264E+00  0.2481E+00
  0.2701E+00  0.2921E+00  0.3140E+00  0.3356E+00  0.3565E+00  0.3766E+00  0.3955E+00
  0.4129E+00  0.4286E+00  0.4425E+00  0.4545E+00  0.4646E+00  0.4730E+00  0.4798E+00
  0.4851E+00  0.4892E+00  0.4923E+00  0.4946E+00  0.4962E+00  0.4974E+00  0.4983E+00
  0.4988E+00  0.4992E+00  0.4995E+00  0.4997E+00  0.4998E+00  0.4999E+00  0.4999E+00
  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5001E+00  0.5002E+00
  0.4993E+00  0.4986E+00  0.4997E+00  0.4870E+00  0.4577E+00  0.4177E+00  0.3451E+00
  0.2443E+00  0.1485E+00  0.7868E-01  0.3717E-01  0.1608E-01  0.6517E-02  0.2512E-02
  0.9293E-03  0.3318E-03  0.1145E-03  0.3809E-04  0.1208E-04  0.3554E-05  0.9159E-06
  0.1828E-06


Final statistics for mf = 11:   24 steps,    38 res,     9 jacobians,
                              rwork size =  7447,  iwork size =   317
Final output is correct to within  0.75E+01  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-02  atol =  0.1E-02   mf = 12



 At time t = 0.10000E+00  current h = 0.22305E-01  current order = 2  current nst =   17

 At time t = 0.20000E+00  current h = 0.41173E-01  current order = 2  current nst =   20

 At time t = 0.30000E+00  current h = 0.51349E-01  current order = 2  current nst =   22

 At time t = 0.40000E+00  current h = 0.59258E-01  current order = 2  current nst =   24

 Values of PDE component i =  1
  0.1712E-02  0.3440E-02  0.5197E-02  0.6997E-02  0.8855E-02  0.1078E-01  0.1279E-01
  0.1488E-01  0.1708E-01  0.1937E-01  0.2178E-01  0.2429E-01  0.2691E-01  0.2963E-01
  0.3247E-01  0.3540E-01  0.3841E-01  0.4153E-01  0.4468E-01  0.4797E-01  0.5119E-01
  0.5464E-01  0.5793E-01  0.6141E-01  0.6483E-01  0.6828E-01  0.7182E-01  0.7520E-01
  0.7881E-01  0.8223E-01  0.8579E-01  0.8925E-01  0.9278E-01  0.9628E-01  0.9978E-01
  0.1033E+00  0.1068E+00  0.1103E+00  0.1138E+00  0.1173E+00  0.1208E+00  0.1243E+00
  0.1278E+00  0.1313E+00  0.1348E+00  0.1383E+00  0.1417E+00  0.1451E+00  0.1485E+00
  0.1517E+00  0.1549E+00  0.1580E+00  0.1611E+00  0.1639E+00  0.1667E+00  0.1693E+00
  0.1717E+00  0.1739E+00  0.1760E+00  0.1778E+00  0.1794E+00  0.1808E+00  0.1819E+00
  0.1828E+00  0.1834E+00  0.1837E+00  0.1836E+00  0.1833E+00  0.1826E+00  0.1814E+00
  0.1800E+00  0.1779E+00  0.1756E+00  0.1726E+00  0.1692E+00  0.1652E+00  0.1605E+00
  0.1553E+00  0.1493E+00  0.1427E+00  0.1352E+00  0.1271E+00  0.1182E+00  0.1088E+00
  0.9905E-01  0.8926E-01  0.7972E-01  0.7062E-01  0.6210E-01  0.5418E-01  0.4687E-01
  0.4013E-01  0.3392E-01  0.2818E-01  0.2285E-01  0.1786E-01  0.1315E-01  0.8653E-02
  0.4292E-02
 Values of PDE component i =  2
  0.7637E-05  0.1793E-04  0.3419E-04  0.6115E-04  0.1058E-03  0.1786E-03  0.2943E-03
  0.4739E-03  0.7456E-03  0.1146E-02  0.1720E-02  0.2523E-02  0.3617E-02  0.5068E-02
  0.6944E-02  0.9312E-02  0.1223E-01  0.1574E-01  0.1988E-01  0.2465E-01  0.3006E-01
  0.3608E-01  0.4269E-01  0.4983E-01  0.5745E-01  0.6552E-01  0.7394E-01  0.8271E-01
  0.9173E-01  0.1010E+00  0.1104E+00  0.1200E+00  0.1297E+00  0.1395E+00  0.1493E+00
  0.1591E+00  0.1690E+00  0.1788E+00  0.1885E+00  0.1982E+00  0.2077E+00  0.2171E+00
  0.2261E+00  0.2349E+00  0.2432E+00  0.2510E+00  0.2583E+00  0.2650E+00  0.2710E+00
  0.2763E+00  0.2808E+00  0.2848E+00  0.2880E+00  0.2907E+00  0.2929E+00  0.2946E+00
  0.2960E+00  0.2970E+00  0.2978E+00  0.2984E+00  0.2989E+00  0.2992E+00  0.2994E+00
  0.2996E+00  0.2997E+00  0.2998E+00  0.2999E+00  0.2999E+00  0.2999E+00  0.2999E+00
  0.2999E+00  0.2998E+00  0.2997E+00  0.2995E+00  0.2990E+00  0.2982E+00  0.2967E+00
  0.2939E+00  0.2892E+00  0.2818E+00  0.2698E+00  0.2505E+00  0.2233E+00  0.1883E+00
  0.1484E+00  0.1090E+00  0.7518E-01  0.4929E-01  0.3109E-01  0.1904E-01  0.1137E-01
  0.6649E-02  0.3811E-02  0.2144E-02  0.1184E-02  0.6409E-03  0.3380E-03  0.1688E-03
  0.6954E-04
 Values of PDE component i =  3
  0.4971E-08  0.1654E-07  0.4938E-07  0.1430E-06  0.4036E-06  0.1108E-05  0.2946E-05
  0.7554E-05  0.1859E-04  0.4372E-04  0.9783E-04  0.2080E-03  0.4199E-03  0.8052E-03
  0.1468E-02  0.2548E-02  0.4220E-02  0.6682E-02  0.1015E-01  0.1482E-01  0.2087E-01
  0.2845E-01  0.3764E-01  0.4847E-01  0.6093E-01  0.7494E-01  0.9039E-01  0.1072E+00
  0.1251E+00  0.1440E+00  0.1638E+00  0.1843E+00  0.2053E+00  0.2268E+00  0.2485E+00
  0.2704E+00  0.2924E+00  0.3143E+00  0.3358E+00  0.3567E+00  0.3768E+00  0.3956E+00
  0.4130E+00  0.4287E+00  0.4426E+00  0.4546E+00  0.4647E+00  0.4731E+00  0.4798E+00
  0.4851E+00  0.4892E+00  0.4923E+00  0.4946E+00  0.4963E+00  0.4974E+00  0.4983E+00
  0.4988E+00  0.4992E+00  0.4995E+00  0.4997E+00  0.4998E+00  0.4999E+00  0.4999E+00
  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5001E+00  0.5002E+00
  0.4993E+00  0.4985E+00  0.4997E+00  0.4872E+00  0.4581E+00  0.4191E+00  0.3472E+00
  0.2462E+00  0.1500E+00  0.7966E-01  0.3769E-01  0.1632E-01  0.6615E-02  0.2551E-02
  0.9450E-03  0.3386E-03  0.1179E-03  0.4004E-04  0.1329E-04  0.4317E-05  0.1366E-05
  0.3917E-06


Final statistics for mf = 12:   24 steps,   128 res,     9 jacobians,
                              rwork size =  7447,  iwork size =   317
Final output is correct to within  0.71E+01  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-02  atol =  0.1E-02   mf = 21



 At time t = 0.10000E+00  current h = 0.17227E-01  current order = 2  current nst =   20

 At time t = 0.20000E+00  current h = 0.33122E-01  current order = 2  current nst =   24

 At time t = 0.30000E+00  current h = 0.33122E-01  current order = 2  current nst =   27

 At time t = 0.40000E+00  current h = 0.33122E-01  current order = 2  current nst =   30

 Values of PDE component i =  1
  0.1714E-02  0.3443E-02  0.5202E-02  0.7005E-02  0.8866E-02  0.1080E-01  0.1281E-01
  0.1491E-01  0.1711E-01  0.1942E-01  0.2183E-01  0.2435E-01  0.2698E-01  0.2971E-01
  0.3255E-01  0.3549E-01  0.3851E-01  0.4162E-01  0.4480E-01  0.4804E-01  0.5134E-01
  0.5468E-01  0.5807E-01  0.6148E-01  0.6492E-01  0.6838E-01  0.7185E-01  0.7533E-01
  0.7882E-01  0.8230E-01  0.8579E-01  0.8928E-01  0.9277E-01  0.9626E-01  0.9975E-01
  0.1032E+00  0.1067E+00  0.1102E+00  0.1137E+00  0.1173E+00  0.1208E+00  0.1243E+00
  0.1278E+00  0.1314E+00  0.1349E+00  0.1384E+00  0.1419E+00  0.1453E+00  0.1487E+00
  0.1520E+00  0.1552E+00  0.1583E+00  0.1613E+00  0.1642E+00  0.1670E+00  0.1696E+00
  0.1720E+00  0.1742E+00  0.1762E+00  0.1780E+00  0.1796E+00  0.1809E+00  0.1820E+00
  0.1829E+00  0.1834E+00  0.1837E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1813E+00
  0.1798E+00  0.1779E+00  0.1755E+00  0.1726E+00  0.1692E+00  0.1652E+00  0.1606E+00
  0.1553E+00  0.1493E+00  0.1426E+00  0.1351E+00  0.1269E+00  0.1180E+00  0.1085E+00
  0.9878E-01  0.8900E-01  0.7947E-01  0.7039E-01  0.6188E-01  0.5397E-01  0.4668E-01
  0.3996E-01  0.3376E-01  0.2804E-01  0.2273E-01  0.1777E-01  0.1308E-01  0.8605E-02
  0.4268E-02
 Values of PDE component i =  2
  0.8292E-05  0.1926E-04  0.3621E-04  0.6381E-04  0.1089E-03  0.1817E-03  0.2968E-03
  0.4749E-03  0.7439E-03  0.1141E-02  0.1711E-02  0.2511E-02  0.3603E-02  0.5056E-02
  0.6939E-02  0.9318E-02  0.1225E-01  0.1578E-01  0.1993E-01  0.2471E-01  0.3012E-01
  0.3614E-01  0.4273E-01  0.4984E-01  0.5743E-01  0.6546E-01  0.7385E-01  0.8257E-01
  0.9156E-01  0.1008E+00  0.1102E+00  0.1197E+00  0.1294E+00  0.1392E+00  0.1491E+00
  0.1591E+00  0.1691E+00  0.1791E+00  0.1891E+00  0.1990E+00  0.2088E+00  0.2183E+00
  0.2274E+00  0.2362E+00  0.2444E+00  0.2520E+00  0.2591E+00  0.2655E+00  0.2712E+00
  0.2762E+00  0.2806E+00  0.2844E+00  0.2875E+00  0.2902E+00  0.2924E+00  0.2941E+00
  0.2955E+00  0.2966E+00  0.2975E+00  0.2981E+00  0.2986E+00  0.2990E+00  0.2993E+00
  0.2995E+00  0.2996E+00  0.2998E+00  0.2998E+00  0.2999E+00  0.2999E+00  0.2999E+00
  0.2999E+00  0.2999E+00  0.2998E+00  0.2995E+00  0.2992E+00  0.2985E+00  0.2972E+00
  0.2948E+00  0.2904E+00  0.2820E+00  0.2686E+00  0.2491E+00  0.2219E+00  0.1869E+00
  0.1474E+00  0.1085E+00  0.7511E-01  0.4945E-01  0.3134E-01  0.1929E-01  0.1159E-01
  0.6825E-02  0.3944E-02  0.2239E-02  0.1249E-02  0.6833E-03  0.3643E-03  0.1837E-03
  0.7624E-04
 Values of PDE component i =  3
  0.6162E-07  0.1377E-06  0.2513E-06  0.4520E-06  0.8610E-06  0.1779E-05  0.3932E-05
  0.8985E-05  0.2057E-04  0.4616E-04  0.1003E-03  0.2094E-03  0.4185E-03  0.7992E-03
  0.1457E-02  0.2534E-02  0.4208E-02  0.6681E-02  0.1016E-01  0.1485E-01  0.2091E-01
  0.2846E-01  0.3758E-01  0.4830E-01  0.6060E-01  0.7441E-01  0.8966E-01  0.1062E+00
  0.1240E+00  0.1428E+00  0.1626E+00  0.1831E+00  0.2045E+00  0.2264E+00  0.2488E+00
  0.2714E+00  0.2941E+00  0.3166E+00  0.3384E+00  0.3594E+00  0.3791E+00  0.3974E+00
  0.4141E+00  0.4290E+00  0.4422E+00  0.4536E+00  0.4633E+00  0.4714E+00  0.4781E+00
  0.4834E+00  0.4877E+00  0.4910E+00  0.4935E+00  0.4954E+00  0.4968E+00  0.4978E+00
  0.4985E+00  0.4990E+00  0.4993E+00  0.4996E+00  0.4997E+00  0.4998E+00  0.4999E+00
  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5001E+00  0.5002E+00  0.5007E+00
  0.5022E+00  0.5016E+00  0.4950E+00  0.4878E+00  0.4672E+00  0.4132E+00  0.3339E+00
  0.2415E+00  0.1512E+00  0.8188E-01  0.3937E-01  0.1726E-01  0.7067E-02  0.2745E-02
  0.1022E-02  0.3663E-03  0.1269E-03  0.4241E-04  0.1362E-04  0.4154E-05  0.1175E-05
  0.2869E-06


Final statistics for mf = 21:   30 steps,    45 res,    11 jacobians,
                              rwork size =  5368,  iwork size =   317
Final output is correct to within  0.14E+02  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-02  atol =  0.1E-02   mf = 22



 At time t = 0.10000E+00  current h = 0.17219E-01  current order = 2  current nst =   20

 At time t = 0.20000E+00  current h = 0.33216E-01  current order = 2  current nst =   24

 At time t = 0.30000E+00  current h = 0.33216E-01  current order = 2  current nst =   27

 At time t = 0.40000E+00  current h = 0.33216E-01  current order = 2  current nst =   30

 Values of PDE component i =  1
  0.1713E-02  0.3441E-02  0.5199E-02  0.7001E-02  0.8861E-02  0.1079E-01  0.1280E-01
  0.1490E-01  0.1710E-01  0.1940E-01  0.2181E-01  0.2433E-01  0.2695E-01  0.2968E-01
  0.3252E-01  0.3545E-01  0.3846E-01  0.4156E-01  0.4473E-01  0.4797E-01  0.5126E-01
  0.5460E-01  0.5797E-01  0.6138E-01  0.6481E-01  0.6826E-01  0.7172E-01  0.7520E-01
  0.7868E-01  0.8216E-01  0.8565E-01  0.8915E-01  0.9264E-01  0.9614E-01  0.9964E-01
  0.1031E+00  0.1066E+00  0.1102E+00  0.1137E+00  0.1172E+00  0.1207E+00  0.1243E+00
  0.1278E+00  0.1314E+00  0.1349E+00  0.1384E+00  0.1419E+00  0.1453E+00  0.1487E+00
  0.1520E+00  0.1552E+00  0.1584E+00  0.1614E+00  0.1643E+00  0.1670E+00  0.1696E+00
  0.1720E+00  0.1742E+00  0.1762E+00  0.1780E+00  0.1796E+00  0.1810E+00  0.1820E+00
  0.1829E+00  0.1834E+00  0.1837E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1814E+00
  0.1798E+00  0.1779E+00  0.1755E+00  0.1726E+00  0.1691E+00  0.1651E+00  0.1605E+00
  0.1552E+00  0.1493E+00  0.1426E+00  0.1351E+00  0.1270E+00  0.1181E+00  0.1087E+00
  0.9897E-01  0.8921E-01  0.7968E-01  0.7060E-01  0.6208E-01  0.5416E-01  0.4685E-01
  0.4011E-01  0.3390E-01  0.2816E-01  0.2283E-01  0.1785E-01  0.1314E-01  0.8644E-02
  0.4287E-02
 Values of PDE component i =  2
  0.8220E-05  0.1912E-04  0.3600E-04  0.6353E-04  0.1086E-03  0.1813E-03  0.2965E-03
  0.4746E-03  0.7437E-03  0.1141E-02  0.1711E-02  0.2511E-02  0.3603E-02  0.5056E-02
  0.6939E-02  0.9317E-02  0.1225E-01  0.1577E-01  0.1992E-01  0.2470E-01  0.3010E-01
  0.3611E-01  0.4268E-01  0.4978E-01  0.5735E-01  0.6535E-01  0.7372E-01  0.8242E-01
  0.9140E-01  0.1006E+00  0.1100E+00  0.1196E+00  0.1293E+00  0.1391E+00  0.1490E+00
  0.1590E+00  0.1691E+00  0.1792E+00  0.1892E+00  0.1991E+00  0.2088E+00  0.2183E+00
  0.2275E+00  0.2362E+00  0.2444E+00  0.2521E+00  0.2591E+00  0.2655E+00  0.2712E+00
  0.2763E+00  0.2806E+00  0.2844E+00  0.2876E+00  0.2902E+00  0.2924E+00  0.2941E+00
  0.2955E+00  0.2966E+00  0.2975E+00  0.2981E+00  0.2986E+00  0.2990E+00  0.2993E+00
  0.2995E+00  0.2997E+00  0.2998E+00  0.2998E+00  0.2999E+00  0.2999E+00  0.2999E+00
  0.2999E+00  0.2999E+00  0.2998E+00  0.2995E+00  0.2992E+00  0.2984E+00  0.2971E+00
  0.2946E+00  0.2900E+00  0.2819E+00  0.2688E+00  0.2494E+00  0.2221E+00  0.1872E+00
  0.1478E+00  0.1089E+00  0.7549E-01  0.4976E-01  0.3156E-01  0.1943E-01  0.1168E-01
  0.6879E-02  0.3975E-02  0.2257E-02  0.1259E-02  0.6888E-03  0.3673E-03  0.1853E-03
  0.7689E-04
 Values of PDE component i =  3
  0.6430E-08  0.2134E-07  0.6312E-07  0.1798E-06  0.4960E-06  0.1323E-05  0.3404E-05
  0.8433E-05  0.2008E-04  0.4585E-04  0.1003E-03  0.2099E-03  0.4196E-03  0.8010E-03
  0.1459E-02  0.2537E-02  0.4212E-02  0.6685E-02  0.1017E-01  0.1486E-01  0.2092E-01
  0.2848E-01  0.3762E-01  0.4836E-01  0.6068E-01  0.7453E-01  0.8981E-01  0.1064E+00
  0.1242E+00  0.1431E+00  0.1628E+00  0.1835E+00  0.2048E+00  0.2267E+00  0.2490E+00
  0.2716E+00  0.2943E+00  0.3167E+00  0.3385E+00  0.3595E+00  0.3792E+00  0.3975E+00
  0.4141E+00  0.4290E+00  0.4422E+00  0.4536E+00  0.4633E+00  0.4714E+00  0.4780E+00
  0.4834E+00  0.4877E+00  0.4910E+00  0.4935E+00  0.4954E+00  0.4968E+00  0.4978E+00
  0.4985E+00  0.4990E+00  0.4993E+00  0.4996E+00  0.4997E+00  0.4998E+00  0.4999E+00
  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5001E+00  0.5002E+00  0.5008E+00
  0.5024E+00  0.5018E+00  0.4949E+00  0.4881E+00  0.4679E+00  0.4136E+00  0.3344E+00
  0.2425E+00  0.1521E+00  0.8249E-01  0.3969E-01  0.1741E-01  0.7129E-02  0.2769E-02
  0.1031E-02  0.3705E-03  0.1289E-03  0.4358E-04  0.1433E-04  0.4592E-05  0.1427E-05
  0.4012E-06


Final statistics for mf = 22:   30 steps,   155 res,    11 jacobians,
                              rwork size =  5368,  iwork size =   317
Final output is correct to within  0.14E+02  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-05  atol =  0.1E-05   mf = 11



 At time t = 0.10000E+00  current h = 0.79980E-02  current order = 4  current nst =   60

 At time t = 0.20000E+00  current h = 0.10116E-01  current order = 4  current nst =   71

 At time t = 0.30000E+00  current h = 0.12090E-01  current order = 4  current nst =   89

 At time t = 0.40000E+00  current h = 0.12812E-01  current order = 4  current nst =  106

 Values of PDE component i =  1
  0.1710E-02  0.3434E-02  0.5188E-02  0.6985E-02  0.8839E-02  0.1076E-01  0.1277E-01
  0.1486E-01  0.1705E-01  0.1934E-01  0.2174E-01  0.2425E-01  0.2687E-01  0.2960E-01
  0.3243E-01  0.3536E-01  0.3838E-01  0.4149E-01  0.4468E-01  0.4793E-01  0.5124E-01
  0.5460E-01  0.5799E-01  0.6143E-01  0.6489E-01  0.6836E-01  0.7186E-01  0.7536E-01
  0.7887E-01  0.8239E-01  0.8591E-01  0.8943E-01  0.9295E-01  0.9647E-01  0.9999E-01
  0.1035E+00  0.1070E+00  0.1106E+00  0.1141E+00  0.1176E+00  0.1211E+00  0.1246E+00
  0.1281E+00  0.1316E+00  0.1350E+00  0.1384E+00  0.1418E+00  0.1452E+00  0.1485E+00
  0.1517E+00  0.1549E+00  0.1580E+00  0.1609E+00  0.1638E+00  0.1665E+00  0.1691E+00
  0.1715E+00  0.1737E+00  0.1758E+00  0.1776E+00  0.1793E+00  0.1806E+00  0.1818E+00
  0.1827E+00  0.1833E+00  0.1836E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1814E+00
  0.1799E+00  0.1780E+00  0.1756E+00  0.1727E+00  0.1693E+00  0.1653E+00  0.1606E+00
  0.1554E+00  0.1494E+00  0.1427E+00  0.1353E+00  0.1271E+00  0.1183E+00  0.1089E+00
  0.9916E-01  0.8935E-01  0.7978E-01  0.7067E-01  0.6212E-01  0.5420E-01  0.4688E-01
  0.4015E-01  0.3394E-01  0.2820E-01  0.2286E-01  0.1788E-01  0.1316E-01  0.8659E-02
  0.4295E-02
 Values of PDE component i =  2
  0.7174E-05  0.1708E-04  0.3312E-04  0.6016E-04  0.1053E-03  0.1792E-03  0.2967E-03
  0.4789E-03  0.7536E-03  0.1157E-02  0.1734E-02  0.2538E-02  0.3631E-02  0.5078E-02
  0.6948E-02  0.9306E-02  0.1221E-01  0.1572E-01  0.1985E-01  0.2462E-01  0.3004E-01
  0.3608E-01  0.4271E-01  0.4988E-01  0.5755E-01  0.6566E-01  0.7415E-01  0.8298E-01
  0.9208E-01  0.1014E+00  0.1109E+00  0.1206E+00  0.1303E+00  0.1401E+00  0.1500E+00
  0.1599E+00  0.1697E+00  0.1794E+00  0.1891E+00  0.1986E+00  0.2079E+00  0.2170E+00
  0.2258E+00  0.2343E+00  0.2425E+00  0.2501E+00  0.2573E+00  0.2639E+00  0.2700E+00
  0.2754E+00  0.2802E+00  0.2843E+00  0.2878E+00  0.2907E+00  0.2931E+00  0.2949E+00
  0.2964E+00  0.2975E+00  0.2983E+00  0.2988E+00  0.2992E+00  0.2995E+00  0.2997E+00
  0.2998E+00  0.2999E+00  0.2999E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.2999E+00  0.2999E+00  0.2997E+00  0.2994E+00  0.2989E+00  0.2980E+00  0.2964E+00
  0.2937E+00  0.2890E+00  0.2814E+00  0.2693E+00  0.2509E+00  0.2245E+00  0.1899E+00
  0.1499E+00  0.1099E+00  0.7540E-01  0.4903E-01  0.3061E-01  0.1852E-01  0.1091E-01
  0.6277E-02  0.3530E-02  0.1941E-02  0.1042E-02  0.5461E-03  0.2774E-03  0.1334E-03
  0.5328E-04
 Values of PDE component i =  3
  0.1970E-09  0.1996E-08  0.1198E-07  0.5563E-07  0.2186E-06  0.7563E-06  0.2357E-05
  0.6709E-05  0.1762E-04  0.4304E-04  0.9826E-04  0.2107E-03  0.4262E-03  0.8157E-03
  0.1482E-02  0.2562E-02  0.4229E-02  0.6681E-02  0.1013E-01  0.1479E-01  0.2084E-01
  0.2843E-01  0.3766E-01  0.4855E-01  0.6109E-01  0.7522E-01  0.9082E-01  0.1078E+00
  0.1259E+00  0.1450E+00  0.1650E+00  0.1857E+00  0.2069E+00  0.2284E+00  0.2500E+00
  0.2717E+00  0.2932E+00  0.3143E+00  0.3350E+00  0.3550E+00  0.3741E+00  0.3922E+00
  0.4092E+00  0.4248E+00  0.4389E+00  0.4514E+00  0.4623E+00  0.4715E+00  0.4791E+00
  0.4852E+00  0.4899E+00  0.4933E+00  0.4958E+00  0.4974E+00  0.4985E+00  0.4992E+00
  0.4996E+00  0.4998E+00  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.4998E+00
  0.4994E+00  0.4982E+00  0.4949E+00  0.4861E+00  0.4652E+00  0.4219E+00  0.3479E+00
  0.2496E+00  0.1516E+00  0.7801E-01  0.3480E-01  0.1387E-01  0.5059E-02  0.1712E-02
  0.5396E-03  0.1582E-03  0.4287E-04  0.1061E-04  0.2351E-05  0.4522E-06  0.7093E-07
  0.7792E-08


Final statistics for mf = 11:  106 steps,   148 res,    23 jacobians,
                              rwork size =  7447,  iwork size =   317
Final output is correct to within  0.67E+01  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-05  atol =  0.1E-05   mf = 12



 At time t = 0.10000E+00  current h = 0.64185E-02  current order = 3  current nst =   61

 At time t = 0.20000E+00  current h = 0.13037E-01  current order = 4  current nst =   71

 At time t = 0.30000E+00  current h = 0.68710E-02  current order = 3  current nst =   89

 At time t = 0.40000E+00  current h = 0.84070E-02  current order = 3  current nst =  100

 Values of PDE component i =  1
  0.1710E-02  0.3434E-02  0.5188E-02  0.6985E-02  0.8839E-02  0.1076E-01  0.1276E-01
  0.1486E-01  0.1705E-01  0.1934E-01  0.2174E-01  0.2425E-01  0.2687E-01  0.2960E-01
  0.3243E-01  0.3536E-01  0.3838E-01  0.4149E-01  0.4467E-01  0.4793E-01  0.5124E-01
  0.5460E-01  0.5800E-01  0.6142E-01  0.6489E-01  0.6836E-01  0.7186E-01  0.7535E-01
  0.7888E-01  0.8238E-01  0.8591E-01  0.8942E-01  0.9295E-01  0.9646E-01  0.1000E+00
  0.1035E+00  0.1070E+00  0.1105E+00  0.1141E+00  0.1176E+00  0.1211E+00  0.1246E+00
  0.1281E+00  0.1316E+00  0.1350E+00  0.1384E+00  0.1418E+00  0.1452E+00  0.1485E+00
  0.1517E+00  0.1549E+00  0.1580E+00  0.1609E+00  0.1638E+00  0.1665E+00  0.1691E+00
  0.1715E+00  0.1737E+00  0.1758E+00  0.1776E+00  0.1793E+00  0.1806E+00  0.1818E+00
  0.1827E+00  0.1833E+00  0.1836E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1814E+00
  0.1799E+00  0.1780E+00  0.1756E+00  0.1727E+00  0.1693E+00  0.1653E+00  0.1606E+00
  0.1554E+00  0.1494E+00  0.1427E+00  0.1353E+00  0.1271E+00  0.1183E+00  0.1089E+00
  0.9916E-01  0.8935E-01  0.7978E-01  0.7067E-01  0.6212E-01  0.5420E-01  0.4688E-01
  0.4015E-01  0.3394E-01  0.2820E-01  0.2286E-01  0.1788E-01  0.1316E-01  0.8659E-02
  0.4295E-02
 Values of PDE component i =  2
  0.7174E-05  0.1708E-04  0.3312E-04  0.6016E-04  0.1053E-03  0.1792E-03  0.2967E-03
  0.4789E-03  0.7536E-03  0.1157E-02  0.1734E-02  0.2538E-02  0.3631E-02  0.5078E-02
  0.6948E-02  0.9306E-02  0.1221E-01  0.1572E-01  0.1985E-01  0.2462E-01  0.3004E-01
  0.3608E-01  0.4271E-01  0.4988E-01  0.5755E-01  0.6566E-01  0.7415E-01  0.8298E-01
  0.9208E-01  0.1014E+00  0.1109E+00  0.1206E+00  0.1303E+00  0.1401E+00  0.1500E+00
  0.1599E+00  0.1697E+00  0.1794E+00  0.1891E+00  0.1986E+00  0.2079E+00  0.2170E+00
  0.2258E+00  0.2343E+00  0.2425E+00  0.2501E+00  0.2573E+00  0.2639E+00  0.2700E+00
  0.2754E+00  0.2802E+00  0.2843E+00  0.2878E+00  0.2907E+00  0.2931E+00  0.2949E+00
  0.2964E+00  0.2975E+00  0.2983E+00  0.2988E+00  0.2992E+00  0.2995E+00  0.2997E+00
  0.2998E+00  0.2999E+00  0.2999E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.2999E+00  0.2999E+00  0.2997E+00  0.2994E+00  0.2989E+00  0.2980E+00  0.2964E+00
  0.2937E+00  0.2890E+00  0.2814E+00  0.2693E+00  0.2509E+00  0.2245E+00  0.1899E+00
  0.1499E+00  0.1099E+00  0.7540E-01  0.4902E-01  0.3061E-01  0.1852E-01  0.1091E-01
  0.6278E-02  0.3530E-02  0.1941E-02  0.1042E-02  0.5461E-03  0.2775E-03  0.1334E-03
  0.5329E-04
 Values of PDE component i =  3
  0.2003E-09  0.2003E-08  0.1198E-07  0.5562E-07  0.2185E-06  0.7561E-06  0.2356E-05
  0.6708E-05  0.1762E-04  0.4303E-04  0.9826E-04  0.2107E-03  0.4262E-03  0.8157E-03
  0.1482E-02  0.2562E-02  0.4229E-02  0.6681E-02  0.1013E-01  0.1479E-01  0.2084E-01
  0.2843E-01  0.3766E-01  0.4855E-01  0.6109E-01  0.7522E-01  0.9082E-01  0.1078E+00
  0.1259E+00  0.1450E+00  0.1650E+00  0.1857E+00  0.2069E+00  0.2284E+00  0.2500E+00
  0.2717E+00  0.2932E+00  0.3143E+00  0.3350E+00  0.3550E+00  0.3741E+00  0.3922E+00
  0.4092E+00  0.4248E+00  0.4389E+00  0.4514E+00  0.4623E+00  0.4715E+00  0.4791E+00
  0.4852E+00  0.4899E+00  0.4933E+00  0.4958E+00  0.4974E+00  0.4985E+00  0.4992E+00
  0.4996E+00  0.4998E+00  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.4998E+00
  0.4994E+00  0.4982E+00  0.4949E+00  0.4861E+00  0.4652E+00  0.4219E+00  0.3479E+00
  0.2496E+00  0.1516E+00  0.7801E-01  0.3480E-01  0.1387E-01  0.5061E-02  0.1713E-02
  0.5401E-03  0.1584E-03  0.4295E-04  0.1063E-04  0.2361E-05  0.4556E-06  0.7202E-07
  0.8103E-08


Final statistics for mf = 12:  100 steps,   353 res,    22 jacobians,
                              rwork size =  7447,  iwork size =   317
Final output is correct to within  0.12E+02  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-05  atol =  0.1E-05   mf = 21



 At time t = 0.10000E+00  current h = 0.57520E-02  current order = 4  current nst =   72

 At time t = 0.20000E+00  current h = 0.92622E-02  current order = 5  current nst =   85

 At time t = 0.30000E+00  current h = 0.12274E-01  current order = 5  current nst =   95

 At time t = 0.40000E+00  current h = 0.12274E-01  current order = 5  current nst =  103

 Values of PDE component i =  1
  0.1710E-02  0.3434E-02  0.5188E-02  0.6985E-02  0.8839E-02  0.1076E-01  0.1277E-01
  0.1486E-01  0.1705E-01  0.1934E-01  0.2174E-01  0.2425E-01  0.2687E-01  0.2960E-01
  0.3243E-01  0.3536E-01  0.3838E-01  0.4149E-01  0.4468E-01  0.4793E-01  0.5124E-01
  0.5460E-01  0.5799E-01  0.6143E-01  0.6489E-01  0.6836E-01  0.7186E-01  0.7536E-01
  0.7887E-01  0.8239E-01  0.8591E-01  0.8943E-01  0.9295E-01  0.9647E-01  0.9999E-01
  0.1035E+00  0.1070E+00  0.1106E+00  0.1141E+00  0.1176E+00  0.1211E+00  0.1246E+00
  0.1281E+00  0.1316E+00  0.1350E+00  0.1384E+00  0.1418E+00  0.1452E+00  0.1485E+00
  0.1517E+00  0.1549E+00  0.1580E+00  0.1609E+00  0.1638E+00  0.1665E+00  0.1691E+00
  0.1715E+00  0.1737E+00  0.1758E+00  0.1776E+00  0.1793E+00  0.1806E+00  0.1818E+00
  0.1827E+00  0.1833E+00  0.1836E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1814E+00
  0.1799E+00  0.1780E+00  0.1756E+00  0.1727E+00  0.1693E+00  0.1653E+00  0.1606E+00
  0.1554E+00  0.1494E+00  0.1427E+00  0.1353E+00  0.1271E+00  0.1183E+00  0.1089E+00
  0.9916E-01  0.8935E-01  0.7978E-01  0.7067E-01  0.6212E-01  0.5420E-01  0.4688E-01
  0.4015E-01  0.3394E-01  0.2820E-01  0.2286E-01  0.1788E-01  0.1316E-01  0.8659E-02
  0.4295E-02
 Values of PDE component i =  2
  0.7174E-05  0.1708E-04  0.3312E-04  0.6016E-04  0.1053E-03  0.1792E-03  0.2967E-03
  0.4789E-03  0.7536E-03  0.1157E-02  0.1734E-02  0.2538E-02  0.3631E-02  0.5078E-02
  0.6948E-02  0.9306E-02  0.1221E-01  0.1572E-01  0.1985E-01  0.2462E-01  0.3004E-01
  0.3608E-01  0.4271E-01  0.4988E-01  0.5755E-01  0.6566E-01  0.7415E-01  0.8298E-01
  0.9208E-01  0.1014E+00  0.1109E+00  0.1206E+00  0.1303E+00  0.1401E+00  0.1500E+00
  0.1599E+00  0.1697E+00  0.1794E+00  0.1891E+00  0.1986E+00  0.2079E+00  0.2170E+00
  0.2258E+00  0.2343E+00  0.2425E+00  0.2501E+00  0.2573E+00  0.2639E+00  0.2700E+00
  0.2754E+00  0.2802E+00  0.2843E+00  0.2878E+00  0.2907E+00  0.2931E+00  0.2949E+00
  0.2964E+00  0.2975E+00  0.2983E+00  0.2988E+00  0.2993E+00  0.2995E+00  0.2997E+00
  0.2998E+00  0.2999E+00  0.2999E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.2999E+00  0.2999E+00  0.2997E+00  0.2994E+00  0.2989E+00  0.2980E+00  0.2964E+00
  0.2937E+00  0.2890E+00  0.2814E+00  0.2693E+00  0.2509E+00  0.2245E+00  0.1899E+00
  0.1499E+00  0.1099E+00  0.7540E-01  0.4903E-01  0.3061E-01  0.1852E-01  0.1091E-01
  0.6277E-02  0.3530E-02  0.1941E-02  0.1042E-02  0.5460E-03  0.2774E-03  0.1333E-03
  0.5327E-04
 Values of PDE component i =  3
  0.1721E-09  0.1938E-08  0.1186E-07  0.5543E-07  0.2183E-06  0.7559E-06  0.2356E-05
  0.6708E-05  0.1762E-04  0.4304E-04  0.9826E-04  0.2107E-03  0.4262E-03  0.8157E-03
  0.1482E-02  0.2562E-02  0.4228E-02  0.6681E-02  0.1013E-01  0.1479E-01  0.2084E-01
  0.2843E-01  0.3766E-01  0.4855E-01  0.6109E-01  0.7522E-01  0.9082E-01  0.1078E+00
  0.1259E+00  0.1450E+00  0.1650E+00  0.1857E+00  0.2069E+00  0.2284E+00  0.2500E+00
  0.2717E+00  0.2932E+00  0.3143E+00  0.3350E+00  0.3550E+00  0.3741E+00  0.3922E+00
  0.4092E+00  0.4248E+00  0.4389E+00  0.4514E+00  0.4623E+00  0.4715E+00  0.4791E+00
  0.4852E+00  0.4899E+00  0.4933E+00  0.4958E+00  0.4974E+00  0.4985E+00  0.4992E+00
  0.4996E+00  0.4998E+00  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.4998E+00
  0.4994E+00  0.4982E+00  0.4949E+00  0.4861E+00  0.4652E+00  0.4218E+00  0.3479E+00
  0.2497E+00  0.1517E+00  0.7801E-01  0.3480E-01  0.1387E-01  0.5057E-02  0.1711E-02
  0.5392E-03  0.1581E-03  0.4283E-04  0.1060E-04  0.2350E-05  0.4530E-06  0.7165E-07
  0.8130E-08


Final statistics for mf = 21:  103 steps,   122 res,    17 jacobians,
                              rwork size =  5368,  iwork size =   317
Final output is correct to within  0.11E+02  times local error tolerance. 




******************************************************************************************

Run with rtol =  0.1E-05  atol =  0.1E-05   mf = 22



 At time t = 0.10000E+00  current h = 0.57538E-02  current order = 4  current nst =   72

 At time t = 0.20000E+00  current h = 0.92575E-02  current order = 5  current nst =   85

 At time t = 0.30000E+00  current h = 0.12294E-01  current order = 5  current nst =   95

 At time t = 0.40000E+00  current h = 0.12294E-01  current order = 5  current nst =  103

 Values of PDE component i =  1
  0.1710E-02  0.3434E-02  0.5188E-02  0.6985E-02  0.8839E-02  0.1076E-01  0.1277E-01
  0.1486E-01  0.1705E-01  0.1934E-01  0.2174E-01  0.2425E-01  0.2687E-01  0.2960E-01
  0.3243E-01  0.3536E-01  0.3838E-01  0.4149E-01  0.4468E-01  0.4793E-01  0.5124E-01
  0.5460E-01  0.5799E-01  0.6143E-01  0.6489E-01  0.6836E-01  0.7186E-01  0.7536E-01
  0.7887E-01  0.8239E-01  0.8591E-01  0.8943E-01  0.9295E-01  0.9647E-01  0.9999E-01
  0.1035E+00  0.1070E+00  0.1106E+00  0.1141E+00  0.1176E+00  0.1211E+00  0.1246E+00
  0.1281E+00  0.1316E+00  0.1350E+00  0.1384E+00  0.1418E+00  0.1452E+00  0.1485E+00
  0.1517E+00  0.1549E+00  0.1580E+00  0.1609E+00  0.1638E+00  0.1665E+00  0.1691E+00
  0.1715E+00  0.1737E+00  0.1758E+00  0.1776E+00  0.1793E+00  0.1806E+00  0.1818E+00
  0.1827E+00  0.1833E+00  0.1836E+00  0.1836E+00  0.1832E+00  0.1825E+00  0.1814E+00
  0.1799E+00  0.1780E+00  0.1756E+00  0.1727E+00  0.1693E+00  0.1653E+00  0.1606E+00
  0.1554E+00  0.1494E+00  0.1427E+00  0.1353E+00  0.1271E+00  0.1183E+00  0.1089E+00
  0.9916E-01  0.8935E-01  0.7978E-01  0.7067E-01  0.6212E-01  0.5420E-01  0.4688E-01
  0.4015E-01  0.3394E-01  0.2820E-01  0.2286E-01  0.1788E-01  0.1316E-01  0.8659E-02
  0.4295E-02
 Values of PDE component i =  2
  0.7174E-05  0.1708E-04  0.3312E-04  0.6016E-04  0.1053E-03  0.1792E-03  0.2967E-03
  0.4789E-03  0.7536E-03  0.1157E-02  0.1734E-02  0.2538E-02  0.3631E-02  0.5078E-02
  0.6948E-02  0.9306E-02  0.1221E-01  0.1572E-01  0.1985E-01  0.2462E-01  0.3004E-01
  0.3608E-01  0.4271E-01  0.4988E-01  0.5755E-01  0.6566E-01  0.7415E-01  0.8298E-01
  0.9208E-01  0.1014E+00  0.1109E+00  0.1206E+00  0.1303E+00  0.1401E+00  0.1500E+00
  0.1599E+00  0.1697E+00  0.1794E+00  0.1891E+00  0.1986E+00  0.2079E+00  0.2170E+00
  0.2258E+00  0.2343E+00  0.2425E+00  0.2501E+00  0.2573E+00  0.2639E+00  0.2700E+00
  0.2754E+00  0.2802E+00  0.2843E+00  0.2878E+00  0.2907E+00  0.2931E+00  0.2949E+00
  0.2964E+00  0.2975E+00  0.2983E+00  0.2988E+00  0.2993E+00  0.2995E+00  0.2997E+00
  0.2998E+00  0.2999E+00  0.2999E+00  0.3000E+00  0.3000E+00  0.3000E+00  0.3000E+00
  0.2999E+00  0.2999E+00  0.2997E+00  0.2994E+00  0.2989E+00  0.2980E+00  0.2964E+00
  0.2937E+00  0.2890E+00  0.2814E+00  0.2693E+00  0.2509E+00  0.2245E+00  0.1899E+00
  0.1499E+00  0.1099E+00  0.7540E-01  0.4903E-01  0.3061E-01  0.1852E-01  0.1091E-01
  0.6277E-02  0.3530E-02  0.1941E-02  0.1042E-02  0.5460E-03  0.2774E-03  0.1333E-03
  0.5327E-04
 Values of PDE component i =  3
  0.1916E-09  0.1973E-08  0.1191E-07  0.5547E-07  0.2183E-06  0.7558E-06  0.2356E-05
  0.6708E-05  0.1762E-04  0.4304E-04  0.9826E-04  0.2107E-03  0.4262E-03  0.8157E-03
  0.1482E-02  0.2562E-02  0.4228E-02  0.6681E-02  0.1013E-01  0.1479E-01  0.2084E-01
  0.2843E-01  0.3766E-01  0.4855E-01  0.6109E-01  0.7522E-01  0.9082E-01  0.1078E+00
  0.1259E+00  0.1450E+00  0.1650E+00  0.1857E+00  0.2069E+00  0.2284E+00  0.2500E+00
  0.2717E+00  0.2932E+00  0.3143E+00  0.3350E+00  0.3550E+00  0.3741E+00  0.3922E+00
  0.4092E+00  0.4248E+00  0.4389E+00  0.4514E+00  0.4623E+00  0.4715E+00  0.4791E+00
  0.4852E+00  0.4899E+00  0.4933E+00  0.4958E+00  0.4974E+00  0.4985E+00  0.4992E+00
  0.4996E+00  0.4998E+00  0.4999E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00
  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.5000E+00  0.4998E+00
  0.4994E+00  0.4982E+00  0.4949E+00  0.4861E+00  0.4652E+00  0.4218E+00  0.3479E+00
  0.2497E+00  0.1517E+00  0.7801E-01  0.3480E-01  0.1387E-01  0.5057E-02  0.1711E-02
  0.5392E-03  0.1581E-03  0.4283E-04  0.1060E-04  0.2350E-05  0.4530E-06  0.7161E-07
  0.8102E-08


Final statistics for mf = 22:  103 steps,   292 res,    17 jacobians,
                              rwork size =  5368,  iwork size =   317
Final output is correct to within  0.12E+02  times local error tolerance. 


******************************************************************************************

Run completed:   0 errors encountered

==========Source and Output for LSODIS Demonstration Program====================

c-----------------------------------------------------------------------
c Demonstration program for the DLSODIS package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c This program solves a semi-discretized form of the Burgers equation,
c
c     u  = -(u*u/2)  + eta * u
c      t           x          xx
c
c for  -1 .le. x .le. 1, t .ge. 0.
c Here eta = 0.05.
c Boundary conditions: u(-1,t) = u(1,t) and du/dx(-1,t) = du/dx(1,t).
c Initial profile: square wave
c     u(0,x) = 0    for 1/2 .lt. abs(x) .le. 1
c     u(0,x) = 1/2  for abs(x) = 1/2
c     u(0,x) = 1    for 0 .le. abs(x) .lt. 1/2
c
c An ODE system is generated by a simplified Galerkin treatment
c of the spatial variable x.
c
c Reference:
c R. C. Y. Chin, G. W. Hedstrom, and K. E. Karlsson,
c A Simplified Galerkin Method for Hyperbolic Equations,
c Math. Comp., vol. 33, no. 146 (April 1979), pp. 647-658.
c
c The problem is run with the DLSODIS package with a 12-node mesh,
c for various appropriate values of the method flag mf.
c Output is on unit lout, set to 6 in a data statement below.
c
c Problem specific data:
c npts  = number of unknowns (npts = 0 mod 4)
c nnz   = number of non-zeros in Jacobian before fill in
c nnza  = number of non-zeros in Jacobian after fill in
c lrwk  = length of real work array (taking into account fill in)
c liwk  = length of integer work array
c ipia  = pointer to ia in iw (ia(j) = iw(ipia+j-1)
c ipja  = pointer to ja in iw (ja(j) = iw(ipja+j-1)
c ipic  = pointer to ic in iw array (ic(j) = iw(ipic+j-1))
c ipjc  = pointer to jc in iw array (jc(j) = iw(ipjc+j-1))
c-----------------------------------------------------------------------
      parameter (npts = 12, nnz = 3*npts, nnza = 5*npts)
      parameter (lrwk = 20+3*nnza+28*npts)
      parameter (liwk = 32+2*nnza+2*npts)
      parameter (ipia = 31, ipja = 31+npts+1)
      parameter (ll = 32+npts+nnz, ipic = ll,ipjc = ll+npts+1)
c
      external res, addasp, jacsp
      integer i, io, istate, itol, iw, j, lout, lrw, liw, lyh,
     1   meth, miter, mf, moss, n, nout, nerr, n14, n34
      double precision eta, delta, r4d, eodsq, a, b, t, tout, tlast,
     1    tinit, errfac
      double precision atol, rtol, rw, y, ydoti, elkup
      double precision zero, fourth, half, one, hun
      dimension y(npts), ydoti(npts), tout(4), atol(2), rtol(2)
      dimension rw(lrwk), iw(liwk)
c Pass problem parameters in the Common block test1.
      common /test1/ r4d, eodsq
c
c Set problem parameters and run parameters
      data eta/0.05d0/, a/-1.0d0/, b/1.0d0/
      data zero/0.0d0/, fourth/0.25d0/, half/0.5d0/, one/1.0d0/,
     1   hun/100.0d0/
      data tinit/0.0d0/, tlast/0.4d0/
      data tout/0.10d0, 0.20d0, 0.30d0, 0.40d0/
      data lout/6/, nout/4/
      data itol/1/, rtol/1.0d-3, 1.0d-6/, atol/1.0d-3, 1.0d-6/
c
      nerr = 0
      lrw = lrwk
      liw = liwk
c
c Compute the mesh width delta and other parameters.
      delta = (b - a)/npts
      r4d = fourth/delta
      eodsq = eta/delta**2
      n14 = npts/4 + 1
      n34 = 3 * (npts/4) + 1
      n = npts
c
c Set the initial profile (for output purposes only).
      do 10 i = 1,n
   10   y(i) = zero
      y(n14) = half
      do 20 i = n14+1,n34-1
   20   y(i) = one
      y(n34) = half
c
      write(lout,1000)
      write(lout,1100) eta,a,b,tinit,tlast,n
      write(lout,1200) (y(i), i=1,n)
c
c Set the initial sparse data structures for coefficient matrix A
c and the Jacobian matrix C
      call struct(iw(ipia), iw(ipja), iw(ipic), iw(ipjc), n)
c
c The j loop is over error tolerances.
      do 200 j = 1,2
c
c This method flag loop is for demonstration only.
      do 200 moss = 0,4
        do 100 meth = 1,2
          do 100 miter = 1,2
   35       mf = 100*moss + 10*meth + miter
c
c Set the initial profile.
            do 40 i = 1,n
   40         y(i) = zero
            y(n14) = half
            do 50 i = n14+1,n34-1
   50         y(i) = one
            y(n34) = half
c
            t = tinit
            istate = 0
c
            write(lout,1500) itol, rtol(j), atol(j), mf
c
c output loop for each case
            do 80 io = 1,nout
c
c Call DLSODIS and print results.
              call dlsodis (res, addasp, jacsp, n, y, ydoti, t,
     1             tout(io), itol, rtol(j), atol(j), 1, istate, 0,
     2             rw, lrw, iw, liw, mf)
              write(lout,2000) t, rw(11), iw(14), (y(i), i=1,n)
c
c If istate is not 2 on return, print message and go to next case.
              if (istate .ne. 2) then
                write(lout,4000) mf, t, istate
                nerr = nerr + 1
                go to 100
                endif
   80       continue
            write(lout,3000) mf, iw(11), iw(12), iw(13),
     1                       iw(17), iw(18), iw(20), iw(21)
c
c Estimate final error and print result.
            lyh = iw(22)
            errfac = elkup(n, y, rw(lyh), itol, rtol(j), atol(j), lout)
            if (errfac .lt. hun) then
              write(lout,5000) errfac
            else
              write(lout,5001) errfac
              nerr = nerr + 1
              endif
  100     continue
  200   continue
c
      write(lout,6000) nerr
      stop
c
 1000 format(20x,' Demonstration Program for DLSODIS' )
 1100 format(//10x,'-- Simplified Galerkin solution of ',
     1       'Burgers equation --'///
     1       13x,'Diffusion coefficient is eta =',d10.2/
     1       13x,'Uniform mesh on interval',d12.3,' to ',d12.3/
     2       13x,'Periodic boundary conditions'/
     2       13x,'Initial data are as follows:'//20x,'t0 = ',d12.5/
     2       20x,'tlast = ',d12.5/20x,'n  = ',i3//)
c
 1200 format(/'Initial profile:',/20(6d12.4/))
c
 1500 format(///85('*')///'Run with itol =',i2,'  rtol =',d12.2,
     1       '  atol =',d12.2,'   mf = ',i3//)
c
 2000 format(' Output for time t =',d12.5,'  current h =',d12.5,
     1       '  current order =',i2/20(6d12.4/))
c
 3000 format(/'Final statistics for mf = ',i3,': ',
     1       i5,' steps,',i6,' res,',i6,' Jacobians,'/
     2       20x,' rw size =',i6,',    iw size =',i6/
     3       20x,i4,' extra res for each jac,',i4,' decomps')
c
 4000 format(/'Final time reached for mf = ',i3,
     1       ' was t = ',d12.5/25x,'at which istate = ',i2//)
 5000 format('Final output is correct to within ',d9.2,
     1       '  times local error tolerance.'/)
 5001 format('Final output is wrong by ',d9.2,
     1       '  times local error tolerance.'/)
 6000 format(///85('*')//
     1       'Run completed: number of errors encountered =',i3)
c
c end of main program for the DLSODIS demonstration program
      end

      subroutine struct(ia, ja, ic, jc, n)
c This subroutine computes the initial sparse data structure of
c the mass (ia,ja) and Jacobian (ic,jc) matrices.
c
      integer ia(*), ja(*), ic(*), jc(*), n,  jj, k, l, m
c
      write(6,1200)
      k = 0
      do 33 l = 1,n
         ia(l) = (l-1)*3+1
         ic(l) = (l-1)*3+1
         do 32 m = l,l+2
            k = k + 1
            ja(k) = m - 1
            jc(k) = m - 1
   32    continue
   33 continue
      ia(n+1) = 3*n + 1
      ic(n+1) = 3*n+1
      ja(1) = n
      jc(1) = n
      ja(k) = 1
      jc(k) = 1
c
      write(6,1300) (ia(jj),jj=1,n+1)
      write(6,1350) (ja(jj),jj=1,k)
      write(6,1400) (ic(jj),jj=1,n+1)
      write(6,1450) (jc(jj),jj=1,k)
      return
1200  format('Initial sparse data structures'/)
1300  format(' ia  ',15i4/10(5x,15i4/))
1350  format(' ja  ',15i4/10(5x,15i4/))
1400  format(' ic  ',15i4/10(5x,15i4/))
1450  format(' jc  ',15i4/10(5x,15i4/))
      end

      subroutine res (n, t, y, v, r, ires)
c This subroutine computes the residual vector
c   r = g(t,y) - A(t,y)*v .
c If ires = -1, only g(t,y) is returned in r, since A(t,y) does
c not depend on y.
c No changes need to be made to this routine if n is changed.
c
      integer n, ires,  i
      double precision t, y(n), v(n), r(n), r4d, eodsq, one, four, six,
     1   fact1, fact4
      common /test1/ r4d, eodsq
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      call gfun (n, t, y, r)
      if (ires .eq. -1) return
c
      fact1 = one/six
      fact4 = four/six
c
      r(1) = r(1) - (fact4*v(1) + fact1*(v(2) + v(n)))
      do 10 i = 2,n-1
         r(i) = r(i) - (fact4*v(i) + fact1*(v(i-1) + v(i+1)))
   10 continue
      r(n) = r(n) - (fact4*v(n) + fact1*(v(1) + v(n-1)))
      return
c end of subroutine res for the DLSODIS demonstration program
      end

      subroutine gfun (n, t, y, g)
c This subroutine computes the right-hand side function g(y,t).
c It uses r4d = 1/(4*delta) and eodsq = eta/delta**2
c from the Common block test1.
c
      integer n, i
      double precision t, y(n), g(n), r4d, eodsq, two
      common /test1/ r4d, eodsq
      data two/2.0d0/
c
      g(1) = r4d*(y(n)**2 - y(2)**2) + eodsq*(y(2) - two*y(1) + y(n))
c
      do 20 i = 2,n-1
        g(i) = r4d*(y(i-1)**2 - y(i+1)**2)
     1        + eodsq*(y(i+1) - two*y(i) + y(i-1))
   20   continue
c
      g(n) = r4d*(y(n-1)**2 - y(1)**2) + eodsq*(y(1)-two*y(n)+y(n-1))
c
      return
c end of subroutine gfun for the DLSODIS demonstration program
      end

      subroutine addasp (n, t, y, j, ip, jp, pa)
c This subroutine computes the sparse matrix A by columns, adds it to
c pa, and returns the sum in pa.
c The matrix A is periodic tridiagonal, of order n, with nonzero elements
c (reading across) of  1/6, 4/6, 1/6, with 1/6 in the lower left and
c upper right corners.
c
      integer n, j, ip(*), jp(*),  jm1, jp1
      double precision t, y(n), pa(n), fact1, fact4, one, four, six
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
c Compute the elements of A.
      fact1 = one/six
      fact4 = four/six
      jm1 = j - 1
      jp1 = j + 1
      if (j .eq. n) jp1 = 1
      if (j .eq. 1) jm1 = n
c
c Add the matrix A to the matrix pa (sparse).
      pa(j) = pa(j) + fact4
      pa(jp1) = pa(jp1) + fact1
      pa(jm1) = pa(jm1) + fact1
      return
c end of subroutine addasp for the DLSODIS demonstration program
      end

      subroutine jacsp (n, t, y, s, j, ip, jp, pdj)
c This subroutine computes the Jacobian dg/dy = d(g-A*s)/dy by
c columns in sparse matrix format.  Only nonzeros are loaded.
c It uses r4d = 1/(4*delta) and eodsq = eta/delta**2 from the Common
c block test1.
c
      integer n, j, ip(*), jp(*),  jm1, jp1
      double precision t, y(n), s(n), pdj(n), r4d, eodsq, two, diag, r2d
      common /test1/ r4d, eodsq
      data two/2.0d0/
c
      diag = -two*eodsq
      r2d = two*r4d
      jm1 = j - 1
      jp1 = j + 1
      if (j .eq. 1) jm1 = n
      if (j .eq. n) jp1 = 1
c
      pdj(jm1) = -r2d*y(j) + eodsq
      pdj(j) = diag
      pdj(jp1) = r2d*y(j) + eodsq
      return
c end of subroutine jacsp for the DLSODIS demonstration program
      end

      double precision function elkup (n, y, ewt, itol, rtol,atol, lout)
c This routine looks up approximately correct values of y at t = 0.4,
c ytrue = y12 or y120 depending on whether n = 12 or 120.
c These were obtained by running with very tight tolerances.
c The returned value is
c   elkup = norm of [ (y - ytrue) / (rtol*abs(ytrue) + atol) ].
c
      integer n, itol, lout, i
      double precision y(n), ewt(n), rtol, atol, y12(12), y120(120),
     1    y120a(16), y120b(16), y120c(16), y120d(16), y120e(16),
     2    y120f(16), y120g(16), y120h(8), dvnorm
      equivalence (y120a(1),y120(1)), (y120b(1),y120(17)),
     1      (y120c(1),y120(33)), (y120d(1),y120(49)),
     1      (y120e(1),y120(65)),
     1      (y120f(1),y120(81)), (y120g(1),y120(97)),
     1      (y120h(1),y120(113))
      data y12/
     1 1.60581860d-02, 3.23063251d-02, 1.21903380d-01, 2.70943828d-01,
     1 4.60951522d-01, 6.57571216d-01, 8.25154453d-01, 9.35644796d-01,
     1 9.90167557d-01, 9.22421221d-01, 5.85764902d-01, 1.81112615d-01/
      data y120a /
     1 1.89009068d-02, 1.63261891d-02, 1.47080563d-02, 1.39263623d-02,
     1 1.38901341d-02, 1.45336989d-02, 1.58129308d-02, 1.77017162d-02,
     1 2.01886844d-02, 2.32742221d-02, 2.69677715d-02, 3.12854037d-02,
     1 3.62476563d-02, 4.18776225d-02, 4.81992825d-02, 5.52360652d-02/
      data y120b /
     1 6.30096338d-02, 7.15388849d-02, 8.08391507d-02, 9.09215944d-02,
     1 1.01792784d-01, 1.13454431d-01, 1.25903273d-01, 1.39131085d-01,
     1 1.53124799d-01, 1.67866712d-01, 1.83334757d-01, 1.99502830d-01,
     1 2.16341144d-01, 2.33816600d-01, 2.51893167d-01, 2.70532241d-01/
      data y120c /
     1 2.89693007d-01, 3.09332757d-01, 3.29407198d-01, 3.49870723d-01,
     1 3.70676646d-01, 3.91777421d-01, 4.13124817d-01, 4.34670077d-01,
     1 4.56364053d-01, 4.78157319d-01, 5.00000270d-01, 5.21843218d-01,
     1 5.43636473d-01, 5.65330432d-01, 5.86875670d-01, 6.08223037d-01/
      data y120d /
     1 6.29323777d-01, 6.50129662d-01, 6.70593142d-01, 6.90667536d-01,
     1 7.10307235d-01, 7.29467947d-01, 7.48106966d-01, 7.66183477d-01,
     1 7.83658878d-01, 8.00497138d-01, 8.16665158d-01, 8.32133153d-01,
     1 8.46875019d-01, 8.60868691d-01, 8.74096465d-01, 8.86545273d-01/
      data y120e /
     1 8.98206892d-01, 9.09078060d-01, 9.19160487d-01, 9.28460742d-01,
     1 9.36989986d-01, 9.44763554d-01, 9.51800339d-01, 9.58122004d-01,
     1 9.63751979d-01, 9.68714242d-01, 9.73031887d-01, 9.76725449d-01,
     1 9.79811001d-01, 9.82297985d-01, 9.84186787d-01, 9.85466039d-01/
      data y120f /
     1 9.86109629d-01, 9.86073433d-01, 9.85291781d-01, 9.83673704d-01,
     1 9.81099057d-01, 9.77414704d-01, 9.72431015d-01, 9.65919133d-01,
     1 9.57609585d-01, 9.47193093d-01, 9.34324619d-01, 9.18631922d-01,
     1 8.99729965d-01, 8.77242371d-01, 8.50830623d-01, 8.20230644d-01/
      data y120g /
     1 7.85294781d-01, 7.46035145d-01, 7.02662039d-01, 6.55609682d-01,
     1 6.05541326d-01, 5.53327950d-01, 4.99999118d-01, 4.46670394d-01,
     1 3.94457322d-01, 3.44389410d-01, 2.97337561d-01, 2.53964948d-01,
     1 2.14705729d-01, 1.79770169d-01, 1.49170367d-01, 1.22758681d-01/
      data y120h /
     1 1.00271052d-01, 8.13689920d-02, 6.56761515d-02, 5.28075160d-02,
     1 4.23908624d-02, 3.40811650d-02, 2.75691506d-02, 2.25853507d-02/
c
      if ((n-12)*(n-120) .ne. 0) go to 300
      if (n .eq. 120) go to 100
c
c Compute local error tolerance using correct y (n = 12).
      call dewset(n, itol, rtol, atol, y12, ewt)
c
c Invert ewt and replace y by the error, y - ytrue.
      do 20  i = 1, 12
        ewt(i) = 1.0d0/ewt(i)
 20     y(i) = y(i) - y12(i)
      go to 200
c
c Compute local error tolerance using correct y (n = 120).
 100  call dewset( n, itol, rtol, atol, y120, ewt )
c
c Invert ewt and replace y by the error, y - ytrue.
      do 120  i = 1, 120
        ewt(i) = 1.0d0/ewt(i)
 120    y(i) = y(i) - y120(i)
c
c Find weighted norm of the error.
 200  elkup = dvnorm (n, y, ewt)
      return
c
c error return
 300  write(lout,400) n
      elkup = 1.0d3
 400  format(/5x,'Illegal use of elkup for n =',i4)
      return
c end of function elkup for the DLSODIS demonstration program
      end

................................................................................

                     Demonstration Program for DLSODIS


          -- Simplified Galerkin solution of Burgers equation --


             Diffusion coefficient is eta =  0.50E-01
             Uniform mesh on interval  -0.100E+01 to    0.100E+01
             Periodic boundary conditions
             Initial data are as follows:

                    t0 =  0.00000E+00
                    tlast =  0.40000E+00
                    n  =  12



Initial profile:
  0.0000E+00  0.0000E+00  0.0000E+00  0.5000E+00  0.1000E+01  0.1000E+01
  0.1000E+01  0.1000E+01  0.1000E+01  0.5000E+00  0.0000E+00  0.0000E+00

Initial sparse data structures

 ia     1   4   7  10  13  16  19  22  25  28  31  34  37
 ja    12   1   2   1   2   3   2   3   4   3   4   5   4   5   6
        5   6   7   6   7   8   7   8   9   8   9  10   9  10  11
       10  11  12  11  12   1
 ic     1   4   7  10  13  16  19  22  25  28  31  34  37
 jc    12   1   2   1   2   3   2   3   4   3   4   5   4   5   6
        5   6   7   6   7   8   7   8   9   8   9  10   9  10  11
       10  11  12  11  12   1



*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf =  11


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf =  11:    10 steps,    14 res,     4 Jacobians,
                     rw size =   421,    iw size =   128
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf =  12


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf =  12:    10 steps,    30 res,     4 Jacobians,
                     rw size =   429,    iw size =   128
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf =  21


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf =  21:    13 steps,    18 res,     4 Jacobians,
                     rw size =   337,    iw size =   128
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf =  22


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf =  22:    13 steps,    34 res,     4 Jacobians,
                     rw size =   345,    iw size =   128
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 111


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 111:    10 steps,    14 res,     4 Jacobians,
                     rw size =   421,    iw size =    30
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 112


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 112:    10 steps,    30 res,     4 Jacobians,
                     rw size =   429,    iw size =    30
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 121


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 121:    13 steps,    18 res,     4 Jacobians,
                     rw size =   337,    iw size =    30
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 122


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 122:    13 steps,    34 res,     4 Jacobians,
                     rw size =   345,    iw size =    30
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 211


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 211:    10 steps,    14 res,     4 Jacobians,
                     rw size =   421,    iw size =    30
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 212


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 212:    10 steps,    30 res,     4 Jacobians,
                     rw size =   429,    iw size =    30
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 221


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 221:    13 steps,    18 res,     4 Jacobians,
                     rw size =   337,    iw size =    30
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 222


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 222:    13 steps,    34 res,     4 Jacobians,
                     rw size =   345,    iw size =    30
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 311


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 311:    10 steps,    14 res,     4 Jacobians,
                     rw size =   421,    iw size =    79
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 312


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 312:    10 steps,    30 res,     4 Jacobians,
                     rw size =   429,    iw size =    79
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 321


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 321:    13 steps,    18 res,     4 Jacobians,
                     rw size =   337,    iw size =    79
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 322


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 322:    13 steps,    34 res,     4 Jacobians,
                     rw size =   345,    iw size =    79
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 411


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 411:    10 steps,    14 res,     4 Jacobians,
                     rw size =   421,    iw size =    79
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 412


 Output for time t = 0.10000E+00  current h = 0.49176E-01  current order = 2
 -0.1463E-02 -0.3833E-02  0.7844E-01  0.3812E+00  0.7671E+00  0.9728E+00
  0.1005E+01  0.1001E+01  0.9973E+00  0.6591E+00  0.1534E+00 -0.9859E-02

 Output for time t = 0.20000E+00  current h = 0.49176E-01  current order = 2
 -0.8416E-02  0.1100E-01  0.1057E+00  0.3273E+00  0.6182E+00  0.8599E+00
  0.9783E+00  0.1002E+01  0.1002E+01  0.7756E+00  0.3045E+00  0.2427E-01

 Output for time t = 0.30000E+00  current h = 0.82814E-01  current order = 2
 -0.3693E-02  0.2323E-01  0.1173E+00  0.2941E+00  0.5236E+00  0.7478E+00
  0.9094E+00  0.9830E+00  0.1004E+01  0.8622E+00  0.4498E+00  0.8959E-01

 Output for time t = 0.40000E+00  current h = 0.82814E-01  current order = 2
  0.1620E-01  0.3267E-01  0.1222E+00  0.2704E+00  0.4595E+00  0.6569E+00
  0.8269E+00  0.9378E+00  0.9903E+00  0.9215E+00  0.5849E+00  0.1806E+00


Final statistics for mf = 412:    10 steps,    30 res,     4 Jacobians,
                     rw size =   429,    iw size =    79
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.60E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 421


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 421:    13 steps,    18 res,     4 Jacobians,
                     rw size =   337,    iw size =    79
                       0 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-02  atol =    0.10E-02   mf = 422


 Output for time t = 0.10000E+00  current h = 0.41065E-01  current order = 3
 -0.1679E-02 -0.3550E-02  0.7815E-01  0.3816E+00  0.7669E+00  0.9726E+00
  0.1006E+01  0.1000E+01  0.9974E+00  0.6590E+00  0.1534E+00 -0.9738E-02

 Output for time t = 0.20000E+00  current h = 0.41065E-01  current order = 3
 -0.9106E-02  0.1195E-01  0.1044E+00  0.3292E+00  0.6187E+00  0.8565E+00
  0.9795E+00  0.1002E+01  0.1002E+01  0.7752E+00  0.3044E+00  0.2472E-01

 Output for time t = 0.30000E+00  current h = 0.41065E-01  current order = 3
 -0.3840E-02  0.2325E-01  0.1168E+00  0.2953E+00  0.5256E+00  0.7457E+00
  0.9071E+00  0.9838E+00  0.1004E+01  0.8621E+00  0.4499E+00  0.8988E-01

 Output for time t = 0.40000E+00  current h = 0.67572E-01  current order = 3
  0.1612E-01  0.3233E-01  0.1218E+00  0.2713E+00  0.4617E+00  0.6568E+00
  0.8235E+00  0.9363E+00  0.9915E+00  0.9222E+00  0.5853E+00  0.1811E+00


Final statistics for mf = 422:    13 steps,    34 res,     4 Jacobians,
                     rw size =   345,    iw size =    79
                       3 extra res for each jac,   4 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf =  11


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf =  11:    28 steps,    38 res,     7 Jacobians,
                     rw size =   421,    iw size =   128
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf =  12


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf =  12:    28 steps,    66 res,     7 Jacobians,
                     rw size =   429,    iw size =   128
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf =  21


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf =  21:    40 steps,    49 res,     7 Jacobians,
                     rw size =   337,    iw size =   128
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf =  22


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf =  22:    40 steps,    77 res,     7 Jacobians,
                     rw size =   345,    iw size =   128
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 111


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 111:    28 steps,    38 res,     7 Jacobians,
                     rw size =   421,    iw size =    30
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 112


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 112:    28 steps,    66 res,     7 Jacobians,
                     rw size =   429,    iw size =    30
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 121


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 121:    40 steps,    49 res,     7 Jacobians,
                     rw size =   337,    iw size =    30
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 122


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 122:    40 steps,    77 res,     7 Jacobians,
                     rw size =   345,    iw size =    30
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 211


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 211:    28 steps,    38 res,     7 Jacobians,
                     rw size =   421,    iw size =    30
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 212


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 212:    28 steps,    66 res,     7 Jacobians,
                     rw size =   429,    iw size =    30
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 221


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 221:    40 steps,    49 res,     7 Jacobians,
                     rw size =   337,    iw size =    30
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 222


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 222:    40 steps,    77 res,     7 Jacobians,
                     rw size =   345,    iw size =    30
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 311


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 311:    28 steps,    38 res,     7 Jacobians,
                     rw size =   421,    iw size =    79
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 312


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 312:    28 steps,    66 res,     7 Jacobians,
                     rw size =   429,    iw size =    79
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 321


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 321:    40 steps,    49 res,     7 Jacobians,
                     rw size =   337,    iw size =    79
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 322


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 322:    40 steps,    77 res,     7 Jacobians,
                     rw size =   345,    iw size =    79
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 411


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 411:    28 steps,    38 res,     7 Jacobians,
                     rw size =   421,    iw size =    79
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 412


 Output for time t = 0.10000E+00  current h = 0.17717E-01  current order = 5
 -0.1666E-02 -0.3562E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.17717E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.24980E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.34040E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 412:    28 steps,    66 res,     7 Jacobians,
                     rw size =   429,    iw size =    79
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.41E+00  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 421


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 421:    40 steps,    49 res,     7 Jacobians,
                     rw size =   337,    iw size =    79
                       0 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************


Run with itol = 1  rtol =    0.10E-05  atol =    0.10E-05   mf = 422


 Output for time t = 0.10000E+00  current h = 0.91135E-02  current order = 4
 -0.1666E-02 -0.3563E-02  0.7808E-01  0.3817E+00  0.7672E+00  0.9718E+00
  0.1006E+01  0.1001E+01  0.9973E+00  0.6590E+00  0.1534E+00 -0.9756E-02

 Output for time t = 0.20000E+00  current h = 0.18660E-01  current order = 5
 -0.8514E-02  0.1107E-01  0.1054E+00  0.3279E+00  0.6190E+00  0.8586E+00
  0.9773E+00  0.1003E+01  0.1002E+01  0.7754E+00  0.3045E+00  0.2439E-01

 Output for time t = 0.30000E+00  current h = 0.18660E-01  current order = 5
 -0.3772E-02  0.2305E-01  0.1171E+00  0.2947E+00  0.5251E+00  0.7475E+00
  0.9073E+00  0.9823E+00  0.1004E+01  0.8626E+00  0.4499E+00  0.8982E-01

 Output for time t = 0.40000E+00  current h = 0.22945E-01  current order = 5
  0.1606E-01  0.3231E-01  0.1219E+00  0.2709E+00  0.4610E+00  0.6576E+00
  0.8252E+00  0.9356E+00  0.9902E+00  0.9224E+00  0.5858E+00  0.1811E+00


Final statistics for mf = 422:    40 steps,    77 res,     7 Jacobians,
                     rw size =   345,    iw size =    79
                       3 extra res for each jac,   7 decomps
Final output is correct to within  0.10E+01  times local error tolerance.




*************************************************************************************

Run completed: number of errors encountered =  0

