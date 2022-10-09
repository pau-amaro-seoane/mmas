c       Copyright (c) 2001, 2002, 2003
c       by James Lombardi, Vassar College, Poughkeepsie, NY
c       and Jessica Sawyer Warren, Rutgers University, NJ.
c    This file is part of Make Me A Star.
c
c    Make Me A Star is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.
c
c    Make Me A Star is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with Make Me A Star; if not, write to the Free Software
c    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
*************************************************************************
* Make Me A Star:
* This file contains the main subroutine makemeastar
* version 1.6
* July 2003
* James Lombardi (lombardi@vassar.edu) and Jessica Sawyer Warren
*************************************************************************
      SUBROUTINE makemeastar(peri,
     $     mProfile1,rProfile1,PProfile1,rhoProfile1,chemicalProfiles1,
     $     mProfile2,rProfile2,PProfile2,rhoProfile2,chemicalProfiles2,
     $     numShells1,numShells2,
     &     verboseness,numberOfChemicalElements,
     $     NumOutLines,Mremnant,Premnant,
     &     Dremnant,Rremnant,ELremnant,jremnant)
**********************************************************************
* Input:
*     periastron = r_p / (R_1+R_2), where r_p is the periastron separation of
*        the (parabolic) orbits, and R_1 and R_2 are the radii of the parent
*        stars when isolated.  Therefore, r_p=0 for a headon collision and
*        increases as the impact parameter increases.
*     verboseness = an integer that sets how talkative the routines should be
*        verboseness=0 : no output
*        verboseness=1 : limited output
*        verboseness=2 : detailed output
*     el = an integer stating how many chemical abundances will be treated
*        (information on the parent abundance profiles must of course be
*         provided in parentfile(1) and parentfile(2))
*     NumOutLines = desired number of rows in the outputted remnant profiles
* Output:
*     Mremnant = an array which gives enclosed mass in the collision product
*     Premnant = an array for pressure profile in the collision product
*     Dremnant = an array for density profile in the collision product
*     Rremnant = an array for radius profile in the collision product
*     ELremnant = an array of size (numberOfShells x numberOfChemicalElements)
*        for chemical composition profiles in the collision product
*     jremnant = an array for the specific angular momentum profile in the
*        collision product
**********************************************************************
* This subroutine reads in data containing enclosed mass M,          *
* pressure P, and density D (in                                      *
* cgs units), and chemical abundances of a particular star.  It will *
* calculate the "entropy" A from P and D, and call                   *
* subroutines SORT, PROFILE, ENERGY, LOSEMASS, and SHOCK.            *
*     ~J.E.S. and J.C.L                                              *
**********************************************************************

      IMPLICIT NONE

* params.h contains declarations and values for N, NROWS, pi, G,
* Msun, & Rsun.  USE CGS UNITS!!!

      INCLUDE 'params.h'
      integer verbosity,NP,NumOutLines,verboseness
      double precision peri
      INTEGER Jpmax(N),I,K,C,Pmax(N),Jr,Jrmax
      DOUBLE PRECISION mProfile1(NROWS),rProfile1(NROWS),
     $     PProfile1(NROWS),rhoProfile1(NROWS)
      DOUBLE PRECISION mProfile2(NROWS),rProfile2(NROWS),
     $     PProfile2(NROWS),rhoProfile2(NROWS)
      DOUBLE PRECISION chemicalProfiles1(NROWS,elmax),
     $     chemicalProfiles2(NROWS,elmax)
      integer numShells1,numShells2
      DOUBLE PRECISION Mp(NROWS,N),Ap(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N),Aptry(N),TotE(N),W,As(NROWS,N)
      DOUBLE PRECISION ELp(NROWS,elmax,N),
     &     ELm(NROWS,elmax,N)
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),ELr(NROWS,elmax),P(NROWS),
     &     D(NROWS),R(NROWS)
      DOUBLE PRECISION Mremnant(NROWS),
     &     ELremnant(elmax,NROWS),
     &     Premnant(NROWS), Dremnant(NROWS),
     &     Rremnant(NROWS),jremnant(NROWS)
      DOUBLE PRECISION INTlow,INThigh,NEWINT,DELNRG,zeroin3,
     &    periastron,TotErem,Jtot,TotJrem,MLsfrac
      COMMON/NRG/TotE,TotErem
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
      COMMON/PERI/periastron
      COMMON/MASSLOSS/MLsfrac
      COMMON/verbiage/verbosity
      EXTERNAL DELNRG
      logical alreadyfudged
      integer numberOfChemicalElements
      integer throwaway,Jrtransition
      common/numel/el
      double precision scaledJtot,Jtotin,unscaledJtotout(NROWS)
      double precision cs,jcs
      double precision jtryin(NROWS),jtryout(NROWS),k2,k3,diff
      character*16 OUTFN
      double precision monm(NROWS),dMp(NROWS),
     $     y2b(NROWS),y2c(NROWS),y2d(NROWS)
      common/delms/monm,dMp,y2b,y2c,y2d
      logical firsttime
      data firsttime/.true./
      
* duplicate (or nearly duplicate) the values of some of the input variables so
* that they can be passed through common blocks:
      NP=NumOutLines-1
      verbosity=verboseness
      el=numberOfChemicalElements

      if(mProfile1(numShells1).ge.mProfile2(numShells2)) then
         Jpmax(1)=numShells1
         DO I=1,Jpmax(1)
            Mp(I,1)=mProfile1(I)
            Rp(I,1)=rProfile1(I)
            Pp(I,1)=PProfile1(I)
            Dp(I,1)=rhoProfile1(I)
            do C=1,el
               Elp(I,C,1)=chemicalProfiles1(I,C)
            enddo
         enddo
         Jpmax(2)=numShells2
         DO I=1,Jpmax(2)
            Mp(I,2)=mProfile2(I)
            Rp(I,2)=rProfile2(I)
            Pp(I,2)=PProfile2(I)
            Dp(I,2)=rhoProfile2(I)
            do C=1,el
               Elp(I,C,2)=chemicalProfiles2(I,C)
            enddo
         enddo
      else
C Switch stars if more massive star was second...
         Jpmax(1)=numShells2
         DO I=1,Jpmax(1)
            Mp(I,1)=mProfile2(I)
            Rp(I,1)=rProfile2(I)
            Pp(I,1)=PProfile2(I)
            Dp(I,1)=rhoProfile2(I)
            do C=1,el
               Elp(I,C,1)=chemicalProfiles2(I,C)
            enddo
         enddo
         Jpmax(2)=numShells1
         DO I=1,Jpmax(2)
            Mp(I,2)=mProfile1(I)
            Rp(I,2)=rProfile1(I)
            Pp(I,2)=PProfile1(I)
            Dp(I,2)=rhoProfile1(I)
            do C=1,el
               Elp(I,C,2)=chemicalProfiles1(I,C)
            enddo
         enddo
      endif

      if(verbosity.ge.2) then
         write(6,*) 'Here are the values of the constants being used:'
         write(6,*) 'c_1=',const1
         write(6,*) 'c_2=',const2
         write(6,*) 'c_3=',const3
         write(6,*) 'c_4=',const4
         write(6,*) 'c_5=',const5
         write(6,*) 'c_6=',const6
         write(6,*) 'c_7=',const7
         write(6,*) 'c_8=',const8
         write(6,*) 'c_9=',const9
         write(6,*) 'c_10=',const10
         write(6,*)
      endif

      if(el.GT.elmax) then
         write(6,*)'ERROR: Increase elmax in params.h to at least',el
         stop
      endif

      if(verbosity.ge.1) write(6,'(1x,a,f7.3)')
     &     'Entering MAKEMEASTAR w/ periastron parameter r_p/(R_1+R_2)='
     &     ,peri

* Now we have to calculate A.

      DO K=1,N
         if(verbosity.ge.2) WRITE(6,*)

         alreadyfudged=.false.
         DO I=1,Jpmax(K)
            Aptry(K)=Pp(I,K)/Dp(I,K)**(5.d0/3.d0)

* We want to fake data for which the "entropy" A does not increase.  Therefore
* the else part of the IF statement fakes bad data.  In particular,
* the "entropy" is forced to increase ever so slightly.
            IF(Aptry(K).GT.Ap(I-1,K) .OR. I.EQ.1) THEN
               Ap(I,K)=Aptry(K)
            else
               Ap(I,K)=1.00001d0*Ap(I-1,K)
               if(.not. alreadyfudged) then
                  if(verbosity.ge.2) then
                     write(6,*)
     &                 'Fudging A profile outside m/M_solar=',
     &                 Mp(I,K)/Msun,' in star',K,' at shell I=',I
                     write(6,*)
                  endif
                  alreadyfudged=.true.
               endif
            endif
         ENDDO

         if(verbosity.ge.2) then
            write(6,*) 'Information about star',K,':'
            WRITE(6,*) 'Jpmax(',K,')=',Jpmax(K)
            WRITE(6,*) 'Mass (in solar units)=',Mp(Jpmax(K),K)/Msun
            WRITE(6,*) 'Radius (in solar units)=',Rp(Jpmax(K),K)/Rsun
         endif

         CALL ENERGY(Mp(1,K),Pp(1,K),Dp(1,K),Rp(1,K),
     &        Jpmax(K),K)
         
      ENDDO

c     First find mass loss for the headon collision case between these two
c     parent stars:
      periastron=0.d0
      CALL LOSEMASS

      if(firsttime) then
         if(verbosity.ge.2) write(6,*) 'Reading mesh spacing file'
         open(30,file='meshspacing')
         read(30,*)
         DO I=1,189
            read(30,*) monm(I), dMp(I)
         ENDDO
         close(30)
         call spline(189,monm,dMp,y2b,y2c,y2d)
         firsttime=.false.
      endif

c Determine the intercept b_1 for the shock heating for the headon case
      INTlow=-0.9d0
      INThigh=0.9d0

      NEWINT=zeroin3(INTlow,INThigh,DELNRG,1.d-17)

      if(verbosity.ge.2) write(6,*)
     &     'The headon collision gave an intercept b_1=',NEWINT

      if(peri.ne.0.d0) then
         periastron=peri
         if(verbosity.ge.2) then
            write(6,*)
            write(6,*)'WE CAN NOW TREAT THE ACTUAL PERIASTRON CASE...'
         endif
         CALL LOSEMASS
         NEWINT=NEWINT-const4*periastron*
     &        LOG10(Mp(Jpmax(1),1)/Mp(Jpmax(2),2))
         if(verbosity.ge.2) then
            write(6,*)'For this periastron separation, b_1=',NEWINT
         endif
      endif

      if(verbosity.ge.2) WRITE(6,*)
      if(verbosity.ge.2) WRITE(6,*) 'Starting final SHOCK'
      CALL SHOCK(NEWINT)
      if(verbosity.ge.2) then
         DO K=1,N
            write(6,*)'central shocked A for star',K,' =',As(1,K)
            write(6,*)'surface shocked A for star',K,' =',As(Jpmax(K),K)
            write(6,*)'shocked A for star',K,' at last bound mass=',
     &           As(Pmax(K),K)
            write(6,*)'mass lost from star',K,'=',Mp(Jpmax(K),K)
     $           -Mp(Pmax(K),K)
         ENDDO
         write(6,*)'fraction of ejecta from star 2=',
     $        (Mp(Jpmax(2),2)-Mp(Pmax(2),2))/
     $        (Mp(Jpmax(1),1)-Mp(Pmax(1),1)
     $        +Mp(Jpmax(2),2)-Mp(Pmax(2),2))
      endif

      if(el.eq.0 .and. .not. studymassloss) then
         if(verbosity.ge.2) WRITE(6,*) 'SKIPPING THE MIXING'
      else
         DO K=1,N
            if(verbosity.ge.2) then
               WRITE(6,*)
               WRITE(6,*) 'Calculating mixing for star',K
            endif
            CALL MIX(K)
         ENDDO
      endif


      if(makeshockfiles) then
         do K=1,2
            WRITE(OUTFN,102) K
 102        FORMAT('shock',I1,'.sbe')
            if(verbosity.ge.2) write(6,*)'WRITING SHOCK FILE...',OUTFN
            open(19,file=OUTFN)
            Do I=1,Jpmax(K)
               write(19,*) Mp(I,K),dlog10(As(I,K)-Ap(I,K)),
     $              dlog10(As(I,K)/Ap(I,K))
            enddo
            close(19)
            if(verbosity.ge.2) write(6,*)'...DONE'
         enddo
      endif
      
* Pmax(K) comes from the LOSEMASS subroutine.

      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) 'Starting final SORT'
      endif
      CALL SORT
      if(verbosity.ge.2) then
         WRITE(6,*) 'surface A for remnant=',Ar(Jrmax)
         WRITE(6,*) 'central A for remnant=',Ar(1)
         WRITE(6,*) 'row max for remnant (Jrmax)=',Jrmax
      endif
      
      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) 'Starting final PROFILE'
      endif
      CALL PROFILE
      
      W=0.d0
      DO K=1,N
         W=TotE(K)+W
      ENDDO

      if(verbosity.ge.2) then      
         WRITE(6,*)
         WRITE(6,*) 'Total energy of parent stars is',W
      endif
      if(verbosity.ge.1) then
         WRITE(6,'(a,f7.3)')' Total mass of remnant (in solar units) is'
     &        ,Mr(Jrmax)/Msun
c         DO K=1,N
c            write(6,*) 'Radius of parent star ',K,' (in solar units)='
c     &           ,Rp(Jpmax(K),K)/Rsun
c            write(6,*) 'Mass of parent star ',K,' (in solar units)='
c     &           ,Mp(Jpmax(K),K)/Msun
c         ENDDO
      endif
      Jtot=Mp(Jpmax(1),1)*Mp(Jpmax(2),2)*
     &     (2*G*periastron*(Rp(Jpmax(1),1)+Rp(Jpmax(2),2))/
     &     (Mp(Jpmax(1),1)+Mp(Jpmax(2),2))
     &     )**0.5d0
      if(verbosity.ge.2) write(6,*)
     &     'Total angular momentum of parent stars is (cgs)',Jtot
      TotJrem=Jtot*(1.d0-const9*MLsfrac)
      if(verbosity.ge.2) then
         write(6,*)'f_L * J_tot=',MLsfrac*Jtot
         write(6,*)'Total angular momentum of remnant is (cgs)',TotJrem
         write(6,*) 'J_remnant/M_remnant = ',TotJrem/Mr(Jrmax)
         write(6,*)'Total system energy E_tot is (cgs)',TotE(1)+TotE(2)
         write(6,*)'Total energy of remnant is (cgs)',TotErem
         write(6,*)'f_L * E_tot=',MLsfrac*(TotE(1)+TotE(2))
      endif

      throwaway=0
      DO Jr=1,Jrmax
         if(P(Jr).lt.0.d0 .and. Jr.lt.Jrmax) then
            throwaway=throwaway+1
         else
            Mremnant(Jr-throwaway)=Mr(Jr)
C            Aremnant(Jr-throwaway)=Ar(Jr)
            Premnant(Jr-throwaway)=dabs(P(Jr))
            Dremnant(Jr-throwaway)=D(Jr)
            Rremnant(Jr-throwaway)=R(Jr)
            do C=1,el
*     Note that the order of the arguments is different in ELremnant than it
*     is in Elr.
               Elremnant(C,Jr-throwaway)=ELr(Jr,C)
            enddo
         endif
      ENDDO

      if(verbosity.ge.2) then
         write(6,*)'Jrmax=',Jrmax
         write(6,*) 'throwing away',throwaway,' lines.'
      endif
      NumOutLines=Jrmax-throwaway

      Do Jr=1,NumOutLines
c jtryout will keep track of j in the outer layers, without yet including
c the contribution due to the k_3 term:
         jtryout(Jr)=(G*Mremnant(Jr)*Rremnant(Jr))**0.5d0*
     $        (Mremnant(Jr)/Mremnant(NumOutLines))**(1.d0/3.d0)

         cs=(5.d0/3.d0*Premnant(Jr)/Dremnant(Jr))**0.5d0
         jcs=cs*Rremnant(Jr)
c jtryin will keep track of j in the inner layers:
         jtryin(Jr)=const10*jcs*
     $        (Mremnant(Jr)/Mremnant(NumOutLines))**(1.d0/3.d0)
      enddo

c unscaledJtotout keeps track of the total angular momentum for the shells
c outside (and including) the shell in question, without the k3 contribution:
      unscaledJtotout(NumOutLines)=jtryout(NumOutLines)*0.5d0*
     $     (Mremnant(NumOutLines)-Mremnant(NumOutLines-1))
      Do Jr=NumOutLines-1,2,-1
         unscaledJtotout(Jr)=unscaledJtotout(Jr+1)+
     $        jtryout(Jr)*0.5d0*(Mremnant(Jr+1)-Mremnant(Jr-1))
      enddo
      unscaledJtotout(1)=unscaledJtotout(2)+
     $     jtryout(1)*0.5d0*(Mremnant(1)+Mremnant(2))
         
c Jtotin keeps track of the total angular momentum for the shells
c inside (and including) the shell in question:
      Jtotin=jtryin(1)*0.5d0*(Mremnant(1)+Mremnant(2))
      do Jr=2,NumOutLines
C get k_2 by matching slope of j...
         k2=(jtryin(Jr)-jtryin(Jr-1))/(jtryout(Jr)-jtryout(Jr-1))

C get k_3 by matching j itself
         k3=0.5d0*(
     $        jtryin(Jr-1)+jtryin(Jr)
     $        -k2*(jtryout(Jr-1)+jtryout(Jr))
     $        )

         scaledJtot=Jtotin+k2*unscaledJtotout(Jr)+
     $        k3*(
     $        Mremnant(NumOutLines)-0.5d0*(Mremnant(Jr-1)+Mremnant(Jr))
     $        )

c         write(6,*) Jr, scaledJtot,TotJrem,
c     $        0.5d0*(Mremnant(Jr)+Mremnant(Jr-1))/Mremnant(NumOutLines)

         if ((scaledJtot-TotJrem)*diff .le. 0.d0 .and. Jr.gt.2) then
            if(dabs(scaledJtot-TotJrem).lt.dabs(diff)) then
               Jrtransition=Jr
            else
               Jrtransition=Jr-1
            endif
            if(verbosity.ge.2) then
               write(6,*) 'Transition occurs near m/M_r=',
     $              0.5d0*
     $              (Mremnant(Jrtransition)+Mremnant(Jrtransition-1))
     $              /Mremnant(NumOutLines)
            endif

C tweak k2 so that get exactly the desired total angular momentum
            k2=(TotJrem-Jtotin-k3*(
     $        Mremnant(NumOutLines)-0.5d0*(Mremnant(Jr-1)+Mremnant(Jr))
     $        ))/unscaledJtotout(Jr)
            goto 234
         endif

C Figure out Jtotin for the next iteration....
         if(Jr.lt.NumOutLines) then
            Jtotin=Jtotin+0.5d0*jtryin(Jr)*
     $           (Mremnant(Jr+1)-Mremnant(Jr-1))
         else
            Jtotin=Jtotin+0.5d0*jtryin(Jr)*
     $           (Mremnant(Jr)-Mremnant(Jr-1))            
         endif

         diff=scaledJtot-TotJrem

      enddo

      if (verbosity.ge.2) then
         write(6,*) 'Using only the form for the outer layer'
      endif
      Jrtransition=1
      k3=0.d0
      k2=TotJrem/unscaledJtotout(1)

 234  continue
      if(verbosity.ge.2) then
         write(6,*)'k2=',k2,' k3=',k3
      endif

      Do Jr=1,Jrtransition-1
         jremnant(Jr)=jtryin(Jr)
      enddo
      Do Jr=Jrtransition,NumOutLines
         jremnant(Jr)=k2*jtryout(Jr)+k3
      enddo
      
      scaledJtot=jremnant(1)*0.5d0*(Mremnant(1)+Mremnant(2))
      do Jr=2, NumOutLines-1
         scaledJtot=scaledJtot+jremnant(Jr)*0.5d0*
     $        (Mremnant(Jr+1)-Mremnant(Jr-1))
      enddo
      scaledJtot=scaledJtot+jremnant(NumOutLines)*0.5d0*
     $     (Mremnant(NumOutLines)-Mremnant(NumOutLines-1))
      if(verbosity.ge.2) then
         write(6,*)'Total J=',scaledJtot,'=',TotJrem
         write(6,*) 'DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif
      END
      
*************************************************************
*************************************************************
*************************************************************
*************************************************************
*************************************************************
*************************************************************
*************************************************************


      FUNCTION DELNRG(INTERCEPT)

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER K,verbosity
      DOUBLE PRECISION W,TotE(N),TotErem,INTERCEPT,DELNRG,MLsfrac
      COMMON/NRG/TotE,TotErem
      COMMON/MASSLOSS/MLsfrac
      COMMON/verbiage/verbosity

      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) '***minimizing to find intercept***'
      endif

      CALL SHOCK(INTERCEPT)
      CALL SORT
      CALL PROFILE

      W=0.d0
      DO K=1,N
         W=TotE(K)+W
      ENDDO

      DELNRG=TotErem/W-1.d0-const7*MLsfrac
      if(verbosity.ge.2)
     & write(6,*)'delta energy=',DELNRG,'(we want this to become small)'

      RETURN
      END
      SUBROUTINE ENERGY(M,P,D,R,Jmax,STAR)

*************************************************************
* This routine will calculate the total energy of a star.   *
*************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER J,Jmax,STAR,T,verbosity
      DOUBLE PRECISION P(NROWS),D(NROWS),M(NROWS),R(NROWS),
     &     Ei,Eg,delM,TotE(N),TotErem
      COMMON/NRG/TotE,TotErem
      COMMON/verbiage/verbosity

* A star's total energy consists of its gravitational energy, Eg,
* and its internal energy, Ei.  These are calculated for the star's
* innermost point and then numerically integrated outward.

 100  Ei=3.d0/2.d0*P(1)/D(1)*M(1)
      Eg=-G*(4.d0*D(1)*pi)**2.d0/15.d0*R(1)**5.d0

      T=2
      IF(R(1).EQ.0.d0) THEN
         Ei=Ei+3.d0/2.d0*P(2)/D(2)*M(2)
         Eg=Eg-G*(4.d0*D(2)*pi)**2.d0/15.d0*R(2)**5.d0
         T=3
      ENDIF

      DO J=T,Jmax
         delM=M(J)-M(J-1)
         Eg=Eg-G*delM*(M(J)/R(J)+M(J-1)/R(J-1))*0.5d0
         Ei=Ei+3.d0/4.d0*(P(J)/max(D(J),1.d-40)+
     $        P(J-1)/max(D(J-1),1.d-40))*delM
      ENDDO

      IF(STAR.GT.N) THEN
         if(verbosity.ge.2)
     &        WRITE(6,*) 'Total energy for remnant =',Eg+Ei,Eg,Ei
         TotErem=Eg+Ei
      ELSE
         if(verbosity.ge.2)
     &        WRITE(6,*) 'Total energy for star',STAR,' =',Eg+Ei
      ENDIF

* If all goes well, the Virial Theorem should hold and 2Ei+Eg=0.

      if(verbosity.ge.2) WRITE(6,*) '2Ei+Eg=',2.d0*Ei+Eg
      
      TotE(STAR)=Eg+Ei

      if(dabs((2.d0*Ei+Eg)/TotE(STAR)).gt.1d-2) then
         write(6,*)
         write(6,*) 'WARNING ... WARNING ... WARNING'
         write(6,*) 'virial/(total energy)=',(2.d0*Ei+Eg)/TotE(STAR)
         write(6,*) 'which seems to be quite far away from 0.'
         write(6,*) 'This ratio would be zero for a star in'
         write(6,*) 'hydrostatic equilibrium.'
         write(6,*) 'Internal energy=',Ei
         write(6,*) 'Gravitational potential energy=',Eg
         if(STAR.GT.N) then
            write(6,*) 'The remnant seems to be far from equilibrium,'
            write(6,*)'which is OK as long as this is not the last call'
            write(6,*) 'to PROFILE....  If it is the last call,'
            write(6,*) 'you should probably request larger arrays'
            write(6,*) 'for your remnant profiles.  If this does not'
            write(6,*) 'help, consider emailing lombardi@vassar.edu'
         else
            write(6,*) 'Parent star',STAR,' seems to be far from'
            write(6,*) 'equilibrium.  Are you sure you the input'
            write(6,*) 'data file for this star has a reasonably'
            write(6,*) 'large number of rows in it?'
c            pause
         endif
      endif
      
      END
      SUBROUTINE LOSEMASS

******************************************************************
* This routine will estimate the mass loss for each parent star. *
* Any data beyond what is "lost" will be ignored.                *
******************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jpmax(N),Pmax(N),I50(N),verbosity,I86(N),I95(N),K
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),
     &     ELp(NROWS,elmax,N),
     &     MLsfrac,Dp(NROWS,N),Ap(NROWS,N),MLoseTotal,
     &     periastron,As(NROWS,N),ELm(NROWS,elmax,N),
     &     Rp50(N),Rp86(N),Rp95(N)
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/MASSLOSS/MLsfrac
      COMMON/PERI/periastron
      COMMON/verbiage/verbosity
      COMMON/index/ I50
      SAVE Rp86, Rp95

      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) 'Entering LOSEMASS FOR PERIASTRON=',periastron
         if(periastron.eq.0.d0)
     &    write(6,*)'EVEN IF YOU DID NOT REQUEST A HEADON COLLISION',
     &        'WE STILL HAVE TO WORK ON THE HEADON CASE FOR A WHILE.'
      endif

      IF(N.NE.2) THEN
         write(6,*) 'WARNING: MASS LOSS DOES NOT WORK UNLESS THERE ARE',
     &            ' TWO STARS'
         STOP
      ENDIF

* There should be no need to find the 50% radius if periastron is non-zero,
* because that should have already been done when treating the periastron=0
* case:
      if(periastron.eq.0.d0) then
         call geti50(I50,I86,I95)
         DO K=1,N
            Rp50(K)=((Mp(I50(K),K)/Mp(Jpmax(K),K)-.5d0)*Rp(I50(K)-1,K)+
     $           (.5d0-Mp(I50(K)-1,K)/Mp(Jpmax(K),K))*Rp(I50(K),K))/
     $           (Mp(I50(K),K)-Mp(I50(K)-1,K))*Mp(Jpmax(K),K)
            Rp86(K)=((Mp(I86(K),K)/Mp(Jpmax(K),K)-.86d0)*Rp(I86(K)-1,K)+
     $           (.86d0-Mp(I86(K)-1,K)/Mp(Jpmax(K),K))*Rp(I86(K),K))/
     $           (Mp(I86(K),K)-Mp(I86(K)-1,K))*Mp(Jpmax(K),K)
            Rp95(K)=((Mp(I95(K),K)/Mp(Jpmax(K),K)-.95d0)*Rp(I95(K)-1,K)+
     $           (.95d0-Mp(I95(K)-1,K)/Mp(Jpmax(K),K))*Rp(I95(K),K))/
     $           (Mp(I95(K),K)-Mp(I95(K)-1,K))*Mp(Jpmax(K),K)
        ENDDO
         if(verbosity.ge.2) then
            write(6,*)'R_{0.5} are',Rp50(1)/Rsun,' and',
     $           Rp50(2)/Rsun,' solar radii *** using these numbers'
            write(6,*)'R_{0.86} are',Rp86(1)/Rsun,' and',
     $           Rp86(2)/Rsun,' solar radii *** using these numbers'
            write(6,*)'R_{0.95} are',Rp95(1)/Rsun,' and',
     $           Rp95(2)/Rsun,' solar radii'
         endif
      endif

* Now we calculate the mass lost.  It depends on the 
* reduced mass, the radius of the star in question, and the radii at
* a mass fractions of 50% and 95% in both stars.  The coefficient
* const1 is determined 
* by calculating this information for a number of cases and averaging.

      MLoseTotal=
     &        const1*Mp(Jpmax(1),1)*Mp(Jpmax(2),2)*
     &        (Rp86(1)+Rp86(2))/
     &        ((Mp(Jpmax(1),1)+Mp(Jpmax(2),2))*
     &        (Rp(I50(1),1)+Rp(I50(2),2) +
     &         const2*periastron*(Rp(Jpmax(1),1)+Rp(Jpmax(2),2)) ))
      MLsfrac=MLoseTotal/(Mp(Jpmax(1),1)+Mp(Jpmax(2),2))

      if(verbosity.ge.1) write(6,'(a,f5.3,a,f5.3)')
     &' Mass loss fraction=',MLsfrac,' for periastron=',periastron
      if(verbosity.ge.2) then
         write(6,*)
     &', corresponding to a remnant mass of',(1.d0-MLsfrac)*
     &        (Mp(Jpmax(1),1)+Mp(Jpmax(2),2))/Msun,' Msun'
c         write(6,*)
c     &', corresponding to a table mass of',(1.d0-MLsfrac)*
c     &        ( NINT(100.d0*Mp(Jpmax(1),1)/Msun)/100.d0
c     $         +NINT(100.d0*Mp(Jpmax(2),2)/Msun)/100.d0 ),' Msun'
      endif
      RETURN
      END

******************************************************************
******************************************************************
******************************************************************
      subroutine geti50(I50,I86,I95)

******************************************************************
* This routine finds the locations in the parent inside of which *
* 50%, 86% and 95% of the total mass is enclosed                 *
******************************************************************
      
      implicit none
      INCLUDE 'params.h'
      INTEGER I,K,Jpmax(N),Pmax(N),I50(N),I86(N),I95(N)
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),
     &     ELp(NROWS,elmax,N),
     &     Dp(NROWS,N),Ap(NROWS,N),
     &     As(NROWS,N),ELm(NROWS,elmax,N)
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As

      DO K=1,N
         DO I=1,NROWS
            IF(Mp(I,K).GE.50.d-2*Mp(Jpmax(K),K)) GOTO 47
         ENDDO
 47      I50(K)=I
         DO I=I50(K)+1,NROWS
            IF(Mp(I,K).GE.86.d-2*Mp(Jpmax(K),K)) GOTO 48
         ENDDO
 48      I86(K)=I
         DO I=I86(K)+1,NROWS
            IF(Mp(I,K).GE.95.d-2*Mp(Jpmax(K),K)) GOTO 49
         ENDDO
 49      I95(K)=I
      ENDDO

      END

      SUBROUTINE MIX(K)

*****************************************************************
* This routine will model the hydrodynamic mixing of the        *
* chemical compositions                                         *
* in the parent stars, after shock heating has occurred and     *
* before mass loss takes place.                                 *
*****************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER C,K,L,I,J,Jpmax(N),Max,Pmax(N),verbosity
c     integer halfway
      DOUBLE PRECISION ELp(NROWS,elmax,N),As(NROWS,N),Mp(NROWS,N),
     &     DENOM(NROWS),CONTRIB,dM(NROWS),alpha,
     &     F(NROWS,NROWS),ELm(NROWS,elmax,N),ELmtot,
     &     ELptot,Dp(NROWS,N),Ap(NROWS,N),Rp(NROWS,N),dlogAdMavg
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/verbiage/verbosity
      common/numel/el
      double precision FBOUND(NROWS)
      character*16 OUTFN

      Max=Jpmax(K)

      DO I=2,Max-1
         dM(I)=0.5d0*(Mp(I+1,K)-Mp(I-1,K))
      ENDDO

      dM(1)=0.5d0*(Mp(2,K)+Mp(1,K))
      dM(Max)=0.5d0*(Mp(Max,K)-Mp(Max-1,K))

      dlogAdMavg=dlog(As(Pmax(K),K))-dlog(As(1,K))

      do I=1,Max
         FBOUND(I)=0.d0
      enddo

      alpha=-const8*dlogAdMavg**2.d0

      if(verbosity.ge.2) then
         write(6,*)'Mp(Max,K)=',Mp(Max,K)
         write(6,*)'dlogAdMavg=',dlogAdMavg,', alpha from paper=',-alpha
c Note that the alpha in the paper has the opposite sign from the alpha in
c this code
      endif

      DO I=1,Max
         DENOM(I)=0.d0
         DO J=1,Max
            F(J,I)=exp(alpha*((Mp(J,K)-Mp(I,K))/Mp(Max,K))**2.d0)
     &            +exp(alpha*((Mp(J,K)+Mp(I,K))/Mp(Max,K))**2.d0)
     &            +exp(alpha*((2.d0*Mp(Max,K)-Mp(J,K)
     &                                 -Mp(I,K))/Mp(Max,K))**2.d0)
            DENOM(I)=DENOM(I)+dM(J)*F(J,I)
         ENDDO
      ENDDO

      DO I=1,Max
         DENOM(I)=dM(I)/DENOM(I)
      ENDDO


      if(studymassloss) then
         WRITE(OUTFN,101) K
 101     FORMAT('fbound',I1,'.sbe')
         open(19,file=OUTFN)
         DO I=1,Max
            DO J=1,Max
C     F(J,I)*DENOM(J)*dM(I) is the mass that went from shell J to shell I
               IF(I.LE.Pmax(K))
     $              FBOUND(J)=FBOUND(J)+F(J,I)*DENOM(J)*dM(I)
            ENDDO
         ENDDO
         
         DO I=1,Max
            FBOUND(I)=FBOUND(I)/dM(I)
            write(19,*) Mp(I,K)/Msun,FBOUND(I)
         enddo
         close(19)
      endif

      DO C=1,el
         ELptot=0.d0
         ELmtot=0.d0
         DO L=1,Max
            CONTRIB=0.d0
            DO I=1,Max
               CONTRIB=CONTRIB+ELp(I,C,K)*F(I,L)*DENOM(I)
C Elp(I,C,K)*F(I,L)*DENOM(I)*dM(L) is the amount of mass of species C that
C went from shell I to shell L in star K
            ENDDO
            ELm(L,C,K)=CONTRIB
            ELptot=ELptot+ELp(L,C,K)*dM(L)
            ELmtot=ELmtot+ELm(L,C,K)*dM(L)
         ENDDO
         if(verbosity.ge.2) WRITE(6,*) 'Element:',C

         IF((ELptot-ELmtot)/ELptot.GT.1.d-6) THEN
            WRITE(6,*) 'WARNING!!  Mass may NOT be conserved!'
            WRITE(6,*) 'ELptot=',ELptot,' ELmtot=',ELmtot
         ENDIF

      ENDDO

      RETURN
      END

      SUBROUTINE PROFILE

********************************************************************
* Given as input the remnant file made by SORT with M, A, chemical *
* data, this subroutine will find values for the remnant's         *
* pressure,  density, and radius.                                  *
********************************************************************

* All units are cgs.

* Define variables: zeroin2 is a root-finding function
* Pclow & Pchigh are lower and upper limits for the guess of the central
* pressure Pcntr

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jrmax,STAR,NP,verbosity
      DOUBLE PRECISION Pcntr,Pclow,Pchigh,Psurf
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),P(NROWS),
     &     D(NROWS),R(NROWS),ELr(NROWS,elmax)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
      COMMON/verbiage/verbosity
      DOUBLE PRECISION PRESSURE,zeroin2
      EXTERNAL PRESSURE
      integer kmax,kount,NMAX
      parameter(NMAX=50)
      double precision dxsav, xp(NROWS),yp(NROWS,NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      double precision
     $     pressure2b(NROWS),pressure2c(NROWS),pressure2d(NROWS),
     $     radius2b(NROWS),radius2c(NROWS),radius2d(NROWS)
      double precision seval

* Define the limits for zeroin2:
* These limits would have to be adjusted if you are creating remnants
* with either exceedingly low or exceedingly high central pressures...
      Pclow=0.0000000000001d0*G*Msun**2.d0/(Rsun**4.d0)
      Pchigh=3000.d0*G*Msun**2.d0/(Rsun**4.d0)

      kmax=0
      Pcntr=zeroin2(Pclow,Pchigh,PRESSURE,1.d-18)
      if(verbosity.ge.2) WRITE(6,*) 'Central pressure=',Pcntr
   
c      stop

*** DO NOT REMOVE THE FOLLOWING ASSIGNMENT STATEMENT - IT IS NEEDED!!***
*** ONE FINAL CALL TO PRESSURE IS NEEDED TO GET CORRECT REMNANT PROFILES!!***
      kmax=NROWS
      dxsav=(Mr(Jrmax)-Mr(1))/NROWS
      Psurf=PRESSURE(Pcntr)

C      print *,'kount=',kount,xp(1),yp(1,1)

      call spline(kount,xp,yp(1,1),pressure2b,pressure2c,pressure2d)
      call spline(kount,xp,yp(1,2),radius2b,radius2c,radius2d)

      P(1)=Pcntr
      D(1)=( P(1)/Ar(1) )**(3.d0/5.d0)
      R(1)=(3.d0*Mr(1)/(4.d0*pi*D(1)))**(1.d0/3.d0)
      
      DO Jr=2,Jrmax

         P(Jr)=max( seval(kount,Mr(Jr),xp,yp(1,1),
     $          pressure2b,pressure2c,pressure2d),    1.d-40)

         R(Jr)=seval(kount,Mr(Jr),xp,yp(1,2),radius2b,radius2c,radius2d)
         D(Jr)=( P(Jr)/Ar(Jr) )**(3.d0/5.d0)

      ENDDO

      if(verbosity.ge.2) WRITE(6,*) 'Surface pressure=',Psurf

      Jr=1
      STAR=N+1
      CALL ENERGY(Mr(Jr),P(Jr),D(Jr),R(Jr),Jrmax,STAR)

      END

**************************************************************
      FUNCTION PRESSURE(Pctry)

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jrmax,T,NP,Jr2
      DOUBLE PRECISION dM,Pctry,Rup
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),P(NROWS),D(NROWS),
     &     R(NROWS),ELr(NROWS,elmax)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
      DOUBLE PRECISION PRESSURE
      double precision Rlow
      EXTERNAL derivs
      double precision vstart(2)
      double precision relerr(2),abserr(2)
      integer iflag
      double precision work(17),ma,mb
      integer iwork(5),ii

      T=2

      P(1)=Pctry
      D(1)=(dabs(Pctry)/Ar(1))**(3.d0/5.d0)
      R(1)=(3.d0*Mr(1)/(4.d0*pi*D(1)))**(1.d0/3.d0)
      
      IF(R(1).EQ.0.d0) THEN
         dM=Mr(2)-Mr(1)
         R(2)=(3.d0*dM/(4.d0*pi*D(1)))**(1.d0/3.d0)
         P(2)=Pctry-2.d0*pi/3.d0*G*R(2)**2.d0*D(1)**2.d0
         D(2)=(dabs(P(2))/Ar(2))**(3.d0/5.d0)
         T=3
         print *,'Check P(2) equation.....'
         stop
      ENDIF

      vstart(1)=P(T-1)
      vstart(2)=R(T-1)

C If you make relerr=1.d-5 or smaller, get spikes in radius vs rp!
C 1.d-5, 1.d-10 gives a spike
C 2.d-5, 1.d-10 looks good
C 2.d-5, 1.d-11 looks good again
C 2.d-5, 1.d-9  best yet
C 2.d-5, 1.d-8 looks good, but 1.d-9 was perhaps a little better
C 3.d-5, 1.d-10 too much noise on top curve

C  7/1/03
C 2.1d-5, 2.1d-5, 5.d-10, 5.d-10 ... many virial warnings...
C 2.1d-5, 3.d-5,  1.d-10, 1.d-10 ... many virial warnings...


      relerr(1)=2.1d-5
      relerr(2)=1.d-4
      abserr(1)=1.d-10
      abserr(2)=1.d-10
      iflag=1
      ma=Mr(T-1)
      mb=Mr(Jrmax)

 22   continue
C      print *,'vstart=',vstart,ma,mb,iflag
      call rkf45(derivs,2,vstart,ma,mb,relerr,abserr,iflag,
     $        work,iwork)

c        print *,'iflag=',iflag
c        stop
         if(iflag.lt.0) then
            iflag=1
            print *,'iflag was negative'
            stop
         endif
         if (iflag.gt.2) then
            print *,'iflag=',iflag
            if(iflag.eq.6) then
               iflag=2
               relerr(1)=relerr(1)*2.d0
               relerr(2)=relerr(2)*2.d0
               print *,'got to',ma,mb,Mr(T-1),Mr(Jrmax),relerr
               goto 22
            endif
            if(iflag.eq.7) then
               iflag=-1
            endif
            if(iflag.eq.8) then
               print *,Jr,Jrmax,Mr(Jr-2),Mr(Jr-1),Mr(Jr),relerr
               print *
               print *,(Mr(ii),ii=1,Jrmax,100)
            endif
            stop
         endif

      PRESSURE=vstart(1)

c      print *, 'PRESSURE=',PRESSURE,'Radius=',vstart(2),Pctry

      RETURN
      END


*****************************************************************
      FUNCTION RADIUS(RJr)

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jrmax,NP
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),P(NROWS),R(NROWS)
      DOUBLE PRECISION PJr,RJr,RADIUS,ELr(NROWS,elmax),D(NROWS)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
c      double precision RADIUSorig

      PJr=P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0+
     &     Mr(Jr)/RJr**4.d0)*(Mr(Jr)-Mr(Jr-1))

cc      RADIUSorig=-RJr+R(Jr-1)+(Mr(Jr)-Mr(Jr-1))/(8.d0*pi)*
cc     &     ((Ar(Jr-1)/P(Jr-1))**(3.d0/5.d0)/
cc     &     (R(Jr-1)**2.d0)+(Ar(Jr)/dabs(PJr))**(3.d0/5.d0)/(RJr**2.d0))
      RADIUS=-RJr+R(Jr-1)+(Mr(Jr)-Mr(Jr-1))/(4.d0*pi)*
     &     ( Ar(Jr-1)*Ar(Jr)/(dabs(P(Jr-1))*dabs(PJr) ))**0.3d0/
     &     (R(Jr-1)*RJr)

C      if(PJr.lt.0.d0) print *,'BAD'

c      print *,'radii=',RADIUSorig,RADIUS,RJr
c     $     -RJr,R(Jr-1),(Mr(Jr)-Mr(Jr-1))/(4.d0*pi),
c     &     ( Ar(Jr-1)*Ar(Jr)/(P(Jr-1)*dabs(P(Jr))))**0.3d0,
c     &     (R(Jr-1)*RJr),
c     &     Ar(Jr-1)*Ar(Jr),(P(Jr-1)*dabs(P(Jr)))

      RETURN
      END






      SUBROUTINE SHOCK(INTERCEPT)

*********************************************************************
* This routine will take the parent star entropy profiles and shock *
* them, creating a new entropy profile for SORT to accept and use   *
* to calculate the remnant's entropy.                               *
*********************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'

* Mto= turn-off mass in solar units
* Rto= turn-off radius in solar units
* (The use of Mto and Rto is simply because of the chosen units for the
* intercept.)
      DOUBLE PRECISION Mto,Rto
      PARAMETER(Mto=0.8d0*Msun, Rto=0.955d0*Rsun)

      INTEGER Jpmax(N),I,Imax,K,Pmax(N),I1,I2,verbosity
      DOUBLE PRECISION Ap(NROWS,N),As(NROWS,N),Dp(NROWS,N),INTERCEPT,
     &     Mp(NROWS,N),Rp(NROWS,N),ELp(NROWS,elmax,N),INTERCEPT2,
     &     MLsfrac,MLoseTotal,massoutsideofI1,massoutsideofI2,
     &     periastron,ELm(NROWS,elmax,N)
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/MASSLOSS/MLsfrac
      COMMON/PERI/periastron
      COMMON/verbiage/verbosity
      double precision RHS

* The equation for the shocked entropy, As, comes from a plot of
* Log(A-Ainit) vs. Log(Arho^5/3).  The value INTERCEPT is the y-intercept
* of the linear function fitted to the plot of the data.  The variable
* "const3" is the slope of that same function.

      DO K=1,N
         Imax=Jpmax(K)
         IF(K.EQ.1) THEN
            DO I=1,Imax
               RHS=INTERCEPT + dlog10(G*Mto**(1.d0/3.d0)*Rto)
     &              +const3*
     $              dlog10(Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))
               As(I,K)=Ap(I,K)+10.d0**RHS
c               As(I,K)=Ap(I,K)+10.d0**INTERCEPT*
c     &          G*Mto**(1.d0/3.d0)*Rto*
c     &          (Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))**const3
            ENDDO
         ELSE
            INTERCEPT2=INTERCEPT+((const4+const5)*periastron-const6)*
     &           LOG10(Mp(Jpmax(1),1)/Mp(Jpmax(K),K))
            DO I=1,Imax
               RHS=INTERCEPT2 + dlog10(G*Mto**(1.d0/3.d0)*Rto)
     &              +const3*
     $              dlog10(Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))
c               As(I,K)=Ap(I,K)+10.d0**INTERCEPT2*
c     &          G*Mto**(1.d0/3.d0)*Rto*
c     &          (Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))**const3
               As(I,K)=Ap(I,K)+10.d0**RHS
            ENDDO
            if(verbosity.ge.2) then
               write(6,*) 'Intercept for star 2=',INTERCEPT2
               write(6,*) 'Intercept for star 1=',INTERCEPT
            endif
         ENDIF
      ENDDO

      MLoseTotal=MLsfrac*(Mp(Jpmax(1),1)+Mp(Jpmax(2),2))
      I1=Jpmax(1)
      I2=Jpmax(2)
      massoutsideofI1=0.d0
      massoutsideofI2=0.d0
      Do while (massoutsideofI1+massoutsideofI2 .lt. MLoseTotal)
         I1=I1-1
         massoutsideofI1=massoutsideofI1 + Mp(I1+1,1) - Mp(I1,1)

         I2=Jpmax(2)
         massoutsideofI2=0.d0
         do while(As(I2,2).gt.As(I1,1))
           I2=I2-1
           massoutsideofI2=massoutsideofI2 + Mp(I2+1,2) - Mp(I2,2)
         enddo
      enddo

      Pmax(1)=I1
      Pmax(2)=I2

      if(verbosity.ge.2) then
         write(6,*)'Pmax(1)=',Pmax(1),' Pmax(2)=',Pmax(2)
         write(6,*)'actual mass loss fraction=',
     &        (massoutsideofI1+massoutsideofI2)/
     &        (Mp(Jpmax(1),1)+Mp(Jpmax(2),2)),
     &        ' for r_p=',periastron
      endif

      RETURN
      END

      
      SUBROUTINE SORT

******************************************************************
* This subroutine will take the parent star data                 *
* of M, A and chemical abundances.  It will create similar arrays*
* for the remnant, sorted by increasing A.  In order             *
* to keep the remnant profile smooth, we'll account for the      *
* profile derivatives.                                           *
******************************************************************

* Define the variables:

* M is the mass enclosed in a star, A is the entropic variable at
* M, and EL is the chemical abundance at M.  Ji is the Jith row in 
* the data file for a star.

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jp(N),Jrmax,K,C,Pmax(N),Jpmax(N),NP
      DOUBLE PRECISION Mp(NROWS,N),ELp(NROWS,elmax,N),As(NROWS,N)
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),ELr(NROWS,elmax),Rp(NROWS,N)
      DOUBLE PRECISION dAr,dAsdMp(N),dArdMr,dArp(N),Mtmp,dAdMtmp,
     &     Armax,MofAr(N),Chemtmp(elmax),dArtmp,Chem(elmax,N),
     &     ELm(NROWS,elmax,N),
     &     Dp(NROWS,N),Ap(NROWS,N)
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      LOGICAL NOTDONE
      common/numel/el
      double precision monm(NROWS),dMp(NROWS),
     $     y2b(NROWS),y2c(NROWS),y2d(NROWS)
      common/delms/monm,dMp,y2b,y2c,y2d
      double precision dmr
      double precision seval
      double precision z2b(NROWS),z2c(NROWS),z2d(NROWS)
      common/splinearrays/z2b,z2c,z2d

* The maximum entropy in the remnant should be the highest in any of
* the parents:

      Armax=0.d0
      DO K=1,N
         Armax=max(Armax,As(Pmax(K),K))
      ENDDO

c      write(6,*)'Armax=',Armax,' Pmax(1)=',Pmax(1),' Pmax(2)=',Pmax(2)

* Initialize the rows:

      DO K=1,N
         Jp(K)=1
      ENDDO
      Jr=1

* Initialize the entropy - the first value of Ar (at the lowest Mr)
* will be the lowest entropy in any of the parent stars.

      Ar(1)=1.d30
      DO K=1,N
         Ar(1)=min(Ar(1),As(1,K))
      ENDDO

C      dlogAr=(DLOG10(Armax)-DLOG10(Ar(1)))/NP

* Define the logical variable that will let the program leave the
* DO loop at the right time, and initialize dAr.

      NOTDONE=.TRUE.

* Begin the DO loop.

      DO WHILE (NOTDONE)

* If Ar gets to be the maximum, then we're done. (hence the IF statement
* below).

c         IF(Ar(Jr).GE.0.999999d0*Armax) NOTDONE=.FALSE.
         IF(Ar(Jr).GE.Armax) NOTDONE=.FALSE.

* We must do a loop over each parent star.

         DO K=1,N
            DO WHILE (As(Jp(K),K).LE.Ar(Jr) .AND. Jp(K).LE.Pmax(K))
               Jp(K)=Jp(K)+1
            ENDDO

            IF(Jp(K).GE.2 .AND. Jp(K).LE.Pmax(K)) THEN
               dAsdMp(K)=(As(Jp(K),K)-As(Jp(K)-1,K))/
     &              (Mp(Jp(K),K)-Mp(Jp(K)-1,K))
               MofAr(K)=(Ar(Jr)-As(Jp(K)-1,K))*Mp(Jp(K),K)/(As(Jp(K),K)-
     &              As(Jp(K)-1,K))+(As(Jp(K),K)-Ar(Jr))*Mp(Jp(K)-1,K)/
     &              (As(Jp(K),K)-As(Jp(K)-1,K))
               dArp(K)=As(Jp(K),K)-As(Jp(K)-1,K)
               DO C=1,el
                  Chem(C,K)=(Ar(Jr)-As(Jp(K)-1,K))*ELm(Jp(K),C,K)/
     &                 (As(Jp(K),K)-As(Jp(K)-1,K))+(As(Jp(K),K)-Ar(Jr))
     &                 *ELm(Jp(K)-1,C,K)/(As(Jp(K),K)-As(Jp(K)-1,K))
               ENDDO
            ELSE IF(Jp(K).EQ.1) THEN
               dAsdMp(K)=1.d30
               DO C=1,el
                  Chem(C,K)=0.d0
               ENDDO
               dArp(K)=1.d30
               MofAr(K)=0.d0
            ELSE IF(Jp(K).GT.Pmax(K)) THEN
               if(Ar(Jr).le.As(Jpmax(K),K)) then
                  dAsdMp(K)=(As(Pmax(K),K)-As(Pmax(K)-1,K))/
     &                 (Mp(Pmax(K),K)-Mp(Pmax(K)-1,K))
               else
                  dAsdMp(K)=1.d30
               endif
               DO C=1,el
                  Chem(C,K)=ELm(Pmax(K),C,K)
               ENDDO
               dArp(K)=1.d30
               MofAr(K)=Mp(Pmax(K),K)
            ENDIF

         ENDDO

         dAdMtmp=0.d0
         DO C=1,el
            Chemtmp(C)=0.d0
         ENDDO
         Mtmp=0.d0
         dArtmp=1.d30
         DO K=1,N
            dAdMtmp=1.d0/dAsdMp(K)+dAdMtmp
            DO C=1,el
               Chemtmp(C)=Chem(C,K)/dAsdMp(K)+Chemtmp(C)
            ENDDO
            Mtmp=MofAr(K)+Mtmp
            dArtmp=min(dArp(K),dArtmp)
         ENDDO

         dArdMr=dAdMtmp**(-1.d0)
         DO C=1,el
            ELr(Jr,C)=dArdMr*Chemtmp(C)
         ENDDO
         Mr(Jr)=Mtmp

         IF(NP.lt.0) then
C N.B.!  The NP here is 1 less than the NP in the driver program...
C Use the smallest reasonable dAr if NP<0 : 
            dAr=dArtmp/abs(NP)

            dmr=seval(189,Mr(Jr)/(Mp(Pmax(1),1)+Mp(Pmax(2),2)),
     $           monm,dMp,y2b,y2c,y2d)
            if(dmr.le.0.d0) then
               write(6,*) 'Oops...'
               stop
            endif
            
            dAr=min(dAr,dArdMr*dmr*(Mp(Pmax(1),1)+Mp(Pmax(2),2)))

         else

            dAr=dArdMr*(Mp(Pmax(1),1)+Mp(Pmax(2),2))/(NP-3.5d0)

C These lines of data will have the *logarithm*
C of the entropic variable A spaced regularly...
c         dAr=10.d0**(DLOG10(Ar(Jr))+dlogAr)-Ar(Jr)

         endif

         IF(Ar(Jr)+dAr.GT.Armax) THEN
            dAr=Armax-Ar(Jr)
         ENDIF       

         Jr=Jr+1
         Ar(Jr)=Ar(Jr-1)+dAr

         IF(Jr.GT.NROWS) THEN
            write(6,*) 'WARNING!  You may need to increase NROWS in'
            write(6,*) 'the header file.  NROWS=',NROWS
c            write(6,*) Ar(Jr-1),dAr,Ar(Jr),Armax,dmr
c            write(6,*) Pmax(1),Pmax(2),Jpmax(1),Jpmax(2)
            STOP
         ENDIF

      ENDDO

      Jrmax=Jr-1

c         call spline(189,monm,dMp,y2b,y2c,y2d)

      call spline(Jrmax,Mr,Ar,z2b,z2c,z2d)


      END

*************************************************************
* zeroin minimizes the radius                               *
*************************************************************
      double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol
C
C This routine is a very slightly modified version of a routine
C from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c


c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = 1.0d0 + eps
c      if (tol1 .gt. 1.0d0) go to 10

      eps=3.d-16

c      PRINT *,'EPS=',eps

c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
      if(fa*fb.gt.0.d0)then
         write(6,*)
         write(6,*)'Maybe for your particular computer you need to'
         write(6,*)'make the parameter eps a little larger in'
         write(6,*)'the zeroin, zeroin2, and/or zeroin3 routines.'
         write(6,*)'Also consider emailing lombardi@vassar.edu'
         print *, a,b,fa,fb
         pause 'root must be bracketed for zeroin'
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end

*************************************************************
* zeroin2 minimizes the pressure                            *
*************************************************************
      double precision function zeroin2(ax,bx,f,tol)
      double precision ax,bx,f,tol
C
C This routine is a very slightly modified version of a routine
C from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin2 abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin2  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign

c
c  compute eps, the relative machine precision
c


c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = 1.0d0 + eps
c      if (tol1 .gt. 1.0d0) go to 10

c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = bx + eps
c      if (tol1 .gt. bx) go to 10

      eps=3.d-16

c      print *,'eps=',eps,bx


c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
      if(fa*fb.gt.0.d0)then
         write(6,*)
         write(6,*)'You must be pushing the limits of this code.'
         write(6,*)'Try expanding the range of Pclow and Pchigh'
         write(6,*)'and see if that helps you get a model...'
         write(6,*)'Also consider emailing lombardi@vassar.edu'
         pause 'root must be bracketed for zeroin2'
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (dabs(xm) .le. tol1) then
         go to 90
      endif

      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c

c   90 zeroin2 = b

 90   if (fb.le.0.d0) then
         zeroin2=c
      else
         zeroin2=b
      endif
      return
      end

*************************************************************
* zeroin3 used to minimize delta Energy                     *
*************************************************************
      double precision function zeroin3(ax,bx,f,tol)
      double precision ax,bx,f,tol
C
C This routine is a very slightly modified version of a routine
C from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin3 abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin3  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c


c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = 1.0d0 + eps
c      if (tol1 .gt. 1.0d0) go to 10

c      print *,'eps=',eps
c      stop
      eps=3.d-16




c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
      if(fa*fb.gt.0.d0)then
         write(6,*)
         write(6,*)
         write(6,*)'You must be pushing the limits of this code.'
         write(6,*)'Try expanding the range of INTlow and INThigh'
         write(6,*)'and see if that helps you get a model...'
         write(6,*)'Also consider emailing lombardi@vassar.edu'
         pause 'root must be bracketed for zeroin3'
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin3 = b
      return
      end

      subroutine spline (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
C This routine is a very slightly modified version of one from the book
C Computer Methods for Mathematical Computations, by George Forsythe, Mike
C Malcolm, and Cleve Moler. The original copy was obtained from the Netlib
C mathematical software file server
C (see http://www2.ucsc.edu/cats/sc/software/netlib/) with the command
C "send spline from fmm".
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.d0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.d0
      c(n) = 0.d0
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.d0*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.d0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.d0*c(i)
   40 continue
      c(n) = 3.d0*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.d0
      d(1) = 0.d0
      b(2) = b(1)
      c(2) = 0.d0
      d(2) = 0.d0
      return
      end

      function seval(n, u, x, y, b, c, d)
      implicit none
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n), seval
C This routine is from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C This copy was obtained from the Netlib mathematical software file server
C (see http://www2.ucsc.edu/cats/sc/software/netlib/) with the command
C "send seval from fmm".
c
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k
      double precision dx
      data i/1/

      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end
      subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c     written by h.a.watts and l.f.shampine
c                   sandia laboratories
c                  albuquerque,new mexico
c
c    rkf45 is primarily designed to solve non-stiff and mildly stiff
c    differential equations when derivative evaluations are inexpensive.
c    rkf45 should generally not be used when the user is demanding
c    high accuracy.
c
c abstract
c
c    subroutine  rkf45  integrates a system of neqn first order
c    ordinary differential equations of the form
c             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
c              where the y(i) are given at t .
c    typically the subroutine is used to integrate from t to tout but it
c    can be used as a one-step integrator to advance the solution a
c    single step in the direction of tout.  on return the parameters in
c    the call list are set for continuing the integration. the user has
c    only to call rkf45 again (and perhaps define a new value for tout).
c    actually, rkf45 is an interfacing routine which calls subroutine
c    rkfs for the solution.  rkfs in turn calls subroutine  fehl which
c    computes an approximate solution over one step.
c
c    rkf45  uses the runge-kutta-fehlberg (4,5)  method described
c    in the reference
c    e.fehlberg , low-order classical runge-kutta formulas with stepsize
c                 control , nasa tr r-315
c
c    the performance of rkf45 is illustrated in the reference
c    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
c                 differential equations-the state of the art ,
c                 sandia laboratories report sand75-0182 ,
c                 to appear in siam review.
c
c
c    the parameters represent-
c      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
c      neqn -- number of equations to be integrated
c      y(*) -- solution vector at t
c      t -- independent variable
c      tout -- output point at which solution is desired
c      relerr,abserr -- relative and absolute error tolerances for local
c            error test. at each step the code requires that
c                 abs(local error) .le. relerr*abs(y) + abserr
c            for each component of the local error and solution vectors
c      iflag -- indicator for status of integration
c      work(*) -- array to hold information internal to rkf45 which is
c            necessary for subsequent calls. must be dimensioned
c            at least  3+6*neqn
c      iwork(*) -- integer array used to hold information internal to
c            rkf45 which is necessary for subsequent calls. must be
c            dimensioned at least  5
c
c
c  first call to rkf45
c
c    the user must provide storage in his calling program for the arrays
c    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
c    in the call list  -      y(neqn) , work(1+8*neqn) , iwork(5)  ,
c    declare f in an external statement, supply subroutine f(t,y,yp) and
c    initialize the following parameters-
c
c      neqn -- number of equations to be integrated.  (neqn .ge. 1)
c      y(*) -- vector of initial conditions
c      t -- starting point of integration , must be a variable
c      tout -- output point at which solution is desired.
c            t=tout is allowed on the first call only, in which case
c            rkf45 returns with iflag=2 if continuation is possible.
c      relerr,abserr -- relative and absolute local error tolerances
c            which must be non-negative. relerr must be a variable while
c            abserr may be a constant. the code should normally not be
c            used with relative error control smaller than about 1.e-8 .
c            to avoid limiting precision difficulties the code requires
c            relerr to be larger than an internally computed relative
c            error parameter which is machine dependent. in particular,
c            pure absolute error is not permitted. if a smaller than
c            allowable value of relerr is attempted, rkf45 increases
c            relerr appropriately and returns control to the user before
c            continuing the integration.
c      iflag -- +1,-1  indicator to initialize the code for each new
c            problem. normal input is +1. the user should set iflag=-1
c            only when one-step integrator control is essential. in this
c            case, rkf45 attempts to advance the solution a single step
c            in the direction of tout each time it is called. since this
c            mode of operation results in extra computing overhead, it
c            should be avoided unless needed.
c
c
c  output from rkf45
c
c      y(*) -- solution at t
c      t -- last point reached in integration.
c      iflag = 2 -- integration reached tout. indicates successful retur
c                   and is the normal mode for continuing integration.
c            =-2 -- a single successful step in the direction of tout
c                   has been taken. normal mode for continuing
c                   integration one step at a time.
c            = 3 -- integration was not completed because relative error
c                   tolerance was too small. relerr has been increased
c                   appropriately for continuing.
c            = 4 -- integration was not completed because more than
c                   3000 derivative evaluations were needed. this
c                   is approximately 500 steps.
c            = 5 -- integration was not completed because solution
c                   vanished making a pure relative error test
c                   impossible. must use non-zero abserr to continue.
c                   using the one-step integration mode for one step
c                   is a good way to proceed.
c            = 6 -- integration was not completed because requested
c                   accuracy could not be achieved using smallest
c                   allowable stepsize. user must increase the error
c                   tolerance before continued integration can be
c                   attempted.
c            = 7 -- it is likely that rkf45 is inefficient for solving
c                   this problem. too much output is restricting the
c                   natural stepsize choice. use the one-step integrator
c                   mode.
c            = 8 -- invalid input parameters
c                   this indicator occurs if any of the following is
c                   satisfied -   neqn .le. 0
c                                 t=tout  and  iflag .ne. +1 or -1
c                                 relerr or abserr .lt. 0.
c                                 iflag .eq. 0  or  .lt. -2  or  .gt. 8
c      work(*),iwork(*) -- information which is usually of no interest
c                   to the user but necessary for subsequent calls.
c                   work(1),...,work(neqn) contain the first derivatives
c                   of the solution vector y at t. work(neqn+1) contains
c                   the stepsize h to be attempted on the next step.
c                   iwork(1) contains the derivative evaluation counter.
c
c
c  subsequent calls to rkf45
c
c    subroutine rkf45 returns with all information needed to continue
c    the integration. if the integration reached tout, the user need onl
c    define a new tout and call rkf45 again. in the one-step integrator
c    mode (iflag=-2) the user must keep in mind that each step taken is
c    in the direction of the current tout. upon reaching tout (indicated
c    by changing iflag to 2),the user must then define a new tout and
c    reset iflag to -2 to continue in the one-step integrator mode.
c
c    if the integration was not completed but the user still wants to
c    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
c    the relerr parameter has been adjusted appropriately for continuing
c    the integration. in the case of iflag=4 the function counter will
c    be reset to 0 and another 3000 function evaluations are allowed.
c
c    however,in the case iflag=5, the user must first alter the error
c    criterion to use a positive value of abserr before integration can
c    proceed. if he does not,execution is terminated.
c
c    also,in the case iflag=6, it is necessary for the user to reset
c    iflag to 2 (or -2 when the one-step integration mode is being used)
c    as well as increasing either abserr,relerr or both before the
c    integration can be continued. if this is not done, execution will
c    be terminated. the occurrence of iflag=6 indicates a trouble spot
c    (solution is changing rapidly,singularity may be present) and it
c    often is inadvisable to continue.
c
c    if iflag=7 is encountered, the user should use the one-step
c    integration mode with the stepsize determined by the code or
c    consider switching to the adams codes de/step,intrp. if the user
c    insists upon continuing the integration with rkf45, he must reset
c    iflag to 2 before calling rkf45 again. otherwise,execution will be
c    terminated.
c
c    if iflag=8 is obtained, integration can not be continued unless
c    the invalid input parameters are corrected.
c
c    it should be noted that the arrays work,iwork contain information
c    required for subsequent integration. accordingly, work and iwork
c    should not be altered.
c
c
      implicit none
      integer neqn,iflag,iwork(5)
      double precision y(neqn),t,tout,relerr(neqn),abserr(neqn),work(1)
c     if compiler checks subscripts, change work(1) to work(3+6*neqn)
c     if compiler checks subscripts, change work(1) to work(1+8*neqn)
c
      external f
c
      integer k1,k2,k3,k4,k5,k6,k1m
      integer k7
c
c
c     compute indices for the splitting of the work array
c
      k1m=neqn+1
      k1=k1m+1
      k2=k1+neqn
      k3=k2+neqn
      k4=k3+neqn
      k5=k4+neqn
      k6=k5+neqn
      k7=k6+neqn
c
c     this interfacing routine merely relieves the user of a long
c     calling list via the splitting apart of two working storage
c     arrays. if this is not compatible with the users compiler,
c     he must use rkfs directly.
c
      call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,work(1),work(k1m),
     1          work(k1),work(k2),work(k3),work(k4),work(k5),work(k6),
     2          work(k7),iwork(1),iwork(2),iwork(3),iwork(4),iwork(5))
c
      return
      end
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,
     1                f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c
c     rkfs integrates a system of first order ordinary differential
c     equations as described in the comments for rkf45 .
c     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
c     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
c     internally by the code and appear in the call list to eliminate
c     local retention of variables between calls. accordingly, they
c     should not be altered. items of possible interest are
c         yp - derivative of solution vector at t
c         h  - an appropriate stepsize to be used for the next step
c         nfe- counter on the number of derivative function evaluations
c
c
      implicit none
      logical hfaild,output
c
      integer  neqn,iflag,nfe,kop,init,jflag,kflag
      double precision  y(neqn),t,tout,relerr(neqn),abserr(neqn),h,
     $     yp(neqn),
     1  f1(neqn),f2(neqn),f3(neqn),f4(neqn),f5(neqn),savre(neqn),
     2  savae(neqn)
c
      external f
c
      double precision  a,ae(neqn),dt,ee(neqn),eeoet,esttol(neqn),
     $     et(neqn),hmin,
     $     remin,
     $     rer,s,
     1  scale(neqn),tol,toln,u26,epsp1,eps,ypk
      double precision esttolmax
c
      integer  k,maxnfe,mflag
c
      double precision  dabs,dmax1,dmin1,dsign
c
c     remin is the minimum acceptable value of relerr.  attempts
c     to obtain higher accuracy with this subroutine are usually
c     very expensive and often unsuccessful.
c
      data remin/1.d-12/
c
c
c     the expense is controlled by restricting the number
c     of function evaluations to be approximately maxnfe.
c     as set, this corresponds to about 500 steps.
c
      data maxnfe/3000/

      integer KMAXX,NMAX,kmax,kount
      parameter(KMAXX=4000,NMAX=50)
      double precision dxsav,xp(KMAXX),ysav(KMAXX,NMAX),xsav
      COMMON /path/ kmax,kount,dxsav,xp,ysav
      integer i

C      print *,'enter y=',y

      kount=0
      if (kmax.gt.0) xsav=t-2.d0*dxsav

c
c
c     check input parameters
c
c
      if (neqn .lt. 1) go to 10
      do k=1,neqn
         if ((relerr(k) .lt. 0.0d0)  .or.  (abserr(k) .lt. 0.0d0))
     $        go to 10
      enddo
      mflag=iabs(iflag)
      if ((mflag .eq. 0) .or. (mflag .gt. 8)) go to 10
      if (mflag .ne. 1) go to 20
c
c     first call, compute machine epsilon
c
      eps = 1.0d0
    5 eps = eps/2.0d0
      epsp1 = eps + 1.0d0
      if (epsp1 .gt. 1.0d0) go to 5

      eps=5.d-19

      u26 = 26.0d0*eps
      go to 50
c
c     invalid input
   10 iflag=8
      return
c
c     check continuation possibilities
c
   20 if ((t .eq. tout) .and. (kflag .ne. 3)) go to 10
      if (mflag .ne. 2) go to 25
c
c     iflag = +2 or -2
      if ((kflag .eq. 3) .or. (init .eq. 0)) go to 45
      if (kflag .eq. 4) go to 40
      do k=1,neqn
         if ((kflag .eq. 5)  .and.  (abserr(k) .eq. 0.0d0)) go to 30
         if ((kflag .eq. 6)  .and.  (relerr(k) .le. savre(k))  .and.
     1        (abserr(k) .le. savae(k))) go to 30
      enddo
      go to 50
c
c     iflag = 3,4,5,6,7 or 8
   25 if (iflag .eq. 3) go to 45
      if (iflag .eq. 4) go to 40
      do k=1,neqn
         if ((iflag .ne. 5) .or. (abserr(k) .le. 0.0d0)) go to 30
      enddo
      goto 45
c
c     integration cannot be continued since user did not respond to
c     the instructions pertaining to iflag=5,6,7 or 8
   30 stop
c
c     reset function evaluation counter
   40 nfe=0
      if (mflag .eq. 2) go to 50
c
c     reset flag value from previous call
   45 iflag=jflag
      if (kflag .eq. 3) mflag=iabs(iflag)
c
c     save input iflag and set continuation flag value for subsequent
c     input checking
   50 jflag=iflag
      kflag=0
c
c     save relerr and abserr for checking input on subsequent calls
c      print *,'mid y=',y
      do k=1,neqn
         savre(k)=relerr(k)
         savae(k)=abserr(k)
      enddo
c
c     restrict relative error tolerance to be at least as large as
c     2*eps+remin to avoid limiting precision difficulties arising
c     from impossible accuracy requests
c
      rer=2.0d0*eps+remin
      do k=1,neqn
         if (relerr(k) .lt. rer) then
c     
c     relative error tolerance too small
            relerr(k)=rer
            iflag=3
            kflag=3
         endif
      enddo
      if(iflag.eq.3 .and. kflag.eq.3) return
c
   55 dt=tout-t
c
      if (mflag .eq. 1) go to 60
      if (init .eq. 0) go to 65
      go to 80
c
c     initialization --
c                       set initialization completion indicator,init
c                       set indicator for too many output points,kop
c                       evaluate initial derivatives
c                       set counter for function evaluations,nfe
c                       estimate starting stepsize
c
   60 init=0
      kop=0
c
      a=t
c      print *,'y=',y
c      stop
      call f(a,y,yp)
c      print *,'yp=',yp
      nfe=1
      if (t .ne. tout) go to 65
      iflag=2
      return
c
c
   65 init=1
      h=dabs(dt)
      toln=0.d0
      do 70 k=1,neqn
        tol=relerr(k)*dabs(y(k))+abserr(k)
        if (tol .le. 0.d0) go to 70
        toln=tol
        ypk=dabs(yp(k))
        if (ypk*h**5 .gt. tol) h=(tol/ypk)**0.2d0
   70 continue

c      print *, 'toln=',toln,tol


      if (toln .le. 0.0d0) h=0.0d0
      h=dmax1(h,u26*dmax1(dabs(t),dabs(dt)))
      jflag=isign(2,iflag)
c
c
c     set stepsize for integration in the direction from t to tout
c
   80 h=dsign(h,dt)
c
c     test to see if rkf45 is being severely impacted by too many
c     output points
c
      if (dabs(h) .ge. 2.0d0*dabs(dt)) kop=kop+1
      if (kop .ne. 100) go to 85
c
c     unnecessary frequency of output
      kop=0
      iflag=7
      return
c
   85 if (dabs(dt) .gt. u26*dabs(t)) go to 95
c
c     if too close to output point,extrapolate and return
c
      do 90 k=1,neqn
   90   y(k)=y(k)+dt*yp(k)
      a=tout

      call f(a,y,yp)

      nfe=nfe+1

      print *,'TOO CLOSE!!!!!!!!!!!!!!!!!'

      go to 300
c
c
c     initialize output point indicator
c
   95 output= .false.
c
c     to avoid premature underflow in the error tolerance function,
c     scale the error tolerances
c
      do k=1,neqn
         scale(k)=2.0d0/relerr(k)
         ae(k)=scale(k)*abserr(k)
      enddo

c      print *,'scale',scale,ae

c
c
c     step by step integration
c
  100 hfaild= .false.
c
c     set smallest allowable stepsize
c
      hmin=u26*dabs(t)
c
c     adjust stepsize if necessary to hit the output point.
c     look ahead two steps to avoid drastic changes in the stepsize and
c     thus lessen the impact of output points on the code.
c
      dt=tout-t
      if (dabs(dt) .ge. 2.0d0*dabs(h)) go to 200
      if (dabs(dt) .gt. dabs(h)) go to 150
c
c     the next successful step will complete the integration to the
c     output point
c
      output= .true.
      h=dt
      go to 200
c
  150 h=0.5d0*dt
c
c
c
c     core integrator for taking a single step
c
c     the tolerances have been scaled to avoid premature underflow in
c     computing the error tolerance function et.
c     to avoid problems with zero crossings,relative error is measured
c     using the average of the magnitudes of the solution at the
c     beginning and end of a step.
c     the error estimate formula has been grouped to control loss of
c     significance.
c     to distinguish the various arguments, h is not permitted
c     to become smaller than 26 units of roundoff in t.
c     practical limits on the change in the stepsize are enforced to
c     smooth the stepsize selection process and to avoid excessive
c     chattering on problems having discontinuities.
c     to prevent unnecessary failures, the code uses 9/10 the stepsize
c     it estimates will succeed.
c     after a step failure, the stepsize is not allowed to increase for
c     the next attempted step. this makes the code more efficient on
c     problems having discontinuities and more effective in general
c     since local extrapolation is being used and extra caution seems
c     warranted.
c
c
c     test number of derivative function evaluations.
c     if okay,try to advance the integration from t to t+h
c
  200 if (nfe .le. maxnfe) go to 220
c
c     too much work
      iflag=4
      kflag=4
      return
c
c     advance an approximate solution over one step of length h
c
  220 call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
      nfe=nfe+5
c
c     compute and test allowable tolerances versus local error estimates
c     and remove scaling of tolerances. note that relative error is
c     measured with respect to the average of the magnitudes of the
c     solution at the beginning and end of the step.
c
      eeoet=0.0d0
      do 250 k=1,neqn
        et(k)=dabs(y(k))+dabs(f1(k))+ae(k)
        if (et(k) .gt. 0.0d0) go to 240
c
c       inappropriate error tolerance
        iflag=5
        return
c
  240   ee(k)=dabs((-2090.0d0*yp(k)+(21970.0d0*f3(k)-15048.0d0*f4(k)))+
     1                        (22528.0d0*f2(k)-27360.0d0*f5(k)))
c        print *,'loop:',k,et(k),ee(k),ee(k)/et(k)
  250   eeoet= dmax1(eeoet,ee(k)/et(k))
c
c        print *,'eeoet=',eeoet

      do k=1,neqn
         esttol(k)=dabs(h)*eeoet*scale(k)/752400.0d0
      enddo

      if (esttol(1) .le. 1.0d0 .and. esttol(2) .le. 1.0d0) go to 260
c
c
c     unsuccessful step
c                       reduce the stepsize , try again
c                       the decrease is limited to a factor of 1/10
c
      hfaild= .true.
      output= .false.

      esttolmax=esttol(1)
      do k=2,neqn
         esttolmax=dmax1(esttolmax,esttol(k))
      enddo

      s=0.1d0
      if (esttolmax .lt. 59049.0d0) s=0.9d0/esttolmax**0.2d0

c      write(6,*) '590... s=',s,esttol(1)
c      do k=1,1
C         write(6,*) eeoet,dabs(h),scale(k),esttol(k)
c      enddo

C      if(esttol(1).lt.30.d0) stop

      h=s*h
      if (dabs(h) .gt. hmin) go to 200
c
c     requested error unattainable at smallest allowable stepsize

      print *,'relerr=',relerr
      print *,'abserr=',abserr
      print *,'esttolmax=',esttolmax
      print *,'esttol=',esttol
      print *, hmin,h

C      stop

      iflag=6
      kflag=6
      return
c
c
c     successful step
c                        store solution at t+h
c                        and evaluate derivatives there
c
  260 t=t+h
      do 270 k=1,neqn
  270   y(k)=f1(k)
      a=t

      if(kmax.gt.0)then
            if(kount.lt.kmax-1)then
               kount=kount+1
               xp(kount)=a
               do 13 i=1,neqn
                  ysav(kount,i)=y(i)
 13            continue
               xsav=a
            else
               print *,'got to kmax...'
           endif
      endif

      call f(a,y,yp)
c      if(y(1).lt.1500)print *, t,y,yp

      nfe=nfe+1
c
c
c                       choose next stepsize
c                       the increase is limited to a factor of 5
c                       if step failure has just occurred, next
c                          stepsize is not allowed to increase
c
      s=5.0d0

      esttolmax=esttol(1)
      do k=2,neqn
         esttolmax=max(esttolmax,esttol(k))
      enddo

      if (esttolmax .gt. 1.889568d-4) s=0.9d0/esttolmax**0.2d0

C      print *, '1.88... s=',s,esttolmax,h

      if (hfaild) s=dmin1(s,1.0d0)
      h=dsign(dmax1(s*dabs(h),hmin),h)
c
c     end of core integrator
c
c
c     should we take another step
c
      if (output) then
C         print *,'NOT TOO CLOSE!!!!!!!!!!!!!!'
         go to 300
      endif
      if (iflag .gt. 0) go to 100
c
c
c     integration successfully completed
c
c     one-step mode
      iflag=-2
      return
c
c     interval mode
  300 t=tout
      iflag=2
      return
c
      end
      subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)
c
c     fehlberg fourth-fifth order runge-kutta method
c
c    fehl integrates a system of neqn first order
c    ordinary differential equations of the form
c             dy(i)/dt=f(t,y(1),---,y(neqn))
c    where the initial values y(i) and the initial derivatives
c    yp(i) are specified at the starting point t. fehl advances
c    the solution over the fixed step h and returns
c    the fifth order (sixth order accurate locally) solution
c    approximation at t+h in array s(i).
c    f1,---,f5 are arrays of dimension neqn which are needed
c    for internal storage.
c    the formulas have been grouped to control loss of significance.
c    fehl should be called with an h not smaller than 13 units of
c    roundoff in t so that the various independent arguments can be
c    distinguished.
c
c
      integer  neqn
      double precision  y(neqn),t,h,yp(neqn),f1(neqn),f2(neqn),
     1  f3(neqn),f4(neqn),f5(neqn),s(neqn)
c
      double precision  ch
      integer  k
c
      ch=h/4.0d0
      do 221 k=1,neqn
  221   f5(k)=y(k)+ch*yp(k)
      call f(t+ch,f5,f1)
c
      ch=3.0d0*h/32.0d0
      do 222 k=1,neqn
  222   f5(k)=y(k)+ch*(yp(k)+3.0d0*f1(k))
      call f(t+3.0d0*h/8.0d0,f5,f2)
c
      ch=h/2197.0d0
      do 223 k=1,neqn
  223   f5(k)=y(k)+ch*(1932.0d0*yp(k)+(7296.0d0*f2(k)-7200.0d0*f1(k)))
      call f(t+12.0d0*h/13.0d0,f5,f3)
c
      ch=h/4104.0d0
      do 224 k=1,neqn
  224   f5(k)=y(k)+ch*((8341.0d0*yp(k)-845.0d0*f3(k))+
     1                            (29440.0d0*f2(k)-32832.0d0*f1(k)))
      call f(t+h,f5,f4)
c
      ch=h/20520.0d0
      do 225 k=1,neqn
  225   f1(k)=y(k)+ch*((-6080.0d0*yp(k)+(9295.0d0*f3(k)-
     1         5643.0d0*f4(k)))+(41040.0d0*f1(k)-28352.0d0*f2(k)))
      call f(t+h/2.0d0,f1,f5)
c
c     compute approximate solution at t+h
c
      ch=h/7618050.0d0
      do 230 k=1,neqn
  230   s(k)=y(k)+ch*((902880.0d0*yp(k)+(3855735.0d0*f3(k)-
     1        1371249.0d0*f4(k)))+(3953664.0d0*f2(k)+
     2        277020.0d0*f5(k)))
c
      return
      end
      subroutine derivs(m,y,dydm)
      implicit none
      include 'params.h'
C     routine to set up differential equations for integrating structure

      double precision m, y(*), dydm(*)
      integer Jrmax,maxJr,Jr,NP
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),ELr(NROWS,elmax)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      double precision z2b(NROWS),z2c(NROWS),z2d(NROWS)
      common/splinearrays/z2b,z2c,z2d
      double precision Aofm,seval

      maxJr=Jrmax
      Aofm=seval(maxJr,m,Mr,Ar,z2b,z2c,z2d)


      dydm(1)=-G*m/(4.d0*pi*y(2)**4)
      dydm(2)=( Aofm/dabs(y(1)) )**0.6d0/(4.d0*pi*y(2)**2)

      return
      end
