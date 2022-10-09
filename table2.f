c       Copyright (c) 2001, 2002, 2003
c       by James Lombardi, Vassar College, Poughkeepsie, NY
c       and Jessica Sawyer Warren, Rutgers University, NJ.
c    This file is distributed with Make Me A Star.
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
* This file contains a sample driver for the makemeastar subroutine
* version 1.6
* July 2003
* James Lombardi (lombardi@vassar.edu) and Jessica Sawyer Warren
*************************************************************************
      program meastar
c Example driver program for the makemeastar subroutine
      implicit none
      include 'params.h'
      integer Jr,C,verbosity,NP
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N), ELp(NROWS,elmax,N)
      integer numShells(N)
      double precision Mproduct(NROWS),Pproduct(NROWS),
     &     Dproduct(NROWS),Rproduct(NROWS),ELproduct(elmax,NROWS),
     &     jproduct(NROWS)
      double precision periastron, R86prod, R50prod, R95prod, R99prod
      double precision R90prod
      double precision periastron2
      integer ii50, ii86, ii90, ii95, ii99
      CHARACTER*64 parentfile(3),outfile,mnmfile,mnmnmfile
      integer imass1,imass2,imass3,icase
      common /productinfo/ Mproduct, Rproduct
      integer elarray(2)

      el = 20
      print *,'Number of chemical species to follow=',el
      verbosity = 1
      print *,'verbosity=',verbosity

c      print *,'Enter case number'

      do icase=1,20
      
         if (icase.eq.1 .or. icase.eq.3 .or. icase.eq.4) then
            imass1=8
            imass2=6
            imass3=6
         else
            imass1=6
            imass2=6
            imass3=8
         endif
         
         print *
         write(6,'(a)')'***********************************************'
         print *,'case and imasses=',icase,imass1,imass2,imass3
         
 202     format('c',i1,'n',i1,'.all')
         write(mnmfile, 202) imass1,imass2
 203     format('c',i1,'n',i1,'n',i1,'.all')
         write(mnmnmfile, 203) imass1,imass2,imass3
      
c     Files to store radii at various enclosed mass fractions
         OPEN(91, FILE=mnmfile)
         OPEN(93, FILE=mnmnmfile)
         
         NP = 0
c         print *,'NP=',NP
         
c         write(6,*) 'Enter rp1/(R1+R2):'
         if(icase.eq.1 .or. icase.eq.2 .or. icase.eq.12) then
            periastron=0.d0
         else if (icase.eq.11 .or. icase.eq.13) then
            periastron=0.125d0
         else
            periastron=0.25d0
         endif
            

 302     format('m0',i1,'c.mdl')
         write(parentfile(1), 302) imass1
         write(parentfile(2), 302) imass2
         print *,'parentfile(1)=',parentfile(1)
         print *,'parentfile(2)=',parentfile(2)

         print *,'r_p1/(R_1+R+2)=',periastron
 102     format('c',i1,'n',i1,f5.3,'.out')
         write(outfile, 102) imass1,imass2,periastron
         print *,'outfile=',outfile
         
         elarray(1)=14
         elarray(2)=14
         call readparentfiles(2,parentfile,verbosity,elarray,
     $        Mp,Rp,Pp,Dp,ELp,numShells)
         
         if(verbosity.ge.2)
     $        print *, 'About to call the makemeastar subroutine...'
         call makemeastar(periastron,
     $        Mp(1,1),Rp(1,1),Pp(1,1),Dp(1,1),ELp(1,1,1),
     $        Mp(1,2),Rp(1,2),Pp(1,2),Dp(1,2),ELp(1,1,2),
     $        numShells(1),numShells(2),
     $        verbosity,
     &        el,NP,Mproduct,Pproduct,Dproduct,Rproduct,
     &        ELproduct,jproduct)
         if(verbosity.ge.2) then
            write(6,*) 'Back from the makemeastar subroutine...'
            write(6,*)
            write(6,*) 'NP has been set to',NP
            write(6,*)
         endif
         
         if(verbosity.ge.1) then
            write(6,'(1x,a,f14.3,a)')
     $           'Collision product radius (before thermal relaxation)='
     $           , Rproduct(NP)/Rsun,' R_sun'
            write(6,*)'... WRITING COLLISION PRODUCT MODEL TO FILE ',
     $           outfile
         endif
         OPEN(50,FILE=outfile)

         do Jr=1,NP
            WRITE(50,54) Mproduct(Jr),Rproduct(Jr),Pproduct(Jr),
     $           Dproduct(Jr),
     &           (real(ELproduct(C,Jr)),C=1,el)            
         enddo
         CLOSE(50)
 54      FORMAT(E17.11,3(1X,E14.8),2F8.5,99(1X,E11.5))

c     Find and write radii at various enclosed mass fractions
         call getindices(NP, ii50, ii86, ii90, ii95, ii99)

         R50prod=
     $        ( (Mproduct(ii50)/Mproduct(NP)-.5d0)*Rproduct(ii50-1)+
     $          (.5d0-Mproduct(ii50-1)/Mproduct(NP))*Rproduct(ii50))/
     $          (Mproduct(ii50)-Mproduct(ii50-1))*Mproduct(NP)
         R86prod=
     $        ( (Mproduct(ii86)/Mproduct(NP)-.86d0)*Rproduct(ii86-1)+
     $        (.86d0-Mproduct(ii86-1)/Mproduct(NP))*Rproduct(ii86))/
     $          (Mproduct(ii86)-Mproduct(ii86-1))*Mproduct(NP)
         R90prod=
     $        ( (Mproduct(ii90)/Mproduct(NP)-.90d0)*Rproduct(ii90-1)+
     $        (.90d0-Mproduct(ii90-1)/Mproduct(NP))*Rproduct(ii90))/
     $        (Mproduct(ii90)-Mproduct(ii90-1))*Mproduct(NP)
         R95prod=
     $        ( (Mproduct(ii95)/Mproduct(NP)-.95d0)*Rproduct(ii95-1)+
     $        (.95d0-Mproduct(ii95-1)/Mproduct(NP))*Rproduct(ii95))/
     $        (Mproduct(ii95)-Mproduct(ii95-1))*Mproduct(NP)
         R99prod=
     $        ( (Mproduct(ii99)/Mproduct(NP)-.99d0)*Rproduct(ii99-1)+
     $          (.99d0-Mproduct(ii99-1)/Mproduct(NP))*Rproduct(ii99))/
     $          (Mproduct(ii99)-Mproduct(ii99-1))*Mproduct(NP)

 544     format(6E15.8)
         write(91, 544) periastron,R50prod/Rsun,R86prod/Rsun,
     $        R95prod/Rsun,R99prod/Rsun,Rproduct(NP)/Rsun

         if(verbosity.ge.2)
     $        write(6,*) icase,
     $        ': in first product, 90 and 95% radii are',
     $        R90prod/Rsun,R95prod/Rsun

         parentfile(1)=outfile
         write(parentfile(2), 302) imass3
         print *,'third parent=',parentfile(2)

            NP = 0
C            print *,'NP=',NP

            elarray(1)=20
            elarray(2)=14
            call readparentfiles(2,parentfile,verbosity,elarray,
     $           Mp,Rp,Pp,Dp,ELp,numShells)

            write(6,*) 'Enter rp2 IN SOLAR RADII:'
C            read(5,*) periastron2
            if(icase.eq.1 .or. icase.eq.2 .or. icase.eq.3 .or.
     $           icase.eq.7 .or. icase.eq.8 .or. icase.eq.9 .or.
     $           icase.eq.10 .or. icase.eq.11) then
               periastron2=0.d0
            else if (icase.eq.4) then
               periastron2=0.0955d0
            else if (icase.eq.5 .or. icase.eq.6) then
               periastron2=0.198d0
            else if (icase.ge.12 .and. icase.le.19) then
               periastron2=0.505d0
            else if (icase.eq.20) then
               periastron2=0.758d0
            endif

            periastron2=periastron2 *6.96d10
     $           /( Rp(numShells(1),1) +Rp(numShells(2),2) )
            print *,'r_p2/(R_1+R+2)=',periastron2


 103        format('c',i1,'n',i1,'n',i1,'_',f5.3,'_',f5.3,'.out')
            write(outfile, 103) imass1,imass2,imass3,periastron,
     $           periastron2
            print *,'outfile=',outfile

            call makemeastar(periastron2,
     $           Mp(1,1),Rp(1,1),Pp(1,1),Dp(1,1),ELp(1,1,1),
     $           Mp(1,2),Rp(1,2),Pp(1,2),Dp(1,2),ELp(1,1,2),
     $           numShells(1),numShells(2),
     $           verbosity,
     &           el,NP,Mproduct,Pproduct,Dproduct,Rproduct,
     &           ELproduct,jproduct)

            OPEN(50,FILE=outfile)
            do Jr=1,NP
               WRITE(50,54) Mproduct(Jr),Rproduct(Jr),Pproduct(Jr),
     $              Dproduct(Jr),
     &              (real(ELproduct(C,Jr)),C=1,el)            
            enddo
            CLOSE(50)

            call getindices(NP, ii50, ii86, ii90, ii95, ii99)
            
            R50prod=
     $           ( (Mproduct(ii50)/Mproduct(NP)-.5d0)*Rproduct(ii50-1)+
     $           (.5d0-Mproduct(ii50-1)/Mproduct(NP))*Rproduct(ii50))/
     $           (Mproduct(ii50)-Mproduct(ii50-1))*Mproduct(NP)
            R86prod=
     $           ( (Mproduct(ii86)/Mproduct(NP)-.86d0)*Rproduct(ii86-1)+
     $           (.86d0-Mproduct(ii86-1)/Mproduct(NP))*Rproduct(ii86))/
     $           (Mproduct(ii86)-Mproduct(ii86-1))*Mproduct(NP)
            R90prod=
     $           ( (Mproduct(ii90)/Mproduct(NP)-.90d0)*Rproduct(ii90-1)+
     $           (.90d0-Mproduct(ii90-1)/Mproduct(NP))*Rproduct(ii90))/
     $           (Mproduct(ii90)-Mproduct(ii90-1))*Mproduct(NP)
            R95prod=
     $           ( (Mproduct(ii95)/Mproduct(NP)-.95d0)*Rproduct(ii95-1)+
     $           (.95d0-Mproduct(ii95-1)/Mproduct(NP))*Rproduct(ii95))/
     $           (Mproduct(ii95)-Mproduct(ii95-1))*Mproduct(NP)
            R99prod=
     $           ( (Mproduct(ii99)/Mproduct(NP)-.99d0)*Rproduct(ii99-1)+
     $           (.99d0-Mproduct(ii99-1)/Mproduct(NP))*Rproduct(ii99))/
     $           (Mproduct(ii99)-Mproduct(ii99-1))*Mproduct(NP)

            if(verbosity.ge.2)
     $           write(6,*) icase,
     $           ': in second product, 90 and 95% radii are',
     $           R90prod/Rsun,R95prod/Rsun
            
 545     format(i3,1x,2G8.2,5E15.8)
            write(93, 545) icase,
     $           periastron,periastron2,
     $           R50prod/Rsun,R86prod/Rsun,
     $           R95prod/Rsun,R99prod/Rsun,Rproduct(NP)/Rsun
            write(6, 545) icase,
     $           periastron,periastron2,
     $           R50prod/Rsun,R86prod/Rsun,
     $           R95prod/Rsun,R99prod/Rsun,Rproduct(NP)/Rsun
            
      close(91)
      close(93)
      
      enddo

      end
**********************************************************************
      SUBROUTINE readparentfiles(N,parentfile,verbosity,el,
     $     Mp,Rp,Pp,Dp,ELp,numShells)
      implicit none
      integer K,N,verbosity,NROWS,JJ,I,el(N),numShells(N),C,elmax
      CHARACTER*64 parentfile(N)
      parameter(NROWS=4000,elmax=20)
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N), ELp(NROWS,elmax,N)
**********************************************************************
* Input:
*     parentfile = ARRAY of character strings for input file names, so ...
*        parentfile(1)=Name of data file for more massive parent star (star 1)
*        parentfile(2)=Name of data file for less massive parent star (star 2)
**********************************************************************

* read in the input files for each parent star K:
      DO K=1,N
         if(verbosity.ge.3) then
            WRITE(6,*)
            WRITE(6,*) 'Input: ',parentfile(K)
         endif

         OPEN(26,FILE=parentfile(K))         

* I refers to rows in the input data file.  JJ refers to rows where the
* enclosed mass M is increasing.   We will throw away any rows where this does
* not occur.  (Some input files may have truncated the value for the mass at
* a small enough number of digits that two rows appear to have the same
* enclosed mass M.)
         JJ=0
         DO I=1,NROWS
            READ(26,*,END=30) Mp(I,K),Rp(I,K),Pp(I,K),
     &           Dp(I,K),(ELp(I,C,K), C=1,el(K))
            do C=el(K)+1,elmax
               ELp(I,C,K)=0.d0
            enddo
            if(K.eq.2 .and. el(1).eq.elmax) then
               Elp(I,el(K)+3,K)=1.d0
               Elp(I,el(K)+6,K)=Elp(I,5,K)
            else if(el(K)+K.le.elmax) then
               ELp(I,el(K)+K,K)=1.d0
               Elp(I,el(K)+K+3,K)=Elp(I,5,K)
            endif
C            print *,K, (Elp(I,C,K),C=5,elmax)
            IF(Mp(I,K).GT.Mp(JJ,K) .OR. JJ.EQ.0) THEN
               JJ=JJ+1
               Mp(JJ,K)=Mp(I,K)
               Pp(JJ,K)=Pp(I,K)
               Dp(JJ,K)=Dp(I,K)
               DO C=1,el(K)
                  ELp(JJ,C,K)=ELp(I,C,K)
               ENDDO
               Rp(JJ,K)=Rp(I,K)
            ENDIF
         ENDDO
            
 30      numShells(K)=JJ
         CLOSE(26)
      ENDDO

      return
      end


      subroutine getindices(NP,I50,I86,I90,I95,I99)

******************************************************************
* This routine finds the locations in the parent inside of which *
* 50%, 86% and 95% of the total mass is enclosed                 *
******************************************************************
      
      implicit none
      INCLUDE 'params.h'
      INTEGER I,I50,I86,I90,I95,I99,NP
      DOUBLE PRECISION Mproduct(NROWS),Rproduct(NROWS)
c     &     ELproduct(NROWS,elmax),
c     &     Dproduct(NROWS),Aproduct(NROWS),
c     &     As(NROWS, N),ELm(NROWS,elmax, N),
c     &     Ap(NROWS,N), Dp(NROWS,N)
      COMMON/productinfo/Mproduct,Rproduct
c      COMMON/NEWENT/Jpmax,Ap,Dp,As

      DO I=1,NP
         IF(Mproduct(I).GE.50.d-2*Mproduct(NP)) GOTO 47
      ENDDO
 47   I50=I
      DO I=I50,NP
         IF(Mproduct(I).GE.86.d-2*Mproduct(NP)) GOTO 48
      ENDDO
 48   I86=I
      DO I=I86,NP
         IF(Mproduct(I).GE.90.d-2*Mproduct(NP)) GOTO 49
      ENDDO
 49   I90=I
      DO I=I90,NP
         IF(Mproduct(I).GE.95.d-2*Mproduct(NP)) GOTO 50
      ENDDO
 50   I95=I
      DO I=I95,NP
         IF(Mproduct(I).GE.99.d-2*Mproduct(NP)) GOTO 51
      ENDDO
 51   I99=I

      END
