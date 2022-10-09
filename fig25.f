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
      double precision periastron2
      integer ii50, ii86, ii95, ii99, iperiastron, iperiastron2
      CHARACTER*64 parentfile(3),outfile,mnmfile,mnmnmfile
      integer ipmax1,ipmax2
      double precision perimax
      integer imass1,imass2,imass3
      common /productinfo/ Mproduct, Rproduct

      el = 0
      print *,'Number of chemical species to follow=',el
      verbosity = 0
      print *,'verbosity=',verbosity
      ipmax1=11
      ipmax2=44
c      ipmax1=11
c      ipmax2=11
      perimax=1.1d0

      
      do imass1=8,4,-2
      do imass2=4,imass1,2
      do imass3=4,8,2
c      do imass1=6,6,-2
c      do imass2=6,6,2
c      do imass3=8,8,2

 202     format('b',i1,'n',i1,'.all')
         write(mnmfile, 202) imass1,imass2
 203     format('b',i1,'n',i1,'n',i1,'.all')
         write(mnmnmfile, 203) imass1,imass2,imass3
      

c     Files to store r50, r86, and r
      OPEN(21, FILE=mnmfile)
      OPEN(23, FILE=mnmnmfile)

C      do iperiastron=0,0
c      do iperiastron=0,ipmax1
      do iperiastron=0,ipmax1,2
c      do iperiastron=4,4
c      do iperiastron=ipmax1-1,ipmax1-1
         NP = 0
         periastron=perimax*iperiastron/ipmax1

 302     format('m0',i1,'c.mdl')
 303     format('m0',i1,'.mdl')
         if(imass1.ne.4) then
            write(parentfile(1), 302) imass1
         else
            write(parentfile(1), 303) imass1
         endif
         if(imass2.ne.4) then
            write(parentfile(2), 302) imass2
         else
            write(parentfile(2), 303) imass2
         endif
            
         print *,'parentfile(1)=',parentfile(1)
         print *,'parentfile(2)=',parentfile(2)

         print *,'periastron/(R_1+R_2)=',periastron
 102     format('b',i1,'n',i1,f5.3,'.out')
         write(outfile, 102) imass1,imass2,periastron
         print *,'outfile=',outfile
         
         call readparentfiles(2,parentfile,verbosity,el,
     $        Mp,Rp,Pp,Dp,ELp,numShells)
         
         if(verbosity.ge.1)
     $        print *, 'About to call the makemeastar subroutine...'
         call makemeastar(periastron,
     $        Mp(1,1),Rp(1,1),Pp(1,1),Dp(1,1),ELp(1,1,1),
     $        Mp(1,2),Rp(1,2),Pp(1,2),Dp(1,2),ELp(1,1,2),
     $        numShells(1),numShells(2),
     $        verbosity,
     &        el,NP,Mproduct,Pproduct,Dproduct,Rproduct,
     &        ELproduct,jproduct)
         if(verbosity.ge.1) then
            write(6,*) 'Back from the makemeastar subroutine...'
            write(6,*)
            write(6,*) 'NP has been set to',NP
            write(6,*)
         endif
         
         if(verbosity.ge.1) then
            write(6,*)'Collision product radius (before thermal
     $           relaxation)
     $           =', Rproduct(NP)/Rsun,' R_sun'
            write(6,*)'... WRITING COLLISION PRODUCT MODEL TO FILE'
         endif
         OPEN(50,FILE=outfile)
         do Jr=1,NP
            WRITE(50,54) Mproduct(Jr),Rproduct(Jr),Pproduct(Jr),
     $           Dproduct(Jr),
     &           (real(ELproduct(C,Jr)),C=1,el)            
         enddo
         CLOSE(50)
 54      FORMAT(E18.12,3(1X,E14.8),2F8.5,99(1X,E11.5))

c     Write rr, rr50, rr86, and rr95 to file
         call getindices(NP, ii50, ii86, ii95, ii99)

         R50prod=
     $        ( (Mproduct(ii50)/Mproduct(NP)-.5d0)*Rproduct(ii50-1)+
     $          (.5d0-Mproduct(ii50-1)/Mproduct(NP))*Rproduct(ii50))/
     $          (Mproduct(ii50)-Mproduct(ii50-1))*Mproduct(NP)
         R86prod=
     $        ( (Mproduct(ii86)/Mproduct(NP)-.86d0)*Rproduct(ii86-1)+
     $        (.86d0-Mproduct(ii86-1)/Mproduct(NP))*Rproduct(ii86))/
     $          (Mproduct(ii86)-Mproduct(ii86-1))*Mproduct(NP)
         R95prod=
     $        ( (Mproduct(ii95)/Mproduct(NP)-.95d0)*Rproduct(ii95-1)+
     $        (.95d0-Mproduct(ii95-1)/Mproduct(NP))*Rproduct(ii95))/
     $        (Mproduct(ii95)-Mproduct(ii95-1))*Mproduct(NP)
         R99prod=
     $        ( (Mproduct(ii99)/Mproduct(NP)-.99d0)*Rproduct(ii99-1)+
     $          (.99d0-Mproduct(ii99-1)/Mproduct(NP))*Rproduct(ii99))/
     $          (Mproduct(ii99)-Mproduct(ii99-1))*Mproduct(NP)

 544     format(6E15.8)
         write(21, 544) periastron,R50prod/Rsun,R86prod/Rsun,
     $        R95prod/Rsun,R99prod/Rsun,Rproduct(NP)/Rsun

         parentfile(1)=outfile
         if(imass3.ne.4) then
            write(parentfile(2), 302) imass3
         else
            write(parentfile(2), 303) imass3
         endif
         print *,'third parent=',parentfile(2)

C         do iperiastron2=6,8
c         do iperiastron2=ipmax+1,ipmax
c         do iperiastron2=ipmax2,0,-1
c         do iperiastron2=0,ipmax2/2
         do iperiastron2=0,ipmax2
c         do iperiastron2=15,15
C         do iperiastron2=0,0

            NP = 0
            periastron2=perimax*iperiastron2/ipmax2

            call readparentfiles(2,parentfile,verbosity,el,
     $           Mp,Rp,Pp,Dp,ELp,numShells)

            call makemeastar(periastron2,
     $           Mp(1,1),Rp(1,1),Pp(1,1),Dp(1,1),ELp(1,1,1),
     $           Mp(1,2),Rp(1,2),Pp(1,2),Dp(1,2),ELp(1,1,2),
     $           numShells(1),numShells(2),
     $           verbosity,
     &           el,NP,Mproduct,Pproduct,Dproduct,Rproduct,
     &           ELproduct,jproduct)

            call getindices(NP, ii50, ii86, ii95, ii99)
            
            R50prod=
     $           ( (Mproduct(ii50)/Mproduct(NP)-.5d0)*Rproduct(ii50-1)+
     $           (.5d0-Mproduct(ii50-1)/Mproduct(NP))*Rproduct(ii50))/
     $           (Mproduct(ii50)-Mproduct(ii50-1))*Mproduct(NP)
            R86prod=
     $           ( (Mproduct(ii86)/Mproduct(NP)-.86d0)*Rproduct(ii86-1)+
     $           (.86d0-Mproduct(ii86-1)/Mproduct(NP))*Rproduct(ii86))/
     $           (Mproduct(ii86)-Mproduct(ii86-1))*Mproduct(NP)
            R95prod=
     $           ( (Mproduct(ii95)/Mproduct(NP)-.95d0)*Rproduct(ii95-1)+
     $           (.95d0-Mproduct(ii95-1)/Mproduct(NP))*Rproduct(ii95))/
     $           (Mproduct(ii95)-Mproduct(ii95-1))*Mproduct(NP)
            R99prod=
     $           ( (Mproduct(ii99)/Mproduct(NP)-.99d0)*Rproduct(ii99-1)+
     $           (.99d0-Mproduct(ii99-1)/Mproduct(NP))*Rproduct(ii99))/
     $           (Mproduct(ii99)-Mproduct(ii99-1))*Mproduct(NP)
            
 545     format(2i3,1x,2G10.4,5E15.8)
            write(23, 545) iperiastron,iperiastron2,
     $           periastron,periastron2,
     $           R50prod/Rsun,R86prod/Rsun,
     $           R95prod/Rsun,R99prod/Rsun,Rproduct(NP)/Rsun
            write(6, 545) iperiastron,iperiastron2,
     $           periastron,periastron2,
     $           R50prod/Rsun,R86prod/Rsun,
     $           R95prod/Rsun,R99prod/Rsun,Rproduct(NP)/Rsun
            
c     Closing the inner periastron loop 
         enddo

      
c     Closing the outer periastron loop 
      enddo
      close(21)
      close(23)

      enddo
      enddo
      enddo

      end
**********************************************************************
      SUBROUTINE readparentfiles(N,parentfile,verbosity,el,
     $     Mp,Rp,Pp,Dp,ELp,numShells)
      implicit none
      integer K,N,verbosity,NROWS,JJ,I,el,numShells(N),C,elmax
      CHARACTER*64 parentfile(N)
      parameter(NROWS=4000,elmax=14)
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
         if(verbosity.ge.2) then
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
     &           Dp(I,K),(ELp(I,C,K), C=1,el)
            IF(Mp(I,K).GT.Mp(JJ,K) .OR. JJ.EQ.0) THEN
               JJ=JJ+1
               Mp(JJ,K)=Mp(I,K)
               Pp(JJ,K)=Pp(I,K)
               Dp(JJ,K)=Dp(I,K)
               DO C=1,el
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


      subroutine getindices(NP,I50,I86,I95,I99)

******************************************************************
* This routine finds the locations in the parent inside of which *
* 50%, 86%, 95% and 99% of the total mass is enclosed                 *
******************************************************************
      
      implicit none
      INCLUDE 'params.h'
      INTEGER I,K,Jpmax(N),Pmax(N),I50,I86,I95,I99,NP
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
      DO I=I50+1,NP
         IF(Mproduct(I).GE.86.d-2*Mproduct(NP)) GOTO 48
      ENDDO
 48   I86=I
      DO I=I86+1,NP
         IF(Mproduct(I).GE.95.d-2*Mproduct(NP)) GOTO 49
      ENDDO
 49   I95=I
      DO I=I95+1,NP
         IF(Mproduct(I).GE.99.d-2*Mproduct(NP)) GOTO 50
      ENDDO
 50   I99=I

      END
