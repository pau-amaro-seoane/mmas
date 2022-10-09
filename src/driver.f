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
      double precision periastron
      CHARACTER*512 parentfile(2),outfile
      
      print *,'SAMPLE DRIVER PROGRAM FOR MAKE ME A STAR!'

***********************************************************************
* Get information from user about what stars she would like to collide:
***********************************************************************
      print *
      print *,'PARENT STARS:'
      print *,'First we need to choose which two stars to collide.'
      print *,'Included with the MMAS software package are the'
      print *,'sample data files m04.mdl, m06c.mdl and m08c.mdl that'
      print *,'model 0.4 solar mass, 0.6 solar mass and 0.8 solar'
      print *,'mass population II main sequence stars, respectively.'
      print *
      print *, 'Enter data file name for the more massive parent star.'
      print *, '(Some systems, including SGIs, require that you type'
      print *, 'single quotes around the filename.)'
      read *, parentfile(1)
      print *, 'Enter data file name for the second parent star.'
      read *, parentfile(2)
c      print *,'parentfiles=',parentfile

      print *
      print *, 'Enter the periastron separation, normalized to the sum'
      print *, 'of the parent star radii (for a headon collision, use'
      print *, 'zero here):'
      read *, periastron
      print *,'periastron/(R_1+R_2)=',periastron
      if(periastron.gt.2.2) then
         print *, 'WARNING: If the stars were to remain perfectly'
         print *, 'spherical as the approached each other, they would'
         print *, 'impact only for normalized periastron separations'
         print *, 'less than 2.  Since the stars do deform, slightly'
         print *, 'larger periastron separations can still lead to'
         print *, 'contact.  Nevertheless, the periastron separation'
         print *, 'requested seems to be rather large... please make'
         print *, 'sure that it is normalized to R_1 + R_2 (or modify'
         print *, 'the code to do that automatically).'
         stop
      endif

      print *,'OUTPUT:'
      print *,'What should the output file be named?'
      read *, outfile
      print *,'outfile=',outfile
      print *
      print *,'Information about the thermodynamic and chemical'
      print *,'composition profiles of the product will be returned'
      print *,'as arrays by the subroutine makemeastar.  Having longer'
      print *,'arrays will give a more accurate model of the product,'
      print *,'but will take longer to calculate.'
      print *
      print *,'How long would you like these arrays to be?'
      print *,'(If you do not know what to enter, something like 600'
      print *,'might be good to try.  Alternatively, if you enter'
      print *,'0 (zero) here, the code will automatically'
      print *,'choose a reasonable number of shells for you.)'
      read *, NP
      print *,'NP=',NP

      print *,'CHEMICAL COMPOSITION ABUNDANCES:'
      print *,'If the data files you are using for the parent stars'
      print *,'contain information about the chemical composition'
      print *,'profiles, then makemeastar can determine how these'
      print *,'elements will be distributed throughout the product.'
      print *,'If you choose not to follow the distribution of'
      print *,'chemical elements, the code will run significantly'
      print *,'faster.  If you choose to treat only some of the'
      print *,'available chemical abundance profiles from your parent'
      print *,'models, then it will be the last columns in the data'
      print *,'file that are ignored.'
      print *
      print *,'So, how many different chemical species do you want to'
      print *,'follow?'
      read *, el
      print *,'Number of chemical species to follow=',el

      print *,'VERBOSITY:'
      print *, 'How verbose would you like the output to be?'
      print *, ' (Enter 0 for no output to screen,'
      print *, '        1 for very limited output,'
      print *, '        2 for detailed output)'
      read *, verbosity
      print *,'verbosity=',verbosity

      call readparentfiles(2,parentfile,verbosity,el,
     $     Mp,Rp,Pp,Dp,ELp,numShells)
      
***********************************************************************
* Now that we know what collision to consider, make the call to the
* main routine:
***********************************************************************
      if(verbosity.ge.1)
     $     print *, 'About to call the makemeastar subroutine...'
      call makemeastar(periastron,
     $     Mp(1,1),Rp(1,1),Pp(1,1),Dp(1,1),ELp(1,1,1),
     $     Mp(1,2),Rp(1,2),Pp(1,2),Dp(1,2),ELp(1,1,2),
     $     numShells(1),numShells(2),
     $     verbosity,
     &     el,NP,Mproduct,Pproduct,Dproduct,Rproduct,
     &     ELproduct,jproduct)
      if(verbosity.ge.2) then
         write(6,*) 'Back from the makemeastar subroutine...'
         write(6,*)
         write(6,*) 'NP has been set to',NP
         write(6,*)
      endif
       

***********************************************************************
* Create a data file for the collision product model:
***********************************************************************
      if(verbosity.ge.1) then
         write(6,*)'Collision product radius (before thermal relaxation)
     $=', Rproduct(NP)/6.9598d10,' R_sun'
         write(6,*)'... WRITING COLLISION PRODUCT MODEL TO FILE ',
     $        outfile
      endif
      OPEN(50,FILE=outfile)
C      write(50,*)'m [M_sun], P [dyne/cm^2], rho [g/cm^3], r [R_sun], '
C     $     ,'X, Z, j [cm^2/s]' 
      do Jr=1,NP

c The following if statement could be used to just print some of the lines...
c         if(Mproduct(Jr).gt.0.95d0*Mproduct(NP).or.mod(Jr,5).eq.1)

c In version 1.4 and earlier, the following was used for the output file:
c         WRITE(50,54) Mproduct(Jr)/1.989d33,Pproduct(Jr),
c     $        Dproduct(Jr),Rproduct(Jr)/6.9598d10,
c     &        (real(ELproduct(C,Jr)),C=1,el),jproduct(Jr)
C 54   FORMAT(E14.8,3(1X,E14.8),2F8.5,99(1X,E10.4))

         WRITE(50,54) Mproduct(Jr),Rproduct(Jr),Pproduct(Jr),
     $        Dproduct(Jr),(real(ELproduct(C,Jr)),C=1,el),jproduct(Jr)     
      enddo
      CLOSE(50)
 54   FORMAT(E17.11,3(1X,E14.8),2F8.5,99(1X,E11.5))

      end
**********************************************************************
      SUBROUTINE readparentfiles(Np,parentfile,verbosity,el,
     $     Mp,Rp,Pp,Dp,ELp,numShells)
      implicit none
      include 'params.h'
      integer K,Np,verbosity,JJ,I,numShells(N),C
      CHARACTER*512 parentfile(N)
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N), ELp(NROWS,elmax,N)
**********************************************************************
* Input:
*     parentfile = ARRAY of character strings for input file names, so ...
*        parentfile(1)=Name of data file for more massive parent star (star 1)
*        parentfile(2)=Name of data file for less massive parent star (star 2)
**********************************************************************

      if(N.ne.Np) then
         write(6,*) 'Expecting ',N,' parents, got',Np
         stop
      endif

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
