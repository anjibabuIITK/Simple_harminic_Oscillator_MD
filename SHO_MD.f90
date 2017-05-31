! A FORTRAN90 program to solve Simple harmonic oscillator problem by using velocity verlet algarithm.
!    
! Authour: ANJI BABU KAPAKAYALA
!          IIT KANPUR, INDIA.
!
   PROGRAM SHO_MD
      IMPLICIT NONE
      REAL*8:: x, v, m, k,dt,a_t,te,a,eke, e, f, v2, a2
!      REAL*8:: force, vel_ver_posi, vel_ver_vel
      INTEGER::i,n  !n=no of steps
      x=1.0 ;v=1.0 ;m=1.0 ;k=1.0                 !mass,force constant ,ntial position and velositiy are taken as 1
      dt=0.01                                    ! Timestep
      n=30000                                    
      OPEN(1,file="sho_md.out",STATUS="NEW")
      OPEN(10,file="vmd_trajectory.xyz",STATUS="NEW")
      CALL force(x,e,f,k)
      a=f/m                                      !acceleration =force/mass
      
      DO i=1,n
        te=0.0                                    !Intilization of energies
         e=0.0
        eke=0.0
        CALL vel_ver_posi(x,dt,v,a)               
        a_t=a
        CALL force(x,e,f,k)
        a=f/m
        CALL vel_ver_vel(a,dt,v,a_t)
!energies      
        eke=0.5*m*(v**2)
        te=eke+e
        WRITE(1,*)i,x,eke,e,te
!=============writing vmd trajectory===============================!
       WRITE(10,*)1  !no of atoms
       WRITE(10,*)
       WRITE(10,*)"A", x,0,0
!===================================================================!
      END DO
      CLOSE(1)
      END PROGRAM SHO_MD
!--------------------FORCE SUBROUTINE-------------------------------!
      SUBROUTINE force(x,e,f,k)
         IMPLICIT NONE
         REAL*8::x, e, f, k, m, a
         e=0.5*k*(x**2)
         f=-k*x
      RETURN
      END SUBROUTINE force
!======================Positions ====================================!
      SUBROUTINE vel_ver_posi(x,dt,v,a) 
        IMPLICIT NONE
        REAL*8:: x, dt, v, a
        x=x+v*dt+0.5*(dt**2)*a
      RETURN
      END SUBROUTINE vel_ver_posi
!=======================================Velocities===================!
      SUBROUTINE vel_ver_vel(a,dt,v,a_t)
       IMPLICIT NONE
       REAL*8::dt, v, a, a_t
       v=v+0.5*dt*(a+a_t)     
      RETURN
      END SUBROUTINE vel_ver_vel
!====================================================================!
