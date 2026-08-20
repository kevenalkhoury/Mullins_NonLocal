************************************************************************
!     Alkhoury, K (2026). A non-local constitutive model for the Mullins effect in
!     filled elastomers. International Journal of Solids and Structures, 340, 114249. 
!     https://doi.org/10.1016/j.ijsolstr.2026.114249
!
!     Keven Alkhoury, April 2026, Abaqus 2024
!     This implementation builds upon the following works:
!
!     Alkhoury, K., Ivko, R., Hossain, M., Nadimpalli, S., & Chester, S. A. (2024).
!     Experiments and modeling of the coupled viscoelasticity and Mullins effect in filled rubber materials.
!     Journal of the Mechanics and Physics of Solids, 188, 105650.
!     https://doi.org/10.1016/j.jmps.2024.105650
!
!     Wang, S., & Chester, S. A. (2018).
!     Experimental characterization and continuum modeling of inelasticity in filled rubber-like materials.
!     International Journal of Solids and Structures, 136-137, 125�136.
!     https://doi.org/10.1016/j.ijsolstr.2017.12.010
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!     statev(1) = LambdaBar ------ Effective stretch
!     statev(2) = nu_nonlocal ------ Nonlocal soft volume fraction
!     statev(3) = Achainmax ------ Maximum amplified stretch
!     statev(4) = numax ------ Maximum local soft volume fraction
!     --------------------------------------------------------------
!     Material Properties Vector
!     --------------------------------------------------------------
!     Mechanical (UMAT):
!     --------------------------------------------------------------
!     Gshear= props(1) ------ Shear moudulus
!     Kbulk = props(2) ------ Bulk modulus
!     lamL = props(3) ------ Locking stretch
!     A = props(4) ------ ``A'' parameter (Mullins model by Qi and Boyce)
!     nu_s0 = props(5) ------ Initial soft volume fraction
!     nu_ss = props(6) ------ Saturation soft volume fraction
!     dE = props(7) ------ Peturbation parameter for numerical jacobian ddsdde
!--------------------------------------------------------------
!     Thermal (UMATHT):
!     --------------------------------------------------------------
!     lc = props(1) ------ Intrinsic characteristic length scale
!     specht = props(2) ------ Specific heat
!     tau = props(3) ------ Viscous regularization parameter
!**********************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      !
      include 'aba_param.inc'
      !
      dimension stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
      !
      character*80 cmname,file1
      !
      integer i,j,ii,jj,iterError
      !
      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),T_tau(3,3),effStr_tau,detF
      real*8 U_t(3,3),U_tau(3,3),U_inv(3,3),R_tau(3,3),stress_power
      real*8 theta_tau,theta_t,detF_tau,numax
      real*8 LambdaBar_t,LambdaBar_tau,Damage_t,Damage_tau
      real*8 nu_t,nu_tau,Achainmax,nu_s0,Histmax_t,Histmax_tau
      real*8 flag,varphi_t,varphi_tau,timeflag,test,test2,Damagemax
      !
      ! Parameters
      !
      real*8 zero,one,two,half,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0)
      !
      ! Read in some material properties
      !
      nu_s0  = props(5)
      !
      ! Identity matrix
      !
      call onem(Iden)
      !
      ! Obtain old and new deformation gradients
      !
      F_t = dfgrd0
      F_tau = dfgrd1
      !
      ! Obtain old and new temperature/soft volume fraction
      !
      theta_t = temp
      theta_tau = temp + dtemp
      !
      ! Obtain state variables from the previous increment
      !
      if((kinc.le.1).and.(kstep.eq.1)) then
            LambdaBar_t = one
            nu_t        = nu_s0
            Achainmax   = one
            numax = nu_s0
      else
            LambdaBar_t = statev(1)
            nu_t        = statev(2)
            Achainmax   = statev(3)
            numax = statev(4)
      endif
      !
      ! Compute the material constitutive response
      !
      flag = one
      !
      call integ(Props,nprops,flag,numax,
     +     F_tau,dtime,timeflag,theta_t,
     +     T_tau,LambdaBar_t,LambdaBar_tau,theta_tau,
     +     nu_t,nu_tau,Achainmax,rpl,drpldt)
      !
      ! Update state variables at the end of the increment
      !
         statev(1)  = LambdaBar_tau
         statev(2)  = nu_tau
         statev(3)  = Achainmax
         statev(4) = numax
      !
      ! Return Abaqus/Standard the Cauchy stress
      !
      if(ntens.eq.6) then
         !
         ! 3D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         stress(5) = T_tau(1,3)
         stress(6) = T_tau(2,3)
      elseif(ntens.eq.4) then
         !
         ! 2D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
      endif
      !
      ! Return Abaqus/Standard the stress-deformation jacobian
      !
      call mechJacobian(props,nprops,F_tau,F_t,
     +           dtime,ntens,T_tau,LambdaBar_t,nu_t,
     +           Achainmax,varphi_t,ddsdde)
      !
      return
      end subroutine umat
****************************************************************************
      !
      ! Constitutive model subroutine
      !
      subroutine integ(props,nprops,flag,numax,
     +     F_tau,dtime,timeflag,theta_t,
     +     T_tau,LambdaBar_t,LambdaBar_tau,theta_tau,
     +     nu_t,nu_tau,Achainmax,rpl,drpldt)
      !
      implicit none
      !
      integer nprops,i,j,k,l,m,n,nargs
      parameter(nargs=8)
      !
      real*8 props(nprops),Gshear,Kbulk,F_tau(3,3),Iden(3,3),T_tau(3,3)
      real*8 theta_tau,SpTanMod(3,3,3,3),detF,trC,TR_tau(3,3),C_tau(3,3)
      real*8 dTRdF(3,3,3,3),Finv(3,3),FinvT(3,3),dtime,Cbar_tau(3,3)
      real*8 Bbar(3,3),trBbar,dTKdF(3,3,3,3),Cinv(3,3),AC(3,3),trAC
      real*8 Teq_tau(3,3),Fbar(3,3),timeflag,detF_tau,C_tau_dis(3,3)
      real*8 F_tau_T(3,3),B_tau(3,3),B_tau_dis(3,3),trB_tau_dis
      real*8 B_tau_dis_0(3,3),C_tau_dis_0(3,3),LambdaBar_t,aux,lamL
      real*8 Amp,nu_s0,nu_ss,LambdaBar_tau,nu_t,nu_tau,Achainmax,X
      real*8 Achain,args(nargs),Teq_tau2(3,3),aux2,lamVL,aux3,dE
      real*8 flag,varphi_t,varphi_tau,Achainold,Xold,AA,lamdot
      real*8 B,D,NN,RR,KK,Achain_t,Achain_tau,trC_tau_dis
      real*8 rate,Q,Rgas,theta,mustar,test,test2,rpl,drpldt
      real*8 theta_t,numax,X_tau,X_t
      !
      ! Parameters
      !
      real*8 zero,one,two,half,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0)
      !
      ! Identity matrix
      !
      call onem(Iden)
      !
      ! Obtain material properties
      !
      Gshear = props(1)
      Kbulk  = props(2)
      lamL   = props(3)
      Amp    = props(4)
      nu_s0  = props(5)
      nu_ss  = props(6)
      dE     = props(7)
      !
      ! Compute the determinant of the deformation gradient
      !
      call mdet(F_tau,detF_tau)
      !
      ! Left Green Cauchy tensor, its distortional part, trace, and deviatoric part
      !
      B_tau = matmul(F_tau,transpose(F_tau))
      B_tau_dis = (detF_tau**(-two/three))* B_tau
      trB_tau_dis = B_tau_dis(1,1) + B_tau_dis(2,2) + B_tau_dis(3,3)
      B_Tau_dis_0 = B_Tau_dis - (one/three)* trB_tau_dis * Iden
      !
      ! Right Green Cauchy tensor, its distortional part, trace
      !
      C_tau = matmul(transpose(F_tau),F_tau)
      C_tau_dis = (detF_tau**(-two/three))* C_tau
      trC_tau_dis = C_tau_dis(1,1) + C_tau_dis(2,2) + C_tau_dis(3,3)
      !
      ! Effective stretch
      !
      LambdaBar_tau = dsqrt(third*trC_tau_dis)
      !
      if (flag.eq.one) then
      !
      ! ``Non-local'' soft volume fraction can only evolve (no Mullins recovery)
      !
      if (theta_tau.ge.numax) then
      nu_tau = theta_tau
      numax = nu_tau
      else
          nu_tau = numax
      endif
      !
      elseif(flag.eq.zero) then
      nu_tau = nu_t
      endif
      !
      ! Update the amplification factor
      !
      X_tau = one+3.5d0*(one-nu_tau)+18.0d0*(one-nu_tau)**(two)
      !
      X_t = one+3.5d0*(one-nu_t)+18.0d0*(one-nu_t)**(two)
      !
      ! Update the amplified stretch
      !
      Achain_tau = dsqrt(X_tau*(LambdaBar_tau*LambdaBar_tau-one)+one)
      !
      ! Compare the updated amplified stretch with the maximum stretch
      ! If the new value exceeds the prior maximum, then update both
      ! soft volume fraction and the maximum stretch
      ! Otherwise, leave them unchanged
      !
      if(Achain_tau.GE.Achainmax.and.flag.eq.one) then
         Achainmax = Achain_tau
      else
          nu_tau = nu_t
          X_tau = one+3.5d0*(one-nu_tau)+18.0d0*(one-nu_tau)**(two)
          LambdaBar_tau = dsqrt(third*trC_tau_dis)
          Achain_tau=dsqrt(X_tau*(LambdaBar_tau*LambdaBar_tau-one)+one)
      endif
      !
      ! Obtain old Achain
      !
      Achain_t = dsqrt(X_t*(LambdaBar_t * LambdaBar_t - one) + one)
      !
      ! Compute the source term and its derivative
      !
      if (Achain_tau.ge.Achainmax) then
          rpl=Amp*(nu_ss-nu_tau)*(lamL-one)*(Achain_tau-Achain_t)
     +    /dtime*(lamL-Achainmax)**two
          drpldt = -Amp*(lamL-one)*(Achain_tau-Achain_t)/
     +    (dtime*(lamL-Achainmax)**two)
      else
      ! ``Non-local'' soft volume fraction can only evolve when stretch exceeds prior maximum stretch
          nu_tau = nu_t
          rpl = numax - theta_tau
          drpldt = -one
      endif
      !
      ! Compute the inverse langevin function
      !
       aux=(three-(Achain_tau/lamL)**two)/(one-(Achain_tau/lamL)**two)
     +      *(Achain_tau/lamL)
      !
      ! Cauchy stress
      !
      T_tau = (one/detF_tau) * third * X_tau * nu_tau * Gshear *
     +    (lamL/Achain_tau) *aux * B_tau_dis_0
     +    + Kbulk * (detF_tau - one)*Iden
      !
      return
      end subroutine integ
***********************************************************************
      !
      ! Subroutine to compute the stress�deformation Jacobian numerically
      ! using a finite difference perturbation of the deformation gradient
      ! Additional implementation details and  background can be
      ! found in:
      ! Chester, S. A. (2011). Mechanics of amorphous polymers and polymer gels.
      ! Doctoral dissertation, Massachusetts Institute of Technology.
      !
      subroutine mechJacobian(props,nprops,F_tau,F_t,
     +           dtime,ntens,T_tau,LambdaBar_t,nu_t,
     +           Achainmax,varphi_t,jac)

      implicit none

      integer i,j,k,nprops,ntens,ii

      real*8 F_tau(3,3),F_t(3,3),F_per(3,3),T_tau(3,3),T_per(3,3),dE,tmp
      real*8 jac(ntens,ntens),E(3,3),props(nprops),dtime,S_t,S_per
      real*8 theta_tau,pnewdt_per,plasticpower_per,eBarP_t,eBarP_per
      real*8 Fp_t(3,3),Fp_per(3,3),eBarPdot_t,eBarPdot_per
      real*8 timeflag,theta_per,rpl_per,drpldt_per
      real*8 LambdaBar_t,nu_t,Achainmax,LambdaBar_per,mu_per
      real*8 Achainmax_per,varphi_t,varphi_per,test_per,test2_per
      real*8 theta_t,numax
      !
      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third,one_half
      real*8 delta
      real*8 flag
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0,
     +     one_half=0.5d0,delta=1.d-5)
      !
      ! Obtain properties
      dE= props(7)
      !
      ! Init
      !
      jac = zero
      !
      do i=1,ntens
      !
         E = zero
         F_per = zero
         T_per = zero
      !
         if(i.eq.1) then
            E(1,1) = dE
         elseif(i.eq.2) then
            E(2,2) = dE
         elseif(i.eq.3) then
            E(3,3) = dE
         elseif(i.eq.4) then
            E(1,2) = dE/two
            E(2,1) = E(1,2)
         elseif(i.eq.5) then
            E(1,3) = dE/two
            E(3,1) = E(1,3)
         else
            E(2,3) = dE/two
            E(3,2) = E(2,3)
         endif
         !
         ! Calculate the perturbed deformation gradient
         !
         F_per = F_tau + matmul(E,F_t)
         !
         ! Check that the perturbed deformation gradient is admissible
         !
         ii = 1
         call mdet(F_per,tmp)
         if(tmp.le.zero) then
            do while(tmp.le.zero)
               !
               ! Change the value of dE in this case
               !
               if(ii.gt.9) then
                  write(80,*) 'this will be a problem for detF_per<0'
                  write(*,*) 'this will be a problem for detF_per<0'
                  call xit
               endif
               E = E/four
               F_per = F_tau + matmul(E,F_t)
               ii = ii + 1
               call mdet(F_per,tmp)
            enddo
               endif
               !
               flag = zero
               !
               call integ(Props,nprops,flag,numax,
     +        F_per,dtime,timeflag,theta_t,
     +        T_per,LambdaBar_t,LambdaBar_per,theta_per,
     +        nu_t,mu_per,Achainmax_per,rpl_per,drpldt_per)
      !
      ! Calculate the tangent using finite difference
      !
      if (ntens .eq. 4) then
          jac(1,i) = (T_per(1,1) - T_tau(1,1)) / dE
          jac(2,i) = (T_per(2,2) - T_tau(2,2)) / dE
          jac(3,i) = (T_per(3,3) - T_tau(3,3)) / dE
          jac(4,i) = (T_per(1,2) - T_tau(1,2)) / dE
      else
          jac(1,i) = (T_per(1,1) - T_tau(1,1)) / dE
          jac(2,i) = (T_per(2,2) - T_tau(2,2)) / dE
          jac(3,i) = (T_per(3,3) - T_tau(3,3)) / dE
          jac(4,i) = (T_per(1,2) - T_tau(1,2)) / dE
          jac(5,i) = (T_per(1,3) - T_tau(1,3)) / dE
          jac(6,i) = (T_per(2,3) - T_tau(2,3)) / dE
      endif
      !
      enddo
      !
      return
      end subroutine mechJacobian
****************************************************************************
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,statev,temp,
     $     dtemp,dtemdx,time,dtime,predef,dpred,cmname,ntgrd,nstatv,
     $     props,nprops,coords,pnewdt,noel,npt,layer,kspt,kstep,kinc)
c
      include 'aba_param.inc'
c
      character*80 cmname
c
      dimension dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     $     dfdg(ntgrd,ntgrd),statev(nstatv),dtemdx(ntgrd),time(2),
     $     predef(1),dpred(1),props(nprops),coords(3)
      real*8 F(3,3),B(3,3),erff,lc,zeta,aux,detF,B_tau(3,3),detF_tau
      real*8 tau
      !
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0,six=6.d0,eight=8.d0)
c
c
      !
      ! Obtain material properties
      !
      lc = props(1)
      specht = props(2)
      tau = props(3)
      !
      aux = (lc**two)/tau
      !
      dudt = specht
      du = dudt*dtemp
      u = u+du
      !
      ! Update conduction term
      !
      DO I=1, NTGRD
         FLUX(I) = -DTEMDX(I)
         DFDG(I,I) = zero
      END DO
      FLUX = FLUX * aux
      DFDG = DFDG * aux
      !
      return
      end
***********************************************************************
c
c
c  The following are all utility routines used in fortran codes
c
c
c
C**********************************************************************
C	THE NEXT SUBROUTINE CALCULATES VARIOUS KINEMATICAL QUANTITIES 
C	ASSOCIATED WITH THE DEFORMATION GRADIENT
C**********************************************************************
	SUBROUTINE SKINEM(F,R,U,E)

C	THIS SUBROUTINE PERFORMS THE RIGHT POLAR DECOMPOSITION
C	[F] = [R][U] OF THE DEFORMATION GRADIENT [F] INTO
C	A ROTATION [R] AND THE RIGHT  STRETCH TENSOR [U].
C	THE EIGENVALUES AND EIGENVECTORS OF [U] AND
C	THE LOGARITHMIC STRAIN [E] = LN [U]
C	ARE ALSO RETURNED.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION F(3,3),FTRANS(3,3), C(3,3), OMEGA(3),
     +	          UEIGVAL(3),EIGVEC(3,3), EIGVECT(3,3), 
     +            U(3,3),E(3,3),UINV(3,3),R(3,3),TEMPM(3,3)
     
	COMMON/ERRORINFO/UMERROR     

C	F(3,3)	-- THE DEFORMATION GRADIENT MATRIX WHOSE
C		   POLAR DECOMPOSITION IS DESIRED.
C	DETF	-- THE DETRMINANT OF [F]; DETF > 0.
C	FTRANS(3,3)	-- THE TRANSPOSE OF [F].
C	R(3,3)	-- THE ROTATION MATRIX; [R]^T [R] = [I];
C		   OUTPUT.
C	U(3,3)	-- THE RIGHT STRETCH TENSOR; SYMMETRIC
C		   AND POSITIVE DEFINITE; OUTPUT.
C	UINV(3,3)	-- THE INVERSE OF [U].
C	C(3,3)	-- THE RIGHT CAUCHY-GREEN TENSOR = [U][U];
C		   SYMMETRIC AND POSITIVE DEFINITE.
C	OMEGA(3)-- THE SQUARES OF THE PRINCIPAL STRETCHES.
C 	UEIGVAL(3)	-- THE PRINCIPAL STRETCHES; OUTPUT.
C	EIGVEC(3,3)	-- MATRIX OF EIGENVECTORS OF [U];OUTPUT.
C	EIGVECT(3,3)    -- TRANSPOSE OF THE ABOVE.
C	E(3,3)	-- THE LOGARITHMIC STRAIN TENSOR, [E]=LN[U];
C		   OUTPUT.
C**********************************************************************

C	STORE THE IDENTITY MATRIX IN  [R], [U], AND [UINV]

	CALL ONEM(R)
	CALL ONEM(U)
	CALL ONEM(UINV)

C	STORE THE ZERO MATRIX IN [E]

	CALL ZEROM(E)

C      	CHECK IF THE DETERMINANT OF [F] IS GREATER THAN ZERO.
C	IF NOT, THEN PRINT DIAGNOSTIC AND STOP.

        CALL MDET(F,DETF)
        IF (DETF .LE. 0.D0) THEN
          WRITE(*,100)
          UMERROR=5.
          RETURN
        ENDIF

C      	CALCULATE THE RIGHT CAUCHY GREEN STRAIN TENSOR [C]

        CALL  MTRANS(F,FTRANS)
        CALL  MPROD(FTRANS,F,C)
 
C	CALCULATE THE EIGENVALUES AND EIGENVECTORS OF  [C]

	CALL SPECTRAL(C,OMEGA,EIGVEC)

C	CALCULATE THE PRINCIPAL VALUES OF [U] AND [E]

	UEIGVAL(1) = DSQRT(OMEGA(1))
	UEIGVAL(2) = DSQRT(OMEGA(2))
	UEIGVAL(3) = DSQRT(OMEGA(3))

	U(1,1) = UEIGVAL(1)
	U(2,2) = UEIGVAL(2)
	U(3,3) = UEIGVAL(3)

	E(1,1) = DLOG( UEIGVAL(1) )
	E(2,2) = DLOG( UEIGVAL(2) )
	E(3,3) = DLOG( UEIGVAL(3) )

C	CALCULATE THE COMPLETE TENSORS [U] AND [E]

	CALL MTRANS(EIGVEC,EIGVECT)
	CALL MPROD(EIGVEC,U,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,U)
	CALL MPROD(EIGVEC,E,TEMPM)
	CALL MPROD(TEMPM,EIGVECT,E)

C	CALCULATE [UINV]

	CALL M3INV(U,UINV)

C	CALCULATE [R]

	CALL MPROD(F,UINV,R)
100     FORMAT(5X,'--ERROR IN KINEMATICS-- THE DETERMINANT OF [F]',
     +         ' IS NOT GREATER THAN 0')

	RETURN
	END

C**********************************************************************
C	THE FOLLOWING SUBROUTINES CALCULATE THE SPECTRAL
C	DECOMPOSITION OF A SYMMETRIC THREE BY THREE MATRIX
C**********************************************************************
	SUBROUTINE SPECTRAL(A,D,V)
C
C	THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C	A SYMMETRIC 3 BY 3 MATRIX [A]. 
C
C	THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C	EIGENVALUES IN ASCENDING ORDER, AND
C	A MATRIX [V] WHOSE COLUMNS CONTAIN THE CORRESPONDING
C	EIGENVECTORS.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER(NP=3)
	DIMENSION D(NP),V(NP,NP)
	DIMENSION A(3,3),E(NP,NP)

	DO 2 I = 1,3
          DO 1 J= 1,3
            E(I,J) = A(I,J)
1	  CONTINUE
2	CONTINUE

	CALL JACOBI(E,3,NP,D,V,NROT)
	CALL EIGSRT(D,V,3,NP)

	RETURN
	END

C**********************************************************************
	SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

C	COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C	MATRIX [A], WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL 
C	NP BY BP ARRAY. ON OUTPUT, ELEMENTS OF [A] ABOVE THE DIAGONAL 
C	ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C	AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C	VECTOR D RETURNS THE EIGENVALUES OF [A] IN ITS FIRST N ELEMENTS.
C	[V] IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C	[A] WHOSE COLUMNS CONTAIN, ON OUTPUT, THE NORMALIZED
C	EIGENVECTORSOF [A]. NROT RETURNS THE NUMBER OF JACOBI ROTATIONS
C	WHICH WERE REQUIRED.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", PAGE 346.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	PARAMETER (NMAX =100)
	DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)

C	INITIALIZE [V] TO THE IDENTITY MATRIX

	DO 12 IP = 1,N	
	  DO 11 IQ = 1,N
	    V(IP,IQ) = 0.D0
11        CONTINUE
          V(IP,IP) = 1.D0
12	CONTINUE

C	INITIALIZE [B] AND [D] TO THE DIAGONAL OF [A], AND Z TO ZERO.
C	THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C	IN EQUATION (11.1.14)

	DO 13 IP = 1,N
	  B(IP) = A(IP,IP)
	  D(IP) = B(IP)
	  Z(IP) = 0.D0
13	CONTINUE
C
	NROT = 0
	DO 24 I = 1,50

C	SUM OFF-DIAGONAL ELEMENTS

          SM = 0.D0
          DO 15 IP = 1, N-1
            DO 14 IQ = IP + 1, N
	      SM = SM + DABS ( A(IP,IQ ))
14          CONTINUE
15        CONTINUE

C	IF SUM = 0., THEN RETURN. THIS IS THE NORMAL RETURN
C	WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE 
C	UNDERFLOW.

          IF ( SM .EQ. 0.D0) RETURN
C
C	  IF ( SM .LT. 1.0D-15) RETURN

C	IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C	|A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE, 
C	SEE EQUATION (11.1.25). THEREAFTER TRESH = 0.

          IF ( I .LT. 4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          ENDIF
C
          DO 22 IP = 1, N-1
            DO 21 IQ = IP+1,N
              G = 100.D0*DABS(A(IP,IQ))

C	AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C	OFF-DIAGONAL ELEMENT IS SMALL.

	      IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP)))
     +            .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G .EQ. DABS(H)) THEN

C	T = 1./(2.*THETA), EQUATION(11.1.10)

	          T =A(IP,IQ)/H
	        ELSE
	          THETA = 0.5D0*H/A(IP,IQ)
	          T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
	          IF (THETA .LT. 0.D0) T = -T
	        ENDIF
	        C = 1.D0/DSQRT(1.D0 + T**2)
	        S = T*C
	        TAU = S/(1.D0 + C)
	        H = T*A(IP,IQ)
	        Z(IP) = Z(IP) - H
	        Z(IQ) = Z(IQ) + H
	        D(IP) = D(IP) - H
	        D(IQ) = D(IQ) + H
	        A(IP,IQ) = 0.D0

C	CASE OF ROTATIONS 1 <= J < P
				
	        DO 16 J = 1, IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
16	        CONTINUE

C	CASE OF ROTATIONS P < J < Q

	        DO 17 J = IP+1, IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
17	        CONTINUE

C	CASE OF ROTATIONS Q < J <= N

	        DO 18 J = IQ+1, N
                  G = A(IP,J)
	          H = A(IQ,J)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(IQ,J) = H + S*(G - H*TAU)
18	        CONTINUE
	        DO 19 J = 1,N
	          G = V(J,IP)
	          H = V(J,IQ)
	          V(J,IP) = G - S*(H + G*TAU)
	          V(J,IQ) = H + S*(G - H*TAU)
19	        CONTINUE
	        NROT = NROT + 1
              ENDIF
21	    CONTINUE
22	  CONTINUE

C	UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z

	  DO 23 IP = 1, N
	    B(IP) = B(IP) + Z(IP)
	    D(IP) = B(IP)
	    Z(IP) = 0.D0
23	  CONTINUE
24	CONTINUE

C	IF THE ALGORITHM HAS REACHED THIS STAGE, THEN
C	THERE ARE TOO MANY SWEEPS, PRINT A DIAGNOSTIC
C	AND STOP.

	WRITE (*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'

	RETURN
	END

C**********************************************************************
	SUBROUTINE EIGSRT(D,V,N,NP)

C	GIVEN THE EIGENVALUES [D] AND EIGENVECTORS [V] AS OUTPUT FROM
C	JACOBI, THIS ROUTINE SORTS THE EIGENVALUES INTO ASCENDING ORDER, 
C	AND REARRANGES THE COLUMNS OF [V] ACCORDINGLY.

C	THIS SUBROUTINE IS TAKEN FROM "NUMERICAL RECIPES", P. 348.
C**********************************************************************

	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION D(NP),V(NP,NP)

	DO 13 I = 1,N-1
	  K = I
	  P = D(I)
	  DO 11 J = I+1,N
	    IF (D(J) .GE. P) THEN
	      K = J
	      P = D(J)
	    END IF
11	  CONTINUE
	  IF (K .NE. I) THEN
	    D(K) = D(I)
	    D(I) = P
	    DO 12 J = 1,N
	      P = V(J,I)
	      V(J,I) = V(J,K)
	      V(J,K) = P
12	    CONTINUE
  	  ENDIF
13	CONTINUE

	RETURN
	END

C**********************************************************************
C	THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************
	SUBROUTINE ZEROV(V,SIZE)

C	THIS SUBROUTINE STORES THE ZERO VECTOR IN A VECTOR V
C	OF SIZE SIZE.
C**********************************************************************

	INTEGER SIZE
	REAL*8 V(0:SIZE-1)

	DO 1 I=0,SIZE
	  V(I) = 0.D0
1	CONTINUE
	
	RETURN
	END

C**********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

        REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
	END

C**********************************************************************
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD(A,B,C)
 
C 	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
	    C(I,J) = 0.D0
	    DO 1 K = 1, 3
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD4(A,B,C)
 
C	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(4,4),B(4,4),C(4,4)

	DO 2 I = 1, 4
   	  DO 2 J = 1, 4
	    C(I,J) = 0.D0
	    DO 1 K = 1, 4
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE DOTPM(A,B,C)

C	THIS SUBROUTINE CALCULATES THE SCALAR PRODUCT OF TWO
C	3 BY 3 MATRICES [A] AND [B] AND STORES THE RESULT IN THE
C	SCALAR C.
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C

	C = 0.D0
	DO 1 I = 1,3
	  DO 1 J = 1,3
            C = C + A(I,J)*B(I,J)
1	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  write(*,10)
	  STOP
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END

C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
        SUBROUTINE INVAR(A,IA,IIA,IIIA)

C	THIS SUBROUTINE CALCULATES THE PRINCIPAL INVARIANTS 
C	IA, IIA, IIIA OF A TENSOR [A].
C**********************************************************************

        REAL*8 A(3,3), AD(3,3),AD2(3,3), DETA, IA,IIA,IIIA

        DO 1 I=1,3
          DO 1 J=1,3
            AD(I,J) = A(I,J)
1       CONTINUE
        IA = AD(1,1) + AD(2,2) + AD(3,3)

C	CALCULATE THE SQUARE OF [AD]

        CALL MPROD(AD,AD,AD2)
        IIA =0.5D0 * ( IA*IA - ( AD2(1,1) + AD2(2,2) + AD2(3,3) ) )

        CALL  MDET(AD,DETA)
        IIIA = DETA

        RETURN
        END

C**********************************************************************
	SUBROUTINE TRACEM(A,TRA)

C	THIS SUBROUTINE CALCULATES THE TRACE OF A 3 BY 3 MATRIX [A]
C	AND STORES THE RESULT IN THE SCALAR TRA
C**********************************************************************

	REAL*8 A(3,3),TRA

	TRA = A(1,1) + A(2,2) + A(3,3)

	RETURN 
	END

C**********************************************************************
	SUBROUTINE DEVM(A,ADEV)

C	THIS SUBROUTINE CALCULATES THE DEVIATORIC PART OF A
C	3 BY 3 MATRIX [A]
C**********************************************************************

	REAL*8 A(3,3),TRA,ADEV(3,3),IDEN(3,3)

	CALL TRACEM(A,TRA)
	CALL ONEM(IDEN)
	CALL ZEROM(ADEV)

	DO 1 I = 1,3
	  DO 1 J = 1,3
	    ADEV(I,J) = A(I,J) - (1.D0/3.D0)*TRA*IDEN(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE EQUIVS(S,SB)

C	THIS SUBROUTINE CALCULATES THE EQUIVALENT SHEAR STRESS SB
C	CORRESPONDING TO A 3 BY 3 STRESS MATRIX [S]
C**********************************************************************

	REAL*8 S(3,3),SDEV(3,3),SDOTS,SB

	SB = 0.D0
	SDOTS = 0.D0

	CALL DEVM(S,SDEV)
	CALL DOTPM(SDEV,SDEV,SDOTS)
	SB = DSQRT(1.5D0*SDOTS)

	RETURN
	END
C **********************************************************************
        SUBROUTINE PRESSURE(A,PRA)
C
C       THIS SUBROUTINE CALCULATES THE MEAN NORMAL PRESSURE
C       OF A 3 BY 3 MATRIX [A]
C       AND STORES THE RESULT IN THE SCALAR PRA
C ----------------------------------------------------------------------
C       VARIABLES
C
        REAL*8 A(3,3),PRA

        PRA = -(1.D0 / 3.D0)*( A(1,1) + A(2,2) + A(3,3) )

        RETURN 
        END
C**********************************************************************
	SUBROUTINE PRTMAT(A,M,N)
C**********************************************************************

	INTEGER M,N
	REAL*8 A(M,N)	  

	DO 10 K=1,M
	  WRITE(80,'(2X,6E12.4,2X)') (A(K,L), L=1,N)
10      CONTINUE

        RETURN
        END

C**********************************************************************
	SUBROUTINE PRTVEC(A,M)
C**********************************************************************

	INTEGER M
	REAL*8 A(M)	  

	WRITE(80,'(2X,6E12.4,2X)') (A(K), K=1,M)

        RETURN
	END
C*************************************************************************	  
c***********************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
C
C	Given an NxN matrix [A], with physical dimension NP, this 
C	routine replaces it by the LU decomposition of a row-wise 
C	permutation of itself. [A] and N are input. [A] is output, 
C	arranged in LU form. INDX is an output vector which records
C	the row permutation effected by the partial pivoting; 
C	D is output as +1 or -1 depending on wheter the nuber of
C	row interchanges was even or odd, respectively. This routine
C	is used in combination with LUBKSB to solve linear equations 
C	or invert a matrix.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	
	PARAMETER (NMAX=100,TINY=1.0E-20)
	DIMENSION A(NP,NP),INDX(N),VV(NMAX)
	D=1.
	DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12	CONTINUE
	DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16	CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19	CONTINUE
	IF(A(N,N).EQ.0.)A(N,N)=TINY
	RETURN
	END
C**********************************************************
       SUBROUTINE LUBKSB(A,N,NP,INDX,B)
C
C	Solves the set of N linear equations [A]{X} = {B}. 
C	Here [A] is input, not as the matrix [A], but as its LU 
C	decomposition, determined by the routine LUDCMP. INDX
C	is input as the permutation vector returned by LUDCMP. {B}
C	is input as the right-hand side vector {B}, and returns
C	with the solution vector {X}. [A], N, NP, INDX are not 
C	modified by this routine, and can be left in place
C	for succesive calls with different right-hand sides {B}.
C	This routine takes into account that {B} will begin with
C	many zero elements, so it is efficient for use in matrix 
C	inversion.
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NP,NP),INDX(N),B(N)
       II=0
       DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12     CONTINUE
       DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14     CONTINUE
       RETURN
       END
c**********************************************************

!
***********************************************************************
