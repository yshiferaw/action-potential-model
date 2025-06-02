c     ================================================================
c     CALCIUM HANDLING AND AUXILIARY SUBROUTINES
c     ================================================================
c     
c     These subroutines handle various aspects of calcium dynamics,
c     stochastic processes, and related calculations in the cardiac
c     myocyte model.
c     ================================================================

c     ================================================================
c     BINOMIAL EVOLUTION FOR STOCHASTIC PROCESSES
c     ================================================================
c     
c     Handles stochastic evolution of calcium spark clusters using 
c     binomial statistics for activation and termination events.
c     ================================================================
	
      subroutine binevol(nt,nx,alpha,beta,dt,ndeltap,ndeltam)                  
      implicit double precision (a-h,o-z)

c     Calculate number of inactive clusters available for activation
      na = nt - nx
      
      if(na.gt.0) then
         xrate = alpha*dt
         call binom(na,xrate,ndeltap)
      else
         ndeltap = 0
      endif 
      
c     Calculate number of active clusters available for termination
200   if(nx.gt.0) then
         xrate = beta*dt
         call binom(nx,xrate,ndeltam)
      else
         ndeltam = 0
      endif 
      
      return
      end

c     ================================================================
c     SARCOPLASMIC RETICULUM CALCIUM UPTAKE (SERCA PUMP)
c     ================================================================
c     
c     Models SR Ca-ATPase pump using Hill kinetics with cooperativity.
c     Equation: J_uptake = Vmax * [Ca]^n / (Km^n + [Ca]^n)
c     ================================================================

      subroutine uptake(ci,vup,xup)
      implicit double precision (a-h,o-z)

      double precision, parameter::Ki=0.30d0      ! half-saturation (μM)
      double precision, parameter::Knsr=800.0d0   ! NSR threshold (μM)
      double precision, parameter::HH=3.00d0      ! Hill coefficient

c     Hill kinetics with cooperativity (HH=3)
      xup = vup*ci**HH/(Ki**HH+ci**HH)
      
      return
      end

c     ================================================================
c     TOTAL CALCIUM CONCENTRATION (FREE + BUFFERED)
c     ================================================================
c     
c     Converts free calcium to total by accounting for binding to
c     calmodulin and SR proteins using rapid equilibrium assumption.
c     ================================================================

      subroutine total(ci,cit)
      implicit double precision (a-h,o-z)

c     Buffer parameters
      bcal = 24.0d0       ! calmodulin total concentration (μM)
      xkcal = 7.0d0       ! calmodulin binding affinity (μM)
      srmax = 47.0d0      ! SR buffer total concentration (μM)
      srkd = 0.6d0        ! SR buffer binding affinity (μM)
       
c     Calculate bound calcium for each buffer system
      bix = bcal*ci/(xkcal+ci)    ! calmodulin-bound calcium
      six = srmax*ci/(srkd+ci)    ! SR protein-bound calcium

c     Total calcium = free + all bound forms
      cit = ci + bix + six

      return
      end 

c     ================================================================
c     FREE CALCIUM CONCENTRATION (INVERSE OF BUFFERING)
c     ================================================================
c     
c     Converts total calcium to free calcium using analytical 
c     approximation to avoid iterative root finding.
c     ================================================================

      subroutine xfree(cit,ci)
      implicit double precision (a-h,o-z)

c     Analytical approximation coefficients
      a = 2.23895
      b = 52.0344
      c = 0.666509

      y = cit

c     Solve quadratic approximation using quadratic formula
      xa = (b+a*c-y)**2 + 4.0*a*c*y
      ci = (-b-a*c+y+dsqrt(xa))/(2.0*a)

      return
      end 

c     ================================================================
c     SODIUM-CALCIUM EXCHANGE CURRENT (NCX)
c     ================================================================
c     
c     Models electrogenic Na/Ca exchanger (3:1 stoichiometry).
c     Can operate forward (Ca efflux) or reverse (Ca influx).
c     ================================================================

      subroutine inaca(v,frt,xnai,xnao,cao,ci,xinacaq)
      implicit double precision (a-h,o-z)

c     Convert calcium units: μM to mM
      cim = ci/1000.0d0

c     Electrochemical driving force calculation
      zw3a = xnai**3*cao*dexp(v*0.35*frt)          ! forward mode
      zw3b = xnao**3*cim*dexp(v*(0.35-1.)*frt)     ! reverse mode
      zw3 = zw3a - zw3b                            ! net driving force
      zw4 = 1.0d0 + 0.2d0*dexp(v*(0.35d0-1.0d0)*frt)  ! voltage factor

c     Calcium-dependent allosteric regulation
      xkdna = 0.3d0    ! regulatory Ca binding affinity (μM)
      aloss = 1.0d0/(1.0d0+(xkdna/ci)**3)

c     Competitive inhibition parameters (all in mM)
      xmcao = 1.3d0      ! extracellular Ca2+ binding
      xmnao = 87.5d0     ! extracellular Na+ binding  
      xmnai = 12.3d0     ! intracellular Na+ binding
      xmcai = 0.0036d0   ! intracellular Ca2+ binding

c     Denominator terms for competitive binding
      yz1 = xmcao*xnai**3 + xmnao**3*cim
      yz2 = xmnai**3*cao*(1.0d0+cim/xmcai)
      yz3 = xmcai*xnao**3*(1.0d0+(xnai/xmnai)**3)
      yz4 = xnai**3*cao + xnao**3*cim
      zw8 = yz1 + yz2 + yz3 + yz4

c     NCX current calculation
      xinacaq = aloss*zw3/(zw4*zw8)

      return
      end 

c     ================================================================
c     L-TYPE CALCIUM CURRENT
c     ================================================================
c     
c     Calculates L-type Ca current using Goldman-Hodgkin-Katz 
c     equation for divalent ions.
c     ================================================================

      subroutine ica(v,frt,cao,ci,pox,rca, xicaq)
      implicit double precision (a-h,o-z)

c     Physical constants and parameters
      xf = 96.485d0        ! Faraday's constant
      pca = 0.00054d0      ! Ca2+ permeability (from Luo-Rudy)
      za = v*2.0d0*frt     ! voltage factor for divalent ion

      factor1 = 4.0d0*pca*xf*frt
      factor = v*factor1

c     Convert calcium units: μM to mM
      cim = ci/1000.0d0

c     Goldman-Hodgkin-Katz equation with special case for V≈0       
      if(dabs(za).lt.0.001d0) then
         rca = factor1*(cim*dexp(za)-0.341d0*(cao))/(2.0d0*frt)
      else   
         rca = factor*(cim*dexp(za)-0.341d0*(cao))/(dexp(za)-1.0d0)         
      endif 

c     L-type calcium current
      xicaq = rca*pox
      
      return
      end 

c     ================================================================
c     L-TYPE CALCIUM CHANNEL MARKOV MODEL
c     ================================================================
c     
c     Detailed 10-state Markov model for L-type Ca channels.
c     Includes boundary-facing and spark-facing channel populations.
c     ================================================================
                   
      subroutine markov(hode,v,ci,c1,c2,xi1,xi2,po,
     + c1s,c2s,xi1s,xi2s,pos,alpha,bts)
      implicit double precision (a-h,o-z)

c     ============================================================
c     CALCIUM-INDEPENDENT RATE CONSTANTS
c     ============================================================      
      a23 = 0.3
      a32 = 3.0
      a42 = 0.00224d0
  	  
c     ============================================================
c     VOLTAGE-DEPENDENT ACTIVATION
c     ============================================================
      vth = 0.0d0
c     vth = -5.0
      s6 = 6.0d0
      poinf = 1.0d0/(1.0d0+dexp(-(v-vth)/s6))
      taupo = 1.0d0   

      a12 = poinf/taupo
      a21 = (1.0d0-poinf)/taupo

c     ============================================================
c     VOLTAGE-DEPENDENT INACTIVATION
c     ============================================================
      vy = -40.0d0
      sy = 4.0d0
      prv = 1.0d0-1.0d0/(1.0d0+dexp(-(v-vy)/sy))

      vyr = -40.0d0
      syr = 10.0d0 
      recov = 10.0d0+4954.0d0*dexp(v/15.6d0)
      recov = 2.0*recov

      tauba = (recov-450.0d0)*prv+450.0d0
      poix = 1.0d0/(1.0d0+dexp(-(v-vyr)/syr))

      a15 = poix/tauba
      a51 = (1.0d0-poix)/tauba

      vx = -40.0d0
      sx = 3.0d0
      tau3 = 3.0d0
      poi = 1.0d0/(1.0d0+dexp(-(v-vx)/sx))
      a45 = (1.0d0-poi)/tau3

c     ============================================================
c     CALCIUM-DEPENDENT RATES
c     ============================================================
      cat = 0.80d0

      zxr = 0.17
      fca = 1.0d0/(1.0d0+(cat/ci)**2)
      a24 = 0.00413d0+zxr*fca
      a34 = 0.00195+zxr*fca

      a43 = a34*(a23/a32)*(a42/a24)
      a54 = a45*(a51/a15)*(a24/a42)*(a12/a21)

      fcax = 1.0d0
      a24s = 0.00413d0+zxr*fcax
      a34s = 0.00195+zxr*fcax

      a43s = a34s*(a23/a32)*(a42/a24s)
      a54s = a45*(a51/a15)*(a24s/a42)*(a12/a21)

c     ============================================================
c     STATE DYNAMICS
c     ============================================================	
      dpo = a23*c1+a43*xi1-(a34+a32)*po-alpha*po+bts*pos
      dc2 = a21*c1+a51*xi2-(a15+a12)*c2+bts*c2s
      dc1 = a12*c2+a42*xi1+a32*po-(a21+a23+a24)*c1+bts*c1s
      dxi1 = a24*c1+a54*xi2+a34*po-(a45+a42+a43)*xi1+bts*xi1s

      dpos = a23*c1s+a43s*xi1s-(a34s+a32)*pos+alpha*po-bts*pos
      dc2s = a21*c1s+a51*xi2s-(a15+a12)*c2s-bts*c2s
      dc1s = a12*c2s+a42*xi1s+a32*pos-(a21+a23+a24s)*c1s-bts*c1s
      dxi1s = a24s*c1s+ a54s*xi2s+a34s*pos-(a45+a42+a43s)*xi1s-bts*xi1s
      dxi2s = a45*xi1s+ a15*c2s -(a51+a54s)*xi2s-bts*xi2s

c     ============================================================
c     STATE VARIABLE UPDATES
c     ============================================================
      po = po+dpo*hode
      c1 = c1+dc1*hode
      c2 = c2+dc2*hode
      xi1 = xi1+dxi1*hode
c     xi2 = xi2+dxi2*hode

      pos = pos+dpos*hode
      c1s = c1s+dc1s*hode
      c2s = c2s+dc2s*hode
      xi1s = xi1s+dxi1s*hode
      xi2s = xi2s+dxi2s*hode	

c     Probability conservation constraint
      xi2 = 1.0d0-c1-c2-po-xi1-pos-c1s-c2s-xi1s-xi2s	
         
      return
      end
