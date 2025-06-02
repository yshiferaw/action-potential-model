c     ================================================================
c     ION CHANNEL SUBROUTINES FOR CARDIAC ELECTROPHYSIOLOGY
c     ================================================================
c     
c     These subroutines calculate ionic currents and update gating
c     variables for various ion channels in cardiac myocytes.
c     Each follows the Hodgkin-Huxley formalism with voltage-dependent
c     rate constants and exponential relaxation to steady states.
c     ================================================================

c     ================================================================
c     FAST SODIUM CURRENT (INa)
c     ================================================================
c     Models the voltage-gated sodium channel responsible for the
c     rapid upstroke of the cardiac action potential. Uses three
c     gating variables: m (activation), h (fast inactivation), 
c     j (slow inactivation).
c     
c     Current: INa = gNa * m³ * h * j * (V - ENa)
c     ================================================================

      subroutine ina(hode,v,frt,xh,xj,xm,xnai,xnao,xina)
      implicit double precision (a-h,o-z)

c     ============================================================
c     SODIUM CHANNEL PARAMETERS
c     ============================================================
      gna = 12.0d0        ! maximum sodium conductance (mS/μF)
      XKMCAM = 0.3d0      ! CaMKII-related parameter (unused here)
      deltax = -0.18d0    ! CaMKII shift parameter (unused here)
       
c     ============================================================
c     SODIUM REVERSAL POTENTIAL
c     ============================================================
      ! Nernst equation: ENa = (RT/F) * ln([Na]o/[Na]i)
      ena = (1.0d0/frt)*dlog(xnao/xnai)

c     ============================================================
c     ACTIVATION GATE (m) - FAST ACTIVATION
c     ============================================================
      ! Controls rapid channel opening during depolarization
      am = 0.32d0*(v+47.13d0)/(1.0d0-dexp(-0.1d0*(v+47.13d0)))
      bm = 0.08d0*dexp(-v/11.0d0)

c     ============================================================
c     CaMKII MODULATION (CURRENTLY DISABLED)
c     ============================================================
      ! These lines would modify channel behavior based on CaMKII activity
      ! camfact = 1.0d0/(1.0d0+(XKMCAM/caM)**4)
      ! vshift = -3.25*camfact
      
      camfact = 0.0       ! no CaMKII effect
      vshift = 0.0        ! no voltage shift
      vx = v - vshift     ! effective voltage

c     ============================================================
c     INACTIVATION GATES (h and j) - VOLTAGE-DEPENDENT KINETICS
c     ============================================================
      if(vx.lt.(-40.0d0)) then
         ! ================================
         ! HYPERPOLARIZED VOLTAGE RANGE (V < -40 mV)
         ! ================================
         
         ! Fast inactivation (h gate) rate constants
         ah = 0.135*dexp((80.0+vx)/(-6.8d0))
         bh = 3.56*dexp(0.079*vx) + 310000.0d0*dexp(0.35d0*vx)

         ! Slow inactivation (j gate) rate constants
         aj1a = -127140.0*dexp(0.2444*vx)
         aj1b = 0.00003474d0*dexp(-0.04391d0*vx)
         aj1c = (vx+37.78)/(1.0d0+dexp(0.311*(vx+79.23)))

         aj = (1.0d0+camfact*deltax)*(aj1a-aj1b)*aj1c
      bj = (0.1212*dexp(-0.01052*vx))/(1.0+dexp(-0.1378d0*(vx+40.14d0)))

      else
         ! ================================
         ! DEPOLARIZED VOLTAGE RANGE (V >= -40 mV)
         ! ================================
         
         ! Fast inactivation (h gate)
         ah = 0.0d0
         bh = 1.0d0/(0.13d0*(1.0d0+dexp((vx+10.66)/(-11.1d0))))
         
         ! Slow inactivation (j gate)
         aj = 0.0d0
         
         bja = 0.3*dexp(-0.0000002535d0*vx)
         bjb = 1.0 + dexp(-0.1d0*(vx+32.0d0))
         bj = bja/bjb
      endif

c     ============================================================
c     TIME CONSTANTS
c     ============================================================
      tauh = 1.0d0/(ah+bh)    ! h gate time constant (ms)
      tauj = 1.0d0/(aj+bj)    ! j gate time constant (ms)
      taum = 1.0d0/(am+bm)    ! m gate time constant (ms)

      ! These assignments don't change the values but are in original
      tauj = tauj
      tauh = tauh

c     ============================================================
c     SODIUM CURRENT CALCULATION
c     ============================================================
      ! Hodgkin-Huxley current: I = g * m³ * h * j * (V - E)
      xina = gna*xh*xj*xm*xm*xm*(v-ena)

c     ============================================================
c     GATING VARIABLE UPDATES
c     ============================================================
      ! Exponential approach to steady state
      xh = ah/(ah+bh) - ((ah/(ah+bh))-xh)*dexp(-hode/tauh)
      xj = aj/(aj+bj) - ((aj/(aj+bj))-xj)*dexp(-hode/tauj)
      xm = am/(am+bm) - ((am/(am+bm))-xm)*dexp(-hode/taum)

      return
      end 

c     ================================================================
c     RAPID DELAYED RECTIFIER POTASSIUM CURRENT (IKr)
c     ================================================================
c     Models the hERG channel responsible for action potential
c     repolarization. Shows inward rectification due to voltage-dependent
c     block. Critical for QT interval and arrhythmia susceptibility.
c     ================================================================

      subroutine ikr(hode,v,frt,xko,xki,xr,xikr)
      implicit double precision (a-h,o-z)

c     ============================================================
c     POTASSIUM REVERSAL POTENTIAL
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)

c     ============================================================
c     EXTRACELLULAR K+ DEPENDENCE AND KINETICS
c     ============================================================
      gss = dsqrt(xko/5.40)    ! conductance scales with sqrt([K]o)
      
      ! Voltage-dependent rate constants
      xkrv1 = 0.00138d0*(v+7.0d0)/(1.0-dexp(-0.123*(v+7.0d0)))
      xkrv2 = 0.00061d0*(v+10.0d0)/(dexp(0.145d0*(v+10.0d0))-1.0d0)
      taukr = 1.0d0/(xkrv1+xkrv2)

c     ============================================================
c     STEADY STATE AND RECTIFICATION
c     ============================================================
      xkrinf = 1.0d0/(1.0d0+dexp(-(v+50.0d0)/7.5d0))  ! activation
      rg = 1.0d0/(1.0d0+dexp((v+33.0d0)/22.4d0))       ! rectification

c     ============================================================
c     IKr CURRENT CALCULATION
c     ============================================================
      gkr = 0.007836d0    ! maximum conductance
      xikr = gkr*gss*xr*rg*(v-ek)

c     ============================================================
c     GATING VARIABLE UPDATE
c     ============================================================
      xr = xkrinf - (xkrinf-xr)*dexp(-hode/taukr)

      return
      end

c     ================================================================
c     SLOW DELAYED RECTIFIER POTASSIUM CURRENT (IKs)
c     ================================================================
c     Models the slow component of delayed rectifier K current.
c     Important for repolarization, especially at fast rates.
c     Shows calcium-dependent regulation.
c     ================================================================

      subroutine iks(hode,v,frt,ci,xnao,xnai,xko,xki,xs1,qks,xiks)
      implicit double precision (a-h,o-z)

c     ============================================================
c     PARAMETERS
c     ============================================================
      prnak = 0.018330d0   ! relative Na+ permeability

c     ============================================================
c     CALCIUM-DEPENDENT MODULATION
c     ============================================================
      ! qks_inf = 0.6d0*ci     ! original Ca-dependent form
      qks_inf = 0.6d0*0.30d0   ! fixed value used in code
      tauqks = 1000.0d0        ! time constant for Ca modulation

c     ============================================================
c     REVERSAL POTENTIAL AND KINETICS
c     ============================================================
      ! Includes Na+ permeability
      eks = (1.0d0/frt)*dlog((xko+prnak*xnao)/(xki+prnak*xnai))
      
      ! Steady state activation
      xs1ss = 1.0/(1.0+dexp(-(v-1.50d0)/16.70d0))
      xs2ss = xs1ss

c     ============================================================
c     TIME CONSTANT CALCULATION
c     ============================================================
      tauxs = 1.0d0/(0.0000719*(v+30.0d0)/(1.0d0-dexp(
     +     -0.148d0*(v+30.0))) + 0.000131d0
     +     *(v+30.0d0)/(dexp(0.0687d0*(v+30.0d0))-1.0d0))

c     ============================================================
c     IKs CURRENT CALCULATION
c     ============================================================
      gksx = 0.200d0    ! maximum conductance
      xiks = gksx*qks*xs1**2*(v-eks)

c     ============================================================
c     GATING VARIABLE UPDATES
c     ============================================================
      xs1 = xs1ss - (xs1ss-xs1)*dexp(-hode/tauxs)
      qks = qks + hode*(qks_inf-qks)/tauqks

      return
      end

c     ================================================================
c     INWARD RECTIFIER POTASSIUM CURRENT (IK1)
c     ================================================================
c     Models the strong inward rectifier that maintains resting
c     potential. Shows voltage-dependent block by Mg2+ and polyamines.
c     ================================================================

      subroutine ik1(hode,v,frt,xki,xko,xik1)
      implicit double precision (a-h,o-z)

c     ============================================================
c     REVERSAL POTENTIAL AND CONDUCTANCE
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)

      gkix = 0.60d0    ! base conductance (reduced in this model)
      gki = gkix*(dsqrt(xko/5.4))    ! [K]o dependence

c     ============================================================
c     RECTIFICATION KINETICS
c     ============================================================
      aki = 1.02/(1.0+dexp(0.2385*(v-ek-59.215)))
      bki = (0.49124*dexp(0.08032*(v-ek+5.476)) + dexp(0.061750
     +     *(v-ek-594.31)))/(1.0+dexp(-0.5143*(v-ek+4.753)))
      
      xkin = aki/(aki+bki)    ! steady state open probability

c     ============================================================
c     IK1 CURRENT CALCULATION
c     ============================================================
      xik1 = gki*xkin*(v-ek)

      return
      end

c     ================================================================
c     ULTRA-RAPID DELAYED RECTIFIER POTASSIUM CURRENT (IKur)
c     ================================================================
c     Models ultra-rapid component, important in atrial cells.
c     Replaces slow Ito in this atrial model.
c     ================================================================

      subroutine ikur(hode,v,frt,xki,xko,xkur,ykur,xikur)
      implicit double precision (a-h,o-z)

c     ============================================================
c     CONDUCTANCE PARAMETERS
c     ============================================================
      ! Note: replaces Ito slow in atrial model
      Gkur = 0.05*0.7
      Gkur = Gkur*1.2

c     ============================================================
c     REVERSAL POTENTIAL
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)

c     ============================================================
c     ACTIVATION KINETICS (xkur gate)
c     ============================================================
      rh1 = (v+6.0)/(-8.6)
      xkurss = 1.0/(1.0+dexp(rh1))
      tauxkur = 9.0d0/(1.0+dexp((v+5.0d0)/12.0)) + 0.5

c     ============================================================
c     INACTIVATION KINETICS (ykur gate)
c     ============================================================
      ykurss = 1.0d0/(1.0d0+dexp((v+7.5)/10.0d0))
      tauykur = 590.0d0/(1.0+dexp((v+60.0)/10.0)) + 3050.0d0

c     ============================================================
c     GATING VARIABLE UPDATES
c     ============================================================
      xkur = xkurss - (xkurss-xkur)*dexp(-hode/tauxkur)
      ykur = ykurss - (ykurss-ykur)*dexp(-hode/tauykur)

c     ============================================================
c     IKur CURRENT CALCULATION
c     ============================================================
      xikur = Gkur*xkur*ykur*(v-ek)

      return
      end

c     ================================================================
c     TRANSIENT OUTWARD POTASSIUM CURRENT (Ito)
c     ================================================================
c     Models fast transient outward current responsible for
c     early repolarization and action potential notch.
c     ================================================================

      subroutine ito(hode,v,frt,xki,xko,xtof,ytof,xito)
      implicit double precision (a-h,o-z)

c     ============================================================
c     REVERSAL POTENTIAL AND CONDUCTANCE
c     ============================================================
      ek = (1.0d0/frt)*dlog(xko/xki)

      gtof = 0.03d0*1.5*0.8
      gtof = gtof*0.9

c     ============================================================
c     ACTIVATION KINETICS (xtof gate)
c     ============================================================
      rg1 = -(v+1.0)/11.0
      xtof_inf = 1.0/(1.0+dexp(rg1))
      
      rg3 = -(v/30.0)**2
      txf = 3.5*dexp(rg3) + 1.5

c     ============================================================
c     INACTIVATION KINETICS (ytof gate)
c     ============================================================
      rg2 = (v+40.5)/11.5
      ytof_inf = 1.0/(1.0+dexp(rg2))

      rg4 = -((v+52.45)/15.8827)**2
      tyf = 25.635*dexp(rg4) + 24.14

c     ============================================================
c     CURRENT CALCULATION AND GATING UPDATES
c     ============================================================
      xitof = gtof*xtof*ytof*(v-ek)
      
      xtof = xtof_inf - (xtof_inf-xtof)*dexp(-hode/txf)
      ytof = ytof_inf - (ytof_inf-ytof)*dexp(-hode/tyf)

      xito = xitof    ! total ito (slow component removed in Grandi)

      return
      end

c     ================================================================
c     SODIUM-POTASSIUM PUMP CURRENT (INaK)
c     ================================================================
c     Models electrogenic Na/K-ATPase that maintains ionic gradients.
c     Exchanges 3 Na+ for 2 K+, creating net outward current.
c     ================================================================

      subroutine inak(v,frt,xko,xnao,xnai,xinak)
      implicit double precision (a-h,o-z)

c     ============================================================
c     PUMP PARAMETERS (Shannon model)
c     ============================================================
      xibarnak = 1.50d0     ! maximum pump current
      xkmko = 1.5d0         ! K+ half-saturation (adjusted for restitution)
      xkmnai = 12.0d0       ! Na+ half-saturation (adjusted for restitution)
      hh = 1.0d0            ! Na+ dependence exponent

c     ============================================================
c     VOLTAGE-DEPENDENT REGULATION
c     ============================================================
      sigma = (dexp(xnao/67.3d0)-1.0d0)/7.0d0
      fnak = 1.0d0/(1.0 + 0.1245*dexp(-0.1*v*frt)
     +      + 0.0365*sigma*dexp(-v*frt))

c     ============================================================
c     PUMP CURRENT CALCULATION
c     ============================================================

      xrv=xko/(xko+xkmko)
      xinak = xibarnak*fnak*(1.0/(1.0+(xkmnai/xnai)**hh))*xrv

      return
      end
