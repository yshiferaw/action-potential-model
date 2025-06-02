c     ================================================================
c     CARDIAC ELECTROPHYSIOLOGY SIMULATION WITH CALCIUM DYNAMICS
c     ================================================================
c     
c     This code simulates electrical activity and calcium handling in
c     cardiac tissue using a 2D spatial grid. It models:
c     - Action potential propagation
c     - Calcium-induced calcium release (CICR)
c     - Ion channel dynamics (Na, K, Ca channels)
c     - Sarcoplasmic reticulum (SR) calcium handling
c     - Calcium sparks and their stochastic properties
c
c     Dependencies: ifx wave2.f atissue-ca3.f currentsv1.f binom.f
c     ================================================================

      implicit double precision (a-h,o-z)
      
c     ================================================================
c     SIMULATION GRID AND TIMING PARAMETERS
c     ================================================================
      parameter (Lx=10,Ly=10)     ! 2D spatial grid dimensions
      parameter(nstim=5)          ! number of stimulus beats
    
c     ================================================================
c     VARIABLE DECLARATIONS
c     ================================================================

      ! Beat cycle length array
      double precision rbcl(nstim)
      
      ! Character arrays (legacy - may be unused)
      character*1 cam(0:9)
      character*15 filenm1(999),filenm2(999), filenm3(999), filenm4(999)
      character*15 filenm5(999),filenm6(999)
      double precision tstim(nstim),ampct(1000)

      ! ============================================================
      ! VOLTAGE AND MEMBRANE POTENTIAL ARRAYS
      ! ============================================================
      double precision v(0:Lx+1,0:Ly+1)    ! membrane voltage (mV)
      double precision vnew(0:Lx+1,0:Ly+1) ! updated voltage after diffusion
      double precision dv(0:Lx,0:Ly)       ! voltage time derivative
      double precision vold(0:Lx,0:Ly)     ! voltage from previous time step
      
      ! ============================================================
      ! ION CHANNEL GATING VARIABLES
      ! ============================================================
      ! These represent the probability of various ion channels being open
      double precision xm(Lx,Ly)    ! INa activation gate
      double precision xh(Lx,Ly)    ! INa fast inactivation gate
      double precision xj(Lx,Ly)    ! INa slow inactivation gate
      double precision xr(Lx,Ly)    ! IKr (rapid delayed rectifier K) gate
      double precision xs1(Lx,Ly)   ! IKs (slow delayed rectifier K) gate
      double precision qks(Lx,Ly)   ! Calcium-dependent gating of IKs
      double precision xkur(Lx,Ly)  ! IKur (ultra-rapid K) activation
      double precision ykur(Lx,Ly)  ! IKur inactivation
      double precision xtof(Lx,Ly)  ! Ito (transient outward K) activation
      double precision ytof(Lx,Ly)  ! Ito inactivation

      ! ============================================================
      ! CALCIUM CONCENTRATION VARIABLES
      ! ============================================================
      ! The model divides each cell into boundary and interior regions
      double precision cb(Lx,Ly)    ! boundary cytosolic Ca concentration (μM)
      double precision ci(Lx,Ly)    ! interior cytosolic Ca concentration (μM)
      double precision csrb(Lx,Ly)  ! boundary SR Ca concentration (μM)
      double precision csri(Lx,Ly)  ! interior SR Ca concentration (μM)
      double precision cnsr(Lx,Ly)  ! network SR Ca concentration (μM)
      double precision cit(Lx,Ly)   ! total interior Ca (free + buffered)
      double precision cbt(Lx,Ly)   ! total boundary Ca (free + buffered)

      ! ============================================================
      ! L-TYPE CALCIUM CHANNEL (LCC) MARKOV MODEL STATES
      ! ============================================================
      ! These track the stochastic states of L-type Ca channels
      double precision po(Lx,Ly)    ! open probability (boundary-facing)
      double precision c1(Lx,Ly)    ! closed state 1 (boundary)
      double precision c2(Lx,Ly)    ! closed state 2 (boundary)
      double precision xi1(Lx,Ly)   ! inactivated state 1 (boundary)
      double precision xi2(Lx,Ly)   ! inactivated state 2 (boundary)
      
      ! LCC states for spark-facing channels
      double precision c1s(Lx,Ly)   ! closed state 1 (spark-facing)
      double precision c2s(Lx,Ly)   ! closed state 2 (spark-facing)
      double precision xi1s(Lx,Ly)  ! inactivated state 1 (spark-facing)
      double precision xi2s(Lx,Ly)  ! inactivated state 2 (spark-facing)
      double precision pos(Lx,Ly)   ! open probability (spark-facing)

      ! ============================================================
      ! CALCIUM SPARK VARIABLES
      ! ============================================================
      ! Sparks are localized Ca release events from the SR
      double precision pi(Lx,Ly)    ! fraction of interior clusters sparking
      double precision pb(Lx,Ly)    ! fraction of boundary clusters sparking
      double precision pox(Lx,Ly)   ! total LCC open probability
      double precision ra(Lx,Ly)    ! spark recruitment variable

      ! Spark cluster counters (stochastic)
      integer nsb(Lx,Ly)   ! number of boundary clusters with sparks
      integer nsi(Lx,Ly)   ! number of interior clusters with sparks
	
      ! ============================================================
      ! ACTION POTENTIAL ANALYSIS ARRAYS
      ! ============================================================
      integer nbeat(0:Lx,0:Ly), nup(0:Lx,0:Ly)
      double precision apa(0:Lx,0:Ly,0:1000)     ! action potential amplitude
      double precision cimax(0:Lx,0:Ly,0:1000)  ! maximum Ca transient
      double precision yap(0:Lx,0:Ly,0:3)       ! AP analysis helper
      double precision yci(0:Lx,0:Ly,0:3)       ! Ca analysis helper
      double precision xamp(0:Lx+1,0:Ly+1)      ! amplitude array
      double precision xint(100,Ly)             ! integration array
      double precision ctmax(1000), peak(300)   ! peak detection
      
      ! Random number generation for stochastic processes
      integer, allocatable :: seed_array(:)
      integer :: i, seed_size, n
      
      ! Action Potential Duration (APD) calculation
      double precision apd(0:Lx,0:Ly,0:1500)  ! APD for each beat
      double precision yapd(0:Lx,0:Ly,0:2)    ! voltage history for APD

c     ================================================================
c     FILE OUTPUT SETUP
c     ================================================================
      ! Open files to record voltage and calcium time series
      open(unit=10,file='vx.dat',status='unknown')   ! voltage
      open(unit=11,file='cbx.dat',status='unknown')  ! boundary Ca
      open(unit=12,file='cix.dat',status='unknown')  ! interior Ca

c     ================================================================
c     PHYSIOLOGICAL CONSTANTS
c     ================================================================

      ! Ion concentrations in millimolar (mM)
      xnao = 136.0d0      ! extracellular sodium concentration
      xki = 140.0d0       ! intracellular potassium concentration  
      xko = 5.40d0        ! extracellular potassium concentration
      cao = 1.8d0         ! extracellular calcium concentration
	
      ! Fundamental physical constants
      temp = 308.0d0      ! temperature in Kelvin (35°C)
      xxr = 8.314d0       ! universal gas constant (J/mol/K)
      xf = 96.485d0       ! Faraday's constant (C/mol)
      frt = xf/(xxr*temp) ! F/RT for Nernst equation calculations

      ! Initialize random number generator
      call seed(1922)

c     ================================================================
c     CURRENT SCALING FACTORS
c     ================================================================
      ! These allow adjustment of current magnitudes for pathological conditions
      gica = 1.0          ! L-type Ca current scaling (1.0 = normal)
      gnaca = 1.0         ! Na-Ca exchange current scaling

c     ================================================================
c     CALCIUM SPARK PARAMETERS
c     ================================================================
      pbx = 0.22          ! threshold for spark-induced spark propagation
      aleak = 0.001       ! baseline SR leak rate
      qq = 0.30           ! spark sensitivity parameter

c     ================================================================
c     PACING PROTOCOL SETUP
c     ================================================================
      rbcl1 = 450.0d0     ! basic cycle length S1 (ms)
      rbcli = 200.0       ! minimum cycle length for restitution
      rbclf = 520.0       ! maximum cycle length for restitution
	
      ! Initialize random seed for reproducible stochastic behavior
      call random_seed(size=n)
      allocate(seed_array(n))
      seed_array = 986565
	
c     ================================================================
c     MAIN SIMULATION LOOPS
c     ================================================================
      
      ! Loop over different restitution curves (currently just 1)
      do kz = 1, 1
         ns1 = 1    ! number of different S2 intervals to test

         ! Loop over different S2 intervals for restitution protocol
         do ju = 1, ns1
            ! Calculate S2 interval (linear interpolation)
            rbcl2 = rbcli + (rbclf-rbcli)*dfloat(ju-1)/dfloat(ns1-1)
            call random_seed(put=seed_array)

c           ========================================================
c           NUMERICAL INTEGRATION PARAMETERS
c           ========================================================
            dt = 0.1d0        ! basic time step (ms)
            mstp0 = 10        ! maximum sub-steps for adaptive integration

            ! Calculate intracellular Na based on cycle length
            ! This empirical relationship mimics Na/K pump adaptation
            xmx = -2.0/250.0
            xnai = xmx*rbcl1 + 16.0   ! intracellular Na (mM)

c           ========================================================
c           STIMULUS PROTOCOL CONFIGURATION
c           ========================================================
            ! Set up the stimulus intervals for restitution protocol
            ! S1-S1-S1-S2-S3 where S1 is basic, S2 is variable, S3 is long
            
            do i = 1, nstim-2
               rbcl(i) = rbcl1     ! S1 intervals (basic cycle length)
            enddo
            
            rbcl(nstim) = 600.0        ! final long interval (S3)
            rbcl(nstim-1) = rbcl2      ! variable S2 interval

c           ========================================================
c           SPATIAL DISCRETIZATION AND DIFFUSION
c           ========================================================
            dx = 0.015d0      ! spatial step in x-direction (cm)
            dy = 0.015d0      ! spatial step in y-direction (cm)
            Dfu = 0.001d0     ! diffusion coefficient (cm²/ms)
            
            ! Calculate diffusion stability parameters
            ! These ensure numerical stability of the diffusion equation
       
       
            slmbdax = Dfu*dt/(4.0d0*dx*dx)  ! x-direction diffusion parameter
            slmbday = Dfu*dt/(4.0d0*dy*dy)  ! y-direction diffusion parameter

c           ========================================================
c           INITIALIZE ALL CELLULAR VARIABLES
c           ========================================================
            do ix = 1, Lx
               do iy = 1, Ly
                  
                  ! ============================================
                  ! CALCIUM INITIAL CONDITIONS
                  ! ============================================
                  ! Start with low cytosolic Ca, loaded SR
                  cb(ix,iy) = 0.3d0      ! boundary cytosolic Ca (μM)
                  ci(ix,iy) = 0.3d0      ! interior cytosolic Ca (μM)

                  cxinit = 800.0         ! initial SR Ca loading (μM)
                  csrb(ix,iy) = cxinit   ! boundary SR Ca
                  csri(ix,iy) = cxinit   ! interior SR Ca  
                  cnsr(ix,iy) = cxinit   ! network SR Ca
	 
                  ! ============================================
                  ! L-TYPE CA CHANNEL STATES (boundary-facing)
                  ! ============================================
                  ! Initialize channels in closed, non-inactivated state
                  po(ix,iy) = 0.0d0      ! open probability
                  c1(ix,iy) = 0.0d0      ! closed state 1
                  c2(ix,iy) = 1.0d0      ! closed state 2 (most channels here)
                  xi1(ix,iy) = 0.0d0     ! inactivated state 1
                  xi2(ix,iy) = 0.0d0     ! inactivated state 2
                  ra(ix,iy) = 0.0d0      ! spark recruitment variable

                  ! ============================================
                  ! L-TYPE CA CHANNEL STATES (spark-facing)
                  ! ============================================
                  pos(ix,iy) = 0.0d0     ! open probability
                  c1s(ix,iy) = 0.0d0     ! closed state 1
                  c2s(ix,iy) = 0.0d0     ! closed state 2
                  xi1s(ix,iy) = 0.0d0    ! inactivated state 1
                  xi2s(ix,iy) = 0.0d0    ! inactivated state 2

                  ! Convert free calcium to total (free + buffered)
                  call total(ci(ix,iy),cit(ix,iy))  ! interior
                  call total(cb(ix,iy),cbt(ix,iy))  ! boundary

                  ! ============================================
                  ! CALCIUM SPARK CLUSTER INITIALIZATION
                  ! ============================================
                  nsb(ix,iy) = 5         ! initial sparking boundary clusters
                  nbt = 1000             ! total boundary clusters
                  nsi(ix,iy) = 5         ! initial sparking interior clusters  
                  nit = 5000             ! total interior clusters
	
                  ! ============================================
                  ! MEMBRANE VOLTAGE AND ION CHANNEL GATES
                  ! ============================================
                  v(ix,iy) = -90.0d0     ! resting potential (mV)
                  
                  ! Sodium channel gates (closed at rest)
                  xm(ix,iy) = 0.001d0    ! activation gate (mostly closed)
                  xh(ix,iy) = 1.0d0      ! fast inactivation (not inactivated)
                  xj(ix,iy) = 1.0d0      ! slow inactivation (not inactivated)
                  
                  ! Potassium channel gates
                  xr(ix,iy) = 0.0d0      ! IKr activation (closed)
                  xs1(ix,iy) = 0.3d0     ! IKs activation
                  qks(ix,iy) = 0.2d0     ! IKs calcium dependence
                  xkur(ix,iy) = 0.01d0   ! IKur activation (mostly closed)
                  ykur(ix,iy) = 1.0d0    ! IKur inactivation (not inactivated)
                  xtof(ix,iy) = 0.02d0   ! Ito activation (mostly closed)
                  ytof(ix,iy) = 0.8d0    ! Ito inactivation

                  ! Override boundary Ca to lower value
                  cb(ix,iy) = 0.1d0
                  vold(ix,iy) = v(ix,iy)
               enddo
            enddo

c           ========================================================
c           CELLULAR VOLUME RATIOS
c           ========================================================
            ! These define the relative sizes of cellular compartments
            vi = 1.0        ! interior cytoplasm volume (normalized)
            vb = 0.3        ! boundary cytoplasm volume
            vbi = vb/vi     ! boundary to interior volume ratio

            vq = 30.0       ! SR volume scaling factor
            visr = vq       ! interior SR volume factor
            vsrin = (1.0/visr)*vi   ! interior SR to cytoplasm ratio

            vbsr = vq       ! boundary SR volume factor  
            vsrbound = (1.0/vbsr)*vb  ! boundary SR to cytoplasm ratio

            vbisr = (vsrbound/vsrin)  ! boundary to interior SR ratio

            vnsr = vq       ! network SR volume factor
            vnsrin = (1.0/vnsr)*vi    ! network SR to cytoplasm ratio

c           ========================================================
c           MAIN TIME INTEGRATION LOOP
c           ========================================================
            t = 0.0    ! total elapsed simulation time (ms)

            ! Loop over each stimulus beat
            do iz = 1, nstim
               nstep = int(rbcl(iz)/dt)  ! time steps in this beat
               
               ! Initialize APD detection variables for this beat
               yapd(ix,iy,1) = 0.0d0
               yapd(ix,iy,2) = 0.0d0

               ! Time integration within each beat
               do ncount = 0, nstep
                  time = dfloat(ncount)*dt  ! time within current beat

                  ! Loop over all spatial grid points
                  do iy = 1, Ly
                     do ix = 1, Lx

c                       ========================================
c                       CALCIUM DYNAMICS CALCULATIONS
c                       ========================================

                        ! Convert total calcium to free calcium concentrations
                        ! This accounts for calcium buffering by proteins
                        call xfree(cit(ix,iy),ci(ix,iy))  ! interior
                        call xfree(cbt(ix,iy),cb(ix,iy))  ! boundary

                        ! Calculate fraction of spark clusters that are active
                        pi(ix,iy) = dfloat(nsi(ix,iy))/dfloat(nit)  ! interior
                        pb(ix,iy) = dfloat(nsb(ix,iy))/dfloat(nbt)  ! boundary

                        ! ====================================
                        ! SR CALCIUM UPTAKE (SERCA PUMP)
                        ! ====================================
                        vupb = 0.1   ! boundary uptake rate
                        vupi = 0.1   ! interior uptake rate
                        call uptake(cb(ix,iy),vupb,xupb)  ! boundary uptake
                        call uptake(ci(ix,iy),vupi,xupi)  ! interior uptake

                        ! ====================================
                        ! Na-Ca EXCHANGE CURRENT
                        ! ====================================
                        ! Calculate Na-Ca exchange at boundary and interior
                        call inaca(v(ix,iy),frt,xnai,xnao,cao,cb(ix,iy),
     +                            xinacaq1)
                        call inaca(v(ix,iy),frt,xnai,xnao,cao,ci(ix,iy),
     +                            xinacaq2)

                        ! Weighted average based on boundary vs interior
                        ru = 1.0  ! weighting factor (1.0 = all boundary)
                        xinacaq = ru*xinacaq1 + (1.0-ru)*xinacaq2
                        xinacaq = gnaca*xinacaq  ! apply scaling factor

                        ! ====================================
                        ! L-TYPE CALCIUM CURRENT
                        ! ====================================
                        pox(ix,iy) = po(ix,iy) + pos(ix,iy)  ! total open prob
                        call ica(v(ix,iy),frt,cao,cb(ix,iy),pox(ix,iy),
     +                           rca,xicaq)
                        xicaq = gica*140.0*xicaq  ! scale and convert units

c                       ========================================
c                       CALCIUM SPARK DYNAMICS - BOUNDARY
c                       ========================================

                        ! Calculate spark initiation rate at boundary
                        ab = 50.0*qq    ! base spark rate parameter
                        csrx = 850.0d0  ! SR Ca threshold for sparks
                        
                        ! SR Ca dependence (sigmoid with steep cutoff)
                        phisr = 1.0/(1.0+(csrx/csrb(ix,iy))**10)

                        ! Total spark initiation rate (LCC-triggered)
                        alphab = ab*dabs(rca)*po(ix,iy)*phisr
                        bts = 1.0/30.0   ! spark termination rate (1/30 ms⁻¹)

                        ! Update L-type Ca channel Markov states
                        call markov(dt,v(ix,iy),cb(ix,iy),c1(ix,iy),
     +                    c2(ix,iy),xi1(ix,iy),xi2(ix,iy),po(ix,iy),
     +                    c1s(ix,iy),c2s(ix,iy),
     +                    xi1s(ix,iy), xi2s(ix,iy),pos(ix,iy),
     +                    alphab,bts)

                        ! Calculate RyR Ca release from boundary SR
                        gsrb = 0.01d0/2.5  ! boundary RyR conductance
                        xryrb = gsrb*csrb(ix,iy)*pb(ix,iy)

c                       ========================================
c                       CALCIUM SPARK DYNAMICS - INTERIOR
c                       ========================================

                        ! Interior spark rate calculation
                        csrxx = 900.0d0  ! SR Ca threshold
                        chx = 850.0d0    ! Ca sensitivity parameter
                        
                        ! SR Ca dependence for interior sparks
                        phi = (csri(ix,iy)/chx)/
     +                       (1.0+(csrxx/csri(ix,iy))**10)
	
                        ! Spark recruitment and propagation mechanism
                        gi = 0.25   ! recruitment strength
                        pra = 0.05  ! recruitment threshold
                        xra = (pra/pi(ix,iy))**6  ! recruitment function
                        ra_inf = phi/(1.0+xra)    ! steady-state recruitment

                        dra = (ra_inf-ra(ix,iy))/60.0  ! recruitment dynamics
	
                        ! Three sources of interior spark initiation:
                        rr = 0.02*3.0*1.5  ! spark propagation rate
                        pbinf = 1.0/(1.0+(pbx/pb(ix,iy))**10)  ! boundary influence

                        ar1 = (aleak/500.0)*(csri(ix,iy)/1000.0)  ! SR leak
                        ar2 = rr*pbinf*phi                        ! propagation  
                        ar3 = gi*ra(ix,iy)                        ! recruitment

                        alphai = ar1 + ar2 + ar3  ! total initiation rate

                        ! Interior RyR calcium release
                        gryri = 0.015  ! interior RyR conductance
                        xryri = gryri*pi(ix,iy)*csri(ix,iy)
                        btsi = 1.0/70.0  ! interior spark termination rate

c                       ========================================
c                       CALCIUM BALANCE EQUATIONS
c                       ========================================

                        ! Net sarcolemmal Ca flux (ICa - INCX)
                        xsarc = -xicaq + xinacaq
	
                        ! SR-cytoplasm Ca diffusion
                        tausri = 50.0  ! time constant for SR equilibration
                        dff = (cnsr(ix,iy)-csri(ix,iy))/tausri

                        ! Intracellular Ca diffusion time constants
                        tau1 = 10.0  ! cytoplasm boundary-interior diffusion
                        tau2 = 10.0  ! SR boundary-NSR diffusion

                        dfbi = 1.0*(cb(ix,iy)-ci(ix,iy))/tau1      ! cyto diff
                        dfbisr = 1.0*(csrb(ix,iy)-cnsr(ix,iy))/tau2  ! SR diff

                        ! Time derivatives for calcium concentrations
                        dcbt = xryrb - xupb + xsarc - dfbi           ! boundary cyto
                        dcsrb = vbsr*(-xryrb+xupb) - dfbisr         ! boundary SR

                        dcit = xryri - xupi + vbi*dfbi              ! interior cyto
                        dcsri = visr*(-xryri) + dff                 ! interior SR
	
                        dnsr = vnsr*(xupi) - dff + vbisr*dfbisr     ! network SR

                        ! Update calcium concentrations
                        cbt(ix,iy) = cbt(ix,iy) + dcbt*dt
                        cit(ix,iy) = cit(ix,iy) + dcit*dt
                        csrb(ix,iy) = csrb(ix,iy) + dcsrb*dt	
                        csri(ix,iy) = csri(ix,iy) + dcsri*dt
                        cnsr(ix,iy) = cnsr(ix,iy) + dnsr*dt
                        ra(ix,iy) = ra(ix,iy) + dra*dt

c                       ========================================
c                       STOCHASTIC SPARK EVOLUTION
c                       ========================================

                        ! Update boundary spark clusters using binomial statistics
                        nsbx = nsb(ix,iy)
                        call binevol(nbt,nsbx,alphab,bts,dt,
     +                              ndeltapx,ndeltamx)

                        ! Ensure spark numbers stay within bounds
                        if(ndeltamx.gt.nsbx.or.ndeltapx.gt.nbt) then
                           nsb(ix,iy) = 0
                        else
                           nsb(ix,iy) = nsb(ix,iy) + ndeltapx - ndeltamx
                        endif 

                        ! Update interior spark clusters
                        nsix = nsi(ix,iy)
                        call binevol(nit,nsix,alphai,btsi,dt,
     +                              ndeltapy,ndeltamy)
	
                        if(ndeltamy.gt.nsix.or.ndeltapy.gt.nit) then
                           nsi(ix,iy) = 0
                        else
                           nsi(ix,iy) = nsi(ix,iy) + ndeltapy - ndeltamy
                        endif 

c                       ========================================
c                       MEMBRANE VOLTAGE DYNAMICS
c                       ========================================

                        ! Convert calcium fluxes to membrane currents
                        wca = 12.0d0                ! conversion factor
                        xinaca = wca*xinacaq        ! Na-Ca exchange current
                        xica = 2.0d0*wca*xicaq      ! L-type Ca current

                        ! Adaptive time stepping based on voltage derivative
                        adq = dabs(dv(ix,iy))
                
                        if(adq.gt.25.0d0) then
                           mstp = 10    ! use fine time steps for fast changes
                        else
                           mstp = 1     ! normal time step
                        endif 
	              
                        hode = dt/dfloat(mstp)  ! actual integration time step

                        ! Sub-stepping loop for voltage integration
                        do iii = 1, mstp  

                           ! ================================
                           ! ION CHANNEL CURRENT CALCULATIONS
                           ! ================================
                           
                           ! Fast sodium current (action potential upstroke)
                           call ina(hode,v(ix,iy),frt,xh(ix,iy),
     &                        xj(ix,iy),xm(ix,iy),xnai,xnao,xina)

                           ! Rapid delayed rectifier K current (repolarization)
                           call ikr(hode,v(ix,iy),frt,xko,xki,
     &                        xr(ix,iy),xikr)

                           ! Slow delayed rectifier K current (late repolarization)
                           call iks(hode,v(ix,iy),frt,cb(ix,iy),xnao,
     &                        xnai,xko,xki,xs1(ix,iy),qks(ix,iy),xiks)

                           ! Inward rectifier K current (maintains resting potential)
                           call ik1(hode,v(ix,iy),frt,xki,xko,xik1)

                           ! Ultra-rapid delayed rectifier K current
                           call ikur(hode,v(ix,iy),frt,xki,xko,
     &                        xkur(ix,iy),ykur(ix,iy),xikur)

                           ! Transient outward K current (early repolarization)
                           call ito(hode,v(ix,iy),frt,xki,xko,
     &                        xtof(ix,iy),ytof(ix,iy),xito)

                           ! Na-K pump current
                           call inak(v(ix,iy),frt,xko,xnao,xnai,xinak)

                           ! ================================
                           ! STIMULUS CURRENT
                           ! ================================
                           if(time.lt.1.0) then
                              stim = 80.0   ! stimulus amplitude (μA/μF)
                           else
                              stim = 0.0    ! no stimulus
                           endif       

                           ! ================================
                           ! TOTAL MEMBRANE CURRENT
                           ! ================================
                           ! Sum all K currents with scaling factors
                           xkzz = xik1 + 0.3*xikr + 0.3*xiks + 
     +                           1.5*xito + 0.3*xikur
                           
                           ! Total current (positive = depolarizing)
                           dvh = -(xina + xkzz + xinaca + xica + xinak) 
     +                          + stim

                           ! Update membrane voltage
                           v(ix,iy) = v(ix,iy) + dvh*hode
                        enddo 
                     enddo  ! end ix loop
                  enddo     ! end iy loop

                  ! ============================================
                  ! SPATIAL DIFFUSION OF VOLTAGE
                  ! ============================================
                  ! Apply diffusion operator twice for numerical stability
                  ! This models electrical coupling between neighboring cells
                  call Euler_Forward(Lx,Ly,v,vnew,slmbdax,slmbday)
                  call Euler_Forward(Lx,Ly,v,vnew,slmbdax,slmbday)

c                 ================================================
c                 ACTION POTENTIAL DURATION (APD) CALCULATION
c                 ================================================
                  ! APD is measured as time from upstroke to repolarization
                  ! to a specific voltage threshold

                  ctoti = 0.0  ! unused variable

                  do iy = 1, Ly
                     do ix = 1, Lx
                        ! Calculate voltage time derivative
                        dv(ix,iy) = (vnew(ix,iy)-vold(ix,iy))/dt
                        vold(ix,iy) = vnew(ix,iy)

                        ! APD measurement using voltage threshold crossing
                        vc = -40.0  ! repolarization threshold (mV)
                        yapd(ix,iy,1) = yapd(ix,iy,2)  ! store previous voltage
                        yapd(ix,iy,2) = v(ix,iy)       ! current voltage
                
                        ! Detect downward crossing of threshold (repolarization)
                        if(v(ix,iy).le.vc.and.yapd(ix,iy,1).gt.vc) then            
                           ! Linear interpolation to find exact crossing time
                           timex = time - dt*(vc-v(ix,iy))/
     +                            (yapd(ix,iy,1)-v(ix,iy))
                           apd(ix,iy,iz) = timex  ! store APD for this beat
                        endif 
                     enddo
                  enddo

c                 ================================================
c                 DATA OUTPUT AND MONITORING
c                 ================================================
                  ! Monitor total spark activity
                  xxns = nsb(2,2) + nsi(2,2)  ! total sparks at center cell
	
                  ! Write time series data every 20 time steps (2 ms intervals)
                  if(mod(ncount,20).eq.0) then
                     write(10,*) t, v(2,2)      ! voltage at center cell
                     write(11,*) t, cb(2,2)     ! boundary Ca at center
                     write(12,*) t, ci(2,2)     ! interior Ca at center
                  endif 

                  t = t + dt  ! increment total simulation time
               enddo ! end time step loop within beat
            enddo ! end beat loop

            ! ====================================================
            ! POST-PROCESSING: RESTITUTION ANALYSIS
            ! ====================================================
            ! Extract APDs from the last two beats for restitution analysis
            xxapd = apd(Lx/2,Ly/2,nstim)      ! APD of final beat
            zapd = apd(Lx/2,Ly/2,nstim-1)     ! APD of penultimate beat
            
            ! Optional: output restitution data
            ! write(20,*) rbcl2-zapd, xxapd  ! (S2-APD1, APD2) data point
            
         enddo ! end rbcl loop (S2 interval variation)
      enddo ! end restitution curve loop

      stop
      end

c     ================================================================
c     SPATIAL DIFFUSION SUBROUTINE
c     ================================================================
c     
c     This subroutine implements the diffusion of membrane voltage
c     between neighboring cells, representing electrical coupling
c     through gap junctions. Uses explicit Euler method with 
c     no-flux boundary conditions.
c     
c     Parameters:
c       LLx, LLy    - grid dimensions
c       v           - voltage array (input/output)
c       vnew        - temporary voltage array
c       slmbdax/y   - diffusion stability parameters
c     ================================================================
      
      subroutine Euler_forward(LLx,LLy,v,vnew,slmbdax,slmbday)
      implicit double precision (a-h,o-z)
        
      double precision  v(0:LLx+1,0:LLy+1),vnew(0:LLx+1,0:LLy+1)
        
c     ============================================================
c     BOUNDARY CONDITIONS: NO-FLUX (SEALED EDGES)
c     ============================================================
c     Set ghost points to enforce zero gradient at boundaries
c     This represents an electrically isolated tissue preparation
      
      ! Corner points (set to nearby interior values)
      v(0,0) = v(2,2)
      v(0,LLy+1) = v(2,LLy-1)
      v(LLx+1,0) = v(LLx-1,2)
      v(LLx+1,LLy+1) = v(LLx-1,LLy-1)

      ! Top and bottom edges (zero gradient in y-direction)
      do ix = 1, LLx
         v(ix,0) = v(ix,2)           ! bottom edge
         v(ix,LLy+1) = v(ix,LLy-1)   ! top edge
      enddo
  
      ! Left and right edges (zero gradient in x-direction)
      do iy = 1, LLy
         v(0,iy) = v(2,iy)           ! left edge
         v(LLx+1,iy) = v(LLx-1,iy)   ! right edge
      enddo
 
c     ============================================================
c     EXPLICIT EULER DIFFUSION STEP
c     ============================================================
c     Solve: ∂V/∂t = D∇²V using finite differences
c     ∇²V ≈ (V[i+1,j] + V[i-1,j] - 2V[i,j])/Δx² + 
c           (V[i,j+1] + V[i,j-1] - 2V[i,j])/Δy²
      
      do ix = 1, LLx
         do iy = 1, LLy
            vnew(ix,iy) = v(ix,iy) + 
     +        slmbdax*(v(ix+1,iy) + v(ix-1,iy) - 2.0d0*v(ix,iy)) +
     +        slmbday*(v(ix,iy+1) + v(ix,iy-1) - 2.0d0*v(ix,iy))
         enddo
      enddo

c     ============================================================
c     UPDATE BOUNDARY CONDITIONS FOR NEW VOLTAGE ARRAY
c     ============================================================
      
      ! Update corner points
      vnew(0,0) = vnew(2,2)
      vnew(0,LLy+1) = vnew(2,LLy-1)
      vnew(LLx+1,0) = vnew(LLx-1,2)
      vnew(LLx+1,LLy+1) = vnew(LLx-1,LLy-1)

      ! Update edges
      do ix = 1, LLx
         vnew(ix,0) = vnew(ix,2)
         vnew(ix,LLy+1) = vnew(ix,LLy-1)
      enddo

      do iy = 1, LLy
         vnew(0,iy) = vnew(2,iy)
         vnew(LLx+1,iy) = vnew(LLx-1,iy)
      enddo

c     ============================================================
c     SECOND DIFFUSION STEP FOR IMPROVED ACCURACY
c     ============================================================
c     Apply diffusion operator again using updated voltage values
c     This helps maintain numerical stability and accuracy
      
      do ix = 1, LLx
         do iy = 1, LLy
            v(ix,iy) = vnew(ix,iy) + 
     +     slmbdax*(vnew(ix+1,iy) + vnew(ix-1,iy) - 2.0d0*vnew(ix,iy)) +
     +     slmbday*(vnew(ix,iy+1) + vnew(ix,iy-1) - 2.0d0*vnew(ix,iy))
         enddo
      enddo

      return
      end

c     ================================================================
c     END OF CARDIAC ELECTROPHYSIOLOGY SIMULATION
c     ================================================================
c
c     SUMMARY OF WHAT THIS CODE DOES:
c     
c     1. INITIALIZATION:
c        - Sets up a 2D grid of cardiac cells (10x10)
c        - Initializes all cellular variables (voltage, Ca, ion channels)
c        - Configures physiological parameters and constants
c     
c     2. STIMULUS PROTOCOL:
c        - Applies a series of electrical stimuli (S1-S1-S2-S3)
c        - Can vary S2 timing for restitution studies
c     
c     3. CELLULAR DYNAMICS (for each cell, each time step):
c        - Calcium handling: uptake, release, diffusion
c        - Stochastic calcium sparks: initiation and termination
c        - Ion channel gating: Na, K, Ca channels
c        - Membrane voltage: integration of all currents
c     
c     4. SPATIAL COUPLING:
c        - Electrical diffusion between neighboring cells
c        - Models gap junction coupling
c     
c     5. MEASUREMENT AND OUTPUT:
c        - Records voltage and calcium time series
c        - Calculates action potential duration (APD)
c        - Can generate restitution curves
c
c     This is a comprehensive model suitable for studying:
c     - Cardiac arrhythmias
c     - Calcium handling disorders  
c     - Drug effects on ion channels
c     - Electrical propagation in tissue
c     - Action potential restitution properties
c     ================================================================
