c
ccc	To compile: "gfortran -ffixed-line-length-132 lepton_event_v17_4.f -o lepton_event_v17_4.exe”
c	To run: “./lepton_event_v17_4.exe”
c
c	version 6 of code weights the distribution of events by xs= acos( costheta ).  This is the most efficient version
c	of the code to run.  The nuclear and atomic form factors are included. 
c	version 8 allows for normal event weighting, or non-standard event weighting, theta=x**2, switchable with ‘standard_phase_space’. 
c		Non-standard is the most efficient. 
c	version 9: histograms J_T phi angle, uses simplified version of W_pol. 
c	version 10: histograms cross section arrays
c	version 11: histogram integrated cross section as a function of x and phi1. 
c	version 12: introduced logicals to do electron, muon, and muon with pion mass hypothesis histogramming. 
c	version 13: added some histogram options and controls
c	version 14: added some more controls, got rid of things in the code that weren’t useful 
c	version 15: brought over some features from the pi+pi- code, fixed some bugs.   Can control shape of the proton form factor. 
c	version 16: consistently use transverse momentum transfer squared everywhere, changed name to lepton_event
c	Uses rejection sampling: 
c	version 17: added feature to run one kinematic point
c	version 17_1: use weighting dcos theta/dx = (1-cos theta)**N , got rid of the option for bending the FF
c	version 17_2: use weighting theta = x**phase_space, with phase_space an integer >1 .   Set phase_space = 0 for dcos theta/dx =1. 
c		      makes a guess for the largest cross section, which seems to be correct
c	version 17_3: add option for log t plot, control proton rms radius independent of dipole FF parameter
c	version 17_4: read brem. file
c
	implicit none
	real*8 E0,theta_min,theta_max
	real*8 pi,cos_max,cos_min,cross_sum,cross_max,cross_test,x_min,x_max,theta1,phi1,costheta1,theta2,phi2,
     &	costheta2,cross_section,total_xscn,xsctn,k1(3),k2(3),m_e,m_part,m_muon,pol,x,q2,failure,
     &	w_mumu,xs1,xs2,xs_max,xs_min,
     &	Atgt(100),tgtlen(100),radlen(100),c_den(100),a_den(100),rho0,hbarc,W,jacobian,
     &	m_pi,phi_JT,delta_w,delta_x,delta_phi,mtgt,ktgt(3),Elab(2),missing_mass,t,delta_t,x_value,
     &	q2_T,W_pol,W_unpol,JS,JT(2),Rexp,xmax,theta1_max,theta2_max,phi1_max,phi2_max,temp,delta_log_t,frac_delta_t,delta_Egamma,
     &	data_array(0:200),error_array(0:200),Egamma_max,E_hi,E_lo,Brem,E_coherent,Egamma
	integer*4 iseed
	integer*8 i,itest(4),nevent,j,nfail,bad_max,i_array,j_array, ztgt,phase_space
	logical hist_w,hist_x,hist_t,hist_phi_JT,hist_Egamma,output_event,hist_log_t,nuc_FF,muon,electron,pion_hypothesis,
     &	brem_init,cobrems
c
c	Standard CPP configuration 
c	data ztgt,E0,pol,theta_min,theta_max /82,5.5,1.0,0.80, 5.3 /   !Standard CPP configuration, min angle to TOF and max angle to MWPC
c
c	Standard GlueX configuration
c	data ztgt,E0,E_coherent,pol,theta_min,theta_max /1,11.0,8.7,1.0,0.90,13.12/	!standard GlueX config., min and max angles in deg. to TOF
	data ztgt,E0,E_coherent,pol,theta_min,theta_max /1,11.0,8.7,1.0,0.090,13.12/	!standard GlueX config., min and max angles in deg. to TOF
c
c	Set tagging interval
	data E_hi,E_lo /8.8,8.6/
c
	data itest(1),itest(2),itest(3),itest(4),nevent 
     &	/100000,1000000,10000000,100000000,10/
	data m_e,m_muon,m_pi,hbarc /0.000511,0.105658,0.139570,0.197/
c	
c  Histogram parameters
	data delta_w,delta_t,delta_x,delta_phi,frac_delta_t,delta_Egamma /0.02,.0002,0.02,5.0,.2,0.05/	!units of GeV, GeV^2, dimensionless, degrees, 
c											fractional bin width in t, GeV
c
c  Target information
	data Atgt(1), tgtlen(1), radlen(1) /1., .0338,63.04/	!tgtlen = # Rad lengths, radlen= rad length of material in g/cm^2
c								this target is based on 30 cm LH2
	data Atgt(6), tgtlen(6), radlen(6), c_den(6), a_den(6)  /12.,  .05, 42.70, 2.45,  .524/	!RL is in g/cm^2
	data Atgt(14),tgtlen(14),radlen(14),c_den(14),a_den(14) /28., .05, 21.82, 3.14,  .537/ !units of c_den and a_den are fm
	data Atgt(20),tgtlen(20),radlen(20),c_den(20),a_den(20) /40., .05, 16.12, 3.51,  .563/ !target length is in units of RL
	data Atgt(26),tgtlen(26),radlen(26),c_den(26),a_den(26) /56., .05, 13.84, 3.971, .5935/ !RL is g/cm^2
	data Atgt(50),tgtlen(50),radlen(50),c_den(50),a_den(50) /116.,.05, 8.82,  5.416, .552/
	data Atgt(82),tgtlen(82),radlen(82) /208.,.05,6.37/  
c
	common /density_rho0/ rho0,c_den,a_den		!the overall normalization factor for nuclear charge densities
c
	double precision ZBQLUAB,zlo,zhi
	data zlo,zhi /0.,1./
c
c      BLOCK DATA ZBQLBD01
*
*       Initializes seed array etc. for random number generator.
*       The values below have themselves been generated using the
*       NAG generator.
*
      DOUBLE PRECISION ZBQLIX(43),B,C
      COMMON /ZBQL0001/ ZBQLIX,B,C
c      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
     +1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     +7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     +2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     +4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     +2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     +1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     +3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     +2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     +3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     +2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     +2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
c
	iseed=0		!set equal to 0 for the program to do a call to the system clock to initialize random number generator
c			!for random number initialization
	call ZBQLINI(ISEED)
c
c  Initializations
c
	temp=-1.
	pi = dacos(temp)
c
c	Find the delta log t step
	delta_log_t=log10(1.+frac_delta_t)

c	Set target mass in GeV
	mtgt=Atgt(ztgt)*.931494
	if (ztgt.eq.1) mtgt=0.93828
c
	do i=0,200
		data_array(i)=0.
		error_array(i)=0.
	end do
c
c	Initialize Brem. distribution: select 1/Egamma or coherent Brems. file
	cobrems=.true.	!set true for reading coherent Brems. file, false for using a 1/Egamma distribution
	if (cobrems) then ! read coherent Brem. file
		brem_init=.true.
		temp=Brem(brem_init,cobrems,E0,Egamma)	!read coherent brems file, then set brem_init = .false.
	endif
c
c	Start logical assignments
c
c   Only one of the histogram logicals can be true, the rest must be false.   The event output can be turned on independent of histograming. 
	hist_w		=.false.
	hist_x		=.false.
	hist_t		=.false.
	hist_phi_JT	=.false.
	hist_log_t	=.false.
	hist_Egamma	=.false.
c
	output_event	=.true.
c
c   Logical assignments
c
	phase_space= 4		!theta=x**phase_space, with integer phase_space >1 . Note: phase_space=4 seems to be fastest. 
c				 !Set phase_space=0 for standard dcos theta/dx =1
	Rexp=real(phase_space)
	nuc_FF=.true.
	muon=.false.	!this is how you change the particle type	
	electron=.not.muon
c
	if (muon) then
		m_part=m_muon
		pion_hypothesis=.false.    !set true if you want calculated invariant masses with pion assumption
	else 
		m_part=m_e
		pion_hypothesis=.false.	!should always have this set false since we can distinguish e from mu
	end if
c
	if (nuc_FF) then	! setup the nuclear form factors
		c_den(ztgt)=c_den(ztgt)/hbarc
		a_den(ztgt)=a_den(ztgt)/hbarc
		call density_init(ztgt) !for initializing the nuclear form factor (rho0 is density const.)
	end if
c
	if (hist_w.or.hist_x.or.hist_t.or.hist_phi_JT.or.hist_log_t.or.hist_Egamma) 
     &	open(unit=2,status='new',file= 'lepton_v17_4_hist.txt')
	if (output_event) 
     &	open(unit=4,status='new',file= 'lepton_v17_4_event.txt')
c
c	End logical assignments
c
c***********************************
c	Evaluate the cross section at one kinematic point at coherent peak
	x=0.5
	theta1=  theta_min*(pi/180.)	
	theta2=  theta_min*(pi/180.)
	phi1=   90.*(pi/180.)
	phi2=   270.*(pi/180.)	
	cross_section=xsctn(E_coherent,ztgt,x,theta1,phi1,theta2,phi2,pol,m_part,nuc_FF,phi_JT,k1,k2)
	print *, ' cross section nb/sr^2 = ', cross_section	
c*************************************
	theta_min=theta_min*pi/180.	!switch to radians
	theta_max=theta_max*pi/180.
	cos_max=cos(theta_min)
	cos_min=cos(theta_max)
c Limits on xs
	if (phase_space.eq.0) then 	 !dcos theta/dx =1
		xs_max=cos_max	
		xs_min=cos_min
	else
		xs_max=theta_max**(1./Rexp)
		xs_min=theta_min**(1./Rexp)
	end if
c
c***************************************************************************
c
	do j=1,4	!loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the maximum 
			!cross section*Brem converges
	cross_max=0.
	do i=1,itest(j)	!find maximum cross section in allowed phase space
		Egamma=ZBQLUAB(E_lo,E_hi)	!get tagged photon energy
		x=0.5		!x_min+(x_max-x_min)*ZBQLUAB(zlo,zhi)	!make a guess for the energy fraction
		phi1=90.*pi/180.	!2.*pi*ZBQLUAB(zlo,zhi)		!make a guess for phi1
		phi2=270.*pi/180.	!2.*pi*ZBQLUAB(zlo,zhi)		!make a guess for phi2
		xs2=ZBQLUAB(xs_min,xs_max)	 
		if (phase_space.eq.0) then	! dcos theta/dx = 1
			theta1=theta_min	!make a guess for theta1
			theta2=acos(xs2)
			jacobian=1.
		else
			theta1=theta_min	!make a guess for theta1
			xs1=theta_min**(1./Rexp)
			theta2=xs2**phase_space
			jacobian=(Rexp*xs1**(phase_space-1)*sin(xs1**phase_space))*(Rexp*xs2**(phase_space-1)*sin(xs2**phase_space))
		end if
		cross_section=xsctn(Egamma,ztgt,x,theta1,phi1,theta2,phi2,pol,m_part,nuc_FF,phi_JT,k1,k2)*jacobian*
     &		Brem(brem_init,cobrems,E0,Egamma)
		if (cross_section.gt.cross_max)  then 
			cross_max=cross_section
			Egamma_max=Egamma
			xmax=x
			theta1_max=theta1*180./pi
			theta2_max=theta2*180./pi
			phi1_max=phi1*180./pi
			phi2_max=phi2*180./pi
		end if
	end do
	print *, 'test events', itest(j), 'maximum xsctn*Brem' , cross_max
	print *, 'Egamma max', Egamma_max, 'x max', xmax, 'theta1 max', theta1_max, 'theta2 max', theta2_max,
     &	'phi1 max', phi1_max, 'phi2 max', phi2_max 
	end do
c
c**********************************************************************************************************
c Loop over 4 samplings of the phase space at coherent peak, each a factor of x10 larger, to see if the integrated cross section converges
c Set limits on energy fraction
	x_max=(E_coherent-m_part)/E_coherent
	x_min=m_part/E_coherent
c
	do j=1,4	!loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the integrated
			!cross section at the coherent peak converges
	cross_sum=0.
	do i=1,itest(j)	
		x=x_min+(x_max-x_min)*ZBQLUAB(zlo,zhi)	!energy fraction
		phi1=2.*pi*ZBQLUAB(zlo,zhi)
		phi2=2.*pi*ZBQLUAB(zlo,zhi)
		xs1=ZBQLUAB(xs_min,xs_max)	
		xs2=ZBQLUAB(xs_min,xs_max)	 
		if (phase_space.eq.0) then	! dcos theta/dx = 1
			theta1=acos(xs1)
			theta2=acos(xs2)
			jacobian=1.
		else
			theta1=xs1**phase_space
			theta2=xs2**phase_space
			jacobian=(Rexp*xs1**(phase_space-1)*sin(xs1**phase_space))*(Rexp*xs2**(phase_space-1)*sin(xs2**phase_space))
		end if
		cross_section=xsctn(E_coherent,ztgt,x,theta1,phi1,theta2,phi2,pol,m_part,nuc_FF,phi_JT,k1,k2)*jacobian
		cross_sum=cross_sum+cross_section
	end do
	total_xscn=cross_sum/float(itest(j))*(xs_max-xs_min)**2*(2.*pi)**2*(x_max-x_min)
	print *, 'test events', itest(j), 'Egamma', E_coherent, 'total cross section nb', total_xscn
	end do
c*************************************************************************
c
c	Start event generation
c
c	Use the widest possible range in x by using the maximum accepted tagged photon energy, then test it
	x_max=(E_hi-m_part)/E_hi	!largest possible x
	x_min=m_part/E_hi		!smallest possible x
c
	nfail=0
	bad_max=0
c
	do i=1,nevent
100		continue
		Egamma=ZBQLUAB(E_lo,E_hi)	!get tagged photon energy
		x=ZBQLUAB(x_min,x_max)		!energy fraction	
c	Test x to make sure it's within the allowed range for the photon energy Egamma
		if (x.ge.((Egamma-m_part)/Egamma).or.x.le.(m_part/Egamma)) go to 100 ! x is out of range, try again
c
		phi1=2.*pi*ZBQLUAB(zlo,zhi)
		phi2=2.*pi*ZBQLUAB(zlo,zhi)
		xs1=ZBQLUAB(xs_min,xs_max)
		xs2=ZBQLUAB(xs_min,xs_max)
c
		if (phase_space.eq.0) then	! dcos theta/dx = 1
			theta1=acos(xs1)
			theta2=acos(xs2)
			jacobian=1.
		else 
			theta1=xs1**phase_space
			theta2=xs2**phase_space
			jacobian=(Rexp*xs1**(phase_space-1)*sin(xs1**phase_space))*(Rexp*xs2**(phase_space-1)*sin(xs2**phase_space))
		end if
c
		cross_section=xsctn(Egamma,ztgt,x,theta1,phi1,theta2,phi2,pol,m_part,nuc_FF,phi_JT,k1,k2)*jacobian*
     &		Brem(brem_init,cobrems,E0,Egamma)
		if (cross_section.gt.cross_max) then
			bad_max=bad_max+1	!an occurrence of cross section larger than cross_max, not supposed to happen
			print *, ' bad max cross section= ', cross_section
		end if
		cross_test=cross_max*ZBQLUAB(zlo,zhi)
		if (cross_test.gt.cross_section) then	!selection fails
			nfail=nfail+1
			go to 100
		end if
c	Event selection succeeds: 
c
		call analysis(Egamma,mtgt,k1,k2,ktgt,w_mumu,t,missing_mass,m_part,pion_hypothesis)	!analyze the event
c							
c		Do the histogramming
c
		if (hist_w) i_array=nint((w_mumu-.200)/delta_w)	!w distribution
		if (hist_x) i_array=nint(x/delta_x)		!	x distribution
		if (hist_t) i_array=nint(t/delta_t)		!	t distribution
		if (hist_phi_JT) i_array=nint(phi_JT*180./pi/delta_phi) !	JT phi distribution in degrees
		if (hist_log_t) i_array=nint((log10(t)+6.)/delta_log_t)	! 10^-6 GeV^2 is bin 0
		if (hist_Egamma) i_array=nint(Egamma/delta_Egamma)	!photon energy distribution
c
		if (i_array.lt.0) i_array=0
		if (i_array.gt.200) i_array=200		
		data_array(i_array)=data_array(i_array)+1.
c
c	3-momentum event output
		if (output_event) write(4,200) Egamma,k1(1),k1(2),k1(3),k2(1),k2(2),k2(3),ktgt(1),ktgt(2),ktgt(3)
200		format(2x,f6.3,1x,9(f10.6,1x))
c
	if (mod(i,100).eq.0) print *,' event # ', i
	end do
c
c   Event generation ends
c	
	failure=float(nfail)/float(nevent)
	print *, 'Failures per event = ', failure , ' Events with cross section exceeding max xsctn = ', bad_max
c

	do i=0,200		!write out the arrays
		if (hist_w) then 
			x_value=float(i)*delta_w+.200
			error_array(i)=sqrt(data_array(i))
			write(2,20) x_value, data_array(i), error_array(i)		
		else if (hist_x) then
			x_value=float(i)*delta_x
			error_array(i)=sqrt(data_array(i))
			write(2,20) x_value, data_array(i), error_array(i)		
		else if (hist_t) then 
			x_value=float(i)*delta_t
			error_array(i)=sqrt(data_array(i))
			write(2,20) x_value, data_array(i), error_array(i)		
		else if (hist_phi_JT) then 
			x_value=float(i)*delta_phi
			if (x_value.eq.0..or.x_value.eq.360.) data_array(i)=2.*data_array(i)  ! this is a binning problem
			error_array(i)=sqrt(data_array(i))
			write(2,20) x_value, data_array(i), error_array(i)		
		else if (hist_log_t) then
			x_value=10.**(float(i)*delta_log_t-6.)
			error_array(i) = sqrt(data_array(i))/(10.**(float(i+1)*delta_log_t-6.)-10.**(float(i-1)*delta_log_t-6.))*2. 
			data_array(i) =        data_array(i)/(10.**(float(i+1)*delta_log_t-6.)-10.**(float(i-1)*delta_log_t-6.))*2. 
			write(2,20) x_value, data_array(i), error_array(i)		
		else if (hist_Egamma) then
			x_value=float(i)*delta_Egamma
			error_array(i)=sqrt(data_array(i))
			write(2,20) x_value,data_array(i),error_array(i)
		end if			
	end do
c
20	format(2x,e10.4,2x,e10.4,2x,e10.4)
c
	stop
	end
c
c*******************************************************************
c
	real*8 function Brem(brem_init,cobrems,E0,Egamma)
	implicit none
	real*8 E0,Egamma,Eg(500),Br(500)
	integer*8 i,imax,ipoint
	logical brem_init,cobrems
	common /Brem_spect/ Eg,Br
c	
	Brem=0.
c
	if (brem_init) then  !open and read coherent Brem file
		open(unit=2,file='CobremsDistribution.dat')
		i=1
10 		continue
		read(2,*,end=20) Eg(i),Br(i) 
c		print *, Eg(i), Br(i)
		i=i+1
		go to 10
20 		continue
		imax=i-1
		close(2)
		brem_init=.false.	!done with initialization
		return
	end if
c
	if(.not.cobrems) then
		Brem=E0/Egamma
		return
	end if
c
	if (cobrems) then !return coherent brems distribution
		ipoint=nint((Egamma +.02)/.04)	
		Brem=Br(ipoint)
		return
	end if
	end
c	
c******************************************************************
c The reference for this is Bakmaev et al., Physics Letters B 660 (2008) 494-500, eqn. 23
c Had to correct a mistake in Eqn. 23 of their paper.   The cos(phi1 + phi2) term should be 
c multiplied by 2.  You can see this by comparing Wp in Eqn. 22 with the vector current part of Eqn. 23
c
	real*8 function xsctn(E0,ztgt,x,theta1,phi1,theta2,phi2,pol,m_part,nuc_FF,phi_JT,k1,k2)	!units of nb/sr^2
	implicit none
	real*8 E0,Z,x,theta1,phi1,theta2,phi2,k1(3),k2(3),pol,W_unpol,W_pol,q2_T
	real*8 m_part,alpha,hbarc
	real*8 pi,E1,E2,k1_mag,k2_mag,p1,p2,q2,c1,c2,JS,JT(2),FF2,FF_nuc,FF_TFM,phi_JT
	integer*8 i,ztgt
	logical nuc_FF
c
	data alpha,hbarc /7.297352e-3,0.197/
c
	pi=acos(-1.)
c
	Z=float(ztgt)
	E1=E0*x
	E2=E0*(1.-x)
	k1_mag=sqrt(E1**2-m_part**2)
	k2_mag=sqrt(E2**2-m_part**2)
	k1(1)=k1_mag*sin(theta1)*cos(phi1)
	k1(2)=k1_mag*sin(theta1)*sin(phi1)
	k1(3)=k1_mag*cos(theta1)
	k2(1)=k2_mag*sin(theta2)*cos(phi2)
	k2(2)=k2_mag*sin(theta2)*sin(phi2)
	k2(3)=k2_mag*cos(theta2)
c
	p1=sqrt(k1(1)**2+k1(2)**2)		!transverse momenta of muon #1, GeV
	p2=sqrt(k2(1)**2+k2(2)**2)		!transverse momenta of muon #2, GeV
	q2_T=(k1(1)+k2(1))**2+(k1(2)+k2(2))**2 ! this is transverse momentum transfer squared
	q2=q2_T+(E0-k1(3)-k2(3))**2  ! this is 3-momentum transfer squared
	c1=p1**2+m_part**2
	c2=p2**2+m_part**2
	JS=1./c1-1./c2	!units of 1/GeV^2	!scalar current, units of GeV^-2
	do i=1,2
		JT(i)=k1(i)/c1+k2(i)/c2	!vector current, units of GeV^-1
	end do
c
	phi_JT=acos(JT(1)/sqrt(JT(1)**2+JT(2)**2))	!phi angle of JT wrt to x axis, radians
	if (JT(2).lt.0.) phi_JT=2.*pi-phi_JT
c
     	W_unpol = m_part**2*JS**2+(x**2+(1.-x)**2)*(JT(1)**2+JT(2)**2)
c     	W_pol = -2.*x*(1.-x)*((p2/c2)**2*cos(2.*phi2)+(p1/c1)**2*cos(2.*phi1)+2.*(p2/c2)*(p1/c1)*cos(phi1+phi2)) !the Bakmaev expression
c	xsctn=2.*alpha**3*Z**2*E0**4*x**2*(1.-x)**2/(pi**2*q2_T**2)*(W_unpol+pol*W_pol) ! note the absence of cos(2phi_JT) in this expression
c     &	*hbarc**2/100.*1.e9*FF2(q2_T,ztgt,nuc_FF) !units of nb/sr^2  The denominator uses the transverse 3-momentum transfer^2, 
	W_pol = -2.*x*(1.-x)*(JT(1)**2+JT(2)**2)	!this is my reduction of the Bakmaev equations
	xsctn=2.*alpha**3*Z**2*E0**4*x**2*(1.-x)**2/(pi**2*q2_T**2)*(W_unpol+pol*cos(2.*phi_JT)*W_pol)  !this contains the cos(2phi_JT) term
     &	*hbarc**2/100.*1.e9*FF2(q2_T,ztgt,nuc_FF) !units of nb/sr^2  The denominator uses the transverse 3-momentum transfer^2, 
c					
	return
	end
c
c****************************************************************************************
c
	real*8 function FF2(q2,ztgt,nuc_FF)
	implicit none
	real*8 q2,hbarc,z,FF_nuc,FF_TFM,alpha(3),b(3),b0,m_e,c,FF
	integer*8 i,ztgt
	logical nuc_FF
	data m_e,hbarc /0.511e-3,0.197/
	data alpha(1),alpha(2),alpha(3),b(1),b(2),b(3) /0.1,0.55,0.35,6.0,1.2,0.3/
c
	z=float(ztgt)
	c=m_e*z**.333
	FF_nuc=FF(q2,ztgt)
c
	FF_TFM=1.
	do i=1,3
		FF_TFM=FF_TFM-alpha(i)*q2/(q2+(b(i)*c)**2)
	end do
c
	if (nuc_FF) then
		FF2=(FF_nuc-FF_TFM)**2
	else
		FF2=1.
	end if
	return
	end	
	
c
c****************************************************************************************
c
	subroutine analysis(E0,mtgt,k1,k2,ktgt,w_mumu,t,missing_mass,m_part,pion_hypothesis)
	implicit none
	logical pion_hypothesis
	real*8 E0,mtgt,k1(3),k2(3),w_mumu,t,m_pi,m_part,m_x,E1,E2,ks(3),pi,ktgt(3),missing_mass
	data m_pi /0.139570/
c
	pi=acos(-1.)
c
	E1=sqrt(k1(1)**2+k1(2)**2+k1(3)**2+m_part**2)	!lepton energies
	E2=sqrt(k2(1)**2+k2(2)**2+k2(3)**2+m_part**2)
	ks(1)=k1(1)+k2(1)	!lepton summed momentum
	ks(2)=k1(2)+k2(2)
	ks(3)=k1(3)+k2(3)
	ktgt(1)=-ks(1)		!target momentum
	ktgt(2)=-ks(2)
	ktgt(3)=E0-ks(3)
	missing_mass=sqrt((E0+mtgt-E1-E2)**2-ktgt(1)**2-ktgt(2)**2-ktgt(3)**2)
	t=ks(1)**2+ks(2)**2+(E0-ks(3))**2-(E0-E1-E2)**2	!4-momentum transfer squared to nucleus, this is positive
c
c		mu mu invariant mass, possibly with pion hypothesis
	m_x=m_part
	if (pion_hypothesis) m_x=m_pi
	E1=sqrt(k1(1)**2+k1(2)**2+k1(3)**2+m_x**2)		!need to put in the mass hypothesis
	E2=sqrt(k2(1)**2+k2(2)**2+k2(3)**2+m_x**2)
	w_mumu=sqrt((E1+E2)**2-ks(1)**2-ks(2)**2-ks(3)**2)
c
	return
	end



c******************************************************************
c
	subroutine density_init(ztgt)
	implicit none
	real*8 pi,c_den(100),a_den(100),rho0,w
	integer*8 i,ztgt
	common /density_rho0/ rho0,c_den,a_den
c
	pi=acos(-1.)
c
	rho0=0.
	if (ztgt.eq.82..or.ztgt.eq.1.) return
c	These equations have to do with Fermi distribution, reference? 
	w=4.*pi*c_den(ztgt)/3.*( (pi*a_den(ztgt))**2+c_den(ztgt)**2 )
	do i=1,10
		w=w+8.*pi*a_den(ztgt)**3*(-1.)**(i-1)*exp(-float(i)*c_den(ztgt)/a_den(ztgt))/float(i)**3
	end do
	rho0=1./w
	return
	end
c
c******************************************************************

	real*8 function FF(Q2,ztgt)
	implicit none
	real*8 Q2, q02,hbarc,Q,gamma,r(12),A(12),rho0,c_den(100),a_den(100),pi,norm,proton_rms
	integer*8 i,ztgt
c
	data q02,proton_rms /0.71,0.879/	!proton dipole form factor parameter GeV^2, proton rms radius fm
	data hbarc /.197/
	common /density_rho0/ rho0,c_den,a_den
c
	data r(1),A(1)   /0.1,0.003845/
	data r(2),A(2)   /0.7,0.009724/
	data r(3),A(3)   /1.6,0.033093/
	data r(4),A(4)   /2.1,0.000120/
	data r(5),A(5)   /2.7,0.083107/
	data r(6),A(6)   /3.5,0.080869/
	data r(7),A(7)   /4.2,0.139957/
	data r(8),A(8)   /5.1,0.260892/
	data r(9),A(9)   /6.0,0.336013/
	data r(10),A(10) /6.6,0.033637/
	data r(11),A(11) /7.6,0.018729/
	data r(12),A(12) /8.7,0.000020/
	data gamma /1.388/
c
c  Select the FF
c
	pi=acos(-1.)
	q=sqrt(Q2)				
c
	if (ztgt.eq.1) then	!proton	
		FF=1./(1.+Q2/q02)**2+2.*Q2/q02-1./6.*Q2*proton_rms**2/hbarc**2		
c
	else if (ztgt.eq.82) then	!lead
		FF=0.
		do i=1,12
		FF=FF + A(i)*(gamma**2*cos(Q*r(i)/hbarc)+2.*r(i)*hbarc/Q*sin(Q*r(i)/hbarc))/(gamma**2+2.*r(i)**2)*
     &		exp(-Q2/4.*gamma**2/hbarc**2)	
		end do
c
	else	!for everything else use 2-parameter fermi model, reference ? 
	FF=4.*pi**2*rho0*a_den(ztgt)**3/( (q*a_den(ztgt))**2*(sinh(pi*q*a_den(ztgt)))**2 )*
     &		( pi*q*a_den(ztgt)*cosh(pi*q*a_den(ztgt))*sin(q*c_den(ztgt))-q*c_den(ztgt)*cos(q*c_den(ztgt))*sinh(pi*q*a_den(ztgt)) )
	do i=1,10
		FF=FF+8.*pi*rho0*a_den(ztgt)**3*(-1.)**(i-1)*float(i)*exp(-float(i)*c_den(ztgt)/a_den(ztgt))/(float(i)**2+(q*a_den(ztgt))**2)**2
	end do
c
	end if
c
	return
	end
c



*******************************************************************
********	FILE: randgen.f				***********
********	AUTHORS: Richard Chandler		***********
********		 (richard@stats.ucl.ac.uk)	***********
********		 Paul Northrop 			***********
********		 (northrop@stats.ox.ac.uk)	***********
********	LAST MODIFIED: 26/8/03			***********
********	See file randgen.txt for details	***********
*******************************************************************


******************************************************************
******************************************************************
      SUBROUTINE ZBQLINI(SEED)
******************************************************************
*       To initialize the random number generator - either
*       repeatably or nonrepeatably. Need double precision
*       variables because integer storage can't handle the
*       numbers involved
******************************************************************
*	ARGUMENTS
*	=========
*	SEED	(integer, input). User-input number which generates
*		elements of the array ZBQLIX, which is subsequently used 
*		in the random number generation algorithm. If SEED=0,
*		the array is seeded using the system clock if the 
*		FORTRAN implementation allows it.
******************************************************************
*	PARAMETERS
*	==========
*	LFLNO	(integer). Number of lowest file handle to try when
*		opening a temporary file to copy the system clock into.
*		Default is 80 to keep out of the way of any existing
*		open files (although the program keeps searching till
*		it finds an available handle). If this causes problems,
*               (which will only happen if handles 80 through 99 are 
*               already in use), decrease the default value.
******************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
******************************************************************
*	VARIABLES
*	=========
*	SEED	See above
*	ZBQLIX	Seed array for the random number generator. Defined
*		in ZBQLBD01
*	B,C	Used in congruential initialisation of ZBQLIX
*	SS,MM,}	System clock secs, mins, hours and days
*	HH,DD }
*	FILNO	File handle used for temporary file
*	INIT	Indicates whether generator has already been initialised
*
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      INTEGER INIT
      DOUBLE PRECISION ZBQLIX(43),B,C
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE INIT

*
*	Ensure we don't call this more than once in a program
*
      IF (INIT.GE.1) THEN
       IF(INIT.EQ.1) THEN
        WRITE(*,1)
        INIT = 2
       ENDIF
       RETURN
      ELSE
       INIT = 1
      ENDIF
*
*       If SEED = 0, cat the contents of the clock into a file
*       and transform to obtain ZQBLIX(1), then use a congr.
*       algorithm to set remaining elements. Otherwise take
*       specified value of SEED.
*
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>	NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
*>>>>>>>	(NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
*>>>>>>>	THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
*>>>>>>>	OF THE FOLLOWING IF BLOCK SHOULD BE	 >>>>>>>
*>>>>>>>	COMMENTED OUT.				 >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (SEED.EQ.0) THEN
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
*>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CALL SYSTEM(' date +%S%M%H%j > zbql1234.tmp')
*
*       Try all file numbers for LFLNO to 999 
*
       FILNO = LFLNO
 10    OPEN(FILNO,FILE='zbql1234.tmp',ERR=11)
       GOTO 12
 11    FILNO = FILNO + 1
       IF (FILNO.GT.999) THEN
        WRITE(*,2)
        RETURN
       ENDIF
       GOTO 10
 12    READ(FILNO,'(3(I2),I3)') SS,MM,HH,DD
       CLOSE(FILNO)
       CALL SYSTEM('rm zbql1234.tmp')
       DSS = DINT((DBLE(SS)/6.0D1) * B)
       DMM = DINT((DBLE(MM)/6.0D1) * B)
       DHH = DINT((DBLE(HH)/2.4D1) * B)
       DDD = DINT((DBLE(DD)/3.65D2) * B)
       TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,B)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<
*<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)       
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
     +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
     +' find an',/5X,
     +'available file number. To rectify the problem, decrease the ',
     +'value of',/5X,
     +'the parameter LFLNO at the start of this routine (in file ',
     +'randgen.f)',/5X,
     +'and recompile. Any number less than 100 should work.')
      END
******************************************************************
      FUNCTION ZBQLU01(DUMMY)
*
*       Returns a uniform random number between 0 & 1, using
*       a Marsaglia-Zaman type subtract-with-borrow generator.
*       Uses double precision, rather than integer, arithmetic 
*       throughout because MZ's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       Ideally, we would explicitly truncate all integer 
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of DNINT to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). This double precision version 
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably OK ...
*
      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
*
*     Update array pointers. Do explicit check for bounds of each to
*     avoid expense of modular arithmetic. If one of them is 0 the others
*     won't be
*
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
*
*     The integer arithmetic there can yield X=0, which can cause 
*     problems in subsequent routines (e.g. ZBQLEXP). The problem
*     is simply that X is discrete whereas U is supposed to 
*     be continuous - hence if X is 0, go back and generate another
*     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
*
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END
******************************************************************
      FUNCTION ZBQLUAB(A,B)
*
*       Returns a random number uniformly distributed on (A,B)
*
      DOUBLE PRECISION A,B,ZBQLU01,ZBQLUAB
      
*
*       Even if A > B, this will work as B-A will then be -ve
*
      IF (A.NE.B) THEN
       ZBQLUAB = A + ( (B-A)*ZBQLU01(0.0D0) )
      ELSE
       ZBQLUAB = A
       WRITE(*,1)
      ENDIF
 1    FORMAT(/5X,'****WARNING**** (function ZBQLUAB) Upper and ',
     +'lower limits on uniform',/5X,'distribution are identical',/)
      END
******************************************************************
      FUNCTION ZBQLEXP(MU)
*
*       Returns a random number exponentially distributed with
*       mean MU
*
      DOUBLE PRECISION MU,ZBQLEXP,ZBQLU01

      ZBQLEXP = 0.0D0

      IF (MU.LT.0.0D0) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      ZBQLEXP = -DLOG(ZBQLU01(0.0D0))*MU

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLEXP',/)

      END
******************************************************************
      FUNCTION ZBQLNOR(MU,SIGMA)
*
*       Returns a random number Normally distributed with mean
*       MU and standard deviation |SIGMA|, using the Box-Muller
*       algorithm
*
      DOUBLE PRECISION THETA,R,ZBQLNOR,ZBQLU01,PI,MU,SIGMA
      DOUBLE PRECISION SPARE
      INTEGER STATUS
      SAVE STATUS,SPARE,PI
      DATA STATUS /-1/

      IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

      IF (STATUS.LE.0) THEN
       THETA = 2.0D0*PI*ZBQLU01(0.0D0)
       R = DSQRT( -2.0D0*DLOG(ZBQLU01(0.0D0)) )
       ZBQLNOR = (R*DCOS(THETA))
       SPARE = (R*DSIN(THETA))
       STATUS = 1
      ELSE
       ZBQLNOR = SPARE
       STATUS = 0
      ENDIF
      
      ZBQLNOR = MU + (SIGMA*ZBQLNOR)

      END
******************************************************************
      FUNCTION ZBQLBIN(N,P)
*
*       Returns a random number binomially distributed (N,P)
*
      DOUBLE PRECISION P,ZBQLBET1
      DOUBLE PRECISION PP,PPP,G,Y,TINY
      INTEGER N,ZBQLBIN,ZBQLGEO,IZ,NN

      TINY = 1.0D-8
      ZBQLBIN = 0

      IF (.NOT.( (P.GE.0.0D0).AND.(P.LE.1.0D0) )) THEN
       WRITE(*,1)
       RETURN
      ELSEIF(N.LE.0) THEN
       WRITE(*,1)
       RETURN
      ENDIF
*
*	First step: if NP > 10, say, things will be expensive, and 
*	we can get into the right ballpark by guessing a value for
*	ZBQLBIN (IZ, say), and simulating Y from a Beta distribution 
*	with parameters IZ and NN-IZ+1 (NN starts off equal to N).
*	If Y is less than PP (which starts off as P) then the IZth order 
*	statistic from NN U(0,1) variates is less than P, and we know 
*	that there are at least IZ successes. In this case we focus on 
*	the remaining (NN-IZ) order statistics and count how many are
*	less than PP, which is binomial (NN-IZ,(PP-Y)/(1-Y)). 
*	Otherwise, if Y is greater than PP there must be less 
*	than IZ successes, so we can count the number of order statistics
*	under PP, which is binomial (IZ-1,P/Y). When we've got NN*PP
*	small enough, we go to the next stage of the algorithm and 
*	generate the final bits directly.
*
      NN = N
      PP = P
 10   IZ = INT(DBLE(NN)*PP) + 1
      IF ( (IZ.GT.10).AND.(IZ.LT.NN-10) ) THEN
       Y = ZBQLBET1(DBLE(IZ),DBLE(NN-IZ+1))
       IF (Y.LT.PP) THEN
        ZBQLBIN = ZBQLBIN + IZ
        NN = NN - IZ
        PP = (PP-Y) / (1.0D0-Y)
       ELSE
        NN = IZ-1
        PP = PP/Y
       ENDIF
       GOTO 10
      ENDIF
*
*	PP is the probability of the binomial we're currently
*	simulating from. For the final part, we simulate either number 
*	of failures or number of success, depending which is cheaper.
*      
 20   IF (PP.GT.0.5) THEN
       PPP = 1.0D0-PP
      ELSE
       PPP = PP
      ENDIF

      G = 0
      IZ = 0
*
*     ZBQLGEO falls over for miniscule values of PPP, so ignore these
*     (tiny probability of any successes in this case, anyway)
* 
      IF (PPP.GT.TINY) THEN
 30    G = G + ZBQLGEO(PPP)
       IF (G.LE.NN) THEN
        IZ = IZ + 1
        GOTO 30
       ENDIF
      ENDIF

      IF (PP.GT.0.5) IZ = NN - IZ
      ZBQLBIN = ZBQLBIN + IZ

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLBIN',/)
      END
******************************************************************
      FUNCTION ZBQLGEO(P)
*
*       Returns a random number geometrically distributed with 
*       parameter P ie. mean 1/P
* 
  
      DOUBLE PRECISION P,ZBQLU01,U,TINY
      INTEGER ZBQLGEO

      TINY = 1.0D-12
      ZBQLGEO = 0
 
      IF (.NOT.( (P.GE.0.0D0).AND.(P.LE.1.0D0) )) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      IF (P.GT.0.9D0) THEN
 10    ZBQLGEO = ZBQLGEO + 1 
       U = ZBQLU01(0.0D0)
       IF (U.GT.P) GOTO 10
      ELSE
       U = ZBQLU01(0.0D0)
*
*	For tiny P, 1-p will be stored inaccurately and log(1-p) may
*	be zero. In this case approximate log(1-p) by -p
*
       IF (P.GT.TINY) THEN
        ZBQLGEO = 1 + INT( DLOG(U)/DLOG(1.0D0-P) )
       ELSE
        ZBQLGEO = 1 + INT(-DLOG(U)/P)
       ENDIF
      ENDIF

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLGEO',/)
      END
******************************************************************
      FUNCTION ZBQLPOI(MU)
*
*       Returns a random number Poisson distributed with mean MU
*
      
      DOUBLE PRECISION ZBQLU01,X,Y,MU,PI
      DOUBLE PRECISION ZBQLLG,ZBQLGAM,MU1,TMP1,TMP2,T
      INTEGER ZBQLPOI,ZBQLBIN,K,INIT
      SAVE INIT,PI
      DATA INIT /0/

      IF (INIT.EQ.0) THEN
       PI = 4.0D0*DATAN(1.0D0)
       INIT = 1
      ENDIF

      ZBQLPOI = 0

      IF (MU.LT.0.0D0) THEN
       WRITE(*,1)
       RETURN
      ENDIF
*
*      For small MU, generate exponentials till their sum exceeds 1
*      (equivalently, uniforms till their product falls below e^-MU)
*
      IF (MU.LE.1.0D3) THEN
       MU1 = MU
*
*     For values of MU less than 1000, use order statistics - the Kth
*     event in a Poisson process of rate MU has a Gamma distribution
*     with parameters (MU,K); if it's greater than 1 we know that there 
*     are less than K events in (0,1) (and the exact number is binomial)
*     and otherwise the remaining number is another Poisson. Choose K so
*     that we'll get pretty close to 1 in the first go but are unlikely
*     to overshoot it.
*
 19    IF (MU1.GT.1.0D1) THEN        
        K = INT(MU1-DSQRT(MU1))
        Y = ZBQLGAM(DBLE(K),MU1)
        IF (Y.GT.1.0D0) THEN
         ZBQLPOI = ZBQLPOI + ZBQLBIN(K-1,(1.0D0/Y))
         RETURN
        ENDIF
        ZBQLPOI = ZBQLPOI + K
        MU1 = MU1  * (1.0D0-Y)
        GOTO 19
       ENDIF
       Y = DEXP(-MU1)
       X = 1.0D0
 20    X = X*ZBQLU01(0.0D0)
       IF (X.GT.Y) THEN
        ZBQLPOI = ZBQLPOI + 1
        GOTO 20
       ENDIF
*
*     For really huge values of MU, use rejection sampling as in 
*     Press et al (1992) - large numbers mean some accuracy may be
*     lost, but it caps the execution time.
*
      ELSE
       TMP1 = DSQRT(2.0D0*MU)
       TMP2 = ZBQLLG(MU+1.0D0)-(MU*DLOG(MU))
 30    Y = DTAN(PI*ZBQLU01(0.0D0))
       ZBQLPOI = INT(MU + (TMP1*Y))
       IF (ZBQLPOI.LT.0) GOTO 30
       X = DBLE(ZBQLPOI)
       T = (X*DLOG(MU)-ZBQLLG(X+1.0D0)) + TMP2
       IF (DABS(T).LT.1.0D2) THEN
        T = 0.9D0*(1.0D0+(Y*Y))*DEXP(T)
        IF (ZBQLU01(0.0D0).GT.T) GOTO 30
       ELSE
        T = DLOG(0.9D0*(1.0D0+(Y*Y))) + T
        IF (DLOG(ZBQLU01(0.0D0)).GT.T) GOTO 30
       ENDIF
      ENDIF 

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLPOI',/)
      END
******************************************************************
      FUNCTION ZBQLGAM(G,H)
*
*       Returns a random number with a gamma distribution with mean
*       G/H and variance G/(H^2). (ie. shape parameter G & scale
*       parameter H)
*
      DOUBLE PRECISION C,D,R,ZBQLGAM,ZBQLU01,G,H,A,z1,z2,B1,B2,M
      DOUBLE PRECISION U1,U2,U,V,TEST,X
      double precision c1,c2,c3,c4,c5,w

      ZBQLGAM = 0.0D0

      IF ( (G.LE.0.0D0).OR.(H.LT.0.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      IF (G.LT.1.0D0) THEN
889    u=ZBQLU01(0.0d0)
       v=ZBQLU01(0.0d0)
       if (u.gt.exp(1.0d0)/(g+exp(1.0d0))) goto 891
       ZBQLGAM=((g+exp(1.0d0))*u/exp(1.0d0))**(1.0d0/g)
       if (v.gt.exp(-ZBQLGAM)) then
	goto 889
       else
	goto 892
       endif
891    ZBQLGAM=-log((g+exp(1.0d0))*(1.0d0-u)/(g*exp(1.0d0)))
       if (v.gt.ZBQLGAM**(g-1.0)) goto 889
892    ZBQLGAM=ZBQLGAM/h
       RETURN
      ELSEIF (G.LT.2.0D0) THEN
       M = 0.0D0
      elseif (g.gt.10.0d0) then
       c1=g-1.0d0
       c2=(g-1.0d0/(6.0d0*g))/c1
       c3=2.0d0/c1
       c4=c3+2.0d0
       c5=1.0d0/sqrt(g)
777    u=ZBQLU01(0.0d0)
       v=ZBQLU01(0.0d0)
       if (g.gt.2.50d0) then
	u=v+c5*(1.0d0-1.860d0*u)
       endif 
       if (u.le.0.0d0.or.u.ge.1.0d0) goto 777 
       w=c2*v/u 
       if (c3*u+w+1.0d0/w.le.c4) goto 778 
       if (c3*log(u)-log(w)+w.ge.1.0d0) goto 777
778    ZBQLGAM=c1*w/h 
       return
      ELSE
       M = -(G-2.0D0) 
      ENDIF
      R = 0.50D0
      a = ((g-1.0d0)/exp(1.0d0))**((g-1.0d0)/(r+1.0d0))
      C = (R*(M+G)+1.0D0)/(2.0D0*R)
      D = M*(R+1.0D0)/R
      z1 = C-DSQRT(C*C-D)
*
*     On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
*     is not exactly zero - this needs trapping if negative.
*
      IF ((Z1-M.LT.0.0D0).AND.(Z1-M.GT.-1.0D-12)) Z1 = M
      z2 = C+DSQRT(C*C-D)
      B1=(z1*(z1-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z1-M)/(R+1.0D0))
      B2=(z2*(z2-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z2-M)/(R+1.0D0))
50    U1=ZBQLU01(0.0D0)
      U2=ZBQLU01(0.0D0)
      U=A*U1
      V=B1+(B2-B1)*U2
      X=V/(U**R)
      IF (X.LE.M) GOTO 50
      TEST = ((X-M)**((G-1)/(R+1)))*EXP(-(X-M)/(R+1.0D0))
      IF (U.LE.TEST) THEN
       ZBQLGAM = (X-M)/H
      ELSE
       GOTO 50
      ENDIF
 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLGAM',/5X, '(both parameters must be positive)',/)

      END
***************************************************************
      FUNCTION ZBQLBET1(NU1,NU2)
*
*       Returns a random number, beta distributed with degrees
*       of freedom NU1 and NU2.
*
      DOUBLE PRECISION NU1,NU2,ZBQLGAM,ZBQLBET1,ZBQLU01,X1,X2

      ZBQLBET1 = 0.0D0

      IF ( (NU1.LE.0.0).OR.(NU2.LE.0.0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF
*      
*       If parameters are too small, gamma subroutine tends to return zero
*       as all the probability goes to the origin and we get rounding
*       errors, even with double precision. In this case, we use Johnk's
*       method, suitably scaled to avoid rounding errors as much as possible.
*
      
      IF ( (NU1.LT.0.9D0).AND.(NU2.LT.0.9D0) ) THEN
 10    X1 = ZBQLU01(0.0D0)
       X2 = ZBQLU01(0.0D0)
       IF ( (X1**(1.0D0/NU1))+(X2**(1.0D0/NU2)).GT.1.0D0) GOTO 10    
       X1 = (DLOG(X2)/NU2) - (DLOG(X1)/NU1)
       ZBQLBET1 = (1.0D0 + DEXP(X1))**(-1)
       IF (ZBQLBET1.GT.1.0D0) GOTO 10    
      ELSE
       X1 = ZBQLGAM(NU1,1.0D0)
       X2 = ZBQLGAM(NU2,1.0D0)
       ZBQLBET1 = X1/(X1+X2)
      ENDIF
       
      RETURN

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLBET1',/5X,
     +'(both degrees of freedom must be positive)',/)

      END
***************************************************************
      FUNCTION ZBQLWEI(A,B)
*
*       Returns a random number, Weibull distributed with shape parameter
*       A and location parameter B, i.e. density is
*	f(x) = ( A/(B**A) ) * x**(A-1) * EXP( -(x/B)**A )
*
      DOUBLE PRECISION A,B,ZBQLU01,ZBQLWEI,U

      ZBQLWEI = 0.0D0

      IF ( (A.LE.0.0).OR.(B.LE.0.0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF
 
      U = ZBQLU01(0.0D0)
      ZBQLWEI = B * ( (-DLOG(U))**(1.0D0/A) )

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLWEI',/5X,
     +'(both parameters must be positive)',/)
      END
***************************************************************
      FUNCTION ZBQLNB(R,P)
*
*       Returns a pseudo-random number according to a Negative
*	Binomial distribution with parameters (R,P). NB these are
*	both DOUBLE - it copes with non-integer R as well. The
*       form of the distribution is *not* the no. of trials to 
*       the Rth success - see documentation for full spec.
*
      DOUBLE PRECISION R,P,ZBQLGAM,Y
      INTEGER ZBQLNB,ZBQLPOI

      ZBQLNB = 0

      IF ( (R.LE.0.0D0).OR.(P.LE.0.0D0).OR.(P.GE.1.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      Y = ZBQLGAM(R,1.0D0)
      Y = Y*P/(1.0D0-P)
      ZBQLNB = ZBQLPOI(Y)

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLNB')
      END
***************************************************************
      FUNCTION ZBQLPAR(A,B)
*
*     Returns a random number, Pareto distributed with parameters
*     A and B. The density is A*(B**A) / (B+X)**(A+1) for X > 0.
*     (this is slightly nonstandard - see documentation in 
*     randgen.txt). The algorithm is straightforward - it uses the
*     inverse CDF method.
*
      DOUBLE PRECISION A,B,ZBQLPAR,ZBQLU01,U

      ZBQLPAR = 0.0D0

      IF ( (A.LE.0.0D0).OR.(B.LE.0.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF
 
      U = ZBQLU01(0.0D0)
      ZBQLPAR = B * (U**(-1.0D0/A)-1.0D0)

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLPAR',/5X,
     +'(both parameters must be positive)',/)
      END
***************************************************************
      FUNCTION ZBQLLG(X)
*
*     Returns log(G(X)) where G is the Gamma function. The algorithm is
*     that given in Press et al (1992), Section 6.1, although this
*     version also allows for arguments less than 1.
*
      DOUBLE PRECISION X,Z,Z2,ZBQLLG,PI,RLN2P,C(0:6),TMP,SUM
      INTEGER INIT,I
      SAVE INIT,C,RLN2P,PI
      DATA INIT /0/
      DATA (C(I),I=0,6) / 
     +              1.000000000190015D0,76.18009172947146D0,
     +              -86.50532032941677D0,24.01409824083091D0,
     +              -1.231739572450155D0,0.1208650973866179D-2,
     +              -0.5395239384953D-5/

      IF (INIT.EQ.0) THEN
        PI = 4.0D0*DATAN(1.0D0)
        RLN2P = 0.5D0*DLOG(2.0D0*PI)
        INIT = 1
      ENDIF
*
*     Compute for x > 1, then use transformation if necessary. Z is
*     our working argument.
*
      IF (X.GE.1.0D0) THEN
       Z = X
      ELSE 
       Z = 2.0D0-X
       Z2 = 1.0D0-X
      ENDIF

      IF (DABS(Z-1.0D0).LT.1.0D-12) THEN
       ZBQLLG = 0.0D0
       RETURN
      ENDIF

      TMP = Z + 4.5D0
      TMP = ( (Z-0.5D0)*DLOG(TMP) ) - TMP + RLN2P

      SUM = C(0)
      DO 50 I=1,6
       SUM = SUM + (C(I)/(Z+DBLE(I-1)))
 50   CONTINUE
      ZBQLLG = TMP + DLOG(SUM)
*
*     Transformation required if X<1
*
      IF (X.LT.1.0D0) THEN
       TMP = PI*Z2
       ZBQLLG = DLOG(TMP/DSIN(TMP)) - ZBQLLG
      ENDIF

      END
	
			



	



				


