! MECH 510, Programming Assignment 3
! Timothy Chui, 37695129
! Program to numerically solve the incompressible energy equation
! Time advance schemes: RK2, IE
! Space discretization scheme: 2nd-Order Centred
! Dec 2016


! =================== modules =========================

MODULE meshmod
    ! Declare variables to pass arround
	
	! Discretization parameters
    INTEGER, PARAMETER :: imax = 400                          ! Size of mesh (x-direction)
	INTEGER, PARAMETER :: jmax = 160                          ! Size of mesh (y-direction)
 
    REAL*8, PARAMETER  :: xmax = 15.0D+0                      ! Maximum x-position
	REAL*8, PARAMETER  :: ymax = 1.0D+0                       ! Maximum y-position
    REAL*8, PARAMETER  :: tmax = 40.0D+0                      ! Maximum time
	
	INTEGER, PARAMETER :: scheme = 2                          ! 1 for RK2, 2 for IE
	
	INTEGER, PARAMETER :: problem = 3                         ! 1 for test (square domain), 2 for channel, 3 for real (long channel)
	
	! Required data
	REAL*8, PARAMETER  :: u0 = 1.0                            ! Constant for u(x,y) (test problems)
	REAL*8, PARAMETER  :: v0 = 1.0                            ! Constant for v(x,y)
    REAL*8, PARAMETER  :: T0 = 1.0                            ! Constant for T(x,y)
	REAL*8, PARAMETER  :: Re = 50.0                           ! Reynolds number
	REAL*8, PARAMETER  :: Pr = 0.7                            ! Prandtl number 
	REAL*8, PARAMETER  :: Ec = 0.001                          ! Eckert number
	REAL*8, PARAMETER  :: ubar = 3.0                          ! Constant for u(x,y) (channel problems)
	
	
    ! Values that do not need to be changed
    INTEGER            :: i,j,n                               ! Space and time indices
    REAL*8, PARAMETER  :: pi = 3.141592654D+0                 ! Value of pi
    REAL*8, PARAMETER  :: dx = xmax/imax                      ! Cell sizes
	REAL*8, PARAMETER  :: dy = ymax/jmax

	REAL*8              :: delta_t                            ! Time step variable
    INTEGER             :: maxn                               ! Max time step (i.e. last n index)
	
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmesh                        ! Mesh of solutions; includes ghost cells
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: umesh,vmesh                  ! Mesh of velocities; includes ghost cells
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FIn,FI_ana                   ! Computed and analytic flux integrals at time level n
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: source_mesh,source_mesh_ana  ! Computed and analytic sources 
	
END MODULE meshmod

! Program structure
! -------------------

! Main Program
! ------------

! Subroutines
!--------------

!!! Initialization subroutines !!!
! init_mesh (Mesh initialization of initial conditions and sources)
! init_mesh_test (Apply test initial conditions)
! init_mesh_channel (Apply real initial conditions)
! init_source (Calculate source terms)

!!! Integration subroutines !!!
! set_bcs (update boundary conditions)
! rk2 (RK2 integrator)
! ie (IE integrator)
! set_lhs (For IE approximate factorization)
! calc_integrals (Flux integral calculator)
! Thomas (Thomas algorithm for IE)

!!! Output subroutines !!!
! output_integrals (output flux integrals and sources for testing)
! output_Tmesh (output solutions)


! Functions
! ----------

! flux (compute flux at faces)
! source (compute source terms)
! analytic (compute analytic fluxes for testing) 
! source_analytic (compute analytic sources for testing)


! -------------------

! =================== main program =========================

PROGRAM mech510_3   

    USE meshmod
    IMPLICIT NONE

    REAL*8   ::   tn = 0.0D+0    ! Time at index n
	
	! Determine time step
	IF (problem == 2) THEN ! Short channel problem
		IF (scheme == 1) THEN ! RK2 (use hardcoded max stable time steps)
			IF (jmax == 10) THEN
				delta_t = 0.8*0.020
			ELSE IF (jmax == 20) THEN
				delta_t = 0.8*0.015
			ELSE IF (jmax == 40) THEN
				delta_t = 0.8*0.0095
			ELSE IF (jmax == 80) THEN
				delta_t = 0.8*0.0025
			ELSE 
				! An automatic case for a non-standard mesh for the channel
				delta_t = (Re*Pr*dy**2)/2. - 0.75*(Re*Pr*dy**2)/2.
			END IF
		ELSE 
			delta_t = 0.01! IE time step
		END IF
			
	ELSE ! Other problems
	
		IF (scheme == 1) THEN ! RK2
			delta_t = (Re*Pr*dy**2)/2. - 0.75*(Re*Pr*dy**2)/2. ! Use automatic case
		ELSE
			delta_t = 0.1     ! IE
		END IF
		
	END IF
	
		
	maxn = nint(tmax/delta_t) ! Determine how many steps to take
	
	
	WRITE(*,*) delta_t
	
	! Initialize mesh
	CALL init_mesh

	! Output mesh at start
	CALL output_Tmesh(tn)

	! Determine which integrator to use
	IF (scheme == 1) THEN 
		CALL rk2(0,maxn-1,tn) ! RK2
	ELSE
		CALL ie(0,maxn,tn)    ! IE
	END IF
	
	! Output final solution mesh
	CALL output_Tmesh(tn)
	
END PROGRAM mech510_3


! =================== subroutines =========================


SUBROUTINE init_mesh           ! Subroutine to initialize T, velocities and sources

	USE meshmod
    IMPLICIT NONE
	
	IF (problem == 1) THEN     ! Testbed
		CALL init_mesh_test
		
	ELSE 
		CALL init_mesh_channel ! Channel problems
		
	END IF
	
END SUBROUTINE init_mesh
	
! =======================

SUBROUTINE init_mesh_test     ! Initialize mesh with CV average values of T(x,y), u(x,y), v(x,y) and source terms for test case
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8 :: xa,xb,yc,yd    ! Positions of cell edges
	
	! Pass through array: initialize Tmesh and vector field
	DO j = 0,jmax+1          ! Loop over rows
	
	    yc = (j-1)*dy        ! Bottom edge
		yd = j*dy            ! Top edge
		
	    DO i = 0, imax+1     ! Loop over columns
		
		    xa = (i-1)*dx
			xb = i*dx
		  	
            ! Compute control volume averages			
		    Tmesh(j,i) = (T0/((pi**2)*dx*dy))*(sin(pi*xb)-sin(pi*xa))*(cos(pi*yc)-cos(pi*yd))
	        umesh(j,i) = (u0/(pi*dx*dy))*(1./2.*(yd**2-yc**2)*(cos(pi*xa)-cos(pi*xb)))
			vmesh(j,i) = (v0/(pi*dx*dy))*(1./2.*(xb**2-xa**2)*(sin(pi*yd)-sin(pi*yc)))
			
			! Testing with other velocities (legacy)
			!umesh(j,i) = 0
			!(j,i) = 0
			!Tmesh(j,i) = T0*cos(pi*(i-1./2.)*dx)*(sin(pi*(j-1./2.)*dy))
	        !umesh(j,i) = u0*((j-1./2.)*dy)*sin(pi*(i-1./2.)*dx)
			!vmesh(j,i) = v0*((i-1./2.)*dx)*cos(pi*(j-1./2.)*dy)
		END DO
	END DO
	
	! Call subroutine to calculate source terms
	CALL init_source


END SUBROUTINE init_mesh_test


! =======================

SUBROUTINE init_mesh_channel    ! Initialize mesh with CV average values of T(x,y), u(x,y), v(x,y) and source terms for channel cases
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8                   :: xa,xb,yc,yd,ymid  ! Positions of cell edges
	REAL*8                   :: dudx,dudy,dvdx,dvdy    ! Flux terms for u and v
	REAL*8                   :: source                 ! Source at face
	REAL*8                   :: source_analytic        ! Analytic source term
	
	
	! First pass through array: initialize Tmesh and vector field
	DO j = 0,jmax+1          ! Loop over rows
	
	    yc = (j-1)*dy        ! Bottom edge
		yd = j*dy            ! Top edge
		
	    DO i = 0, imax+1     ! Loop over columns
		
		    xa = (i-1)*dx
			xb = i*dx
		  	
            ! Compute control volume averages			
		    Tmesh(j,i) = 1./(2.*dy*dx)*(yd**2-yc**2)*(xb-xa)
	        !umesh(j,i) = 1./(2.*dy*dx)*(3.*ubar*(yd**2-yc**2) - 2.*ubar*(yd**3 - yc**3))*(xb-xa)
			
			ymid = (j-1./2.)*dy
			umesh(j,i) = 6.*ubar*ymid*(1.-ymid)
			vmesh(j,i) = 0.
			
		END DO
	END DO
	
	! Call subroutine to calculate source terms
	CALL init_source
	
END SUBROUTINE init_mesh_channel

! =======================

SUBROUTINE init_source  ! Subroutine to initialize source term

	USE meshmod
    IMPLICIT NONE
	
	REAL*8                    :: dudx,dudy,dvdx,dvdy    ! Flux terms for u and v
	REAL*8                    :: source                 ! Source at face
	REAL*8                    :: source_analytic        ! Analytic source term
	
	! Second pass through array; calculate source terms
	DO j = 1,jmax          ! Loop over rows

		DO i = 1,imax     ! Loop over columns

			! Compute source terms along x and y directions
			dudx = (source(umesh(j,i),umesh(j,i+1),dx) + source(umesh(j,i-1),umesh(j,i),dx))/2.
			dvdy = (source(vmesh(j,i),vmesh(j+1,i),dy) + source(vmesh(j-1,i),vmesh(j,i),dy))/2.
			
			dudy = (source(umesh(j,i),umesh(j+1,i),dy) + source(umesh(j-1,i),umesh(j,i),dy))/2.
			dvdx = (source(vmesh(j,i),vmesh(j,i+1),dx) + source(vmesh(j,i-1),vmesh(j,i),dx))/2.
			
			! Assemble source term
			source_mesh(j,i) = Ec/Re*(2.*dudx**2 + 2.*dvdy**2 + (dvdx + dudy)**2)
			
			! Calculate analytic source terms at CV centroid
			source_mesh_ana(j,i) = source_analytic((i-1./2.)*dx, (j-1./2.)*dy)


		END DO
		
	END DO
	
END SUBROUTINE init_source


! =======================

SUBROUTINE set_bcs(tn, Tmeshn)   ! Subroutine to update boundary conditions

    USE meshmod
    IMPLICIT NONE
    
    REAL*8                                :: tn                 ! Time at index n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1)  :: Tmeshn             ! Mesh at time level n (could be interim mesh)
	REAL*8                                :: Tw,Ttop,Tbot       ! Wall values
    REAL*8                                :: yc,yd,xmid         ! Cell positions
	
	IF (problem == 1) THEN   ! Handle test problem
		DO j = 1, jmax
			Tmeshn(j,0) = Tmeshn(j,1)
			!Tmeshn(j,imax+1) = 2.*Tw-Tmeshn(j,imax)
			Tmeshn(j,imax+1) = Tmeshn(j,imax)
		END DO
		
		DO i = 1, imax
			Tmeshn(0,i) = Tmeshn(1,i)
			!Tmeshn(jmax+1,i) = 2.*Tw-Tmeshn(jmax,i)
			Tmeshn(jmax+1,i) = Tmeshn(jmax,i)
		END DO
	
	ELSE IF (problem == 2) THEN  ! Handle channel problem (Neumann outflow)
		Ttop = 1
		Tbot = 0
		
		DO j = 1, jmax
			yc = (j-1)*dy
			yd = j*dy
			
			! Left boundary
			Tw = ((-2.*yc**2 + 3*Ec*Pr*ubar**2*(1./10.*((2*yc-1)**5-(2*yd-1)**5)-yc+yd) &
		    		+2.*yd**2)/(4.*dy))		
					
			Tmeshn(j,0) = 2.*Tw - Tmeshn(j,1) ! Dirichlet (left)

			Tmeshn(j,imax+1) = Tmeshn(j,imax) ! Neumann (right)
		END DO
		
		DO i = 1, imax
			Tmeshn(0,i) = 2.*Tbot - Tmeshn(1,i)
			Tmeshn(jmax+1,i) = 2.*Ttop - Tmeshn(jmax,i)
		END DO
	
	ELSE  ! Handle real problem
		Tbot = 0
		
		DO j = 1, jmax
			yc = (j-1)*dy
			yd = j*dy
			Tw = ((-2.*yc**2 + 3*Ec*Pr*ubar**2*(1./10.*((2*yc-1)**5-(2*yd-1)**5)-yc+yd) &
		    		+2.*yd**2)/(4.*dy))		
			
			Tmeshn(j,0) = 2.*Tw - Tmeshn(j,1)

			Tmeshn(j,imax+1) = Tmeshn(j,imax)
		END DO
		
		DO i = 1, imax
		
			xmid = (i-1./2.)*dx
			
			! Handle space-varying top wall BC
			IF (xmid >= 1. .AND. xmid <= 3.) THEN
				Ttop = 2. - cos(pi*(xmid-1.))
			ELSE
				Ttop = 1.
			END IF
				
				
			Tmeshn(0,i) = 2.*Tbot - Tmeshn(1,i) ! Dirichlet (bottom)
			Tmeshn(jmax+1,i) = 2.*Ttop - Tmeshn(jmax,i) ! Dirichlet (top)
		END DO
	
	END IF

END SUBROUTINE set_bcs

! =======================

SUBROUTINE rk2(nmin, nmax, tn)                        ! Subroutine for time integration (RK2)

    USE meshmod
    IMPLICIT NONE
    
    INTEGER                              :: nmin         ! Initial time index
    INTEGER                              :: nmax         ! Final time index
    REAL*8                               :: Tw           ! Wall value
    REAL*8                               :: tn           ! Time at index n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: mesh_1       ! Interim mesh for time advance scheme
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI_1         ! Interim flux integrals for time advance scheme
	
	REAL*8, DIMENSION(1:jmax,1:imax)     :: rhsn           ! Right-hand-side of energy equation at time level n
	REAL*8                               :: rhs1           ! Interim right-hand-side of energy equation at time level n

    DO n = nmin,nmax  ! Want to obtain solution at maxn, hence last loop is maxn-1
        
        ! ------Calculate T(1)------
        tn = delta_t*n                     ! Current time
        
        CALL set_bcs(tn, Tmesh)         ! Update boundary conditions of solution mesh
        !CALL calc_integral(mesh, FIn, Tw)  ! Evaluate flux integrals at time level n
				
		! Calculate flux integrals
        CALL calc_integrals(Tmesh,FIn)
        
        ! Go from left to right across mesh and calculate interim values
        DO j = 1,jmax
		    DO i = 1,imax
			    rhsn(j,i) = source_mesh(j,i)-FIn(j,i)
                mesh_1(j,i) = Tmesh(j,i) + delta_t*rhsn(j,i)
			END DO
        END DO
        

        ! ------Calculate T(n+1)------
        tn = delta_t*(n+1)                   ! Time n+1
		
        CALL set_bcs(tn, mesh_1)         ! Update boundary conditions of interim mesh
        
		CALL calc_integrals(mesh_1, FI_1) ! Evaluate interim flux integrals at time level n+1
        
        ! Go from left to right across mesh and calculate new solution
        DO j = 1,jmax
		    DO i = 1,imax
			    rhs1 = source_mesh(j,i)-FI_1(j,i)
                Tmesh(j,i) = Tmesh(j,i) + delta_t*(rhs1 + rhsn(j,i))/2.  ! Overwrite old solution array
		    END DO
        END DO

    END DO
	

    
END SUBROUTINE rk2


! =======================

SUBROUTINE ie(nmin, nmax, tn)                        ! Subroutine for time integration (IE); rows first

    USE meshmod
    IMPLICIT NONE
    
    INTEGER                              :: nmin         ! Initial time index
    INTEGER                              :: nmax         ! Final time index
    REAL*8                               :: Tw           ! Wall value
    REAL*8                               :: tn           ! Time at index n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: mesh_1       ! Interim mesh for time advance scheme
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI_1         ! Interim flux integrals for time advance scheme
	
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: dT_tilda     ! Interim dT (rows first)
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: dT           ! T(n+1) = T(n) + dT
	
	REAL*8, DIMENSION(0:imax+1)          :: rhs = 0.0D+0         ! Right-hand-side of energy equation at time level n
	REAL*8                               :: rhs1           ! Interim right-hand-side of energy equation at time level n
	
	REAL*8                               :: alpha
    REAL*8, DIMENSION(0:imax+1,3)        :: lhsx
	REAL*8, DIMENSION(0:jmax+1,3)        :: lhsy
	REAL*8                               :: beta
	
	
	DO n = nmin,nmax  ! Want to obtain solution at maxn, hence last loop is maxn-1
        
        ! ------Calculate T(1)------
        tn = delta_t*n                     ! Current time
        
        CALL set_bcs(tn, Tmesh)         ! Update boundary conditions of solution mesh
        CALL calc_integrals(Tmesh,FIn)      ! Calculate flux integrals
		
		
        ! Go across columns from row to row; first part of approx factorization
		alpha = delta_t/(Re*Pr*dx**2)
        DO j = 1,jmax
		    !WRITE(*,*) j
            rhs = (source_mesh(j,:) - FIn(j,:))*delta_t
			
			! Set LHS for approximate factorization (first part)
			lhsx = 0.0D+0
			beta = umesh(j,1)*delta_t/(2.*dx)
   
			CALL set_lhs(lhsx,imax,alpha,beta,1,1,1,-1) ! Dirichlet BCS on left, Neumann BCs on right
			CALL Thomas(lhsx,rhs,imax)
			
			dT_tilda(j,:) = rhs

        END DO
		
		! Go across rows from column to column; second part of approx factorization
		alpha = delta_t/(Re*Pr*dy**2)
		DO i = 1,imax
			! Set LHS for approximate factorization (second part)
			lhsy = 0.0D+0
			beta = vmesh(1,i)*delta_t/(2.*dy)
			dT_tilda(0,i) = 0.0D+0
			dT_tilda(jmax+1,i) = 0.0D+0
			CALL set_lhs(lhsy,jmax,alpha,beta,1,1,1,1) ! Dirichlet BCs on top and bottom
			CALL Thomas(lhsy,dT_tilda(:,i),jmax)
			
			dT(:,i) = dT_tilda(:,i)
        END DO
		
		! Update solution
		DO j = 1,jmax
			DO i = 1,imax
				Tmesh(j,i) = Tmesh(j,i) + dT(j,i)
			END DO
		END DO 
		
    END DO
	
END SUBROUTINE ie


! =======================

SUBROUTINE set_lhs(lhs,vsize,alpha,betas,bca,bcb,bcc,bcd)     ! Initialize LHS for approximate factorization, for IE time advance
	
	USE meshmod
    IMPLICIT NONE
	

	INTEGER                              :: vsize             ! Size of column (in x, first call) or row (in y, second call)
	REAL*8, DIMENSION(0:vsize+1,3)       :: lhs 
	REAL*8                               :: alpha
	REAL*8                               :: betas
	INTEGER                              :: bca,bcb,bcc,bcd
	INTEGER                              :: jind,iind
	

	DO jind = 0,vsize+1
	
		iind = jind 
		
		IF (jind == 0) THEN     ! Boundary conditions (left and bottom)
		    lhs(jind,1) = 0
			lhs(jind,2) = bca   
			lhs(jind,3) = bcb
					
		ELSE IF (jind  == vsize+1) THEN  ! Boundary conditions (right and top)
		    lhs(jind,1) = bcc 
			lhs(jind,2) = bcd
			lhs(jind,3) = 0
		
		ELSE
			
			lhs(jind,1) = -alpha-betas ! Main diagonals
			lhs(jind,2) = 1.+2.*alpha
			lhs(jind,3) = betas-alpha
			
		END IF
	END DO
		    
		

END SUBROUTINE set_lhs

! =======================

SUBROUTINE calc_integrals(Tmeshn,FI)  ! Subroutine to evaluate flux integrals at time level n, with sources

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmeshn       ! Mesh at time level n (could be interim mesh)
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI      ! Mesh at time level n (could be interim mesh)
	
    REAL*8                               :: FIx, FIy     ! Fluxes in the x and y directions
	REAL*8                               :: flux         ! Flux at face
	REAL*8                               :: analytic     ! Analytic value of flux integral at cell centroids
    REAL*8                               :: analytic_cv  ! CV average analytic flux integral

	! Calculate flux integral, i.e. (right face) - (left face)
    DO j = 1,jmax      ! Loop over rows (interior CVs)
		
		DO i = 1,imax  ! Loop over rows (interior CVs)
 		
		    ! Compute flux integrals along x and y directions
			FIx = (1./dx*(flux(Tmeshn(j,i),Tmeshn(j,i+1),umesh(j,i),umesh(j,i+1),dx) & 
			       - flux(Tmeshn(j,i-1),Tmeshn(j,i),umesh(j,i-1),umesh(j,i),dx)))
			FIy = (1./dy*(flux(Tmeshn(j,i),Tmeshn(j+1,i),vmesh(j,i),vmesh(j+1,i),dy) &
			       - flux(Tmeshn(j-1,i),Tmeshn(j,i),vmesh(j-1,i),vmesh(j,i),dy)))
			
			! Sum them to get full flux integral
			FI(j,i) = FIx + FIy
			
			! Calculate analytic flux integral at CV centroid
			FI_ana(j,i) = analytic((i-1./2.)*dx, (j-1./2.)*dy)
			
		END DO
		
	END DO

END SUBROUTINE calc_integrals

! =======================

subroutine Thomas(LHS, RHS, Size) ! Written by Carl Olliver-Gooch
	integer Size
	double precision LHS(0:Size+1,3), RHS(0:Size+1)

	integer i
  
	LHS(0,1) = 0
	LHS(Size+1,3) = 0
    !Forward elimination 
	do i = 0, Size
		LHS(i,3) = LHS(i,3) / LHS(i,2)
		RHS(i)   = RHS(i)   / LHS(i,2)

		LHS(i+1,2) = LHS(i+1,2) - LHS(i,3)*LHS(i+1,1)
		RHS(i+1)   = RHS(i+1)   - RHS(i)  *LHS(i+1,1)
	enddo

   ! Last line of elimination
	RHS(Size+1) = RHS(Size+1) / LHS(Size+1,2)

    !Back-substitution 
	do i = Size, 0, -1
		RHS(i) = RHS(i) - RHS(i+1)*LHS(i,3)
	enddo
	
	return
 end


! =======================

SUBROUTINE output_integrals ! Output fluxes and sources for testing
    
    USE meshmod
    IMPLICIT NONE

    INTEGER :: ierror                                                  ! For error catching
    CHARACTER(len=100) :: fname1,fname2,fname3,fname4,fname5,fname6    ! Output filenames
    CHARACTER(len=100) :: fmt                                          ! Output format
        
	! Set filename format
	fmt = "(A8,I0,A6,I0,A6,I0,A4)"  
	
    ! Set filenames
    WRITE(fname1, fmt) "FI_comp_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),".txt"
    WRITE(*,*) "Writing output to ", fname1
    OPEN(UNIT=20, FILE=fname1,STATUS="REPLACE", IOSTAT=ierror)
	
	WRITE(fname2, fmt) "FI_anal_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),".txt"
    WRITE(*,*) "Writing output to ", fname2
    OPEN(UNIT=30, FILE=fname2,STATUS="REPLACE", IOSTAT=ierror)
	
	WRITE(fname3, fmt) "FI_erro_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),".txt"
    WRITE(*,*) "Writing output to ", fname3
    OPEN(UNIT=40, FILE=fname3,STATUS="REPLACE", IOSTAT=ierror)
	

    WRITE(fname4, fmt) "sourc_c_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),".txt"
    WRITE(*,*) "Writing output to ", fname4
    OPEN(UNIT=50, FILE=fname4,STATUS="REPLACE", IOSTAT=ierror)
	
	WRITE(fname5, fmt) "sourc_a_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),".txt"
    WRITE(*,*) "Writing output to ", fname5
    OPEN(UNIT=60, FILE=fname5,STATUS="REPLACE", IOSTAT=ierror)
	
	WRITE(fname6, fmt) "sourc_e_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),".txt"
    WRITE(*,*) "Writing output to ", fname6
    OPEN(UNIT=70, FILE=fname6,STATUS="REPLACE", IOSTAT=ierror)
	
  
	DO j=1,jmax   ! Loop over each row
        WRITE(20,"(999F15.10)") (FIn(j,i), i=1,imax)             ! Implied loop to write out every column (computed flux integral)
		WRITE(30,"(999F15.10)") (FI_ana(j,i), i=1,imax)         ! Analytic flux integral
		WRITE(40,"(999F15.10)") (FI_ana(j,i)-FIn(j,i), i=1,imax) ! Flux integral error
		
		WRITE(50,"(999F15.10)") (source_mesh(j,i), i=1,imax)             ! Implied loop to write out every column (computed flux integral)
		WRITE(60,"(999F15.10)") (source_mesh_ana(j,i), i=1,imax)         ! Analytic flux integral
		WRITE(70,"(999F15.10)") (source_mesh_ana(j,i)-source_mesh(j,i), i=1,imax) ! Flux integral error
    ENDDO


END SUBROUTINE output_integrals

! =======================

SUBROUTINE output_Tmesh(tn) ! Output final solutions for plotting
    
    USE meshmod
    IMPLICIT NONE
    
	REAL*8   ::   tn     ! Time at index n
    INTEGER :: ierror                                                  ! For error catching
    CHARACTER(len=100) :: fname                                        ! Output filenames
    CHARACTER(len=100) :: fmt                                          ! Output format
        
	! Set filename format
	fmt = "(A8,I0,A6,I0,A6,I0,A3,I0,A8,I0,A3,I0,A4)"  
	
    ! Set filenames
    WRITE(fname, fmt) "Tmesh_cn",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_t_",nint(tn),"_scheme_",scheme,"_p_",problem,".txt"
    WRITE(*,*) "Writing output to ", fname
    OPEN(UNIT=90, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)
  
	DO j=1,jmax   ! Loop over each row
        WRITE(90,"(999F15.10)") (Tmesh(j,i), i=1,imax)             ! Implied loop to write out every column (computed flux integral)
    ENDDO

	CLOSE(UNIT=90)
END SUBROUTINE output_Tmesh

! =================== functions ==========================

FUNCTION flux(Ta, Tb, ua, ub, ds)   ! Function to calculate flux at face s+1/2,r (direction independent function)

    USE meshmod
    IMPLICIT NONE

    REAL*8                        :: Ta, Tb          ! Temperatures at adjacent cells along direction s
	REAL*8                        :: ua, ub          ! Velocities at adjacent cells along direction s
    REAL*8                        :: ds              ! dx or dy
	
	REAL*8                        :: term1, term2    ! Terms in the flux calculation
    REAL*8                        :: flux            ! Flux at face s+1/2

    
    ! Compute (uT - 1/(Re*Pr)*dT/dx) = (term1 - term2)
	
	term1 = (ub*Tb + ua*Ta)/2.      ! uT
	term2 = 1./(Re*Pr)*(Tb - Ta)/ds ! 1/(Re*Pr)*dT/dx
	
	flux = term1 - term2            ! Full flux at face
	
	RETURN

END FUNCTION flux


! =================


FUNCTION source(ua, ub, ds)   ! Function to calculate source term fluxes at face s+1/2,r (direction independent function)

    USE meshmod
    IMPLICIT NONE


	REAL*8                        :: ua, ub          ! Velocities at adjacent cells along direction s
    REAL*8                        :: ds              ! dx or dy
	
    REAL*8                        :: source          ! Flux at face s+1/2

    source = (ub - ua)/(ds)
	
	RETURN

END FUNCTION source


! =================

FUNCTION analytic(xpos, ypos)    ! Function to calculate analytic flux integral at a CV centroid

    USE meshmod
    IMPLICIT NONE

    REAL*8     :: xpos,ypos             ! Position of the cell centre
	REAL*8     :: first,second,third    ! Terms in the flux integral
	
    REAL*8     :: analytic              ! Analytic flux integral
    
    ! Compute analytic flux integral
	first = u0*T0*pi*cos(2.*pi*xpos)*ypos*sin(pi*ypos)
	second = v0*T0*pi*xpos*cos(pi*xpos)*cos(2.*pi*ypos)
	third = (2.*T0*(pi**2)*cos(pi*xpos)*sin(pi*ypos))/(Re*Pr)
	
	analytic = first + second + third
	
	RETURN

END FUNCTION analytic

! =================

FUNCTION source_analytic(xpos, ypos)    ! Function to calculate analytic source integral at a CV centroid

    USE meshmod
    IMPLICIT NONE

    REAL*8     :: xpos,ypos             ! Position of the cell centre
	REAL*8     :: first,second,third    ! Terms in the source integral
	
    REAL*8     :: source_analytic       ! Analytic source integral
    
    ! Compute analytic flux integral
	first = 2.*(u0*pi*cos(pi*xpos)*ypos)**2
	second = 2.*(v0*pi*xpos*sin(pi*ypos))**2
	third = (u0*sin(pi*xpos) + v0*cos(pi*ypos))**2
	
	source_analytic = Ec/Re*(first + second + third)
	
	RETURN

END FUNCTION source_analytic









