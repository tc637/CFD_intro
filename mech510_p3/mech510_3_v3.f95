! MECH 510, Programming Assignment 3
! Timothy Chui, 37695129
! Program to numerically solve the incompressible energy equation
! Time advance scheme: RK2
! Space discretization scheme: 2nd-Order Upwind (default), 2nd-Order Centred, 1st-Order Upwind
! Nov 2016


! =================== modules =========================

MODULE meshmod
    ! Declare variables to pass arround
	
	! Discretization parameters
    INTEGER, PARAMETER :: imax = 10                        ! Size of mesh (x-direction)
	INTEGER, PARAMETER :: jmax = 10                       ! Size of mesh (y-direction)
 
    REAL*8, PARAMETER  :: xmax = 1.0D+0                    ! Maximum x-position
	REAL*8, PARAMETER  :: ymax = 1.0D+0                    ! Maximum y-position
    REAL*8, PARAMETER  :: tmax = 1.0D+0                    ! Maximum time
	
	REAL*8, PARAMETER  :: delta_t = 0.001              ! Time step from CFL condition
    INTEGER, PARAMETER :: maxn = nint(tmax/delta_t)    ! Max time step (i.e. last n index)

	! Required data
	REAL*8, PARAMETER  :: u0 = 1.0                         ! Constant for u(x,y)
	REAL*8, PARAMETER  :: v0 = 1.0                         ! Constant for v(x,y)
    REAL*8, PARAMETER  :: T0 = 1.0                         ! Constant for T(x,y)
	REAL*8, PARAMETER  :: Re = 50.0                        ! Reynolds number
	REAL*8, PARAMETER  :: Pr = 0.7                         ! Prandtl number 
	REAL*8, PARAMETER  :: Ec = 0.1                         ! Eckert number
	
    
    ! Values that do not need to be changed
    INTEGER            :: i,j,n                            ! Space and time indices
    REAL*8, PARAMETER  :: pi = 3.141592654D+0              ! Value of pi
    REAL*8, PARAMETER  :: dx = xmax/imax                   ! Grid size
	REAL*8, PARAMETER  :: dy = ymax/jmax

    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmesh          ! Mesh of solutions; includes ghost cells
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: umesh,vmesh    ! Mesh of velocities; includes ghost cells
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI_ana          ! Analytic flux integrals at time level n
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: source_mesh_ana  ! Analytic sources at time level n 
	
END MODULE meshmod

! Program structure
! -------------------

! Main Program


! -------------------

! =================== main program =========================

PROGRAM mech510_3  

    USE meshmod
    IMPLICIT NONE
	
    REAL*8   ::   tn     ! Time at index n

    
	! Initialize mesh
	CALL init_mesh
	
	! Calculate flux integrals
    !CALL calc_integrals(Tmesh)
	
	CALL rk2(0, maxn-1)
	
	! Output computed and analytic flux integrals
	!CALL output

END PROGRAM mech510_3


! =================== subroutines =========================

SUBROUTINE init_mesh         ! Initialize mesh with CV average values of T(x,y), u(x,y) and v(x,y)
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8 :: xa,xb,yc,yd    ! Positions of cell edges
	
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
			
			!Tmesh(j,i) = T0*cos(pi*(i-1./2.)*dx)*(sin(pi*(j-1./2.)*dy))
	        !umesh(j,i) = u0*((j-1./2.)*dy)*sin(pi*(i-1./2.)*dx)
			!vmesh(j,i) = v0*((i-1./2.)*dx)*cos(pi*(j-1./2.)*dy)
		END DO
	END DO
	
END SUBROUTINE init_mesh

! =======================

SUBROUTINE rk2(nmin, nmax)                        ! Subroutine for time integration (RK2)

    USE meshmod
    IMPLICIT NONE
    
    INTEGER                              :: nmin         ! Initial time index
    INTEGER                              :: nmax         ! Final time index
    !REAL*8                        :: Tw           ! Wall value
    REAL*8                               :: tn           ! Time at index n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: mesh_1       ! Interim mesh for time advance scheme
    
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FIn          ! Flux integrals of domain at time level n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI_1         ! Interim flux integrals for time advance scheme
	
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: source_mesh_n      ! Flux integrals of domain at time level n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: source_mesh_1      ! Interim flux integrals for time advance scheme
	
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: rhsn           ! Right-hand-side of energy equation at time level n
	REAL*8                               :: rhs1           ! Interim right-hand-side of energy equation at time level n

    DO n = nmin,nmax  ! Want to obtain solution at maxn, hence last loop is maxn-1
        
        ! ------Calculate T(1)------
        tn = delta_t*n                     ! Current time
        
        !CALL set_bcs(Tw, tn, mesh)         ! Update boundary conditions of solution mesh
        !CALL calc_integral(mesh, FIn, Tw)  ! Evaluate flux integrals at time level n
		
		! Calculate flux integrals
        CALL calc_integrals(Tmesh,FIn,source_mesh_n)
        
        ! Go from left to right across mesh and calculate interim values
        DO j = 1,jmax
		    DO i = 1,imax
			    rhsn(j,i) = source_mesh_n(j,i)-FIn(j,i)
                mesh_1(j,i) = Tmesh(j,i) + delta_t*rhsn(j,i)
			END DO
        END DO
        
        ! ------Calculate T(n+1)------
        tn = delta_t*(n+1)                   ! Time n+1
        
        !CALL set_bcs(Tw, tn, mesh_1)         ! Update boundary conditions of interim mesh
        
		CALL calc_integrals(mesh_1, FI_1, source_mesh_1) ! Evaluate interim flux integrals at time level n+1
        
        ! Go from left to right across mesh and calculate new solution
        DO j = 1,jmax
		    DO i = 1,imax
			    rhs1 = source_mesh_1(j,i) - FI_1(j,i)
                Tmesh(j,i) = Tmesh(j,i) + delta_t*(rhs1 + rhsn(j,i))/2.  ! Overwrite old solution array
		    END DO
        END DO
        
    END DO
    
END SUBROUTINE rk2


! =======================

SUBROUTINE calc_integrals(Tmeshn,FI,source_mesh)  ! Subroutine to evaluate flux integrals at time level n, with sources

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmeshn       ! Mesh at time level n (could be interim mesh)
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI      ! Mesh at time level n (could be interim mesh)
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: source_mesh       ! Mesh at time level n (could be interim mesh)
	
    REAL*8                               :: FIx, FIy     ! Fluxes in the x and y directions
	REAL*8                               :: flux         ! Flux at face
	REAL*8                               :: analytic     ! Analytic value of flux integral at cell centroids
    REAL*8                               :: analytic_cv  ! CV average analytic flux integral
	
	REAL*8                               :: dudx,dudy,dvdx,dvdy    ! Flux terms for u and v
	REAL*8                               :: source                 ! Source at face
	REAL*8                               :: source_analytic        ! Analytic source term
    
    
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
			!FI_ana(j,i) = analytic_cv((i-1)*dx, i*dx, (j-1)*dy, j*dy)
			
		    ! Compute source fluxes along x and y directions
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

END SUBROUTINE calc_integrals



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

! =================

FUNCTION analytic_cv(xa,xb,yc,yd)  ! Function to calculate CV-averaged analytic flux integral         

    USE meshmod
    IMPLICIT NONE

    REAL*8     :: xa,xb,yc,yd           ! Position of the cell edges
	REAL*8     :: first,second,third    ! Terms in the flux integral
	
    REAL*8     :: analytic_cv           ! Analytic flux integral
    
    ! Compute CV-averaged analytic flux integral
	first = ((-u0*T0/(2.*pi**2*dx*dy))*((sin(2.*pi*xa)-sin(2.*pi*xb)) &
	         *(-sin(pi*yc)+pi*yc*cos(pi*yc)+sin(pi*yd)-pi*yd*cos(pi*yd))))
			 
	second = ((v0*T0/(2.*pi**2*dx*dy))*((sin(2.*pi*yc)-sin(2.*pi*yd)) & 
	         *(pi*xa*sin(pi*xa)+cos(pi*xa)-pi*xb*sin(pi*xb)-cos(pi*xb))))
			 
	third = (2.*T0/(Re*Pr*dx*dy))*(sin(pi*xa)-sin(pi*xb))*(cos(pi*yd)-cos(pi*yc))
	
	analytic_cv = first + second + third
	
	RETURN

END FUNCTION analytic_cv





