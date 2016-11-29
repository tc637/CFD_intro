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
	INTEGER, PARAMETER :: jmax = 10                        ! Size of mesh (y-direction)
 
    REAL*8, PARAMETER  :: xmax = 1.0D+0                    ! Maximum x-position
	REAL*8, PARAMETER  :: ymax = 1.0D+0                    ! Maximum y-position
    REAL*8, PARAMETER  :: tmax = 1.0D+0                    ! Maximum time

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
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FI, FI_ana     ! Computed and analytic integrals at time level n
	
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
    CALL calc_integral(Tmesh)
	
	! Output computed and analytic flux integrals
	CALL output

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

SUBROUTINE calc_integral(Tmeshn)    ! Subroutine to evaluate flux integral at time level n 

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmeshn       ! Mesh at time level n (could be interim mesh)
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
			!FI_ana(j,i) = analytic_cv((i-1)*dx, i*dx, (j-1)*dy, j*dy)
		END DO
	END DO

END SUBROUTINE calc_integral


! =======================

SUBROUTINE output ! Output final solutions for plotting
    
    USE meshmod
    IMPLICIT NONE
    
    INTEGER :: ierror                                  ! For error catching
    CHARACTER(len=100) :: fname1,fname2,fname3         ! Output filenames
    CHARACTER(len=100) :: fmt                          ! Output format
        
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
	
  
	DO j=1,jmax   ! Loop over each row
        WRITE(20,"(999F15.10)") (FI(j,i), i=1,imax)             ! Implied loop to write out every column (computed flux integral)
		WRITE(30,"(999F15.10)") (FI_ana(j,i), i=1,imax)         ! Analytic flux integral
		WRITE(40,"(999F15.10)") (FI_ana(j,i)-FI(j,i), i=1,imax) ! Flux integral error
    ENDDO


END SUBROUTINE output

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


