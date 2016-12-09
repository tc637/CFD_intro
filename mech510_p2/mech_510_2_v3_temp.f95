! MECH 510, Programming Assignment 2
! Timothy Chui, 37695129
! Program to numerically solve the one-dimensional wave equation
! Time advance scheme: RK2
! Space discretization scheme: 2nd-Order Upwind
! Nov 2016


! =================== modules =========================

MODULE meshmod
    ! Declare variables to pass arround
    ! Parameters to change
    INTEGER, PARAMETER :: meshsize = MESHSIZE          ! Size of mesh
    REAL*8, PARAMETER  :: cfl = 0.4                    ! CFL number (u*dt/dx)
    REAL*8, PARAMETER  :: xmax = 1.0D+0                ! Maximum x-position
    REAL*8, PARAMETER  :: tmax = 1.0D+0                ! Maximum time
    REAL*8, PARAMETER  :: u = 2.0D+0                   ! Wave speed
	INTEGER, PARAMETER :: ic = 1                       ! Switch to decide which IC to use (1=sine, 2=lin)
	INTEGER, PARAMETER :: bc = BC                      ! Switch to decide which 3 to use at i=3/2 (1=2nd-up,2=2nd-cen,3=1st-up)
    
    ! Values that do not need to be changed
    INTEGER            :: i, n                         ! Space and time indices
    REAL*8, PARAMETER  :: pi = 3.141592654D+0          ! Value of pi
    REAL*8, PARAMETER  :: delta_x = xmax/meshsize      ! Grid size
    REAL*8, PARAMETER  :: delta_t = cfl*delta_x/u      ! Time step from CFL condition
    INTEGER, PARAMETER :: maxn = nint(tmax/delta_t)    ! Max time step (i.e. last n index)

    REAL*8, DIMENSION(0:meshsize) :: mesh              ! Mesh of solutions; include ghost cell at 0

END MODULE meshmod


! =================== main program =========================

PROGRAM mech510_2  

    USE meshmod
    IMPLICIT NONE
    
    REAL*8                        :: tn           ! Time at index n
    REAL*8, DIMENSION(0:meshsize) :: mesh_1       ! Interim mesh for time advance scheme
    REAL*8, DIMENSION(0:meshsize) :: FIn          ! Flux integrals of domain at time level n
    REAL*8, DIMENSION(0:meshsize) :: FI_1         ! Interim flux integrals for time advance scheme
    REAL*8                        :: Tw           ! Wall value
    
    
    ! Initialize mesh
    CALL set_ics
    
    ! Loop through time and calculate flux integrals (using RK2)
    DO n = 0,maxn-1  ! Want to obtain solution at maxn, hence last loop is maxn-1
        !WRITE(*,"(999F15.10)") mesh
        
        ! ------Calculate T(1)------
        tn = delta_t*n                 ! Current time
        
        CALL set_bcs(Tw, tn, mesh)         ! Update boundary conditions of solution mesh
        CALL calc_integral(mesh, FIn, Tw)  ! Evaluate flux integrals at time level n
        
        ! Go from left to right across mesh and calculate interim values
        DO i = 1,meshsize
            mesh_1(i) = mesh(i) - delta_t*FIn(i)
        END DO
        
        ! ------Calculate T(n+1)------
        tn = delta_t*(n+1)               ! Time n+1
        
        CALL set_bcs(Tw, tn, mesh_1)         ! Update boundary conditions of interim mesh
        CALL calc_integral(mesh_1, FI_1, Tw) ! Evaluate interim flux integrals at time level n+1
        
        ! Go from left to right across mesh and calculate new solution
        DO i = 1,meshsize
            mesh(i) = mesh(i) - delta_t*(FIn(i) + FI_1(i))/2.  ! Overwrite old solution array
        END DO
        
    END DO
    
	IF (ic == 1) THEN ! Validation for sinusoidal initial condition
		tn = tmax
		CALL calc_analytic(tn)  ! Calculate and output analytic solution at time tn
	END IF
    
    CALL output             ! Output computed solution

END PROGRAM mech510_2

! =================== subroutines =========================

SUBROUTINE set_ics            ! Subroutine to set initial conditions (CV-average values)

    USE meshmod
    IMPLICIT NONE
    
    REAL*8 :: xa, xb          ! Limits of integration for cell boundaries
    
	IF (ic == 1) THEN ! Sinusoidal initial condition
		DO i = 1,meshsize         ! Initialize interior values
			xb = i*delta_x
			xa = (i-1)*delta_x
			
			! CV-averages of T(x,0) 
			mesh(i) = 1./(2.*pi*delta_x)*(cos(2.*pi*xb) - cos(2.*pi*xa))
		END DO
		
	ELSE              ! Linear initial condition
		mesh = 0.0D+0
		
		DO i = 1,nint(1./delta_x)       ! Initialize interior values
			xb = i*delta_x
			xa = (i-1)*delta_x
			
			! CV-averages of T(x,0) 
			mesh(i) = -1./delta_x*(1./2.*xb**2 - 1./2.*xa**2)
		END DO
	END IF
		
END SUBROUTINE set_ics

! ===================


SUBROUTINE set_bcs(Tw, tn, meshn)  ! Subroutine to update boundary conditions

    USE meshmod
    IMPLICIT NONE
    
    REAL*8                        :: Tw       ! Wall value
    REAL*8                        :: tn       ! Time at index n
    REAL*8, DIMENSION(0:meshsize) :: meshn    ! Mesh at time level n (could be interim mesh)
    
    ! Set wall value at time tn
    Tw = sin(4.*pi*tn)

    ! Set ghost-cell value by extrapolation
    meshn(0) = 2.0*Tw - meshn(1)

END SUBROUTINE set_bcs

! ===================

SUBROUTINE calc_integral(meshn, FI, Tw)  ! Subroutine to evaluate flux integral at time level n 

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:meshsize) :: meshn    ! Mesh at time level n (could be interim mesh)
    REAL*8, DIMENSION(0:meshsize) :: FI       ! Flux integrals of domain at time level n (could be interim integrals)
    REAL*8                        :: Tw       ! Wall value
    REAL*8                        :: flux     ! Flux at face
    
    ! Calculate flux integral, i.e. (right face) - (left face)
    DO i = 1,meshsize
        FI(i) = u*(flux(i,meshn,Tw) - flux(i-1,meshn,Tw))/delta_x ! Call function flux to evaluate i+1/2, i-1/2
    END DO

END SUBROUTINE calc_integral


! ===================

SUBROUTINE calc_analytic(tn)          ! Calculate and output analytic solution at time tn for validation

    USE meshmod
    IMPLICIT NONE
    
    REAL*8             :: tn          ! Time at index n
    REAL*8             :: analytic    ! Analytic (CV-averaged) solution for validation
    CHARACTER(len=100) :: fmt         ! Output format
    CHARACTER(len=100) :: fname       ! Output filename
    REAL*8             :: xa, xb      ! Limits of integration for cell boundaries
    INTEGER            :: ierror      ! Error catching
 
    ! No need to use array for analytic solution since no time-stepping is needed
    ! Just output cell values from left to right directly to file
    
    ! Set filename format
	fmt = "(A9,I0,A6,I0,A6,I0,A5,I0,A4,I1,A4,I1,A4)"  
    
    ! Set filename
    WRITE(fname, fmt) "analytic_",meshsize,"_xmax_",int(xmax),"_tmax_",int(tmax), &
						"_cfl_",int(cfl*1000.),"_ic_",int(ic),"_bc_",int(bc),".txt"
	
	WRITE(*,*) "Writing analytic solution to ", fname
    
    OPEN(UNIT=10, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)
    
    DO i = 1,meshsize
        xb = i*delta_x
        xa = (i-1)*delta_x
        
        ! Calculate CV-average analytic solutions and output to file
        WRITE(10,"(999F15.10)") analytic(xa,xb,tn)  ! Use function analytic
    END DO

END SUBROUTINE calc_analytic

! ===================


SUBROUTINE output  ! Output final solution for plotting
    
    USE meshmod
    IMPLICIT NONE
    
    INTEGER :: ierror            ! For error catching
    CHARACTER(len=100) :: fname  ! Output filename
    CHARACTER(len=100) :: fmt    ! Output format
        
	! Set filename format
	fmt = "(A8,I0,A6,I0,A6,I0,A5,I0,A4,I1,A4,I1,A4)"  
	
    ! Set filename
    WRITE(fname, fmt) "compute_",meshsize,"_xmax_",nint(xmax),"_tmax_",nint(tmax), &
						"_cfl_",nint(cfl*1000.),"_ic_",int(ic),"_bc_",int(bc),".txt"
    
    WRITE(*,*) "Writing output to ", fname
    
    OPEN(UNIT=20, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)
    
    DO i=1,meshsize                              ! Loop over each row
        WRITE(20,"(999F15.10)") mesh(i) ! Implied loop to write out every column
    ENDDO
    
END SUBROUTINE output

! =================== functions =========================

FUNCTION analytic(xa, xb, tn)    ! Calculate analytic solution of cell at time index n
    
    USE meshmod
    IMPLICIT NONE

    REAL*8 :: xa, xb      ! Positions of cell edges
    REAL*8 :: tn          ! Time at index n
    REAL*8 :: analytic    ! CV average of analytic solution
     
	analytic = 1./(2.*pi*delta_x)*(cos(2.*pi*(2.*tn - xb)) - cos(2.*pi*(2.*tn - xa)))

    RETURN
    
END FUNCTION analytic

! ===================

FUNCTION flux(ind, meshn, Tw) ! Subroutine to calculate at face i+1/2

    USE meshmod
    IMPLICIT NONE

    INTEGER                       :: ind      ! Index of attached CV (i) to face (i+1/2)
    REAL*8, DIMENSION(0:meshsize) :: meshn    ! Mesh at time level n (could be interim mesh)
    REAL*8                        :: Tw       ! Wall value
    REAL*8                        :: flux     ! Flux at face i+1/2
    
    IF (ind == 0) THEN
        flux = Tw   ! Flux at i = 1/2
    ELSE IF (ind == 1 .AND. bc == 2) THEN     ! 2nd-centred at i=3/2 
		flux = (meshn(ind+1) + meshn(ind))/2
	ELSE IF (ind == 1 .AND. bc == 3) THEN     ! 1st-upwind at i=3/2
		flux = (meshn(ind))     
	ELSE                                      ! 2nd-upwind for all other cases
        flux = (3.*meshn(ind) - meshn(ind-1))/2. 
    END IF
    
    RETURN

END FUNCTION flux
