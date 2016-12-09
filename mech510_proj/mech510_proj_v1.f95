! MECH 510, Final Project
! Timothy Chui, 37695129
! Program to numerically solve the 2D incompressible Navier-Stokes equations (with artificial compressibility)
! Time advance scheme: IE
! Space discretization scheme: 2nd-Order Centred
! Dec 2016


! =================== modules =========================

MODULE meshmod
    ! Declare variables to pass arround
	
	! Discretization parameters
    INTEGER, PARAMETER :: imax = 40                          ! Size of mesh (x-direction)
	INTEGER, PARAMETER :: jmax = 40                          ! Size of mesh (y-direction)
 
    REAL*8, PARAMETER  :: xmax = 1.0D+0                      ! Maximum x-position
	REAL*8, PARAMETER  :: ymax = 1.0D+0                       ! Maximum y-position
    REAL*8, PARAMETER  :: tmax = 40.0D+0                      ! Maximum time
		
	INTEGER, PARAMETER :: problem = 1                         ! 1 for test (square domain)
	
	! Required data
	REAL*8, PARAMETER  :: u0 = 1.0                            ! Constant for u(x,y) (test problems)
	REAL*8, PARAMETER  :: v0 = 1.0                            ! Constant for v(x,y)
    REAL*8, PARAMETER  :: P0 = 1.0                            ! Constant for P(x,y)
	REAL*8, PARAMETER  :: Re = 10.0                           ! Reynolds number
    REAL*8, PARAMETER  :: beta = 1.0                          ! Parameter for artificial compressibility
    
    ! Trig factors
	REAL*8             :: Cx,C2x,Sx,S2x                       ! cos(pix),cos(2pix),sin(pix),sin(2pix)
    REAL*8             :: Cy,C2y,Sy,S2y                       ! cos(piy),cos(2piy),sin(piy),sin(2piy)
	
	
    ! Values that do not need to be changed
    INTEGER            :: i,j,n                               ! Space and time indices
    REAL*8, PARAMETER  :: pi = 3.141592654D+0                 ! Value of pi
    REAL*8, PARAMETER  :: dx = xmax/imax                      ! Cell sizes
	REAL*8, PARAMETER  :: dy = ymax/jmax

	REAL*8              :: delta_t                            ! Time step variable
    INTEGER             :: maxn                               ! Max time step (i.e. last n index)
	
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: solmesh_n         ! Mesh of solutions at n; includes ghost cells (3rd dim: 1 = p, 2 = u, 3 = v)
	REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: FI_n,FI_ana       ! Computed and analytic flux integrals at time level n
	
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

! =================== main program =========================

PROGRAM mech510_proj   

    USE meshmod
    IMPLICIT NONE

    REAL*8   ::   tn = 0.0D+0    ! Time at index n
	
    CALL init_mesh
	CALL calc_integrals(solmesh_n, FI_n)
    CALL output_integrals
    
END PROGRAM mech510_proj


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

SUBROUTINE init_mesh_test ! Initialize mesh with centroid values of p(x,y,0), u(x,y,0), v(x,y,0) for the test case
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8 :: xpos,ypos   ! Positions of cell centroids
	
	! Pass through array: initialize Tmesh and vector field
	DO j = 0,jmax+1                 ! Loop over rows
	
	    ypos = (j - 1./2.)*dy       ! Centroid
		
	    DO i = 0, imax+1            ! Loop over columns
        
            xpos = (i - 1./2.)*dx   ! Centroid position
        
            CALL update_trig(xpos,ypos)    ! Set trig factors
            
            solmesh_n(j,i,1) = P0*Cx*Cy    ! Initialize P(x,y)
            solmesh_n(j,i,2) = u0*Sx*S2y   ! Initialize u(x,y)
            solmesh_n(j,i,3) = v0*S2x*Sy   ! Initialize v(x,y)
		
		END DO
        
	END DO
	
END SUBROUTINE init_mesh_test

! =======================

SUBROUTINE update_trig(xpos,ypos) ! Subroutine to update trig factors
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8 :: xpos,ypos   ! Positions of cell centroids
	
    Cx = cos(pi*xpos)
    Sx = sin(pi*xpos)
    
    Cy = cos(pi*ypos)
    Sy = sin(pi*ypos)
    
    C2x = cos(2.*pi*xpos)
    S2x = sin(2.*pi*xpos)
    
    C2y = cos(2.*pi*ypos)
    S2y = sin(2.*pi*ypos)
	
END SUBROUTINE update_trig

! =======================

SUBROUTINE calc_integrals(solmesh,FI)  ! Subroutine to evaluate flux integrals at time level n, with sources

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: solmesh      ! Mesh at time level n (could be interim mesh)
	REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: FI           ! Flux integral at time level n (could be interim mesh)
	
	REAL*8, DIMENSION(3)                   :: fluxa, fluxb ! Fluxes at left, right (bottom, top) faces of CV (P, u, v)
    REAL*8, DIMENSION(3)                   :: fflux, gflux ! Flux integrals along x- and y-directions (P, u, v)
    
    REAL*8                                 :: xpos,ypos    ! Positions of cell centroids


	! Calculate flux integral, i.e. (right face) - (left face)
    DO j = 1,jmax      ! Loop over rows (interior CVs)
    
        ypos = (j - 1./2.)*dy     ! Centroid
		
		DO i = 1,imax  ! Loop over columns (interior CVs)
        
            xpos = (i - 1./2.)*dx ! Centroid
            
            
            ! Calculate fluxes in x-direction
            CALL calc_fflux(i,i+1,solmesh,fluxb)
            CALL calc_fflux(i-1,i,solmesh,fluxa)
            
            ! "F" flux integral
            fflux = (fluxb-fluxa)/dx
            
            ! Calculate fluxes in y-direction
            CALL calc_gflux(j,j+1,solmesh,fluxb)
            CALL calc_gflux(j-1,j,solmesh,fluxa)
            
            ! "G" flux integral
            gflux = (fluxb-fluxa)/dy
            
            ! Sum directional integrals and switch signs
            FI(j,i,:) = -(fflux + gflux)
            
            ! For FI validation
            IF (problem == 1) THEN  
                ! Calculate analytic flux integral at CV centroid
                CALL calc_analytic(xpos,ypos)
            END IF
            
		END DO
        
	END DO

END SUBROUTINE calc_integrals

! =======================


SUBROUTINE calc_fflux(iinda,iindb,solmesh,flux)    ! Subroutine to calculate x-direction fluxes at a face

    USE meshmod
    IMPLICIT NONE
    
    INTEGER                                :: iinda, iindb     ! Indices of CVs adjacent to face along x-direction
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: solmesh          ! Mesh at time level n (could be interim mesh)
    REAL*8, DIMENSION(3)                   :: flux             ! Flux at face (P, u, v)                      

    ! P-flux along x-direction
    flux(1) = (solmesh(j,iindb,2)+solmesh(j,iinda,2))/(2.*beta)
    
    ! u-flux along x-direction
    flux(2) = (((solmesh(j,iindb,2)+solmesh(j,iinda,2))/2.)**2 + (solmesh(j,iindb,1)+solmesh(j,iinda,1))/2. &
                -1./Re*((solmesh(j,iindb,2)-solmesh(j,iinda,2))/dx))
                
    ! v-flux along x-direction
    flux(3) = (((solmesh(j,iindb,2)+solmesh(j,iinda,2))/2.)*((solmesh(j,iindb,3)+solmesh(j,iinda,3))/2.) &
               -1./Re*((solmesh(j,iindb,3)-solmesh(j,iinda,3))/dx))

END SUBROUTINE calc_fflux

! =======================


SUBROUTINE calc_gflux(jinda,jindb,solmesh,flux)    ! Subroutine to calculate y-direction fluxes at a face

    USE meshmod
    IMPLICIT NONE
    
    INTEGER                                :: jinda, jindb    ! Indices of CVs adjacent to face along y-direction
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: solmesh         ! Mesh at time level n (could be interim mesh)
    REAL*8, DIMENSION(3)                   :: flux            ! Flux at face (P, u, v)                          

    ! P-flux along y-direction
    flux(1) = (solmesh(jindb,i,3)+solmesh(jinda,i,3))/(2.*beta)
    
    ! u-flux along y-direction
    flux(2) = (((solmesh(jindb,i,2)+solmesh(jinda,i,2))/2.)*((solmesh(jindb,i,3)+solmesh(jinda,i,3))/2.) &
           -1./Re*((solmesh(jindb,i,2)-solmesh(jinda,i,2))/dy))
      
    ! v-flux along y-direction
    flux(3) = (((solmesh(jindb,i,3)+solmesh(jinda,i,3))/2.)**2 + (solmesh(jindb,i,1)+solmesh(jinda,i,1))/2. &
            -1./Re*((solmesh(jindb,i,3)-solmesh(jinda,i,3))/dy))

END SUBROUTINE calc_gflux

! =======================

SUBROUTINE calc_analytic(xpos, ypos)    ! Function to calculate analytic flux integral at a CV centroid

    USE meshmod
    IMPLICIT NONE

    REAL*8   :: xpos,ypos                     ! Position of the cell centre
	REAL*8   :: t1,t2,t3,t4,t5,t6,t7,t8,t9    ! Terms in the flux integral (across all 3 components)

    ! Update trig factors
    CALL update_trig(xpos, ypos)
    
    ! Compute analytic flux integral for P
	t1 = -pi/beta*(u0*Cx*S2y + v0*S2x*Cy)
    
    FI_ana(j,i,1) = t1
    
    ! Compute analytic flux integral for u
    t2 = P0*pi*Sx*Cy
    t3 = -u0**2*pi*S2x*S2y**2
    t4 = -u0*v0*pi*Sx*S2x*(Cy*S2y + 2.*C2y*Sy)
    t5 = -u0*(5.*(pi**2)*Sx*S2y/Re)
    
    FI_ana(j,i,2) = t2 + t3 + t4 + t5
    
    ! Compute analytic flux integral for v
    t6 = P0*pi*Cx*Sy
    t7 = -v0**2*pi*(S2x**2)*S2y
    t8 = -u0*v0*pi*Sy*S2y*(Cx*S2x + 2.*C2x*Sx)
    t9 = -v0*(5.*(pi**2)*S2x*Sy/Re)

    FI_ana(j,i,3) = t6 + t7 + t8 + t9

END SUBROUTINE calc_analytic

! =======================

SUBROUTINE output_integrals ! Output computed and analytic flux integrals for testing
    
    USE meshmod
    IMPLICIT NONE

    INTEGER :: ierror                            ! For error catching
    CHARACTER(len=100), DIMENSION(6) :: fnames   ! Output filenames
    CHARACTER(len=100) :: fmt                    ! Output format
    INTEGER :: ind,varind                        ! Filename and variable indexing
        
	! Set filename format
	fmt = "(A8,I0,A6,I0,A6,I0,A5,I0,A4)"  
	
    ! (1 -> 3, analytic), (4 -> 6, computed)
    DO ind = 1,6
        
        ! Analytic flux integrals
        IF (ind <= 3) THEN
            varind = ind   ! Variable index
            
            ! Set filename
            WRITE(fnames(ind), fmt) "FI_anal_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_var_",varind,".txt"
            
            ! Open and write to file
            OPEN(UNIT=10*ind, FILE=fnames(ind), STATUS="REPLACE", IOSTAT=ierror)
            DO j=1,jmax   ! Loop over each row
                WRITE(10*ind,"(999F15.10)") (FI_ana(j,i,varind), i=1,imax) 
            END DO
            CLOSE(UNIT=10*ind)
            
        ! Compute flux integrals
        ELSE
            varind = ind-3 ! Variable index (subtract by 3 to get 1, 2, 3)
            
            ! Set filename
            WRITE(fnames(ind), fmt) "FI_comp_",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_var_",varind,".txt"
            
            ! Open and write to file
            OPEN(UNIT=10*ind, FILE=fnames(ind), STATUS="REPLACE", IOSTAT=ierror)
            DO j=1,jmax   ! Loop over each row
                WRITE(10*ind,"(999F15.10)") (FI_n(j,i,varind), i=1,imax) 
            END DO
            CLOSE(UNIT=10*ind)
            
        END IF
        
        ! Report output filename
        WRITE(*,*) "Writing to ", fnames(ind)
         
    END DO

END SUBROUTINE output_integrals






