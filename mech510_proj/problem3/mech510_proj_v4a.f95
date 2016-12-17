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
    INTEGER, PARAMETER :: imax = 20                            ! Size of mesh (x-direction)
	INTEGER, PARAMETER :: jmax = 20                            ! Size of mesh (y-direction)
 
    REAL*8, PARAMETER  :: xmax = 1.0D+0                        ! Maximum x-position
	REAL*8, PARAMETER  :: ymax = 1.0D+0                        ! Maximum y-position
    
	INTEGER, PARAMETER :: problem = 3                          ! 1 for FI/Jacobian tests (square domain), 2 for IE tests, 3 for real problem
	
	! Required data
	REAL*8, PARAMETER  :: u0 = 1.0                             ! Constant for u(x,y) (test problems)
	REAL*8, PARAMETER  :: v0 = 1.0                             ! Constant for v(x,y)
    REAL*8, PARAMETER  :: P0 = 1.0                             ! Constant for P(x,y)
	REAL*8, PARAMETER  :: Re = 100.0                           ! Reynolds number (test data uses Re = 100 !!!)
    REAL*8, PARAMETER  :: beta = 1.0                           ! Parameter for artificial compressibility
    REAL*8, PARAMETER  :: apres = 0.0                          ! Parameter for pressure source term to handle oscillations
    
    REAL*8, PARAMETER  :: Tw = 1.0                             ! Top wall velocity
    
    REAL*8, PARAMETER  :: over = 1.0                           ! Overrelaxation parameter
    
    ! Trig factors
	REAL*8             :: Cx,C2x,Sx,S2x                        ! cos(pix),cos(2pix),sin(pix),sin(2pix)
    REAL*8             :: Cy,C2y,Sy,S2y                        ! cos(piy),cos(2piy),sin(piy),sin(2piy)
	
	
    ! Values that do not need to be changed
    INTEGER            :: i,j,n                                ! Space and time indices
    REAL*8, PARAMETER  :: pi = 3.141592654D+0                  ! Value of pi
    REAL*8, PARAMETER  :: dx = xmax/imax                       ! Cell sizes
	REAL*8, PARAMETER  :: dy = ymax/jmax

    ! Time variables
	REAL*8, PARAMETER   :: dt = 0.05                           ! Time step variable
    REAL*8, PARAMETER   :: tol = 1.0D-6                      ! Tolerance for solution change before confirming steady-state
    INTEGER             :: maxn                               ! Max time step (i.e. last n index)

	
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)      :: solmesh_n   ! Mesh of solutions at n; includes ghost cells (3rd dim: 1 = p, 2 = u, 3 = v)
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)      :: d_solmesh   ! Change in solution of mesh
	REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)      :: FI_n,FI_ana ! Computed and analytic flux integrals at time level n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3,3,6)  :: jacobs_n    ! Scaled Jacobians at time level n (Ax,Bx,Cx,Ay,By,Cy)
	
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
	
    ! Test flux integrals and Jacobian calculations
    IF (problem == 1) THEN
        CALL calc_integrals(solmesh_n, FI_n, jacobs_n)
        CALL output_integrals
        CALL test_jacobian
     
    ! Stationary walls tests  
    ELSE IF (problem == 2) THEN
        CALL set_bcs(solmesh_n)
        CALL ie(tn)
        CALL output_mesh(tn)
    
    ! Moving top wall problem   
    ELSE
        
        CALL set_bcs(solmesh_n)
        !CALL calc_integrals(solmesh_n, FI_n, jacobs_n)
        !CALL output_integrals
        !CALL output_mesh(tn)
        CALL ie(tn)
        CALL output_mesh(tn)
        
    END IF
    
END PROGRAM mech510_proj


! =================== subroutines =========================

SUBROUTINE init_mesh ! Initialize mesh with centroid values of p(x,y,0), u(x,y,0), v(x,y,0) f
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8 :: xpos,ypos   ! Positions of cell centroids
	
	! Pass through array: initialize Tmesh and vector field
	DO j = 0,jmax+1                 ! Loop over rows
	
        IF (problem == 3) THEN
            ypos = -(j - 1./2.)*dy      ! Top-left corner is (0,0)
        ELSE
            ypos = (j - 1./2.)*dy       ! Bottom-left corner is (0,0)
        END IF 
		
	    DO i = 0, imax+1            ! Loop over columns
        
            xpos = (i - 1./2.)*dx   ! Centroid position
        
            CALL update_trig(xpos,ypos)    ! Set trig factors
            
            solmesh_n(j,i,1) = P0*Cx*Cy    ! Initialize P(x,y)
            solmesh_n(j,i,2) = u0*Sx*S2y   ! Initialize u(x,y)
            solmesh_n(j,i,3) = v0*S2x*Sy   ! Initialize v(x,y)
		
		END DO
        
	END DO
	
END SUBROUTINE init_mesh

! =======================

SUBROUTINE set_bcs(solmesh)  ! Subroutine to set ghost cells to enforce boundary conditions
    
    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3) :: solmesh         ! Mesh at time level n (could be interim mesh)
    
    ! Stationary walls
    
    ! Assume Neumann for pressure, Dirichlet for u and v (no-slip)
    DO j=0,jmax+1  ! Update left and right BCs
        ! Left wall (stationary)
        solmesh(j,0,1) = solmesh(j,1,1)  ! P (Neumann)
        solmesh(j,0,2) = -solmesh(j,1,2) ! u (Dirichlet)
        solmesh(j,0,3) = -solmesh(j,1,3) ! v (Dirchlet)
        
        ! Right wall (stationary)
        solmesh(j,imax+1,1) = solmesh(j,imax,1)  ! P (Neumann)
        solmesh(j,imax+1,2) = -solmesh(j,imax,2) ! u (Dirichlet)
        solmesh(j,imax+1,3) = -solmesh(j,imax,3) ! v (Dirchlet)
    END DO

    DO i=0,imax+1  ! Update bottom and top BCs
    
        IF (problem == 3) THEN ! Upside down box, with (0,0) on top-left corner
        
            ! Top wall (moving)
            solmesh(0,i,1) = solmesh(1,i,1)  ! P (Neumann)
            solmesh(0,i,2) = 2.0*Tw-solmesh(1,i,2) ! u (Dirichlet)
            solmesh(0,i,3) = -solmesh(1,i,3) ! v (Dirchlet)
            
            ! Bottom wall(stationary)
            solmesh(jmax+1,i,1) = solmesh(jmax,i,1)  ! P (Neumann)
            solmesh(jmax+1,i,2) = -solmesh(jmax,i,2) ! u (Dirichlet)
            solmesh(jmax+1,i,3) = -solmesh(jmax,i,3) ! v (Dirchlet)

        
        ELSE
            ! Bottom wall (stationary)
            solmesh(0,i,1) = solmesh(1,i,1)  ! P (Neumann)
            solmesh(0,i,2) = -solmesh(1,i,2) ! u (Dirichlet)
            solmesh(0,i,3) = -solmesh(1,i,3) ! v (Dirchlet)
            
            ! Top wall (stationary)
            solmesh(jmax+1,i,1) = solmesh(jmax,i,1)  ! P (Neumann)
            solmesh(jmax+1,i,2) = -solmesh(jmax,i,2) ! u (Dirichlet) 
            solmesh(jmax+1,i,3) = -solmesh(jmax,i,3) ! v (Dirchlet)
        
        END IF
        
    END DO

END SUBROUTINE set_bcs    
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

SUBROUTINE ie(tn)              ! Subroutine for time integration (IE); rows first

    USE meshmod
    IMPLICIT NONE
    
    REAL*8                                :: tn            ! Time at index n

    INTEGER, PARAMETER                    :: maxsize = 200 ! Max possible size of array
    	
	REAL*8, DIMENSION(3,0:maxsize,0:maxsize) :: d_soltilda ! Interim dU (rows first) (3 X j X i)
	
    REAL*8, DIMENSION(3,0:maxsize) :: rhs                  ! RHS of approximate factorization equation
	
    REAL*8, DIMENSION(3,3,3,0:maxsize)      :: lhsx        ! LHS of first equation in approx fact
	REAL*8, DIMENSION(3,3,3,0:maxsize)      :: lhsy        ! LHS of second equation in approx fact
    INTEGER                                 :: iss         ! is steady state? (determine from tolerance; 0 = No, 1 = Yes)
    
    INTEGER                                 :: jstart, jend, jinc ! j-indices; need to reverse order if top-left is (0,0)
	
	n = 0
    iss = 0
    
    IF (problem == 3) THEN
        CALL output_converge(1)
        jstart = jmax
        jend = 1
        jinc = -1
    ELSE
        jstart = 1
        jend = jmax
        jinc = 1
    END IF
    
    WRITE(*,*) jinc
	DO WHILE (iss == 0)  ! Want to obtain solution at maxn, hence last loop is maxn-1
        
        ! ------Calculate T(1)------
        n = n + 1
        
        tn = dt*n                       ! Current time
        
        
        CALL set_bcs(solmesh_n)                             ! Update boundary conditions of solution mesh
        CALL calc_integrals(solmesh_n, FI_n, jacobs_n)      ! Calculate flux integrals
		
        
        ! Go across columns from row to row; first part of approx factorization
        DO j = jstart,jend,jinc
            !WRITE(*,*) j
        
            ! Set RHS with appropriate BCs and flux integrals
            DO i = 0,imax+1
            
                IF (i == 0 .AND. i == imax+1) THEN
                    rhs(:,i) = (/0.0,0.0,0.0/)         ! BCs
                ELSE
                    rhs(:,i) = FI_n(j,i,:)*dt          ! Flux integrals with time step
                END IF
            END DO
                    
            ! Set LHS for approximate factorization (first part)
            lhsx = 0.0D+0

            CALL set_lhs(lhsx,j,imax,maxsize,1)             ! Along constant rows (j)
            CALL SolveBlockTri(lhsx, rhs, imax+2, maxsize)  ! Call subroutine to solve block tridiag problem
            
            d_soltilda(:,j,:) = rhs                         ! Assign solution to new array for clarity
            
            !WRITE(*,"(999F20.15)") (d_soltilda(:,j,i),i=1,imax)
            
        END DO
        
        
		
		! Go across rows from column to column; second part of approx factorization

		DO i = 1,imax
        
            ! Set RHS with appropriate BCs
			d_soltilda(:,0,i) = (/0.0,0.0,0.0/)
			d_soltilda(:,jmax+1,i) = (/0.0,0.0,0.0/)
            
            ! Set LHS for approximate factorization (second part)
	    	lhsy = 0.0D+0
            
			CALL set_lhs(lhsy,i,jmax,maxsize,2)             ! Along constant columns (i)
			CALL SolveBlockTri(lhsy, d_soltilda(:,:,i), jmax+2, maxsize)
			
            ! Write to delta_U array
            DO j = jstart,jend,jinc
                d_solmesh(j,i,:) = (d_soltilda(:,j,i))
            END DO
            
        END DO
		

		! Update solution
        solmesh_n = solmesh_n + over*d_solmesh

        ! Check if solution has reached steady-state
        CALL reach_ss(iss)
        		
    END DO
    
    
    IF (problem == 3) THEN
        CALL output_converge(2)
    END IF
    
    WRITE(*,*) "Beta = ",beta," A = ",apres,"Over = ",over
    WRITE(*,*) "Iteration took ", n, " time steps to reach steady-state." 
	WRITE(*,*)
    
END SUBROUTINE ie

! =======================

SUBROUTINE reach_ss(iss) ! Subroutine to check if solution has converged to steady-state
    USE meshmod
    IMPLICIT NONE
    
    INTEGER              :: iss     ! is steady state? (determine from tolerance; 0 = No, 1 = Yes)
    
    REAL*8               :: maxdP, maxdu, maxdv  ! abs(max(dP)), abs(max(du)), abs(max(dv)) (Linf-norms)
    REAL*8               :: l2dP, l2du, l2dv     ! L2-norms of dP, du and dv
    REAL*8, DIMENSION(3) :: max_array            ! Array of maximum changes (in either L2 or Linf)

    IF (problem == 3) THEN ! Use L2 norm to determine convergence
        l2dP = (sum(d_solmesh(:,:,1)**2)/(float(imax)*float(jmax)))**(1./2.)
        l2du = (sum(d_solmesh(:,:,2)**2)/(float(imax)*float(jmax)))**(1./2.)
        l2dv = (sum(d_solmesh(:,:,3)**2)/(float(imax)*float(jmax)))**(1./2.)
        
        
        WRITE(11*1,"(999F20.18)") l2dP
        WRITE(11*2,"(999F20.18)") l2du
        WRITE(11*3,"(999F20.18)") l2dv
        
        max_array = (/l2dP,l2du,l2dv/) 
        
    ELSE ! Use Linf norm to determine convergence
        maxdP = maxval(abs(d_solmesh(:,:,1)))
        maxdu = maxval(abs(d_solmesh(:,:,2)))
        maxdv = maxval(abs(d_solmesh(:,:,3)))
        max_array = (/maxdP,maxdu,maxdv/) 
    END IF
    
    
    ! Quit out of ie loop if max change is smaller than tolerance
    IF (maxval(max_array) < tol) THEN
        iss = 1
    ELSE ! Otherwise continue with iteration
        iss = 0  
    END IF

END SUBROUTINE reach_ss

! =======================

SUBROUTINE set_lhs(lhs,nd,vsize,maxsize,dir)     ! Initialize LHS for approximate factorization, for IE time advance
	
	USE meshmod
    IMPLICIT NONE
	

	INTEGER                              :: vsize    ! Size of column (in x, first call) or row (in y, second call)
    INTEGER                              :: maxsize  ! Maximum array size allocation
    REAL*8, DIMENSION(3,3,3,0:maxsize)   :: lhs      ! LHS matrix
    INTEGER                              :: nd       ! Positional index (block index)
    INTEGER                              :: dir      ! Direction, either 1 (x) or 2 (y)
    
	INTEGER                              :: knd      ! Index for going down tridiagonal matrix
    REAL*8, DIMENSION(3,3)               :: bc1 = 0.0D+0, bc2 = 0.0D+0, bc3 = 0.0D+0 ! Boundary condition arrays

    ! Blocks on main diagonal always the identity matrix
    bc2(1,1) = 1.0
    bc2(2,2) = 1.0
    bc2(3,3) = 1.0

    DO knd = 0,vsize+1      ! Go down block tridiagonal matrix
    
        IF (dir == 1) THEN  ! Along constant rows
        
            IF (knd == 0) THEN   ! Left wall
                
                bc3(1,1) = -1.0  ! Neumann for pressure
                bc3(2,2) = 1.0   ! Dirichlet for u
                bc3(3,3) = 1.0   ! Dirichlet for v
                
                lhs(:,:,2,knd) = bc2
                lhs(:,:,3,knd) = bc3
                
            ELSE IF (knd == vsize+1) THEN   ! Right wall
            
                bc1(1,1) = -1.0  ! Neumann for pressure
                bc1(2,2) = 1.0   ! Dirichlet for u
                bc1(3,3) = 1.0   ! Dirichlet for v
                
                lhs(:,:,1,knd) = bc1
                lhs(:,:,2,knd) = bc2

            ELSE
                ! Interior tridiagonals
                lhs(:,:,1,knd) = jacobs_n(nd,knd,:,:,1)*dt        ! dt*Ax
                lhs(:,:,2,knd) = bc2 + jacobs_n(nd,knd,:,:,2)*dt  ! I + dt*Bx 
                lhs(:,:,3,knd) = jacobs_n(nd,knd,:,:,3)*dt        ! dt*Cx
                
            END IF
            
        ELSE     ! Along constant columns
        
            IF (knd == 0) THEN   ! Bottom wall
                
                bc3(1,1) = -1.0  ! Neumann for pressure
                bc3(2,2) = 1.0   ! Dirichlet for u
                bc3(3,3) = 1.0   ! Dirichlet for v
                
                lhs(:,:,2,knd) = bc2
                lhs(:,:,3,knd) = bc3
                
            ELSE IF (knd == vsize+1) THEN  ! Top wall
            
                bc1(1,1) = -1.0  ! Neumann for pressure
                bc1(2,2) = 1.0   ! Dirichlet for u
                bc1(3,3) = 1.0   ! Dirichlet for v
                
                lhs(:,:,1,knd) = bc1
                lhs(:,:,2,knd) = bc2

            ELSE
                ! Interior tridiagonals                   
                lhs(:,:,1,knd) = jacobs_n(knd,nd,:,:,4)*dt        ! dt*Ay
                lhs(:,:,2,knd) = bc2 + jacobs_n(knd,nd,:,:,5)*dt  ! I + dt*By
                lhs(:,:,3,knd) = jacobs_n(knd,nd,:,:,6)*dt        ! dt*Cy
                
            END IF
            
        END IF
        
    END DO
    
END SUBROUTINE set_lhs
! =======================

SUBROUTINE calc_integrals(solmesh,FI,jacobs)  ! Subroutine to evaluate flux integrals and Jacobians at time level n

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)       :: solmesh      ! Mesh at time level n (could be interim mesh)
	REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)       :: FI           ! Flux integral at time level n (could be interim mesh)
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3,3,6)   :: jacobs       ! Scaled Jacobians at time level n and position i,j (Ax,Bx,Cx,Ay,By,Cy)
	
	REAL*8, DIMENSION(3)                   :: fluxa, fluxb ! Fluxes at left, right (bottom, top) faces of CV (P, u, v)
    REAL*8, DIMENSION(3)                   :: fflux, gflux ! Flux integrals along x- and y-directions (P, u, v)
    
    ! Jacobian matrices (for F and G)
    REAL*8, DIMENSION(3,3)                 :: jacobpa,jacobpb,jacobma,jacobmb ! (dF_i+0.5/dU_i,dF_i+0.5/dU_i+1,dF_i-0.5/dU_i-1,dF_i-0.5/dU_i)
    
    REAL*8                                 :: xpos,ypos    ! Positions of cell centroids
    
    INTEGER                                :: jstart, jend, jinc ! j-indices; need to reverse order if top-left is (0,0)


    IF (problem == 3) THEN
        jstart = jmax
        jend = 1
        jinc = -1
    ELSE
        jstart = 1
        jend = jmax
        jinc = 1
    END IF
    
    
	! Calculate flux integral, i.e. (right face) - (left face)
    DO j = jstart,jend,jinc     ! Loop over rows (interior CVs)
    
        ypos = (j - 1./2.)*dy     ! Centroid
		
		DO i = 1,imax  ! Loop over columns (interior CVs)
        
            xpos = (i - 1./2.)*dx ! Centroid
            
            ! ------ Handle flux integrals ------
            
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
            
            ! ------ Handle Jacobians ------
        
            ! Calculate F Jacobians
            CALL calc_fjacob(i,i+1,solmesh,jacobpa,jacobpb)  ! Indices i+0.5 (i,i+1)
            CALL calc_fjacob(i-1,i,solmesh,jacobma,jacobmb)  ! Indices i-0.5 (i-1,i)
            
            jacobs(j,i,:,:,1) = -1./dx*(jacobma)                 ! Ax
            jacobs(j,i,:,:,2) = 1./dx*(jacobpa - jacobmb)        ! Bx
            jacobs(j,i,:,:,3) = 1./dx*(jacobpb)                  ! Cx
            
            
            ! Calculate G Jacobians
            
            IF (problem == 3) THEN
                CALL calc_gjacob(j,j-1,solmesh,jacobpa,jacobpb)  ! Indices j+0.5 (j,j+1)
                CALL calc_gjacob(j+1,j,solmesh,jacobma,jacobmb)  ! Indices j-0.5 (j-1,j)
            ELSE
                CALL calc_gjacob(j,j+1,solmesh,jacobpa,jacobpb)  ! Indices j+0.5 (j,j+1)
                CALL calc_gjacob(j-1,j,solmesh,jacobma,jacobmb)  ! Indices j-0.5 (j-1,j)
            END IF
            
            jacobs(j,i,:,:,4) = -1./dy*(jacobma)                 ! Ay
            jacobs(j,i,:,:,5) = 1./dy*(jacobpa - jacobmb)        ! By
            jacobs(j,i,:,:,6) = 1./dy*(jacobpb)                  ! Cy
            
            ! For FI validation
            IF (problem == 1) THEN  
                ! Calculate analytic flux integral at CV centroids
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
        
        ! Output analytic flux integrals
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
            
        ! Output computed flux integrals
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



! =======================

SUBROUTINE output_converge(oc)              ! Output convergence info for problem 3
    
    USE meshmod
    IMPLICIT NONE

    INTEGER :: oc                                ! Open (oc = 1) or close (oc = 2)
    INTEGER :: ierror                            ! For error catching
    CHARACTER(len=100), DIMENSION(6) :: fnames   ! Output filenames
    CHARACTER(len=100) :: fmt                    ! Output format
    INTEGER :: ind,varind                        ! Filename and variable indexing
        
	! Set filename format
	fmt = "(A8,I0,A6,I0,A6,I0,A5,I0,A4,I0,A4,I0,A4)"  
    
    IF (oc == 1) THEN  ! If open
        DO ind = 1,3
            
            ! Output L2 norms
            varind = ind   ! Variable index
            
            ! Set filename
            WRITE(fnames(ind), fmt) "converg_",jmax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_var_",varind, &
                "_dt_",nint(dt*100),"_Tw_",nint(Tw*100),".txt"
            
            ! Open and write to file
            OPEN(UNIT=11*ind, FILE=fnames(ind), STATUS="REPLACE", IOSTAT=ierror)
    
            ! Report output filename
            WRITE(*,*) "Opening convergence file ", fnames(ind)
             
        END DO
        
    ELSE
        DO ind = 1,3
        
            ! Close file
            CLOSE(UNIT=11*ind)
             
        END DO
    
    END IF

END SUBROUTINE output_converge


! =======================

SUBROUTINE output_mesh(tn) ! Output final solutions for plotting

    USE meshmod
    IMPLICIT NONE

    REAL*8   ::   tn                                                   ! Time at index n
    INTEGER :: ierror                                                  ! For error catching
    CHARACTER(len=200) :: fname                                        ! Output filenames
    CHARACTER(len=100) :: fmt                                          ! Output format

    ! Set filename format
    fmt = "(A8,I0,A6,I0,A6,I0,A3,I0,A3,I0,A1,I0,A1,I0,A1,I0,A1,I0,A4)"

    ! Set filenames (pressure)
    WRITE(fname, fmt) "pmesh_cn",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_t_",nint(tn),"_p_",problem,"_", &
       nint(beta*100),"_",nint(apres*100),"_",nint(over*100),"_",nint(Tw*100),".txt"
    WRITE(*,*) "Writing output to ", fname
    OPEN(UNIT=90, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)
    
   ! Set filenames (u velocity)
    WRITE(fname, fmt) "umesh_cn",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_t_",nint(tn),"_p_",problem,"_", &
       nint(beta*100),"_",nint(apres*100),"_",nint(over*100),"_",nint(Tw*100),".txt"
    WRITE(*,*) "Writing output to ", fname
    OPEN(UNIT=91, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)
    
   ! Set filenames (v velocity)
    WRITE(fname, fmt) "vmesh_cn",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_t_",nint(tn),"_p_",problem,"_", &
       nint(beta*100),"_",nint(apres*100),"_",nint(over*100),"_",nint(Tw*100),".txt"
    WRITE(*,*) "Writing output to ", fname
    OPEN(UNIT=92, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)

   DO j=1,jmax   ! Loop over each row
        WRITE(90,"(999F20.15)") (solmesh_n(j,i,1), i=1,imax)     ! Implied loop to write out every column
        WRITE(91,"(999F20.15)") (solmesh_n(j,i,2), i=1,imax)     ! Implied loop to write out every column
        WRITE(92,"(999F20.15)") (solmesh_n(j,i,3), i=1,imax)     ! Implied loop to write out every column 
   ENDDO

   CLOSE(UNIT=90)
   CLOSE(UNIT=91)
   CLOSE(UNIT=92)
    
END SUBROUTINE output_mesh

! =======================

SUBROUTINE calc_fjacob(iinda, iindb, solmesh, fui, fuip1) ! Subroutine to compute F Jacobians
    
    USE meshmod
    IMPLICIT NONE
    
    INTEGER                                  :: iinda, iindb     ! Indices of CVs adjacent to face along x-direction
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)   :: solmesh          ! Mesh at time level n (could be interim mesh)
    REAL*8, DIMENSION(3,3)                   :: fui, fuip1       ! dF/dU_i,j, dF/dU_i+1,j                  

    ! Initialize Jacobians with 0s
    fui = 0.0D+0
    fuip1 = 0.0D+0
    
    ! Compute non-zero derivatives for dF/dU_i,j
    fui(1,2) = 1./(2.*beta)
    fui(2,1) = 1./2.
    fui(2,2) = (solmesh(j,iindb,2) + solmesh(j,iinda,2))/2. + 1./(Re*dx)
    fui(3,2) = (solmesh(j,iindb,3) + solmesh(j,iinda,3))/4.
    fui(3,3) = (solmesh(j,iindb,2) + solmesh(j,iinda,2))/4. + 1./(Re*dx)
    
    ! Compute non-zero derivatives for dF/dU_i+1,j
    fuip1(1,2) = 1./(2.*beta)
    fuip1(2,1) = 1./2.
    fuip1(2,2) = (solmesh(j,iindb,2) + solmesh(j,iinda,2))/2. - 1./(Re*dx)
    fuip1(3,2) = (solmesh(j,iindb,3) + solmesh(j,iinda,3))/4.
    fuip1(3,3) = (solmesh(j,iindb,2) + solmesh(j,iinda,2))/4. - 1./(Re*dx)
    
END SUBROUTINE calc_fjacob

! =======================

SUBROUTINE calc_gjacob(jinda, jindb, solmesh, fuj, fujp1) ! Subroutine to compute G Jacobians
    
    USE meshmod
    IMPLICIT NONE
    
    INTEGER                                  :: jinda, jindb     ! Indices of CVs adjacent to face along y-direction
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)   :: solmesh          ! Mesh at time level n (could be interim mesh)
    REAL*8, DIMENSION(3,3)                   :: fuj, fujp1       ! dG/dU_i,j, dG/dU_i,j+1                  

    ! Initialize Jacobians with 0s
    fuj = 0.0D+0
    fujp1 = 0.0D+0
    
    ! Compute non-zero derivatives for dG/dU_i,j
    fuj(1,3) = 1./(2.*beta)
    fuj(2,2) = (solmesh(jindb,i,3)+solmesh(jinda,i,3))/4. + 1./(Re*dy)
    fuj(2,3) = (solmesh(jindb,i,2)+solmesh(jinda,i,2))/4.
    fuj(3,1) = 1./2.
    fuj(3,3) = (solmesh(jindb,i,3)+solmesh(jinda,i,3))/2. + 1./(Re*dy)
     
    ! Compute non-zero derivatives for dG/dU_i,j+1
    fujp1(1,3) = 1./(2.*beta)
    fujp1(2,2) = (solmesh(jindb,i,3)+solmesh(jinda,i,3))/4. - 1./(Re*dy)
    fujp1(2,3) = (solmesh(jindb,i,2)+solmesh(jinda,i,2))/4.
    fujp1(3,1) = 1./2.
    fujp1(3,3) = (solmesh(jindb,i,3)+solmesh(jinda,i,3))/2. - 1./(Re*dy)
    
END SUBROUTINE calc_gjacob



! =======================


SUBROUTINE test_jacobian ! Subroutine to test correctness of flux Jacobian

    USE meshmod
    IMPLICIT NONE
    
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)      :: solmesh_old     ! Copy of old solution mesh (at time 0)
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)      :: FI_old          ! Copy of old flux integral (at time 1)
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3,3,6)  :: jacobs_old      ! Scaled Jacobians at time level n (Ax,Bx,Cx,Ay,By,Cy)
    REAL*8, DIMENSION(0:jmax+1,0:imax+1,3)      :: lhs,rhs         ! LHS and RHS of equation (1)
    INTEGER*8                                   :: varind          ! Variable index
    
    ! Set delta_U (make 0 except for centre of array)
    d_solmesh = 0.0D+0
    d_solmesh(int(jmax/2),int(imax/2),:) = 1.0D-6
    
    solmesh_old = solmesh_n  ! Retain copy of solution at time n = 0
    FI_old = FI_n            ! Retain copy of flux integral at time n = 0
    jacobs_old = jacobs_n    ! Retain Jacobians at time n = 0
    
    solmesh_n = solmesh_old(:,:,:) + d_solmesh(:,:,:)  ! Update solution to time n = 1

    CALL calc_integrals(solmesh_n, FI_n, jacobs_n)     ! Compute new flux integral
     
    lhs = FI_n - FI_old                                ! LHS of equation (1)

    ! Compute RHS of equation (1), using Jacobians at time n = 0
    DO j = 1,jmax      ! Loop over rows (interior CVs)
        DO i = 1,imax  ! Loop over columns (interior CVs)
            DO varind = 1,3  ! Loop over variables
            
                ! Matrix multiply for each term; flip sign to match flux integral
                rhs(j,i,varind) = -(sum(jacobs_old(j,i,varind,:,1)*d_solmesh(j,i-1,:)) &
                                    + sum(jacobs_old(j,i,varind,:,2)*d_solmesh(j,i,:)) &
                                    + sum(jacobs_old(j,i,varind,:,3)*d_solmesh(j,i+1,:)) &
                                    + sum(jacobs_old(j,i,varind,:,4)*d_solmesh(j-1,i,:)) &
                                    + sum(jacobs_old(j,i,varind,:,5)*d_solmesh(j,i,:)) &
                                    + sum(jacobs_old(j,i,varind,:,6)*d_solmesh(j+1,i,:)))
                                    
            END DO
        END DO
	END DO
        
    ! Output error between LHS and RHS for cells with changed flux integrals
    WRITE(*,*) "Error for P"
    
    DO j=int(jmax/2-1),int(jmax/2+1)
        WRITE(*,"(9999F30.20)") (lhs(j,i,1)-rhs(j,i,1), i=int(imax/2-1),int(imax/2+1))
    END DO
    
    WRITE(*,*)
    WRITE(*,*) "Error for u"
    
    DO j=int(jmax/2-1),int(jmax/2+1)
        WRITE(*,"(9999F30.20)") (lhs(j,i,2)-rhs(j,i,2), i=int(imax/2-1),int(imax/2+1))
    END DO
    
    WRITE(*,*)
    WRITE(*,*) "Error for v"
    
    DO j=int(jmax/2-1),int(jmax/2+1)
        WRITE(*,"(9999F30.20)") (lhs(j,i,3)-rhs(j,i,3), i=int(imax/2-1),int(imax/2+1))
    END DO
    
    WRITE(*,*)
    

END SUBROUTINE test_jacobian



! =======================

! Code for block-tridiagonal matrix problems, written by Carl Ollivier-Gooch

!     Code for solving block-tridiagonal matrix problems, using the
!     Thomas algorithm.  The only subroutine in here that you'll -need-
!     to call is SolveThomas, although things like Add3x3 or AddVec
!     might be useful, too.

!     LHS array is sized as (3,3,3,*).  The first two indices identify
!     the element within a block; the first index is the row and the
!     second is the column of the Jacobian matrix.  The third index
!     tells which block it is: 1 is below the main diagonal, 2 is on
!     the main diagonal, 3 is above the main diagonal.  The last index
!     tells which block row you're looking at (the i or j index from
!     the discretization).

!     RHS array is (3, *).  The first index tells which element of the
!     solution vector, and the second is the block row.

!     Before linking this with your own code, you'll want to remove the
!     main program included here as a test.

      subroutine SpewMatrix(Source)
      double precision Source(3,3)
      write(6,*) Source(1,1), Source(2,1), Source(3,1)
      write(6,*) Source(1,2), Source(2,2), Source(3,2)
      write(6,*) Source(1,3), Source(2,3), Source(3,3)
      return
      end

      subroutine SpewVector(Source)
      double precision Source(3)
      write(6,*) Source(1), Source(2), Source(3)
      return
      end

      subroutine CopyVec(Source, Target)
      double precision Source(3), Target(3)
      Target(1) = Source(1)
      Target(2) = Source(2)
      Target(3) = Source(3)
      return
      end

      subroutine Copy3x3(Source, Target)
      double precision Source(3,3), Target(3,3)

      Target(1,1) = Source(1,1)
      Target(1,2) = Source(1,2)
      Target(1,3) = Source(1,3)
      
      Target(2,1) = Source(2,1)
      Target(2,2) = Source(2,2)
      Target(2,3) = Source(2,3)
      
      Target(3,1) = Source(3,1)
      Target(3,2) = Source(3,2)
      Target(3,3) = Source(3,3)
      return
      end

      subroutine Mult3x3(A, B, C)
      double precision A(3,3), B(3,3), C(3,3)

      C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1) 
      C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2) 
      C(1,3) = A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3) 

      C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1)
      C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
      C(2,3) = A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)

      C(3,1) = A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1)
      C(3,2) = A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2)
      C(3,3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
      return
      end

      subroutine MultVec(A, Vec, Result)
      double precision A(3,3), Vec(3), Result(3)

      Result(1) = A(1,1)*Vec(1) + A(1,2)*Vec(2) + A(1,3)*Vec(3) 
      Result(2) = A(2,1)*Vec(1) + A(2,2)*Vec(2) + A(2,3)*Vec(3) 
      Result(3) = A(3,1)*Vec(1) + A(3,2)*Vec(2) + A(3,3)*Vec(3)
      return
      end

      subroutine Add3x3(A, B, Factor, C)
      double precision A(3,3), B(3,3), Factor, C(3,3)

      C(1,1) = A(1,1) + Factor * B(1,1)
      C(1,2) = A(1,2) + Factor * B(1,2)
      C(1,3) = A(1,3) + Factor * B(1,3)
      
      C(2,1) = A(2,1) + Factor * B(2,1)
      C(2,2) = A(2,2) + Factor * B(2,2)
      C(2,3) = A(2,3) + Factor * B(2,3)
      
      C(3,1) = A(3,1) + Factor * B(3,1)
      C(3,2) = A(3,2) + Factor * B(3,2)
      C(3,3) = A(3,3) + Factor * B(3,3)
      return
      end

      subroutine AddVec(A, B, Factor, C)
      double precision A(3), B(3), Factor, C(3)
      C(1) = A(1) + Factor * B(1)
      C(2) = A(2) + Factor * B(2)
      C(3) = A(3) + Factor * B(3)
      return
      end

      subroutine Invert3x3(Block, Inverse)
      double precision Block(3,3), Inverse(3,3)
      double precision DetInv

      DetInv = 1. / (+ Block(1,1)*Block(2,2)*Block(3,3) &
           + Block(1,2)*Block(2,3)*Block(3,1) &
           + Block(1,3)*Block(2,1)*Block(3,2) & 
           - Block(1,3)*Block(2,2)*Block(3,1) &
           - Block(1,2)*Block(2,1)*Block(3,3) &
           - Block(1,1)*Block(2,3)*Block(3,2))

!     Expand by minors to compute the inverse 
      Inverse(1,1) = + DetInv * (Block(2,2)*Block(3,3)-Block(3,2)*Block(2,3)) 
      Inverse(2,1) = - DetInv * (Block(2,1)*Block(3,3)-Block(3,1)*Block(2,3)) 
      Inverse(3,1) = + DetInv * (Block(2,1)*Block(3,2)-Block(3,1)*Block(2,2)) 
      Inverse(1,2) = - DetInv * (Block(1,2)*Block(3,3)-Block(3,2)*Block(1,3)) 
      Inverse(2,2) = + DetInv * (Block(1,1)*Block(3,3)-Block(3,1)*Block(1,3)) 
      Inverse(3,2) = - DetInv * (Block(1,1)*Block(3,2)-Block(3,1)*Block(1,2)) 
      Inverse(1,3) = + DetInv * (Block(1,2)*Block(2,3)-Block(2,2)*Block(1,3)) 
      Inverse(2,3) = - DetInv * (Block(1,1)*Block(2,3)-Block(2,1)*Block(1,3)) 
      Inverse(3,3) = + DetInv * (Block(1,1)*Block(2,2)-Block(2,1)*Block(1,2)) 

      return
      end

      subroutine SolveBlockTri(LHS, RHS, NRows, MaxSize)
      double precision LHS(3, 3, 3, MaxSize), RHS(3, MaxSize)
      integer NRows, MaxSize

      integer j
      double precision Inv(3,3), TempMat(3,3), TempVec(3)
      double precision TempMat2(3,3), TVec2(3)

      do 10 j = 1, NRows-1
!     Compute the inverse of the main block diagonal.
         call Invert3x3(LHS(1,1,2,j), Inv)
!     Scale the right-most block diagonal by the inverse.
         call Mult3x3(Inv, LHS(1,1,3,j), TempMat)
         call Copy3x3(TempMat, LHS(1,1,3,j))
!     Scale the right-hand side by the inverse.
         call MultVec(Inv, RHS(1,j), TempVec)
         call CopyVec(TempVec, RHS(1,j))
         
!     Left-multiply the jth row by the sub-diagonal on the j+1st row and
!     subtract from the j+1st row.  This involves the super-diagonal
!     term and the RHS of the jth row.
         
!     First the LHS manipulation
         call Mult3x3(LHS(1,1,1,j+1), LHS(1,1,3,j), TempMat)
         call Add3x3(LHS(1,1,2,j+1), TempMat, -1.d0, TempMat2)
         call Copy3x3(TempMat2, LHS(1,1,2,j+1))

!     Now the RHS manipulation 
         call MultVec(LHS(1,1,1,j+1), RHS(1,j), TempVec)
         call AddVec(RHS(1,j+1), TempVec, -1.d0, TVec2)
         call CopyVec(TVec2, RHS(1,j+1))
 10   end do
      
!     Done with forward elimination loop 

!     Compute the inverse of the last main block diagonal.
      j = NRows
      
      call Invert3x3(LHS(1,1,2,j), Inv)

!     Scale the right-hand side by the inverse. 
      
      call MultVec(Inv, RHS(1,j), TempVec)
      call CopyVec(TempVec, RHS(1,j))
      
!     Now do the back-substitution. 
      do 20 j = NRows - 1, 1, -1
!     Matrix-vector multiply and subtract.
         RHS(1,j) = RHS(1,j) - (LHS(1,1,3,j)*RHS(1,j+1) &
             + LHS(1,2,3,j)*RHS(2,j+1) &
             + LHS(1,3,3,j)*RHS(3,j+1))
         RHS(2,j) = RHS(2,j) - (LHS(2,1,3,j)*RHS(1,j+1) &
             + LHS(2,2,3,j)*RHS(2,j+1) &
             + LHS(2,3,3,j)*RHS(3,j+1))
         RHS(3,j) = RHS(3,j) - (LHS(3,1,3,j)*RHS(1,j+1) &
             + LHS(3,2,3,j)*RHS(2,j+1) &
             + LHS(3,3,3,j)*RHS(3,j+1))
 20   end do
      return
      end

      subroutine InitLHS(LHS, NRows)
      double precision LHS(3,3,3,*)
      integer NRows, i

      do 10 i = 1, NRows
         LHS(1,1,1,i) = 1.-2.
         LHS(1,2,1,i) = 2.
         LHS(1,3,1,i) = 3.
         LHS(2,1,1,i) = 4.
         LHS(2,2,1,i) = 5.-2.
         LHS(2,3,1,i) = 6.
         LHS(3,1,1,i) = 7.
         LHS(3,2,1,i) = 8.
         LHS(3,3,1,i) = 0.-2.
         
         LHS(1,1,2,i) = 1.
         LHS(1,2,2,i) = 2.
         LHS(1,3,2,i) = 3.
         LHS(2,1,2,i) = 4.
         LHS(2,2,2,i) = 5.
         LHS(2,3,2,i) = 6.
         LHS(3,1,2,i) = 7.
         LHS(3,2,2,i) = 8.
         LHS(3,3,2,i) = 0.
         
         LHS(1,1,3,i) = 1.-3.
         LHS(1,2,3,i) = 2.
         LHS(1,3,3,i) = 3.
         LHS(2,1,3,i) = 4.
         LHS(2,2,3,i) = 5.-3.
         LHS(2,3,3,i) = 6.
         LHS(3,1,3,i) = 7.
         LHS(3,2,3,i) = 8.
         LHS(3,3,3,i) = 0.-3.
 10   end do
      return
      end

      

