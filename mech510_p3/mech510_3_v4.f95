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
	
	INTEGER, PARAMETER :: scheme = 2                        ! 1 for RK2, 2 for IE
	

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

	!REAL*8, PARAMETER  :: delta_t = (Re*Pr*dx**2)/4. ! Time step from stability considerations
	REAL*8, PARAMETER  :: delta_t = 0.05 ! Time step from stability considerations
    INTEGER, PARAMETER :: maxn = nint(tmax/delta_t)    ! Max time step (i.e. last n index)
	
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmesh          ! Mesh of solutions; includes ghost cells
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: umesh,vmesh    ! Mesh of velocities; includes ghost cells
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: FIn,FI_ana          ! Computed and analytic flux integrals at time level n
	REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: source_mesh,source_mesh_ana  ! Computed and analytic sources 
	
END MODULE meshmod

! Program structure
! -------------------

! Main Program


! -------------------

! =================== main program =========================

PROGRAM mech510_3  

    USE meshmod
    IMPLICIT NONE

    REAL*8   ::   tn = 0.0D+0    ! Time at index n
	
	WRITE(*,*) maxn
    WRITE(*,*) delta_t
	! Initialize mesh
	CALL init_mesh
	CALL output_Tmesh(tn)
	! Calculate flux integrals
    !CALL calc_integrals(Tmesh,FIn)
	
	!CALL rk2(0,maxn-1,tn)
	
	CALL ie(0,maxn-1,tn)
	
	! Output computed and analytic flux integrals
	!CALL output_integrals
	
	
	CALL output_Tmesh(tn)
	
END PROGRAM mech510_3


! =================== subroutines =========================

SUBROUTINE init_mesh     ! Initialize mesh with CV average values of T(x,y), u(x,y), v(x,y) and source terms
	
	USE meshmod
    IMPLICIT NONE
	
	REAL*8 :: xa,xb,yc,yd    ! Positions of cell edges
	REAL*8                               :: dudx,dudy,dvdx,dvdy    ! Flux terms for u and v
	REAL*8                               :: source                 ! Source at face
	REAL*8                               :: source_analytic        ! Analytic source term
	
	
	! First pass through array: initialize Tmesh and vector field
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
			umesh(j,i) = 0
			vmesH(j,i) = 0
			!Tmesh(j,i) = T0*cos(pi*(i-1./2.)*dx)*(sin(pi*(j-1./2.)*dy))
	        !umesh(j,i) = u0*((j-1./2.)*dy)*sin(pi*(i-1./2.)*dx)
			!vmesh(j,i) = v0*((i-1./2.)*dx)*cos(pi*(j-1./2.)*dy)
			
		END DO
	END DO
	
	
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

END SUBROUTINE init_mesh


! =======================

SUBROUTINE set_bcs(Tw, tn, Tmeshn)             ! Subroutine to update boundary conditions

    USE meshmod
    IMPLICIT NONE
    
    REAL*8                        :: Tw       ! Wall value
    REAL*8                        :: tn       ! Time at index n
    REAL*8, DIMENSION(0:jmax+1,0:imax+1) :: Tmeshn    ! Mesh at time level n (could be interim mesh)
    
    ! Set wall value at time tn
    Tw = sin(4.*pi*tn)

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
	
	REAL*8, DIMENSION(1:jmax,1:imax)   :: rhsn           ! Right-hand-side of energy equation at time level n
	REAL*8                               :: rhs1           ! Interim right-hand-side of energy equation at time level n

    DO n = nmin,nmax  ! Want to obtain solution at maxn, hence last loop is maxn-1
        
        ! ------Calculate T(1)------
        tn = delta_t*n                     ! Current time
        
        CALL set_bcs(Tw, tn, Tmesh)         ! Update boundary conditions of solution mesh
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
		
        CALL set_bcs(Tw, tn, mesh_1)         ! Update boundary conditions of interim mesh
        
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


SUBROUTINE set_lhs(lhs,vsize,alpha,betas,bca,bcb,bcc,bcd)     ! Initialize LHS for approximate factorization, for IE time advance
	
	USE meshmod
    IMPLICIT NONE
	

	INTEGER                              :: vsize             ! Size of column (in x, first call) or row (in y, second call)
	REAL*8, DIMENSION(0:vsize+1,3)       :: lhs 
	REAL*8                               :: alpha
	REAL*8, DIMENSION(1:vsize)           :: betas
	INTEGER                              :: bca,bcb,bcc,bcd
	INTEGER                              :: jind,iind
	

	DO jind = 0,vsize+1
	    !WRITE(*,*) jind
		iind = jind 
		
		IF (jind == 0) THEN     ! Boundary conditions (left and bottom)
		    lhs(jind ,1) = 0
			lhs(jind ,2) = bca   
			lhs(jind ,3) = bcb
					
		ELSE IF (jind  == vsize+1) THEN  ! Boundary conditions (right and top)
		    lhs(jind ,1) = bcc 
			lhs(jind ,2) = bcd
			lhs(jind ,3) = 0
		
		ELSE
			
			lhs(jind ,1) = -alpha-betas(iind-1)
			lhs(jind ,2) = 1.+2.*alpha
			lhs(jind ,3) = betas(iind+1)-alpha
			
		END IF
	END DO
		    
		

END SUBROUTINE set_lhs

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
	REAL*8, DIMENSION(0:imax+1)          :: betasx
	REAL*8, DIMENSION(0:jmax+1)          :: betasy
	
	
	DO n = nmin,nmax  ! Want to obtain solution at maxn, hence last loop is maxn-1
        
        ! ------Calculate T(1)------
        tn = delta_t*n                     ! Current time
        
		
        CALL set_bcs(Tw, tn, Tmesh)         ! Update boundary conditions of solution mesh
        CALL calc_integrals(Tmesh,FIn)      ! Calculate flux integrals
		
        ! Go across columns from row to row; first part of approx factorization
		alpha = delta_t/(Re*Pr*dx**2)
        DO j = 1,jmax
		    !WRITE(*,*) j
            rhs = (source_mesh(j,:) - FIn(j,:))*delta_t
			
			! Set LHS for approximate factorization (first part)
			lhsx = 0.0D+0
			betasx = umesh(j,:)*delta_t/(2.*dx)

			CALL set_lhs(lhsx,imax,alpha,betasx,-1,1,-1,1)
			
			CALL Thomas(lhsx,rhs,imax)
			
			dT_tilda(j,:) = rhs

        END DO
		

		! Go across rows from column to column; second part of approx factorization
		!alpha = delta_t/(Re*Pr*dy**2)
		DO i = 1,imax
			! Set LHS for approximate factorization (second part)
			lhsy = 0.0D+0
			betasy = vmesh(:,i)*delta_t/(2.*dy)

			CALL set_lhs(lhsy,jmax,alpha,betasy,-1,1,-1,1)
			
			CALL Thomas(lhsy,dT_tilda(:,i),jmax)
			
			dT(:,i) = dT_tilda(:,i)
        END DO
		
		
		DO j = 1,jmax
			DO i = 1,imax
				Tmesh(j,i) = Tmesh(j,i) + dT(j,i)
			END DO
		END DO 
		
    END DO
	
END SUBROUTINE ie
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
			!FI_ana(j,i) = analytic_cv((i-1)*dx, i*dx, (j-1)*dy, j*dy)
			
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


SUBROUTINE output_Tmesh(tn) ! Output final solutions for plotting
    
    USE meshmod
    IMPLICIT NONE
    
	REAL*8   ::   tn     ! Time at index n
    INTEGER :: ierror                                                  ! For error catching
    CHARACTER(len=100) :: fname                                        ! Output filenames
    CHARACTER(len=100) :: fmt                                          ! Output format
        
	! Set filename format
	fmt = "(A8,I0,A6,I0,A6,I0,A3,I0,A8,I0,A4)"  
	
    ! Set filenames
    WRITE(fname, fmt) "Tmesh_cn",imax,"_xmax_",nint(xmax),"_ymax_",nint(ymax),"_t_",nint(tn),"_scheme_",scheme,".txt"
    WRITE(*,*) "Writing output to ", fname
    OPEN(UNIT=90, FILE=fname,STATUS="REPLACE", IOSTAT=ierror)
  
	DO j=1,jmax   ! Loop over each row
        WRITE(90,"(999F15.10)") (Tmesh(j,i), i=1,imax)             ! Implied loop to write out every column (computed flux integral)
    ENDDO

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





