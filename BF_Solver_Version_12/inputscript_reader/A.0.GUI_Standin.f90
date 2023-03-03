PROGRAM GUI_STANDIN
	IMPLICIT NONE
	INTEGER :: I

	! DOMAIN SETUP
	REAL(8) :: Lx, Ly, Lz
	INTEGER :: Nx, Ny, Nz
	INTEGER :: Grid_type
	REAL(8) :: Ug_bottom, Vg_bottom
	REAL(8) :: Ug_top, Vg_top
	REAL(8) :: Ug_left, Vg_left
	REAL(8) :: Ug_right, Vg_right
	
	
	! GAS PHASE VARIABLES
	INTEGER :: Condition
	REAL(8) :: gas_flowrate
	REAL(8) :: rho_g
	REAL(8) :: mu_g
	
	! SOLID PHASE VARIABLES
	INTEGER :: N_Particle_type
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: Num, Insertion_rate, Deletion_rate
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Insertion_min_x,Insertion_min_y
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Insertion_max_x,Insertion_max_y
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Deletion_min_x,Deletion_min_y
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Deletion_max_x,Deletion_max_y
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: dp
	
	REAL(8) :: phi_s
	REAL(8) :: rho_s
	REAL(8) :: coeff_friction, coeff_restitution

	INTEGER :: N_Pressure_Ports
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Port_X, Port_Y
	
	! REAL(8), ALLOCATABLE, DIMENSION(:)	:: phi_s
	! REAL(8), ALLOCATABLE, DIMENSION(:)	:: rho_s
	! REAL(8), ALLOCATABLE, DIMENSION(:)	:: coeff_friction, coeff_restitution
	
	! ALLOCATE (Num(1:N_Particle_type))
	! ALLOCATE (Insertion_min_x(1:N_Particle_type), Insertion_min_y(1:N_Particle_type))
	! ALLOCATE (Insertion_max_x(1:N_Particle_type), Insertion_max_y(1:N_Particle_type))
	! ALLOCATE (dp(1:N_Particle_type), phi_s(1:N_Particle_type), rho_s(1:N_Particle_type))
	! ALLOCATE (coeff_friction(1:N_Particle_type), coeff_restitution(1:N_Particle_type))
	
	
	INTEGER :: Particle_droprate
	! FINES PHASE VARIABLES
	REAL(8) :: fine_mass_flowrate
	REAL(8) :: dpf
	REAL(8) :: phi_f
	REAL(8) :: rho_f
	REAL(8) :: mu_f

	! LIQUID PHASE VARIABLES
	REAL(8)	:: Liquid_flow_rate
	REAL(8) :: rho_l
	REAL(8)	:: mu_l
	REAL(8)	:: contact_angle
	REAL(8)	:: surface_tension
	
	
	REAL(8) :: Raceway_Voidage
	REAL(8) :: Domain_Voidage
	
	INTEGER :: Flag_Raceway_Computation

	REAL(8)	:: r, D_T, h
	
	REAL(8) :: mark
	REAL(8) :: time_step
	
	REAL(8) :: Area_particle, Area_cell, PI
	REAL(8) :: Kn, sim_time, time_ste
	
	! CONVERGENCE PARAMETERS
	REAL(8) :: tol
	REAL(8) :: mass_bal
	REAL(8) :: relaxu
	REAL(8) :: relaxv
	REAL(8) :: relaxp
	REAL(8) :: relaxke
	REAL(8) :: relaxde
	REAL(8) :: relaxepfd
	REAL(8) :: Omega_u
	REAL(8) :: Omega_v
	REAL(8) :: Omega_p
	REAL(8) :: Omega_ke
	REAL(8) :: Omega_de
	REAL(8) :: Omega_uf
	REAL(8) :: Omega_vf
	
	! TURBULENCE MODEL	
	REAL(8) :: c_mu
	REAL(8) :: c_1
	REAL(8) :: c_2
	REAL(8) :: sigma_k
	REAL(8) :: sigma_e
	REAL(8) :: lm 
	

	! ROTAMETER DETAILS
	INTEGER :: N_Rotameters
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Rotameter_distance
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Rotameter_height
	REAL(8), ALLOCATABLE, DIMENSION(:)	:: Rotameter_opening

	
	
	CHARACTER(50) :: Structure_Name, Filename
	INTEGER, ALLOCATABLE, DIMENSION(:) ::Structure_Type
	INTEGER :: N_tuyeres, N_Structures
	
	REAL(8), ALLOCATABLE, DIMENSION(:) ::A_Xc, B_Xc, C_Xc
	REAL(8), ALLOCATABLE, DIMENSION(:) ::A_Yc, B_Yc, C_Yc
	REAL(8), ALLOCATABLE, DIMENSION(:) ::radius
	
	REAL(8), ALLOCATABLE, DIMENSION(:) ::Block_length, Block_width, Voidage_Geom
	REAL(8), ALLOCATABLE, DIMENSION(:) ::Tuyere_height, Tuyere_opening,Tuyere_protrusion,Tuyere_angle,Gas_Inlet
	
	PI = 3.142
!-----------------------------------------------------------------------------
	OPEN(11, FILE = "1.0.Inputscript_USER.inp", FORM="FORMATTED")
		READ(11,*) !		PURE GAS SOLVER
		READ(11,*) !--------------------------------------
		READ(11,*) Lx, Ly, Lz
		READ(11,*) Nx, Ny, Nz
		READ(11,*) Grid_type
		READ(11,*) Ug_top, Vg_top
		READ(11,*) Ug_bottom, Vg_bottom
		READ(11,*) Ug_left, Vg_left
		READ(11,*) Ug_right, Vg_right	
		READ(11,*) Domain_Voidage
		READ(11,*) Raceway_Voidage
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	GAS DETAILS
		READ(11,*) !--------------------------------------
		READ(11,*) Condition
		READ(11,*) gas_flowrate
		READ(11,*) rho_g
		READ(11,*) mu_g
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	SOLID DETAILS
		READ(11,*) !--------------------------------------	
		READ(11,*) Flag_Raceway_Computation
		READ(11,*) N_Particle_type
		
		ALLOCATE ( Num(1:N_Particle_type) )
		ALLOCATE ( Insertion_min_x(1:N_Particle_type), Insertion_min_y(1:N_Particle_type) )
		ALLOCATE ( Insertion_max_x(1:N_Particle_type), Insertion_max_y(1:N_Particle_type) )
		ALLOCATE ( Insertion_rate(1:N_Particle_type) )
		ALLOCATE ( Deletion_min_x(1:N_Particle_type), Deletion_min_y(1:N_Particle_type) )
		ALLOCATE ( Deletion_max_x(1:N_Particle_type), Deletion_max_y(1:N_Particle_type) )
		ALLOCATE ( Deletion_rate(1:N_Particle_type) )
		ALLOCATE ( dp(1:N_Particle_type) )
		DO i = 1, N_Particle_type
			READ(11,*) Num(I)
			READ(11,*) Insertion_min_x(I), Insertion_min_y(I)
			READ(11,*) Insertion_max_x(I), Insertion_max_y(I)
			READ(11,*) Insertion_rate(I)
			READ(11,*) Deletion_min_x(I), Deletion_min_y(I)
			READ(11,*) Deletion_max_x(I), Deletion_max_y(I)
			READ(11,*) Deletion_rate(I)
			READ(11,*) dp(I)
		END DO
		READ(11,*) phi_s
		READ(11,*) rho_s
		READ(11,*) coeff_friction
		READ(11,*) coeff_restitution
		READ(11,*) Kn
		READ(11,*) sim_time	
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	LIQUID DETAILS
		READ(11,*) !--------------------------------------
		READ(11,*) Liquid_flow_rate
		READ(11,*) rho_l
		READ(11,*) mu_l
		READ(11,*) contact_angle
		READ(11,*) surface_tension
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	FINE DETAILS
		READ(11,*) !--------------------------------------
		READ(11,*) fine_mass_flowrate
		READ(11,*) dpf
		READ(11,*) phi_f
		READ(11,*) rho_f
		READ(11,*) mu_f
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	TUYERE DETAILS
		READ(11,*) !--------------------------------------
		READ(11,*) N_Tuyeres
		IF(N_Tuyeres.GT.0) THEN
			ALLOCATE (Tuyere_height(1:N_Tuyeres), Tuyere_opening(1:N_Tuyeres), Tuyere_protrusion(1:N_Tuyeres),Tuyere_angle(1:N_Tuyeres))
			ALLOCATE (Gas_Inlet(1:N_Tuyeres))
			DO i = 1, N_Tuyeres
				READ(11,*) ! Tuyere Label
				READ(11,*) Tuyere_height(I)
				READ(11,*) Tuyere_opening(I)
				READ(11,*) Tuyere_protrusion(I)
				READ(11,*) Tuyere_angle(I)
			END DO
		END IF
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	ROTAMETERS DETAILS
		READ(11,*) !--------------------------------------
		READ(11,*) N_Rotameters
		IF (N_Rotameters.GT.0) THEN
			ALLOCATE (Rotameter_distance(1:N_Rotameters), Rotameter_height(1:N_Rotameters), Rotameter_opening(1:N_Rotameters))
			DO i = 1, N_Rotameters
				READ(11,*) ! Rotameter Label
				READ(11,*) Rotameter_distance(I)
				READ(11,*) Rotameter_height(I)
				READ(11,*) Rotameter_opening(I)
			END DO
		END IF
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	INTERNAL STRUCTURES DETAILS
		READ(11,*) !--------------------------------------
		READ(11,*) N_Structures		
		IF(N_Structures.GT.0) THEN
			ALLOCATE (A_Xc(1:N_Structures), B_Xc(1:N_Structures), C_Xc(1:N_Structures))
			ALLOCATE (A_Yc(1:N_Structures), B_Yc(1:N_Structures), C_Yc(1:N_Structures))
			ALLOCATE (Block_length(1:N_Structures), Block_width(1:N_Structures), radius(1:N_Structures))
			ALLOCATE (Structure_Type(1:N_Structures))
			ALLOCATE (Voidage_Geom(1:N_Structures))
			DO I = 1, N_Structures
				READ(11,*) Structure_Name
				IF(Structure_Name .EQ. "TRIANGLE") THEN
					Structure_Type(i) = 1
					READ(11,*) C_Xc(i), C_Yc(i)
					READ(11,*) B_Xc(i), B_Yc(i)
					READ(11,*) A_Xc(i), A_Yc(i)
					READ(11,*) Voidage_Geom(i)
				ELSE IF(Structure_Name .EQ. "CIRCLE") THEN
					Structure_Type(i) = 2
					READ(11,*) B_Xc(i), B_Yc(i)
					READ(11,*) radius(i)
					READ(11,*) Voidage_Geom(i)
				ELSE IF(Structure_Name .EQ. "BLOCK") THEN
					Structure_Type(i) = 3
					READ(11,*) Block_length(i), Block_width(i)
					READ(11,*) B_Xc(i), B_Yc(i)
					READ(11,*) Voidage_Geom(i)
				END IF
			END DO
		ELSE IF(N_Structures.EQ.-1) THEN
			READ(11,*) Filename
		END IF
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	CONVERGENCE PARAMETERS
		READ(11,*) !--------------------------------------
		READ(11,*) tol
		READ(11,*) mass_bal
		READ(11,*) relaxu
		READ(11,*) relaxv
		READ(11,*) relaxp
		READ(11,*) relaxke
		READ(11,*) relaxde
		READ(11,*) relaxepfd
		READ(11,*) Omega_u
		READ(11,*) Omega_v
		READ(11,*) Omega_p
		READ(11,*) Omega_ke
		READ(11,*) Omega_de
		READ(11,*) Omega_uf
		READ(11,*) Omega_vf
		READ(11,*) mark
		READ(11,*) time_step
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	PRESSURE PORTS
		READ(11,*) !--------------------------------------		
		READ(11,*) N_Pressure_Ports
		IF(N_Pressure_Ports.GT.0) THEN
			ALLOCATE (Port_X(1:N_Pressure_Ports), Port_Y(1:N_Pressure_Ports))
			DO I = 1, N_Pressure_Ports
				READ(11,*) Port_X(I), Port_Y(I)
			END DO
		END IF
		READ(11,*) !--------------------------------------
		READ(11,*) !-->	TURBULENCE MODELLING PARAMETERS
		READ(11,*) !--------------------------------------
		READ(11,*) c_1
		READ(11,*) c_2
		READ(11,*) c_mu
		READ(11,*) sigma_k
		READ(11,*) sigma_e
	CLOSE(11)
!--------------------------------------------------------------------
!--------------------------------------------------------------------	
	OPEN(21, FILE = "0.3.Input_BEAST.tmp", FORM="FORMATTED")
		WRITE(21,*) Grid_type
		WRITE(21,*) Lx
		WRITE(21,*) Ly
		WRITE(21,*) Lz

		WRITE(21,*) Nx
		WRITE(21,*) Ny
		WRITE(21,*) Nz

		WRITE(21,*) Ug_bottom
		WRITE(21,*) Ug_top
		WRITE(21,*) Ug_left
		WRITE(21,*) Ug_right
		
		WRITE(21,*) Vg_bottom
		WRITE(21,*) Vg_top
		WRITE(21,*) Vg_left
		WRITE(21,*) Vg_right

		WRITE(21,*) Condition
		WRITE(21,*) gas_flowrate
		WRITE(21,*) rho_g
		WRITE(21,*) mu_g
		
		WRITE(21,*) Liquid_flow_rate
		WRITE(21,*) rho_l
		WRITE(21,*) mu_l
		WRITE(21,*) contact_angle
		WRITE(21,*) surface_tension
		
		WRITE(21,*) fine_mass_flowrate
		WRITE(21,*) dpf
		WRITE(21,*) phi_f
		WRITE(21,*) rho_f
		WRITE(21,*) mu_f
				
		WRITE(21,*) tol
		WRITE(21,*) mass_bal
		
		WRITE(21,*) relaxu
		WRITE(21,*) relaxv
		WRITE(21,*) relaxp
		WRITE(21,*) relaxke
		WRITE(21,*) relaxde
		WRITE(21,*) relaxepfd

		WRITE(21,*) Omega_u
		WRITE(21,*) Omega_v
		WRITE(21,*) Omega_p
		WRITE(21,*) Omega_ke
		WRITE(21,*) Omega_de
		WRITE(21,*) Omega_uf
		WRITE(21,*) Omega_vf
		
		WRITE(21,*) c_1
		WRITE(21,*) c_2
		WRITE(21,*) c_mu
		WRITE(21,*) sigma_k
		WRITE(21,*) sigma_e
			
		WRITE(21,*) N_Pressure_Ports
		DO I = 1, N_Pressure_Ports
			WRITE(21,*) Port_X(I)
			WRITE(21,*) Port_Y(I)
		END DO
		
		WRITE(21,*) Flag_Raceway_Computation
		WRITE(21,*) N_Particle_type
		DO I = 1, N_Particle_type
			WRITE(21,*) Num(I)
			WRITE(21,*) dp(I)
			WRITE(21,*) Insertion_min_x(I)
			WRITE(21,*) Insertion_min_y(I)
			WRITE(21,*) Insertion_max_x(I)
			WRITE(21,*) Insertion_max_y(I)
			WRITE(21,*) Insertion_rate(I)
			WRITE(21,*) Deletion_min_x(I)
			WRITE(21,*) Deletion_min_y(I)
			WRITE(21,*) Deletion_max_x(I)
			WRITE(21,*) Deletion_max_y(I)
			WRITE(21,*) Deletion_rate(I)
		END DO
		WRITE(21,*) phi_s
		WRITE(21,*) rho_s
		WRITE(21,*) coeff_friction
		WRITE(21,*) coeff_restitution
		WRITE(21,*) Kn
		WRITE(21,*) sim_time
		WRITE(21,*) mark
		WRITE(21,*) time_step
		
		WRITE(21,*) Raceway_Voidage
		WRITE(21,*) Domain_Voidage

	CLOSE(21)
	
	
	OPEN(22, FILE = "0.4.Input_Geometry.tmp", FORM="FORMATTED")
		WRITE(22,*) N_tuyeres
		IF(N_tuyeres.GT.0) THEN
			DO I =1, N_tuyeres
				WRITE(22,*) Tuyere_height(I)
				WRITE(22,*) Tuyere_opening(I)
				WRITE(22,*) Tuyere_protrusion(I)
				WRITE(22,*) Tuyere_angle(I)
				WRITE(22,*) Gas_Inlet(I)
			END DO
		END IF
		WRITE(22,*) N_Rotameters
		DO i = 1, N_Rotameters
			WRITE(22,*) Rotameter_distance(I)
			WRITE(22,*) Rotameter_height(I)
			WRITE(22,*) Rotameter_opening(I)
		END DO
		WRITE(22,*) N_Structures
		IF(N_Structures.GT.0) THEN
			DO I = 1, N_Structures
				IF(Structure_Type(I) .EQ. 1) THEN
					WRITE(22,*)Structure_Type (I)
					WRITE(22,*) C_Xc(i), C_Yc(i)
					WRITE(22,*) B_Xc(i), B_Yc(i)
					WRITE(22,*) A_Xc(i), A_Yc(i)
					WRITE(22,*) Voidage_Geom(i)
				ELSE IF(Structure_Type (I).EQ. 2) THEN
					WRITE(22,*)Structure_Type(I)
					WRITE(22,*) B_Xc(i), B_Yc(i)
					WRITE(22,*) radius(i)
					WRITE(22,*) Voidage_Geom(i)
				ELSE IF(Structure_Type (I).EQ. 3) THEN
					WRITE(22,*)Structure_Type(I)
					WRITE(22,*) Block_length(i), Block_width(i)
					WRITE(22,*) B_Xc(i), B_Yc(i)
					WRITE(22,*) Voidage_Geom(i)
				END IF
			END DO
		ELSE IF(N_Structures.EQ.-1) THEN
			WRITE(22,*) Filename
		END IF
	CLOSE(22)	
END PROGRAM GUI_STANDIN

