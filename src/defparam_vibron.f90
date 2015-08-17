MODULE vibron_u4
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  ! Type definitions
  TYPE :: so4_bas
     INTEGER(KIND = I4B) :: N_U4_val
     INTEGER(KIND = I4B) :: omega_SO4_val
     INTEGER(KIND = I4B) :: L_val
  END TYPE so4_bas
  !
  TYPE :: u3_bas
     INTEGER(KIND = I4B) :: N_U4_val
     INTEGER(KIND = I4B) :: np_U3_val
     INTEGER(KIND = I4B) :: L_val
  END TYPE u3_bas
  !
  TYPE :: exp_level
     REAL(KIND = DP) :: exp_term_value
     REAL(KIND = DP) :: exp_err_value
     INTEGER(KIND = I4B) :: exp_v
     INTEGER(KIND = I4B) :: exp_L
     CHARACTER(LEN = 1) :: spin_symm
     CHARACTER(LEN = 3) :: reference
  END TYPE exp_level
  !
  TYPE term_values
     TYPE(exp_level) :: exper_level
     TYPE(term_values), POINTER :: next
  END TYPE term_values
  !
  ! Parameters
  !
  TYPE(u3_bas), DIMENSION(:), ALLOCATABLE :: U3_Basis !  U(4) > U(3) > SO(3) basis
  !
  TYPE(so4_bas), DIMENSION(:), ALLOCATABLE :: SO4_Basis !  U(4) > SO(4) > SO(3) basis
  !
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Ham_U4_mat ! Hamiltonian matrix
  !
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: eigenval_vec ! Hamiltonian Eigenvalues
  !
  INTEGER(KIND = I4B) :: ierr
  !
  INTEGER(KIND = I4B) :: N_val ! U(4) N value
  ! 
  INTEGER(KIND = I4B) :: L_min ! Minimum SO(3) angular momentum 
  !
  INTEGER(KIND = I4B) :: L_max ! Minimum SO(3) angular momentum 
  !
  INTEGER(KIND = I4B), PARAMETER :: nHamoper = 13_I4B ! Number of Hamiltonian operators
  !
  REAL(KIND = DP), DIMENSION(1:nHamoper) :: Ham_parameter = 0.0_DP ! Hamiltonian parameters 
  ! 
  LOGICAL :: eigenvec   ! If .T. compute eigenvalues and eigenvectors
  LOGICAL :: excitation ! If .T. compute excitation energy with respect to L = 0 g.s.
  LOGICAL :: save_avec  ! If .T. save eigenvector components.
  LOGICAL :: save_ham   ! If .T. save hamiltonian matrix
  !
  REAL(KIND = DP) :: GS_energy ! L = 0 Ground state energy 
  !
  !
  REAL(KIND = DP), PARAMETER :: Zero_Parameter = 1.0E-20_DP ! Zero limit to discard a parameter evaluation
  !
  REAL(KIND = SP) :: time_check, time_check_ref
  !
  TYPE(exp_level), DIMENSION(:), ALLOCATABLE :: exp_data
  !
  ! Control verbosity
  INTEGER(KIND = I4B) :: Iprint = 1_I4B
  !
  !
CONTAINS
  !
  !
  FUNCTION DIM_L_BLOCK(N_val, L_val)
    !
    IMPLICIT NONE
    !
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    !
    INTEGER(KIND = I4B) :: DIM_L_BLOCK 
    !
    !
    DIM_L_BLOCK = (N_val - L_val + 2_I4B - MOD(N_val + L_val,2_I4B))/2_I4B
    !
    !
    !
  END FUNCTION DIM_L_BLOCK
  !
  SUBROUTINE SO4_BASIS_VIBRON(N_val, L_val) 
    !
    ! Subroutine to build the U(4) > SO(4) > SO(3) basis
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index, omega
    !
    !
    !
    IF (MOD(N_val + L_val, 2) == 0) THEN
       !
       omega = L_val ! same parity case
       !
    ELSE
       !
       omega = L_val + 1 ! different parity case
       !
    ENDIF
    !
    index = 1_I4B
    !
    DO WHILE (omega <= N_val) 
       !
       SO4_Basis(index)%omega_SO4_val = omega
       !
       index = index + 1_I4B
       omega = omega + 2_I4B
       !
    ENDDO
    !
    SO4_Basis(1:Index-1)%N_U4_val = N_val
    SO4_Basis(1:Index-1)%L_val = L_val
    !
    !
  END SUBROUTINE SO4_BASIS_VIBRON
  !
  !
  SUBROUTINE U3_BASIS_VIBRON(N_val, L_val) 
    !
    ! Subroutine to build the U(4) > U(3) > SO(3) basis
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    !
    ! Local Variables
    INTEGER(KIND = I4B) :: index, np
    !
    !
    !
    !
    index = 1_I4B
    !
    np = L_val
    !
    DO WHILE (np <= N_val) 
       !
       U3_Basis(index)%np_U3_val = np
       !
       index = index + 1_I4B
       np = np + 2_I4B
       !
    ENDDO
    !
    U3_Basis(1:Index-1)%N_U4_val = N_val
    U3_Basis(1:Index-1)%L_val = L_val
    !
    !
  END SUBROUTINE U3_BASIS_VIBRON
  !
  !
  SUBROUTINE U3_HAMILTONIAN_VIBRON(N_val, L_val, dim_block) 
    !
    ! Subroutine to build the U(4) > U(3) > SO(3) Hamiltonian
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    INTEGER(KIND = I4B) :: np, oper, U3_state
    !
    !
    !
    ! Build Hamiltonian
    operator : DO oper = 1, nHamoper
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(Ham_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
       CASE (1) ! n ---> DIAGONAL
          !
          DO U3_state = 1, dim_block
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             !
             Ham_U4_mat(U3_state, U3_state) = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(1)*REAL(np,DP)
             !
          ENDDO
          !
          !
       CASE (2) ! n^2 ---> DIAGONAL
          !
          DO U3_state = 1, dim_block
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state) = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(2)*REAL(np**2,DP)
             !
          ENDDO
          !
          !
       CASE (3) ! (L^2 + D^2)/N 
          !
          DO U3_state = 1, dim_block
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state) = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(3)*f_D2L2_diag(N_val, np, L_val)/REAL(N_val,DP)
          ENDDO
          !
          !
          ! Non-Diagonal contribution  np, np + 2
          DO U3_state = 1, dim_block - 1
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state + 1) = Ham_U4_mat(U3_state, U3_state + 1) + &
                  Ham_parameter(3)*f_D2L2_nondiag(N_val, np, L_val)/REAL(N_val,DP)
             !
          ENDDO
          !
       CASE (4) ! L^2 ---> DIAGONAL
          !
          DO U3_state = 1, dim_block
             Ham_U4_mat(U3_state, U3_state)  = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(4)*REAL(L_val*(L_val + 1_I4B),DP)
          ENDDO
          !
          !
       CASE (5) ! L^4 ---> DIAGONAL
           !
          DO U3_state = 1, dim_block
             Ham_U4_mat(U3_state, U3_state)  = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(5)*REAL(L_val*(L_val + 1_I4B)*L_val*(L_val + 1_I4B))
          ENDDO
          !
          !
       CASE (6) ! [L^2(L^2 + D^2)]/N 
          !
          DO U3_state = 1, dim_block
             np = U3_Basis(U3_state)%np_U3_val
             Ham_U4_mat(U3_state, U3_state) = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(6)* &
                  REAL(L_val*(L_val + 1_I4B),DP)*f_D2L2_diag(N_val, np, L_val)/REAL(N_val,DP)
          ENDDO
          !
          !
          ! Non-Diagonal contribution  np, np + 2
          DO U3_state = 1, dim_block - 1
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state + 1) = Ham_U4_mat(U3_state, U3_state + 1) + &
                  Ham_parameter(6)*REAL(L_val*(L_val + 1_I4B),DP)*f_D2L2_nondiag(N_val, np, L_val)/REAL(N_val,DP)
             !
          ENDDO
          !
          !
       CASE (7) ! (L^2 + D^2)^2/N^2
          !
          DO U3_state = 1, dim_block
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state) = Ham_U4_mat(U3_state, U3_state) + &
                  Ham_parameter(7)*( &
                  f_D2L2_diag(N_val, np, L_val)*f_D2L2_diag(N_val, np, L_val) + & 
                  f_D2L2_nondiag(N_val, np, L_val)**2 + f_D2L2_nondiag(N_val, np-2, L_val)**2 & 
                  )/REAL(N_val*N_val,DP)
          ENDDO
          !
          ! Non-Diagonal contribution  np, np + 2
          DO U3_state = 1, dim_block - 1
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state + 1) = Ham_U4_mat(U3_state, U3_state + 1) + &
                  Ham_parameter(7)*( &
                  f_D2L2_diag(N_val, np, L_val)*f_D2L2_nondiag(N_val, np, L_val) + &
                  f_D2L2_diag(N_val, np+2, L_val)*f_D2L2_nondiag(N_val, np, L_val) &
                  )/REAL(N_val*N_val,DP)
             !
          ENDDO
          !
          ! Non-Diagonal contribution  np, np + 4
          DO U3_state = 1, dim_block - 2
             !
             np = U3_Basis(U3_state)%np_U3_val
             !
             Ham_U4_mat(U3_state, U3_state + 2) = Ham_U4_mat(U3_state, U3_state + 2) + &
                  Ham_parameter(7)*( &
                  f_D2L2_nondiag(N_val, np, L_val)*f_D2L2_nondiag(N_val, np + 2, L_val) &
                  )/REAL(N_val*N_val,DP)
             !
          ENDDO
          !
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       END SELECT
       !
    ENDDO operator
    !
    !
    !
  END SUBROUTINE U3_HAMILTONIAN_VIBRON
  !
  !
  FUNCTION f_D2L2_diag(N_val, np, L_val)
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, np, L_val
    !
    REAL(KIND = DP) ::  f_D2L2_diag
    !
    f_D2L2_diag =  REAL(N_val*(2.0_DP*np+3_I4B) - 2_I4B*np*(np+1_I4B) + L_val*(L_val + 1_I4B),DP) 
    !
    !
    !
  END FUNCTION f_D2L2_diag
  !
  !
  FUNCTION f_D2L2_nondiag(N_val, np, L_val)
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, np, L_val
    !
    REAL(KIND = DP) ::  f_D2L2_nondiag
    !
    IF (np < 0) THEN
       f_D2L2_nondiag = 0.0_DP
    ELSE
       f_D2L2_nondiag =  SQRT( &
            REAL((N_val - np - 1),DP) * REAL(N_val - np,DP) * REAL(np - L_val + 2_I4B,DP) * REAL(np + L_val + 3_I4B,DP) &
            )
    ENDIF
    !
  END FUNCTION f_D2L2_nondiag
  !  
  !
  SUBROUTINE SO4_HAMILTONIAN_VIBRON(N_val, L_val, dim_block) 
    !
    ! Subroutine to build the U(4) > SO(4) > SO(3) Hamiltonian
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    !
    INTEGER(KIND = I4B) :: omega, oper, SO4_state
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Build Hamiltonian
    operator : DO oper = 1, nHamoper
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(Ham_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
          !
       CASE (1) ! n
          !
          ! Diagonal contribution 
          DO SO4_state = 1, dim_block
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(1)*f_n_diag(N_val, omega, L_val)
             !
          ENDDO
          !
          ! Non-Diagonal contribution  omega, omega + 2
          DO SO4_state = 1, dim_block - 1
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 1) = Ham_U4_mat(SO4_state, SO4_state + 1) + &
                  Ham_parameter(1)*f_n_nondiag(N_val, omega, L_val)
             !
          ENDDO
          !
          !
       CASE (2) ! n^2
          !
          !
          ! Diagonal contribution 
          DO SO4_state = 1, dim_block
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(2)*( &
                  f_n_diag(N_val, omega, L_val)**2 + &
                  f_n_nondiag(N_val, omega, L_val)**2 + &
                  f_n_nondiag(N_val, omega - 2, L_val)**2 &
                  )
             !
          ENDDO
          !
          ! Non-Diagonal contribution omega, omega + 2
          DO SO4_state = 1, dim_block - 1
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 1) = Ham_U4_mat(SO4_state, SO4_state + 1) + &
                  Ham_parameter(2)*( &
                  f_n_diag(N_val, omega, L_val)*f_n_nondiag(N_val, omega, L_val) + &
                  f_n_diag(N_val, omega + 2, L_val)*f_n_nondiag(N_val, omega, L_val) &
                  )
             !
          ENDDO
          !
          ! Non-Diagonal contribution omega, omega + 4
          DO SO4_state = 1, dim_block - 2
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 2) = Ham_U4_mat(SO4_state, SO4_state + 2) + &
                  Ham_parameter(2)*( &
                  f_n_nondiag(N_val, omega, L_val)*f_n_nondiag(N_val, omega + 2, L_val) &
                  )
             !
          ENDDO
          !
       CASE (3) ! (L^2 + D^2)/N ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(3)*REAL(omega*(omega + 2_I4B), DP)/REAL(N_val,DP)
          ENDDO
          !
          !
       CASE (4) ! L^2 ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(4)*REAL(L_val*L_val + L_val)
          ENDDO
          !
          !
       CASE (5) ! L^4 ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(5)*REAL(L_val*L_val*(L_val + 1_I4B)*(L_val + 1_I4B))
          ENDDO
          !
          !
       CASE (6) ! [L^2(L^2 + D^2)]_+/N ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(6)*REAL(L_val*(L_val + 1_I4B)*omega*(omega + 2_I4B), DP)/REAL(N_val,DP)
          ENDDO
          !
          !
       CASE (7) ! (L^2 + D^2)^2/N^2 ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(7)*REAL(omega*(omega + 2_I4B), DP)**2/REAL(N_val,DP)**2
          ENDDO
          !
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       END SELECT
       !
    ENDDO operator
    !
    !
    !
    !
  END SUBROUTINE SO4_HAMILTONIAN_VIBRON
  !
  !
  FUNCTION f_n_diag(N_val, omega, L_val)
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, omega, L_val
    !
    REAL(KIND = DP) ::  f_n_diag
    !
    !write(*,*) "f_n_d"!, N_val, omega, L_val!, f_n_diag
    IF (L_val == 0 .AND. omega /= 0) THEN
       f_n_diag =  REAL(N_val-1_I4B,DP)/2.0_DP
    ELSE IF (omega == 0) THEN
       f_n_diag =  REAL(3_I4B*N_val,DP)/4.0_DP
    ELSE
       f_n_diag =  REAL(N_val-1_I4B,DP)/2.0_DP + &
            REAL((N_val+2_I4B)*L_val*(L_val + 1_I4B),DP)/REAL(2_I4B*omega*(omega + 2_I4B), DP)
    ENDIF
    !
    !
    !
  END FUNCTION f_n_diag
  !
  FUNCTION f_n_nondiag(N_val, omega, L_val)
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val, omega, L_val
    !
    REAL(KIND = DP) ::  f_n_nondiag
    !
    !
    IF (omega < 0) THEN
       f_n_nondiag = 0.0_DP
    ELSE
       f_n_nondiag =  SQRT( &
            REAL(N_val-omega,DP)*POCCHAMMER_S(omega-L_val+1_I4B,2_I4B)* &
            REAL(N_val+omega+4_I4B,DP)*POCCHAMMER_S(omega+L_val+2_I4B,2_I4B)/&
            (REAL(16_I4B*(omega + 2_I4B),DP)*POCCHAMMER_S(omega + 1_I4B, 3_I4B)) &
            )
    ENDIF
    !
    !
  END FUNCTION f_n_nondiag
  !
  !
  FUNCTION POCCHAMMER_S(a, b)
    !
    ! Pocchammer symbol as defined in Frank and Isacker 5.176 p.162
    ! (rising factorial)
    !
    ! by Currix TM
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=I4B), INTENT(IN) :: a, b
    REAL(KIND=DP) :: POCCHAMMER_S
    !
    INTEGER(KIND=I4B) :: index  
    !
    POCCHAMMER_S = REAL(a,DP)
    !
    DO index = 1, b - 1
       POCCHAMMER_S = POCCHAMMER_S * REAL(a + index,DP)
    ENDDO
    !
    !
    !
  END FUNCTION POCCHAMMER_S
  !
  !
  SUBROUTINE read_experimental_data(exp_data_filename, l_exp_min, l_exp_max, v_exp_min, v_exp_max, num_exp_data)
    ! 
    ! Subroutine to read experimental rovibrational data for diatomic molecules
    ! FORMAT:
    !   v   K   o/p   wavenumber    flag    ref
    !
    ! Lines starting with '#' are comments and should be skipped
    !
    !
    CHARACTER(LEN = 64), INTENT(IN) :: exp_data_filename
    INTEGER(KIND = I4B), INTENT(IN) :: l_exp_min, l_exp_max, v_exp_min, v_exp_max
    INTEGER(KIND = I4B), INTENT(OUT) :: num_exp_data
    !
    !
    ! LOCAL VARIABLES
    TYPE(term_values), TARGET  :: first_term_val  
    TYPE(term_values), POINTER :: current_term_val, temp_term_val
    !
    REAL(KIND = DP) :: term_value_READ, error_value_READ
    INTEGER(KIND = I4B) ::  L_value_READ, v_value_READ
    CHARACTER(LEN = 1) :: spin_symmetry_READ
    CHARACTER(LEN = 3) :: reference_READ
    CHARACTER(LEN = 300) :: line_READ
    INTEGER(KIND = I4B) :: index
    !
    !
    num_exp_data = 0
    !
    ! Read All Experimental Data
    !
    OPEN(UNIT=20, FILE=exp_data_filename, STATUS='OLD')
    !
    ! Define data head 
    DO 
       READ(UNIT=20, FMT=*)   line_READ
       !
       IF (line_READ(1:1) /= '#') THEN ! First line that is not a comment
          !
          BACKSPACE(UNIT = 20)
          READ(UNIT=20, FMT=*)   v_value_READ, L_value_READ, term_value_READ, error_value_READ,  &
               spin_symmetry_READ, reference_READ
          !
          IF (Iprint > 1) print*, "exp head ", v_value_READ, L_value_READ, term_value_READ, error_value_READ,  &
               spin_symmetry_READ, reference_READ
          !
          EXIT
       ENDIF
    ENDDO
    !
    first_term_val%exper_level%exp_term_value = term_value_READ
    first_term_val%exper_level%exp_err_value = error_value_READ
    first_term_val%exper_level%exp_v = v_value_READ
    first_term_val%exper_level%exp_L = L_value_READ
    first_term_val%exper_level%spin_symm = spin_symmetry_READ
    first_term_val%exper_level%reference = reference_READ
    !
    NULLIFY(first_term_val%next)
    current_term_val => first_term_val
    !
    num_exp_data = 1_I4B
    !
    ! Fill the pointer with experimental data
    DO
       !
       READ(UNIT=20, FMT=*, END=33)   line_READ
       !
       IF (line_READ(1:1) /= '#') THEN ! This line is not a comment
          !
          BACKSPACE(UNIT = 20)
          READ(UNIT=20, FMT=*)   v_value_READ, L_value_READ, term_value_READ, error_value_READ,  &
               spin_symmetry_READ, reference_READ
          !
          IF (Iprint > 1) print*, "exp data ", v_value_READ, L_value_READ, term_value_READ, error_value_READ,  &
               spin_symmetry_READ, reference_READ
          !
          !
          ALLOCATE(temp_term_val, STAT = ierr)
          IF (IERR /= 0) STOP 'STOP :: POINTER ALLOCATION FAILED. Sayonara baby.' 
          !
          temp_term_val%exper_level%exp_term_value = term_value_READ
          temp_term_val%exper_level%exp_err_value = error_value_READ
          temp_term_val%exper_level%exp_v = v_value_READ
          temp_term_val%exper_level%exp_L = L_value_READ
          temp_term_val%exper_level%spin_symm = spin_symmetry_READ
          temp_term_val%exper_level%reference = reference_READ
          !
          NULLIFY(temp_term_val%next)
          current_term_val%next => temp_term_val
          current_term_val => temp_term_val
          !
          num_exp_data = num_exp_data + 1
          !
       ENDIF
       !
       CYCLE
       !
33     EXIT
       !
    ENDDO
    !
    CLOSE(UNIT = 20)
    !
    ! Save experimental data
    !
    ! ALLOCATE ARRAY EXP_DATA
    ALLOCATE(exp_data(1:num_exp_data), STAT = IERR)    
    IF (IERR /= 0) THEN
       WRITE(UNIT = *, FMT = *) "exp_data allocation request denied."
       STOP
    ENDIF
    !
    ! Read pointer
    current_term_val => first_term_val
    !
    index = 1
    DO
       IF (.NOT.ASSOCIATED(current_term_val)) EXIT
       !
       ! Check limits
       IF (current_term_val%exper_level%exp_v >= v_exp_min .AND. &
            current_term_val%exper_level%exp_v <= v_exp_max .AND. &
            current_term_val%exper_level%exp_L >= l_exp_min .AND. &
            current_term_val%exper_level%exp_L <= l_exp_max) THEN
          exp_data(index)%exp_term_value = current_term_val%exper_level%exp_term_value
          exp_data(index)%exp_err_value = current_term_val%exper_level%exp_err_value
          exp_data(index)%exp_v = current_term_val%exper_level%exp_v
          exp_data(index)%exp_L = current_term_val%exper_level%exp_L
          exp_data(index)%spin_symm = current_term_val%exper_level%spin_symm
          exp_data(index)%reference = current_term_val%exper_level%reference
          !
          index = index + 1
          !
       ENDIF
       !
       current_term_val => current_term_val%NEXT
       !
    ENDDO
    !
    num_exp_data = index - 1
    !
    !
    !
  END SUBROUTINE read_experimental_data
  !
  FUNCTION Inv_Part_Ratio(eigenvector)
    !
    ! Subroutine to compute the IPR for a 2DVM eigenstate
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: eigenvector ! Eigenvector components (U(2) basis)
    !
    REAL(KIND = DP) :: Inv_Part_Ratio
    !
    !
    Inv_Part_Ratio = 1.0_DP/SUM(eigenvector**4)
    !
  END FUNCTION Inv_Part_Ratio
  !
  !
  SUBROUTINE SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, Basis, Ham_mat)
    !
    ! Subroutine to save the components of the eigenvectors 
    ! of the Vibron Model Hamiltonian
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum    
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    CHARACTER(LEN=*), INTENT(IN) :: Basis
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: Ham_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: basis_index
    !
    CHARACTER(LEN=65) :: output_filename
    !
    !
    ! Build filename
    IF ( N_val < 10) THEN !to avoid spaces
       WRITE(output_filename, '("eigvec_",A,"_N",I1,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( N_val < 100) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I2,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( N_val < 1000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I3,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( N_val < 10000) THEN 
       WRITE(output_filename, '("eigvec_",A,"_N",I4,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE
       WRITE(output_filename, '("eigvec_",A,"_N",I6,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ENDIF
    !
    OPEN(UNIT = 76, FILE = output_filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    !
    WRITE(UNIT=76, FMT=*) "# N = ", N_val, "; L = ", L_val," ; ", Basis,  " basis"
    WRITE(UNIT=76, FMT=*) "# epsilon = ",  Ham_parameter(1), "; alpha = ", Ham_parameter(2), ";  beta = ", Ham_parameter(3), &
         ";  gamma = ", Ham_parameter(4)
    WRITE(UNIT=76, FMT=*) "# gamma_2 = ", Ham_parameter(5), "; kappa = ", Ham_parameter(6), ";  beta_2 = ", Ham_parameter(7)
    WRITE(UNIT=76, FMT=*) "# gamma_3 = ", Ham_parameter(8), "; kappa_2 = ", Ham_parameter(9), ";  beta_3 = ", Ham_parameter(10)
    WRITE(UNIT=76, FMT=*) "# delta = ", Ham_parameter(11), "; phi = ", Ham_parameter(12), "; beta_4 = ", Ham_parameter(13)
    !
    DO basis_index = 1, dim_block
       WRITE(UNIT=76, FMT=*) Ham_mat(basis_index, 1:dim_block)
    ENDDO
    !
  END SUBROUTINE SAVE_EIGENV_COMPONENTS
  !
  SUBROUTINE SAVE_HAM_MATRIX(N_val, L_val, dim_block, Basis, Ham_mat)
    !
    ! Subroutine to save the Vibron Model Hamiltonian matrix
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum    
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    CHARACTER(LEN=*), INTENT(IN) :: Basis
    !
    REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: Ham_mat ! Hamiltonian matrix
    !
    INTEGER(KIND = I4B) :: basis_index
    !
    CHARACTER(LEN=65) :: output_filename
    !
    !
    ! Build filename
    IF ( N_val < 10) THEN !to avoid spaces
       WRITE(output_filename, '("ham_matrix_",A,"_N",I1,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( dim_block < 100) THEN 
       WRITE(output_filename, '("ham_matrix_",A,"_N",I2,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( dim_block < 1000) THEN 
       WRITE(output_filename, '("ham_matrix_",A,"_N",I3,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE IF ( dim_block < 10000) THEN 
       WRITE(output_filename, '("ham_matrix_",A,"_N",I4,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ELSE
       WRITE(output_filename, '("ham_matrix_",A,"_N",I6,"_L",I1,".dat")') TRIM(Basis), N_val, L_val
    ENDIF
    !
    OPEN(UNIT = 76, FILE = output_filename, STATUS = "UNKNOWN", ACTION = "WRITE")
    !
    WRITE(UNIT=76, FMT=*) "# N = ", N_val, "; L = ", L_val," ; dim = ", dim_block, " ; ", Basis,  " basis"
    WRITE(UNIT=76, FMT=*) "# epsilon = ",  Ham_parameter(1), "; alpha = ", Ham_parameter(2), ";  beta = ", Ham_parameter(3), &
         ";  gamma = ", Ham_parameter(4)
    WRITE(UNIT=76, FMT=*) "# gamma_2 = ", Ham_parameter(5), "; kappa = ", Ham_parameter(6), ";  beta_2 = ", Ham_parameter(7)
    WRITE(UNIT=76, FMT=*) "# gamma_3 = ", Ham_parameter(8), "; kappa_2 = ", Ham_parameter(9), ";  beta_3 = ", Ham_parameter(10)
    WRITE(UNIT=76, FMT=*) "# delta = ", Ham_parameter(11), "; phi = ", Ham_parameter(12), "; beta_4 = ", Ham_parameter(13)
    !
    DO basis_index = 1, dim_block
       WRITE(UNIT=76, FMT=*) Ham_mat(basis_index, 1:dim_block)
    ENDDO
    !
  END SUBROUTINE SAVE_HAM_MATRIX
  !
  !
  SUBROUTINE SO4_HAMILTONIAN_VIBRON_FIT(N_val, L_val, dim_block) 
    !
    ! Subroutine to build the U(4) > SO(4) > SO(3) Hamiltonian
    !
    !  by Currix TM.
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = I4B), INTENT(IN) :: N_val ! U(4) [N]
    !
    INTEGER(KIND = I4B), INTENT(IN) :: L_val ! Angular momentum
    !
    INTEGER(KIND = I4B), INTENT(IN) :: dim_block ! Angular momentum L_val block dimension
    !
    !
    INTEGER(KIND = I4B) :: omega, oper, SO4_state
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Build Hamiltonian
    operator : DO oper = 1, nHamoper
       !
       IF (Iprint > 1) WRITE(*,*) "Operator number ", oper
       !
       IF (ABS(Ham_parameter(oper)) < Zero_Parameter) CYCLE ! Skip not needed operators
       !
       SELECT CASE (oper)
          !
          !
       CASE (1) ! n  --> \epsilon
          !
          ! Diagonal contribution 
          DO SO4_state = 1, dim_block
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(1)*f_n_diag(N_val, omega, L_val)
             !
          ENDDO
          !
          ! Non-Diagonal contribution  omega, omega + 2
          DO SO4_state = 1, dim_block - 1
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 1) = Ham_U4_mat(SO4_state, SO4_state + 1) + &
                  Ham_parameter(1)*f_n_nondiag(N_val, omega, L_val)
             !
          ENDDO
          !
          !
       CASE (2) ! n^2  --> \alpha
          !
          !
          ! Diagonal contribution 
          DO SO4_state = 1, dim_block
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(2)*( &
                  f_n_diag(N_val, omega, L_val)**2 + &
                  f_n_nondiag(N_val, omega, L_val)**2 + &
                  f_n_nondiag(N_val, omega - 2, L_val)**2 &
                  )
             !
          ENDDO
          !
          ! Non-Diagonal contribution omega, omega + 2
          DO SO4_state = 1, dim_block - 1
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 1) = Ham_U4_mat(SO4_state, SO4_state + 1) + &
                  Ham_parameter(2)*( &
                  f_n_diag(N_val, omega, L_val)*f_n_nondiag(N_val, omega, L_val) + &
                  f_n_diag(N_val, omega + 2, L_val)*f_n_nondiag(N_val, omega, L_val) &
                  )
             !
          ENDDO
          !
          ! Non-Diagonal contribution omega, omega + 4
          DO SO4_state = 1, dim_block - 2
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 2) = Ham_U4_mat(SO4_state, SO4_state + 2) + &
                  Ham_parameter(2)*( &
                  f_n_nondiag(N_val, omega, L_val)*f_n_nondiag(N_val, omega + 2, L_val) &
                  )
             !
          ENDDO
          !
       CASE (3) ! (L^2 + D^2)/N   --> \beta ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(3)*REAL(omega*(omega + 2_I4B), DP)/REAL(N_val,DP)
          ENDDO
          !
          !
       CASE (4) ! L^2   --> \gamma ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(4)*REAL(L_val*L_val + L_val)
          ENDDO
          !
          !
       CASE (5) ! L^4   --> \gamma_2 ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(5)*REAL(L_val*L_val*(L_val + 1_I4B)*(L_val + 1_I4B))
          ENDDO
          !
          !
       CASE (6) ! [L^2(L^2 + D^2)]_+/N   --> \kappa ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(6)*REAL(L_val*(L_val + 1_I4B)*omega*(omega + 2_I4B), DP)/REAL(N_val,DP)
          ENDDO
          !
          !
       CASE (7) ! (L^2 + D^2)^2/N^2   --> \beta_2  ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(7)*REAL(omega*(omega + 2_I4B), DP)**2/REAL(N_val,DP)**2
          ENDDO
          !
       CASE (8) ! L^6   --> \gamma_3 ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(8)*REAL(L_val*L_val*L_val*(L_val + 1_I4B)*(L_val + 1_I4B)*(L_val + 1_I4B))
          ENDDO
          !
          !
       CASE (9) ! {[L^2(L^2 + D^2)]_+/N}^2   --> \kappa_2 ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state)  = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(9)*REAL(L_val*(L_val + 1_I4B)*omega*(omega + 2_I4B), DP)**2/REAL(N_val,DP)**2
          ENDDO
          !
          !
       CASE (10) ! (L^2 + D^2)^3/N^3   --> \beta_3  ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(10)*REAL(omega*(omega + 2_I4B), DP)**3/REAL(N_val,DP)**3
          ENDDO
          !
          !
       CASE (11) ! n L^2  --> \delta
          !
          ! Diagonal contribution 
          DO SO4_state = 1, dim_block
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(11)* & 
                  f_n_diag(N_val, omega, L_val)*REAL(L_val*L_val + L_val)
             !
          ENDDO
          !
          ! Non-Diagonal contribution  omega, omega + 2
          DO SO4_state = 1, dim_block - 1
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 1) = Ham_U4_mat(SO4_state, SO4_state + 1) + &
                  Ham_parameter(11)*f_n_nondiag(N_val, omega, L_val)*REAL(L_val*L_val + L_val)
             !
          ENDDO
          !
          !
       CASE (12) ! [n(L^2 + D^2) + (L^2 + D^2)n]/N  --> \phi
          !
          ! Diagonal contribution 
          DO SO4_state = 1, dim_block
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(12)* & 
                  f_n_diag(N_val, omega, L_val)*REAL(omega*(omega + 2_I4B), DP)/REAL(N_val,DP)
             !
          ENDDO
          !
          ! Non-Diagonal contribution  omega, omega + 2
          DO SO4_state = 1, dim_block - 1
             !
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             !
             Ham_U4_mat(SO4_state, SO4_state + 1) = Ham_U4_mat(SO4_state, SO4_state + 1) + &
                  Ham_parameter(12)* &
                  ( &
                  f_n_nondiag(N_val, omega, L_val)*REAL(omega*(omega + 2_I4B), DP)/REAL(N_val,DP) + & 
                  f_n_nondiag(N_val, omega, L_val)*REAL((omega + 2_I4B)*(omega + 4_I4B), DP)/REAL(N_val,DP) &
                  )
             !
          ENDDO
          !
          !
       CASE DEFAULT
          !
          STOP 'You should not be here. Invalid nonzero parameter number. Sayonara baby.'
          !
       CASE (13) ! (L^2 + D^2)^4/N^4   --> \beta_4  ---> DIAGONAL
          !
          DO SO4_state = 1, dim_block
             omega = SO4_Basis(SO4_state)%omega_SO4_val
             Ham_U4_mat(SO4_state, SO4_state) = Ham_U4_mat(SO4_state, SO4_state) + &
                  Ham_parameter(13)*REAL(omega*(omega + 2_I4B), DP)**4/REAL(N_val,DP)**4
          ENDDO
          !
          !
       END SELECT
       !
    ENDDO operator
    !
    !
    !
    !
  END SUBROUTINE SO4_HAMILTONIAN_VIBRON_FIT
  !
END MODULE vibron_u4
