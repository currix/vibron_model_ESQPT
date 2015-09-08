PROGRAM vibron_eigensystem_u3_ipr
  !
  ! Eigenvalues and eigenstates of a U(4) Vibron Model Hamiltonian in the U(4) > U(3) > SO(3) basis
  !
  ! H = epsilon n + alpha n^2 + beta (L^2 + D^2)/N + gamma L^2 + eta L^4 + kappa [L^2(L^2 + D^2)]_+/N
  !     + beta2 [(L^2 + D^2)/N]^2
  !
  USE nrtype
  !
  USE vibron_u4
  !
  ! Lapack 95
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVD
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B) :: L_val
  !
  INTEGER(KIND = I4B) :: dim_block ! L Block dimension
  !
  REAL(KIND = DP) :: epsilon, alpha, beta, gamma, eta, kappa, beta2 ! Hamiltonian parameters
  !
  INTEGER(KIND = I4B) :: U3_state, np 
  !
  !
  ! NAMELISTS
  NAMELIST/par_aux/ Iprint, eigenvec, excitation, save_avec, save_ham
  NAMELIST/par_0/ N_val, L_min, L_max
  NAMELIST/par_1/ epsilon, alpha, beta, gamma, eta, kappa, beta2
  !
  ! 
  ! Initialize time
  CALL CPU_TIME(time_check_ref)
  !
  !
  ! Read parameters
  !
  READ(UNIT = *, NML = par_aux)
  !
  IF (Iprint > 1) PRINT 10
  READ(UNIT = *, NML = par_0)
  !
  IF (Iprint > 1) PRINT 20
  READ(UNIT = *, NML = par_1)
  !
  !
  IF (Iprint > 1) THENºx
     WRITE(UNIT = *, FMT = 5) Iprint, eigenvec, excitation, save_avec, save_ham 
     WRITE(UNIT = *, FMT = 15) N_val, L_min, L_max
     WRITE(UNIT = *, FMT = 25) epsilon, alpha, beta, gamma, eta, kappa, beta2
  ENDIF
  !
  ! TESTS
  IF (L_min > L_max) STOP 'ERROR :: L_MIN > L_MAX. SAYONARA BABY'
  IF (N_val < L_min .OR. N_val < L_max) STOP 'ERROR :: N_VAL < L_VAL. SAYONARA BABY'
  !
  ! HAMILTONIAN PARAMETERS
  !ham_param(1) => epsilon    => n
  Ham_parameter(1) = epsilon
  !ham_param(2) => alpha      => n^2
  Ham_parameter(2) = alpha
  !ham_param(3) => beta       => (D^2 + L^2)/N
  Ham_parameter(3) = beta
  !ham_param(4) => gamma      => L^2
  Ham_parameter(4) = gamma
  !ham_param(5) => eta        => L^4
  Ham_parameter(5) = eta
  !ham_param(6) => kappa      => [L^2(D^2 + L^2)]_+/N
  Ham_parameter(6) = kappa
  !ham_param(7) => beta2      => [(D^2 + L^2)/N]^2
  Ham_parameter(7) = beta2
  !
  !
  ! Compute gs L = 0 energy if excitation && L_min NOT ZERO
  IF (L_min /= 0 .AND. excitation) THEN
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "Computing L = 0 g.s. energy"
     !
     dim_block = DIM_L_BLOCK(N_val, 0_I4B)
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "Block dimension = ", dim_block
     !
     ! ALLOCATE BASIS
     ALLOCATE(U3_Basis(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U3_Basis allocation request denied."
        STOP
     ENDIF
     !
!!!CALL U3_BASIS_VIBRON(N_val, 0_I4B, U3_Basis) ! Build U3 basis
     !
     !
     ! ALLOCATE HAMILTONIAN
     ALLOCATE(Ham_U4_mat(1:dim_block,1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_U4_mat allocation request denied."
        STOP
     ENDIF
     !
     Ham_U4_mat = 0.0_DP
     !
     !
     ! ALLOCATE EIGENVALUES VECTOR
     ALLOCATE(eigenval_vec(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "eigenval_vec allocation request denied."
        STOP
     ENDIF
     !
     eigenval_vec = 0.0_DP
     !
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "local basis dimension -> ", size(U3_Basis)
     !
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = 0 Time (1) :: U3 Basis Building ", time_check - time_check_ref
     time_check_ref = time_check
     !
     ! BUILD HAMILTONIAN
!!!!CALL U3_HAMILTONIAN_VIBRON(N_val, 0_I4B, dim_block, U3_Basis, Ham_U4_mat)
     !
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = 0 ", &
          " Time (2) :: U3 Hamiltonian Building ", time_check - time_check_ref
     time_check_ref = time_check
     !
     ! Diagonalize Hamiltonian matrix (LAPACK95)
     !CALL LA_SYEVD(A=Ham_U4_mat, W=eigenval_vec, JOBZ='N', UPLO='U', IL=1, IU=1)
     !
     GS_energy = eigenval_vec(1)
     !
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "GS_energy = ", GS_energy
     !      
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = 0 ", &
          " Time (3) :: U3 Hamiltonian diagonalization ", time_check - time_check_ref
     time_check_ref = time_check
     !
     !
     ! DEALLOCATE BASIS
     DEALLOCATE(U3_Basis, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U3_Basis deallocation request denied."
        STOP
     ENDIF
     !
     ! DEALLOCATE HAMILTONIAN
     DEALLOCATE(Ham_U4_mat, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_U4_mat deallocation request denied."
        STOP
     ENDIF
     !
     ! DEALLOCATE EIGENVALUES VECTOR
     DEALLOCATE(eigenval_vec, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "eigenval_vec deallocation request denied."
        STOP
     ENDIF
     ! 
  ENDIF
  !
  ! MAIN LOOP
  ang_momentum : DO L_val = L_min, L_max
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "Angular Momentum L = ", L_val
     !
     dim_block = DIM_L_BLOCK(N_val, L_val)
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "Block dimension = ", dim_block
     !
     ! ALLOCATE BASIS
     ALLOCATE(U3_Basis(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U3_Basis allocation request denied."
        STOP
     ENDIF
     !
     CALL U3_BASIS_VIBRON(N_val, L_val) ! Build U3 basis
     !
     !
     ! ALLOCATE HAMILTONIAN
     ALLOCATE(Ham_U4_mat(1:dim_block,1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_U4_mat allocation request denied."
        STOP
     ENDIF
     !
     Ham_U4_mat = 0.0_DP
     !
     !
     ! ALLOCATE EIGENVALUES VECTOR
     ALLOCATE(eigenval_vec(1:dim_block), STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "eigenval_vec allocation request denied."
        STOP
     ENDIF
     !
     eigenval_vec = 0.0_DP
     !
     !
     IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "local basis dimension -> ", size(U3_Basis)
     !
     IF (Iprint > 1) THEN
        DO U3_state = 1, dim_block
           WRITE(UNIT = *, FMT = *) U3_Basis(U3_state) 
        ENDDO
     ENDIF
     !
     !
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, " Time (1) :: U3 Basis Building ", time_check - time_check_ref
     time_check_ref = time_check
     !
     ! BUILD HAMILTONIAN
     CALL U3_HAMILTONIAN_VIBRON(N_val, L_val, dim_block)
     !      
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
          " Time (2) :: U3 Hamiltonian Building ", time_check - time_check_ref
     time_check_ref = time_check
     !
     ! Diagonalize Hamiltonian matrix (LAPACK95)
     !
     IF (Iprint > 0) WRITE(*,*) "Hamiltonian Matrix Diagonalization"
     !
     IF (eigenvec .OR. save_avec) THEN
        IF (Iprint > 0) WRITE(*,*) "Hello (diagonalization proceeding)...."
        CALL LA_SYEVD(A=Ham_U4_mat, W=eigenval_vec, JOBZ='V', UPLO='U')
        IF (Iprint > 0) WRITE(*,*) "(Diagonalization finished) Bye bye...."
     ELSE
        CALL LA_SYEVD(A=Ham_U4_mat, W=eigenval_vec, JOBZ='N', UPLO='U')
     ENDIF
     !
     IF (excitation) THEN
        !
        IF (L_val == 0_I4B) THEN
           GS_energy = eigenval_vec(1)
           IF (Iprint > 0) WRITE(UNIT = *, FMT = *) "GS_energy = ", GS_energy
        ENDIF
        !
        eigenval_vec = eigenval_vec - GS_energy
        !
     ENDIF
     !
     !      
     ! Check time
     CALL CPU_TIME(time_check)
     IF (Iprint == -2) WRITE(UNIT = *, FMT = *) "L_val = ", L_val, &
          " Time (3) :: U3 Hamiltonian diagonalization ", time_check - time_check_ref
     time_check_ref = time_check
     !
     !
     !
     IF (Iprint > 0) THEN
        WRITE(UNIT = *, FMT = *) "L_val = ", L_val
        WRITE(*,*) "index n E_n IPR_n C_n C_max ind_max"
     ENDIF
     !
     DO U3_state = 1, dim_block
        !
        np = U3_Basis(U3_state)%np_U3_val
        !
        ! Output index n energy/N ipr/N C_n_i, i=0,1,2,3 Cmaxj nmaxj j = 1,2,3
        WRITE(UNIT = *, FMT = *)  U3_state, np, eigenval_vec(U3_state)/N_val, &
             Inv_Part_Ratio(Ham_U4_mat(1:dim_block, U3_state))/N_val, &
             Ham_U4_mat(1:4, U3_state), &
             Sorted_Components(Ham_U4_mat(1:dim_block, U3_state),dim_block,4) 
        !
     ENDDO
     !
     ! Save eigenvector components
     IF (save_avec) &
          CALL SAVE_EIGENV_COMPONENTS(N_val, L_val, dim_block, "u3", Ham_U4_mat)
     !
     ! DEALLOCATE BASIS
     DEALLOCATE(U3_Basis, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "U3_Basis deallocation request denied."
        STOP
     ENDIF
     !
     ! DEALLOCATE HAMILTONIAN
     DEALLOCATE(Ham_U4_mat, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "Ham_U4_mat deallocation request denied."
        STOP
     ENDIF
     !
     ! DEALLOCATE EIGENVALUES VECTOR
     DEALLOCATE(eigenval_vec, STAT = IERR)    
     IF (IERR /= 0) THEN
        WRITE(UNIT = *, FMT = *) "eigenval_vec deallocation request denied."
        STOP
     ENDIF
     !
  ENDDO ang_momentum
  !
5 FORMAT(1X, " Iprint = ", I2, "; eigenvec = ", L2, "; excitation = ", L2, "; save_avec = ", L2)
10 FORMAT(1X, "Reading  N_val, L_min, L_max")
15 FORMAT(1X, "N_val = ", I6, "; L_min = ", I6,"; L_max = ", I6)
20 FORMAT(1X, "Reading  epsilon, alpha, beta, gamma, eta, kappa")
25 FORMAT(1X, "epsilon = ", ES14.7, "; alpha = ", ES14.7, /, &
        "beta = ", ES14.7, "; gamma = ", ES14.7, /, &
        "eta = ", ES14.7, "; kappa = ", ES14.7, "; beta2 = ", ES14.7)
  !
  !
  !
END PROGRAM vibron_eigensystem_u3_ipr
