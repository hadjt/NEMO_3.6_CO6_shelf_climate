!----------------------------------------------------------------------
! NEMO system team, System and Interface for oceanic RElocable Nesting
!----------------------------------------------------------------------
!
! MODULE: mpp
!
! DESCRIPTION:
!> @brief
!> This module manage massively parallel processing.
!
!> @details
!> define type TMPP:<br/>
!> @code
!> TYPE(TMPP) :: tl_mpp
!> @endcode
!>
!>    to initialise a mpp structure:<br/>
!> @code
!>    tl_mpp=mpp_init( cd_file, id_mask, 
!>                       [id_niproc,] [id_njproc,] [id_nproc,]
!>                       [id_preci,] [id_precj,] 
!>                       [cd_type,] [id_ew])
!> @endcode
!> or
!> @code
!>    tl_mpp=mpp_init( cd_file, td_var, 
!>                      [id_niproc,] [id_njproc,] [id_nproc,]
!>                      [id_preci,] [id_precj,]
!>                      [cd_type] )
!> @endcode
!> or
!> @code
!>    tl_mpp=mpp_init( td_file [,id_ew] )
!> @endcode
!>       - cd_file is the filename of the global domain file, in which 
!>         MPP will be done (example: Bathymetry) 
!>       - td_file is the file structure of one processor file composing an MPP
!>       - id_mask is the 2D mask of global domain [optional]
!>       - td_var is a variable structure (on T-point) from global domain file.
!>         mask of the domain will be computed using FillValue [optional]
!>       - id_niproc is the number of processor following I-direction to be used
!>         [optional]
!>       - id_njproc is the number of processor following J-direction to be used 
!>         [optional]
!>       - id_nproc is the total number of processor to be used [optional]
!>       - id_preci is the size of the overlap region following I-direction [optional]
!>       - id_precj is the size of the overlap region following J-direction [optional]
!>       - cd_type is the type of files composing MPP [optional]
!>       - id_ew is east-west overlap [optional]<br/> 
!>  
!>    to get mpp name:<br/>
!>    - tl_mpp\%c_name
!>
!>    to get the total number of processor:<br/>
!>    - tl_mpp\%i_nproc
!>
!>    to get the number of processor following I-direction:<br/>
!>    - tl_mpp\%i_niproc
!>
!>    to get the number of processor following J-direction:<br/>
!>    - tl_mpp\%i_njproc
!>
!>    to get the length of the overlap region following I-direction:<br/>
!>    - tl_mpp\%i_preci
!>
!>    to get the length of the overlap region following J-direction:<br/>
!>    - tl_mpp\%i_precj
!>
!>    to get the type of files composing mpp structure:<br/>
!>    - tl_mpp\%c_type
!>
!>    to get the type of the global domain:<br/>
!>    - tl_mpp\%c_dom
!>
!>    MPP dimensions (global domain)<br/>
!>    to get the number of dimensions to be used in mpp strcuture:<br/>
!>    - tl_mpp\%i_ndim
!>
!>    to get the array of dimension structure (4 elts) associated to the
!>    mpp structure:<br/>
!>    - tl_mpp\%t_dim(:)
!>
!>    MPP processor (files composing domain)<br/>
!>    - tl_mpp\%t_proc(:)
!>
!>    to clean a mpp structure:<br/>
!> @code
!>    CALL mpp_clean(tl_mpp)
!> @endcode
!>
!>    to print information about mpp:<br/>
!> @code
!>    CALL mpp_print(tl_mpp)
!> @endcode
!>
!>    to add variable to mpp:<br/>
!> @code
!>    CALL mpp_add_var(td_mpp, td_var)
!> @endcode
!>       - td_var is a variable structure
!>
!>    to add dimension to mpp:<br/>
!> @code
!>    CALL mpp_add_dim(td_mpp, td_dim)
!> @endcode
!>       - td_dim is a dimension structure
!>
!>    to add attribute to mpp:<br/>
!> @code
!>    CALL mpp_add_att(td_mpp, td_att)
!> @endcode
!>       - td_att is a attribute structure
!>
!>    to delete variable from mpp:<br/>
!> @code
!>    CALL mpp_del_var(td_mpp, td_var)
!> @endcode
!>    or
!> @code
!>    CALL mpp_del_var(td_mpp, cd_name)
!> @endcode
!>       - td_var is a variable structure
!>       - cd_name is variable name or standard name
!>
!>    to delete dimension from mpp:<br/>
!> @code
!>    CALL mpp_del_dim(td_mpp, td_dim)
!> @endcode
!>       - td_dim is a dimension structure
!>
!>    to delete attribute from mpp:<br/>
!> @code
!>    CALL mpp_del_att(td_mpp, td_att)
!> @endcode
!>    or
!> @code
!>    CALL mpp_del_att(td_mpp, cd_name)
!> @endcode
!>       - td_att is a attribute structure
!>       - cd_name is attribute name
!>
!>    to overwrite variable to mpp:<br/>
!> @code
!>    CALL mpp_move_var(td_mpp, td_var)
!> @endcode
!>       - td_var is a variable structure
!>
!>    to overwrite dimension to mpp:<br/>
!> @code
!>    CALL mpp_move_dim(td_mpp, td_dim)
!> @endcode
!>       - td_dim is a dimension structure
!>
!>    to overwrite attribute to mpp:<br/>
!> @code
!>    CALL mpp_move_att(td_mpp, td_att)
!> @endcode
!>       - td_att is a attribute structure
!>
!>    to determine domain decomposition type:<br/>
!> @code
!>    CALL mpp_get_dom(td_mpp)
!> @endcode
!>
!>    to get processors to be used:<br/>
!> @code
!>    CALL mpp_get_use( td_mpp, id_imin, id_imax, & 
!>    &                         id_jmin, id_jmax )
!> @endcode
!>       - id_imin 
!>       - id_imax 
!>       - id_jmin 
!>       - id_jmax 
!>
!>    to get sub domains which form global domain contour:<br/>
!> @code
!>    CALL mpp_get_contour( td_mpp ) 
!> @endcode
!>
!>    to get global domain indices of one processor:<br/>
!> @code
!>    il_ind(1:4)=mpp_get_proc_index( td_mpp, id_procid )
!> @endcode
!>       - il_ind(1:4) are global domain indices (i1,i2,j1,j2)
!>       - id_procid is the processor id
!>
!>    to get the processor domain size:<br/>
!> @code
!>    il_size(1:2)=mpp_get_proc_size( td_mpp, id_procid )
!> @endcode
!>       - il_size(1:2) are the size of domain following I and J
!>       - id_procid is the processor id
!>
!> @author
!>  J.Paul
! REVISION HISTORY:
!> @date November, 2013 - Initial Version
!> @date November, 2014 - Fix memory leaks bug
!
!> @note Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!----------------------------------------------------------------------
MODULE mpp
   USE global                          ! global parameter
   USE kind                            ! F90 kind parameter
   USE logger                          ! log file manager
   USE fct                             ! basic useful function
   USE dim                             ! dimension manager
   USE att                             ! attribute manager
   USE var                             ! variable manager
   USE file                            ! file manager
   USE iom                             ! I/O manager
   IMPLICIT NONE
   ! NOTE_avoid_public_variables_if_possible

   ! type and variable
   PUBLIC :: TMPP       !< mpp structure

   ! function and subroutine
   PUBLIC :: mpp_copy           !< copy mpp structure
   PUBLIC :: mpp_init           !< initialise mpp structure
   PUBLIC :: mpp_clean          !< clean mpp strcuture
   PUBLIC :: mpp_print          !< print information about mpp structure
   PUBLIC :: mpp_add_var        !< split/add one variable strucutre in mpp structure
   PUBLIC :: mpp_add_dim        !< add one dimension to mpp structure
   PUBLIC :: mpp_add_att        !< add one attribute strucutre in mpp structure
   PUBLIC :: mpp_del_var        !< delete one variable strucutre in mpp structure
   PUBLIC :: mpp_del_dim        !< delete one dimension strucutre in mpp structure
   PUBLIC :: mpp_del_att        !< delete one attribute strucutre in mpp structure
   PUBLIC :: mpp_move_var       !< overwrite variable structure in mpp structure
   PUBLIC :: mpp_move_dim       !< overwrite one dimension strucutre in mpp structure
   PUBLIC :: mpp_move_att       !< overwrite one attribute strucutre in mpp structure
   PUBLIC :: mpp_recombine_var  !< recombine variable from mpp structure
   PUBLIC :: mpp_get_index      !< return index of mpp 

   PUBLIC :: mpp_get_dom        !< determine domain decomposition type (full, overlap, noverlap)
   PUBLIC :: mpp_get_use        !< get sub domains to be used (which cover "zoom domain")
   PUBLIC :: mpp_get_contour    !< get sub domains which form global domain contour
   PUBLIC :: mpp_get_proc_index !< get processor domain indices
   PUBLIC :: mpp_get_proc_size  !< get processor domain size

   PRIVATE :: mpp__add_proc            ! add one proc strucutre in mpp structure
   PRIVATE :: mpp__del_proc            ! delete one proc strucutre in mpp structure
   PRIVATE :: mpp__del_proc_id         ! delete one proc strucutre in mpp structure, given procesor id
   PRIVATE :: mpp__del_proc_str        ! delete one proc strucutre in mpp structure, given procesor file structure 
   PRIVATE :: mpp__move_proc           ! overwrite proc strucutre in mpp structure
   PRIVATE :: mpp__compute             ! compute domain decomposition
   PRIVATE :: mpp__del_land            ! remove land sub domain from domain decomposition
   PRIVATE :: mpp__optimiz             ! compute optimum domain decomposition
   PRIVATE :: mpp__land_proc           ! check if processor is a land processor
   PRIVATE :: mpp__check_dim           ! check mpp structure dimension with proc or variable dimension
   PRIVATE :: mpp__check_proc_dim      ! check if processor and mpp structure use same dimension
   PRIVATE :: mpp__check_var_dim       ! check if variable  and mpp structure use same dimension
   PRIVATE :: mpp__del_var_name        ! delete variable in mpp structure, given variable name
   PRIVATE :: mpp__del_var_mpp         ! delete all variable in mpp structure
   PRIVATE :: mpp__del_var_str         ! delete variable in mpp structure, given variable structure
   PRIVATE :: mpp__del_att_name        ! delete variable in mpp structure, given variable name
   PRIVATE :: mpp__del_att_str         ! delete variable in mpp structure, given variable structure
   PRIVATE :: mpp__split_var           ! extract variable part that will be written in processor 
   PRIVATE :: mpp__copy_unit           ! copy mpp structure
   PRIVATE :: mpp__copy_arr            ! copy array of mpp structure
   PRIVATE :: mpp__get_use_unit        ! get sub domains to be used (which cover "zoom domain")
   PRIVATE :: mpp__init_mask           ! initialise mpp structure, given file name
   PRIVATE :: mpp__init_var            ! initialise mpp structure, given variable strcuture
   PRIVATE :: mpp__init_file           ! initialise a mpp structure, given file structure 
   PRIVATE :: mpp__init_file_cdf       ! initialise a mpp structure with cdf file
   PRIVATE :: mpp__init_file_rstdimg   ! initialise a mpp structure with rstdimg file
   PRIVATE :: mpp__clean_unit          ! clean mpp strcuture
   PRIVATE :: mpp__clean_arr           ! clean array of mpp strcuture

   TYPE TMPP !< mpp structure

      ! general 
      CHARACTER(LEN=lc)                  :: c_name = ''   !< base name 
      INTEGER(i4)                        :: i_id   = 0    !< mpp id

      INTEGER(i4)                        :: i_niproc = 0  !< number of processors following i
      INTEGER(i4)                        :: i_njproc = 0  !< number of processors following j
      INTEGER(i4)                        :: i_nproc  = 0  !< total number of proccessors used
      INTEGER(i4)                        :: i_preci = 1   !< i-direction overlap region length
      INTEGER(i4)                        :: i_precj = 1   !< j-direction overlap region length
      INTEGER(i4)                        :: i_ew    = -1  !< east-west overlap
      INTEGER(i4)                        :: i_perio = -1  !< NEMO periodicity index
      INTEGER(i4)                        :: i_pivot = -1  !< NEMO pivot point index F(0),T(1)

      CHARACTER(LEN=lc)                  :: c_type = ''   !< type of the files (cdf, cdf4, dimg)
      CHARACTER(LEN=lc)                  :: c_dom  = ''   !< type of domain (full, overlap, nooverlap)

      INTEGER(i4)                        :: i_ndim = 0    !< number of dimensions used in mpp
      TYPE(TDIM),  DIMENSION(ip_maxdim)  :: t_dim         !< global domain dimension

      TYPE(TFILE), DIMENSION(:), POINTER :: t_proc => NULL()     !< files/processors composing mpp

   END TYPE

   INTERFACE mpp_get_use
      MODULE PROCEDURE mpp__get_use_unit 
   END INTERFACE mpp_get_use

   INTERFACE mpp_clean
      MODULE PROCEDURE mpp__clean_unit 
      MODULE PROCEDURE mpp__clean_arr   
   END INTERFACE mpp_clean

   INTERFACE mpp__check_dim
      MODULE PROCEDURE mpp__check_proc_dim !< check if processor and mpp structure use same dimension
      MODULE PROCEDURE mpp__check_var_dim  !< check if variable  and mpp structure use same dimension
   END INTERFACE mpp__check_dim

   INTERFACE mpp__del_proc
      MODULE PROCEDURE mpp__del_proc_id
      MODULE PROCEDURE mpp__del_proc_str
   END INTERFACE mpp__del_proc

   INTERFACE mpp_del_var
      MODULE PROCEDURE mpp__del_var_name
      MODULE PROCEDURE mpp__del_var_str
      MODULE PROCEDURE mpp__del_var_mpp
   END INTERFACE mpp_del_var

   INTERFACE mpp_del_att
      MODULE PROCEDURE mpp__del_att_name
      MODULE PROCEDURE mpp__del_att_str
   END INTERFACE mpp_del_att

   INTERFACE mpp_init
      MODULE PROCEDURE mpp__init_mask
      MODULE PROCEDURE mpp__init_var
      MODULE PROCEDURE mpp__init_file
   END INTERFACE mpp_init

   INTERFACE mpp_copy
      MODULE PROCEDURE mpp__copy_unit  ! copy mpp structure
      MODULE PROCEDURE mpp__copy_arr   ! copy array of mpp structure
   END INTERFACE

CONTAINS
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine copy mpp structure in another one
   !> @details 
   !> mpp file are copied in a temporary array, 
   !> so input and output mpp structure do not point on the same 
   !> "memory cell", and so on are independant. 
   !>
   !> @warning do not use on the output of a function who create or read an
   !> structure (ex: tl_file=file_copy(file_init()) is forbidden).
   !> This will create memory leaks.
   !> @warning to avoid infinite loop, do not use any function inside 
   !> this subroutine
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !> @date November, 2014
   !>    - use function instead of overload assignment operator 
   !> (to avoid memory leak)
   !
   !> @param[in] td_mpp   mpp structure
   !> @return copy of input mpp structure
   !-------------------------------------------------------------------
   FUNCTION mpp__copy_unit( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(IN)  :: td_mpp
      ! function
      TYPE(TMPP) :: mpp__copy_unit

      ! local variable
      TYPE(TFILE) :: tl_file

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      CALL logger_trace("MPP COPY: "//TRIM(td_mpp%c_name)//" in "//&
      &  TRIM(mpp__copy_unit%c_name))

      ! copy mpp variable
      mpp__copy_unit%c_name     = TRIM(td_mpp%c_name)
      mpp__copy_unit%i_id       = td_mpp%i_id
      mpp__copy_unit%i_niproc   = td_mpp%i_niproc
      mpp__copy_unit%i_njproc   = td_mpp%i_njproc
      mpp__copy_unit%i_nproc    = td_mpp%i_nproc
      mpp__copy_unit%i_preci    = td_mpp%i_preci
      mpp__copy_unit%i_precj    = td_mpp%i_precj
      mpp__copy_unit%c_type     = TRIM(td_mpp%c_type)
      mpp__copy_unit%c_dom      = TRIM(td_mpp%c_dom)
      mpp__copy_unit%i_ndim     = td_mpp%i_ndim
      mpp__copy_unit%i_ew       = td_mpp%i_ew
      mpp__copy_unit%i_perio    = td_mpp%i_perio
      mpp__copy_unit%i_pivot    = td_mpp%i_pivot

      ! copy dimension
      mpp__copy_unit%t_dim(:) = dim_copy(td_mpp%t_dim(:))
      
      ! copy file structure
      IF( ASSOCIATED(mpp__copy_unit%t_proc) )THEN
         CALL file_clean(mpp__copy_unit%t_proc(:))
         DEALLOCATE(mpp__copy_unit%t_proc)
      ENDIF
      IF( ASSOCIATED(td_mpp%t_proc) .AND. mpp__copy_unit%i_nproc > 0 )THEN
         ALLOCATE( mpp__copy_unit%t_proc(mpp__copy_unit%i_nproc) )
         DO ji=1,mpp__copy_unit%i_nproc
            tl_file = file_copy(td_mpp%t_proc(ji))
            mpp__copy_unit%t_proc(ji) = file_copy(tl_file)
         ENDDO
         ! clean
         CALL file_clean(tl_file)
      ENDIF

   END FUNCTION mpp__copy_unit
   !-------------------------------------------------------------------
   !> @brief
   !> This subroutine copy an array of mpp structure in another one
   !> @details 
   !> mpp file are copied in a temporary array, 
   !> so input and output mpp structure do not point on the same 
   !> "memory cell", and so on are independant. 
   !>
   !> @warning do not use on the output of a function who create or read an
   !> structure (ex: tl_file=file_copy(file_init()) is forbidden).
   !> This will create memory leaks.
   !> @warning to avoid infinite loop, do not use any function inside 
   !> this subroutine
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !> @date November, 2014
   !>    - use function instead of overload assignment operator 
   !> (to avoid memory leak)
   !>
   !> @param[in] td_mpp   mpp structure
   !> @return copy of input array of mpp structure
   !-------------------------------------------------------------------
   FUNCTION mpp__copy_arr( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), DIMENSION(:), INTENT(IN)  :: td_mpp
      ! function
      TYPE(TMPP), DIMENSION(SIZE(td_mpp(:))) :: mpp__copy_arr

      ! local variable
      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      DO ji=1,SIZE(td_mpp(:))
         mpp__copy_arr(ji)=mpp_copy(td_mpp(ji))
      ENDDO

   END FUNCTION mpp__copy_arr
   !-------------------------------------------------------------------
   !> @brief This subroutine print some information about mpp strucutre.
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !
   !> @param[in] td_mpp mpp structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp_print(td_mpp)
      IMPLICIT NONE

      ! Argument      
      TYPE(TMPP), INTENT(IN) :: td_mpp

      ! local variable
      INTEGER(i4), PARAMETER :: il_freq = 4

      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_proc
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_lci
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_lcj

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk
      INTEGER(i4) :: jl
      INTEGER(i4) :: jm
      !----------------------------------------------------------------

      WRITE(*,'((a,a),2(/3x,a,a),9(/3x,a,i0))')&
      &  "MPP : ",TRIM(td_mpp%c_name), &
      &  " type   : ",TRIM(td_mpp%c_type), &
      &  " dom    : ",TRIM(td_mpp%c_dom), &
      &  " nproc  : ",td_mpp%i_nproc, &
      &  " niproc : ",td_mpp%i_niproc, &
      &  " njproc : ",td_mpp%i_njproc, &
      &  " preci  : ",td_mpp%i_preci, &
      &  " precj  : ",td_mpp%i_precj, &
      &  " ndim   : ",td_mpp%i_ndim,  &
      &  " overlap: ",td_mpp%i_ew, &
      &  " perio  : ",td_mpp%i_perio, &
      &  " pivot  : ",td_mpp%i_pivot

      ! print dimension
      IF(  td_mpp%i_ndim /= 0 )THEN
         WRITE(*,'(/a)') " MPP dimension"
         DO ji=1,ip_maxdim
            IF( td_mpp%t_dim(ji)%l_use )THEN
               CALL dim_print(td_mpp%t_dim(ji))
            ENDIF
         ENDDO
      ENDIF

      ! print file
      IF( td_mpp%i_nproc /= 0 .AND. ASSOCIATED(td_mpp%t_proc) )THEN
         IF( ALL( td_mpp%t_proc(:)%i_iind==0 ) .OR. &
         &   ALL( td_mpp%t_proc(:)%i_jind==0 ) )THEN

            DO ji=1,td_mpp%i_nproc
               CALL file_print(td_mpp%t_proc(ji))
               WRITE(*,'((a),(/3x,a,i0),2(/3x,a,a),4(/3x,a,i0,a,i0)/)')&
               &  " Domain decomposition : ", &
               &  " id          : ",td_mpp%t_proc(ji)%i_pid, &
               &  " used        : ",TRIM(fct_str(td_mpp%t_proc(ji)%l_use)), &
               &  " contour     : ",TRIM(fct_str(td_mpp%t_proc(ji)%l_ctr)), &
               &  " left-bottom : ",td_mpp%t_proc(ji)%i_impp,', ',&
               &  td_mpp%t_proc(ji)%i_jmpp, &
               &  " dimension   : ",td_mpp%t_proc(ji)%i_lci,' x ',&
               &  td_mpp%t_proc(ji)%i_lcj, &
               &  " first indoor indices : ",td_mpp%t_proc(ji)%i_ldi,', ',&
               &  td_mpp%t_proc(ji)%i_ldj, &
               &  " last  indoor indices : ",td_mpp%t_proc(ji)%i_lei,', ',&
               &  td_mpp%t_proc(ji)%i_lej

            ENDDO

            IF( td_mpp%t_proc(1)%i_nvar > 0 )THEN
               WRITE(*,'(/a)') " Variable(s) used : "
               DO ji=1,td_mpp%t_proc(1)%i_nvar
                  WRITE(*,'(3x,a)') TRIM(td_mpp%t_proc(1)%t_var(ji)%c_name) 
               ENDDO
            ENDIF

         ELSE

            DO ji=1,td_mpp%i_nproc
               WRITE(*,'((a, a),(/3x,a,i0),(/3x,a,a),4(/3x,a,i0,a,i0)/)')&
               &  " Domain decomposition : ",TRIM(td_mpp%t_proc(ji)%c_name),&
               &  " id          : ",td_mpp%t_proc(ji)%i_pid, &
               &  " used        : ",TRIM(fct_str(td_mpp%t_proc(ji)%l_use)),&
               &  " left-bottom : ",td_mpp%t_proc(ji)%i_impp,', ',&
               &  td_mpp%t_proc(ji)%i_jmpp, &
               &  " dimension   : ",td_mpp%t_proc(ji)%i_lci,' x ',&
               &  td_mpp%t_proc(ji)%i_lcj, &
               &  " first indoor indices : ",td_mpp%t_proc(ji)%i_ldi,',',&
               &  td_mpp%t_proc(ji)%i_ldj, &
               &  " last  indoor indices : ",td_mpp%t_proc(ji)%i_lei,', ',&
               &  td_mpp%t_proc(ji)%i_lej

            ENDDO
            
            IF( td_mpp%t_proc(1)%i_nvar > 0 )THEN
               WRITE(*,'(/a)') " Variable(s) used : "
               DO ji=1,td_mpp%t_proc(1)%i_nvar
                  WRITE(*,'(3x,a)') TRIM(td_mpp%t_proc(1)%t_var(ji)%c_name) 
               ENDDO
            ENDIF

            ALLOCATE( il_proc(td_mpp%i_niproc,td_mpp%i_njproc) )
            ALLOCATE( il_lci(td_mpp%i_niproc,td_mpp%i_njproc) )
            ALLOCATE( il_lcj(td_mpp%i_niproc,td_mpp%i_njproc) )

            DO jk=1,td_mpp%i_nproc
               ji=td_mpp%t_proc(jk)%i_iind
               jj=td_mpp%t_proc(jk)%i_jind
               il_proc(ji,jj)=jk
               il_lci(ji,jj)=td_mpp%t_proc(jk)%i_lci
               il_lcj(ji,jj)=td_mpp%t_proc(jk)%i_lcj
            ENDDO

            jl = 1
            DO jk = 1,(td_mpp%i_niproc-1)/il_freq+1
               jm = MIN(td_mpp%i_niproc, jl+il_freq-1)
               WRITE(*,*)
               WRITE(*,9401) (ji, ji = jl,jm)
               WRITE(*,9400) ('***', ji = jl,jm-1)
               DO jj = 1, td_mpp%i_njproc
                  WRITE(*,9403) ('   ', ji = jl,jm-1)
                  WRITE(*,9402) jj, ( il_lci(ji,jj), il_lcj(ji,jj), ji = jl,jm)
                  WRITE(*,9404) (il_proc(ji,jj), ji= jl,jm)
                  WRITE(*,9403) ('   ', ji = jl,jm-1)
                  WRITE(*,9400) ('***', ji = jl,jm-1)
               ENDDO
               jl = jl+il_freq
            ENDDO
         
            DEALLOCATE( il_proc )
            DEALLOCATE( il_lci )
            DEALLOCATE( il_lcj )

         ENDIF
      ELSE
         WRITE(*,'(/a)') " Domain decomposition : none"
      ENDIF


9400   FORMAT('     ***',20('*************',a3))
9403   FORMAT('     *     ',20('         *   ',a3))
9401   FORMAT('        ',20('   ',i3,'          '))
9402   FORMAT(' ',i3,' *  ',20(i0,'  x',i0,'   *   '))
9404   FORMAT('     *  ',20('      ',i3,'   *   '))

   END SUBROUTINE mpp_print
   !-------------------------------------------------------------------
   !> @brief
   !> This function initialise mpp structure, given file name, 
   !> and optionaly mask and number of processor following I and J
   !> @detail
   !> - If no total number of processor is defined (id_nproc), optimize 
   !> the domain decomposition (look for the domain decomposition with 
   !> the most land processor to remove)
   !> - length of the overlap region (id_preci, id_precj) could be specify
   !> in I and J direction (default value is 1)
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[in] cd_file   file name of one file composing mpp domain
   !> @param[in] id_mask   domain mask
   !> @param[in] id_niproc number of processors following i
   !> @param[in] id_njproc number of processors following j
   !> @param[in] id_nproc  total number of processors
   !> @param[in] id_preci  i-direction overlap region
   !> @param[in] id_precj  j-direction overlap region
   !> @param[in] cd_type   type of the files (cdf, cdf4, dimg)
   !> @param[in] id_ew     east-west overlap
   !> @param[in] id_perio  NEMO periodicity index
   !> @param[in] id_pivot  NEMO pivot point index F(0),T(1)
   !> @return mpp structure
   !-------------------------------------------------------------------
   TYPE(TMPP) FUNCTION mpp__init_mask(cd_file, id_mask,              &
   &                                  id_niproc, id_njproc, id_nproc,&
   &                                  id_preci, id_precj,            &
                                      cd_type, id_ew, id_perio, id_pivot)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*),            INTENT(IN) :: cd_file
      INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: id_mask
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_niproc
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_njproc
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_nproc
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_preci
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_precj
      CHARACTER(LEN=*),            INTENT(IN), OPTIONAL :: cd_type
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_ew
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_perio
      INTEGER(i4),                 INTENT(IN), OPTIONAL :: id_pivot

      ! local variable
      CHARACTER(LEN=lc)                :: cl_type

      INTEGER(i4)      , DIMENSION(2) :: il_shape

      TYPE(TDIM)                      :: tl_dim

      TYPE(TATT)                      :: tl_att
      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      ! clean mpp
      CALL mpp_clean(mpp__init_mask)

      ! check type
      cl_type=''
      IF( PRESENT(cd_type) ) cl_type=TRIM(ADJUSTL(cd_type))

      IF( TRIM(cl_type) /= '' )THEN
         SELECT CASE(TRIM(cd_type))
            CASE('cdf')
               mpp__init_mask%c_type='cdf'
            CASE('dimg')
               mpp__init_mask%c_type='dimg'
            CASE DEFAULT
               CALL logger_warn( "MPP INIT: type "//TRIM(cd_type)//&
               & " unknown. type dimg will be used for mpp "//&
               &  TRIM(mpp__init_mask%c_name) )
               mpp__init_mask%c_type='dimg'
         END SELECT
      ELSE
         mpp__init_mask%c_type=TRIM(file_get_type(cd_file))
      ENDIF

      ! get mpp name
      mpp__init_mask%c_name=TRIM(file_rename(cd_file))

      ! get global domain dimension
      il_shape(:)=SHAPE(id_mask)

      tl_dim=dim_init('X',il_shape(1))
      CALL mpp_add_dim(mpp__init_mask, tl_dim)

      tl_dim=dim_init('Y',il_shape(2))
      CALL mpp_add_dim(mpp__init_mask, tl_dim)

      ! clean
      CALL dim_clean(tl_dim)

      IF( (       PRESENT(id_niproc)  .AND. (.NOT. PRESENT(id_njproc))) .OR. &
          ((.NOT. PRESENT(id_niproc)) .AND.        PRESENT(id_njproc) ) )THEN
          CALL logger_warn( "MPP INIT: number of processors following I and J "//&
          & "should be both specified")
      ELSE
         ! get number of processors following I and J
         IF( PRESENT(id_niproc) ) mpp__init_mask%i_niproc=id_niproc
         IF( PRESENT(id_njproc) ) mpp__init_mask%i_njproc=id_njproc
      ENDIF

      ! get maximum number of processors to be used
      IF( PRESENT(id_nproc) ) mpp__init_mask%i_nproc = id_nproc

      ! get overlap region length
      IF( PRESENT(id_preci) ) mpp__init_mask%i_preci= id_preci
      IF( PRESENT(id_precj) ) mpp__init_mask%i_precj= id_precj

      ! east-west overlap
      IF( PRESENT(id_ew) ) mpp__init_mask%i_ew= id_ew
      ! NEMO periodicity
      IF( PRESENT(id_perio) ) mpp__init_mask%i_perio= id_perio
      IF( PRESENT(id_pivot) ) mpp__init_mask%i_pivot= id_pivot

      IF( mpp__init_mask%i_nproc  /= 0 .AND. &
      &   mpp__init_mask%i_niproc /= 0 .AND. &
      &   mpp__init_mask%i_njproc /= 0 .AND. &
      &   mpp__init_mask%i_nproc > &
      &   mpp__init_mask%i_niproc * mpp__init_mask%i_njproc )THEN

         CALL logger_error("MPP INIT: invalid domain decomposition ")
         CALL logger_debug("MPP INIT: "// &
         & TRIM(fct_str(mpp__init_mask%i_nproc))//" > "//&
         & TRIM(fct_str(mpp__init_mask%i_niproc))//" x "//&
         & TRIM(fct_str(mpp__init_mask%i_njproc)) )

      ELSE

         IF( mpp__init_mask%i_niproc /= 0 .AND. &
         &   mpp__init_mask%i_njproc /= 0 )THEN
            ! compute domain decomposition
            CALL mpp__compute( mpp__init_mask )
            ! remove land sub domain
            CALL mpp__del_land( mpp__init_mask, id_mask )
         ELSEIF( mpp__init_mask%i_nproc  /= 0 )THEN
            ! optimiz
            CALL mpp__optimiz( mpp__init_mask, id_mask )

         ELSE
            CALL logger_warn("MPP INIT: number of processor to be used "//&
            &                "not specify. force to one.")
            mpp__init_mask%i_nproc  = 1
            ! optimiz
            CALL mpp__optimiz( mpp__init_mask, id_mask )
         ENDIF
         CALL logger_info("MPP INIT: domain decoposition : "//&
         &  'niproc('//TRIM(fct_str(mpp__init_mask%i_niproc))//') * '//&
         &  'njproc('//TRIM(fct_str(mpp__init_mask%i_njproc))//') = '//&
         &  'nproc('//TRIM(fct_str(mpp__init_mask%i_nproc))//')' )

         ! get domain type
         CALL mpp_get_dom( mpp__init_mask )

         DO ji=1,mpp__init_mask%i_nproc

            ! get processor size
            il_shape(:)=mpp_get_proc_size( mpp__init_mask, ji )

            tl_dim=dim_init('X',il_shape(1))
            CALL file_move_dim(mpp__init_mask%t_proc(ji), tl_dim)

            tl_dim=dim_init('Y',il_shape(2))
            CALL file_move_dim(mpp__init_mask%t_proc(ji), tl_dim)            

            ! add type
            mpp__init_mask%t_proc(ji)%c_type=TRIM(mpp__init_mask%c_type)

            ! clean
            CALL dim_clean(tl_dim)
         ENDDO

         ! add global attribute
         tl_att=att_init("DOMAIN_number_total",mpp__init_mask%i_nproc)
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init("DOMAIN_I_number_total",mpp__init_mask%i_niproc)
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init("DOMAIN_J_number_total",mpp__init_mask%i_njproc)
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init("DOMAIN_size_global",mpp__init_mask%t_dim(1:2)%i_len)
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_I_position_first", &
         &                mpp__init_mask%t_proc(:)%i_impp )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_J_position_first", &
         &                mpp__init_mask%t_proc(:)%i_jmpp )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_I_position_last", &
         &                mpp__init_mask%t_proc(:)%i_lci )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_J_position_last", &
         &                mpp__init_mask%t_proc(:)%i_lcj )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_I_halo_size_start", &
         &                mpp__init_mask%t_proc(:)%i_ldi )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_J_halo_size_start", &
         &                mpp__init_mask%t_proc(:)%i_ldj )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_I_halo_size_end", &
         &                mpp__init_mask%t_proc(:)%i_lei )
         CALL mpp_add_att(mpp__init_mask, tl_att)

         tl_att=att_init( "DOMAIN_J_halo_size_end", &
         &                mpp__init_mask%t_proc(:)%i_lej )
         CALL mpp_add_att(mpp__init_mask, tl_att)         

         ! clean
         CALL att_clean(tl_att)
      ENDIF

   END FUNCTION mpp__init_mask
   !-------------------------------------------------------------------
   !> @brief
   !> This function initialise mpp structure, given variable strcuture 
   !> and optionaly number of processor following I and J
   !> @detail
   !> - If no total number of processor is defined (id_nproc), optimize 
   !> the domain decomposition (look for the domain decomposition with 
   !> the most land processor to remove)
   !> - length of the overlap region (id_preci, id_precj) could be specify
   !> in I and J direction (default value is 1)
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[in] cd_file   file name of one file composing mpp domain
   !> @param[in] td_var    variable structure
   !> @param[in] id_niproc number of processors following i
   !> @param[in] id_njproc number of processors following j
   !> @param[in] id_nproc  total number of processors
   !> @param[in] id_preci  i-direction overlap region
   !> @param[in] id_precj  j-direction overlap region
   !> @param[in] cd_type   type of the files (cdf, cdf4, dimg)
   !> @param[in] id_perio  NEMO periodicity index
   !> @param[in] id_pivot  NEMO pivot point index F(0),T(1)
   !> @return mpp structure
   !-------------------------------------------------------------------
   TYPE(TMPP) FUNCTION mpp__init_var( cd_file, td_var,               &
   &                                  id_niproc, id_njproc, id_nproc,&
   &                                  id_preci, id_precj, cd_type,   &
   &                                  id_perio, id_pivot )
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN) :: cd_file
      TYPE(TVAR),       INTENT(IN) :: td_var
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_niproc
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_njproc
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_nproc
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_preci
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_precj
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cd_type
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_perio
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_pivot

      ! local variable
      INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: il_mask
      !----------------------------------------------------------------

      IF( ASSOCIATED(td_var%d_value) )THEN
         ALLOCATE( il_mask(td_var%t_dim(1)%i_len, &
         &                 td_var%t_dim(2)%i_len, &
         &                 td_var%t_dim(3)%i_len) )
         il_mask(:,:,:)=var_get_mask(td_var)
         
         mpp__init_var=mpp_init( cd_file, il_mask(:,:,1),       &
         &                       id_niproc, id_njproc, id_nproc,&
         &                       id_preci, id_precj, cd_type,   &
         &                       id_ew=td_var%i_ew, &
         &                       id_perio=id_perio, id_pivot=id_pivot)

         DEALLOCATE(il_mask)
      ELSE
         CALL logger_error("MPP INIT: variable value not define.")
      ENDIF

   END FUNCTION mpp__init_var
   !-------------------------------------------------------------------
   !> @brief This function initalise a mpp structure given file structure. 
   !> @details 
   !> It reads restart dimg files, or some netcdf files.
   !>
   !> @warning 
   !>  netcdf file must contains some attributes:
   !>    - DOMAIN_number_total 
   !>    - DOMAIN_size_global
   !>    - DOMAIN_number
   !>    - DOMAIN_position_first
   !>    - DOMAIN_position_last
   !>    - DOMAIN_halo_size_start
   !>    - DOMAIN_halo_size_end
   !>  or the file is assume to be no mpp file.
   !>  
   !> 
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !
   !> @param[in] td_file   file strcuture
   !> @param[in] id_ew     east-west overlap
   !> @param[in] id_perio  NEMO periodicity index
   !> @param[in] id_pivot  NEMO pivot point index F(0),T(1)
   !> @return mpp structure
   !-------------------------------------------------------------------
   TYPE(TMPP) FUNCTION mpp__init_file( td_file, id_ew, id_perio, id_pivot )
      IMPLICIT NONE

      ! Argument
      TYPE(TFILE), INTENT(IN) :: td_file
      INTEGER(i4), INTENT(IN), OPTIONAL :: id_ew
      INTEGER(i4), INTENT(IN), OPTIONAL :: id_perio
      INTEGER(i4), INTENT(IN), OPTIONAL :: id_pivot

      ! local variable
      TYPE(TMPP)  :: tl_mpp
      
      TYPE(TFILE) :: tl_file
      
      TYPE(TDIM)  :: tl_dim

      TYPE(TATT)  :: tl_att

      INTEGER(i4) :: il_nproc
      INTEGER(i4) :: il_attid

      INTEGER(i4), DIMENSION(2) :: il_shape
      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      ! clean mpp
      CALL mpp_clean(mpp__init_file)

      ! check file type
      SELECT CASE( TRIM(td_file%c_type) )
         CASE('cdf')
            ! need to read all file to get domain decomposition
            tl_file=file_copy(td_file)

            ! open file
            CALL iom_open(tl_file)

            ! read first file domain decomposition
            tl_mpp=mpp__init_file_cdf(tl_file)

            ! get number of processor/file to be read
            il_nproc = 1
            il_attid = 0

            IF( ASSOCIATED(tl_file%t_att) )THEN
               il_attid=att_get_id( tl_file%t_att, "DOMAIN_number_total" )
            ENDIF
            IF( il_attid /= 0 )THEN
               il_nproc = INT(tl_file%t_att(il_attid)%d_value(1))
            ENDIF

            ! close file
            CALL iom_close(tl_file)

            IF( il_nproc /= 1 )THEN
               DO ji=1,il_nproc

                  ! clean mpp strcuture
                  CALL mpp_clean(tl_mpp)
 
                  ! get filename 
                  tl_file=file_rename(td_file,ji)
 
                  ! open file
                  CALL iom_open(tl_file)

                  ! read domain decomposition
                  tl_mpp = mpp__init_file_cdf(tl_file)
                  IF( ji == 1 )THEN
                     mpp__init_file=mpp_copy(tl_mpp)
                  ELSE
                     IF( ANY( mpp__init_file%t_dim(1:2)%i_len /= &
                                      tl_mpp%t_dim(1:2)%i_len) )THEN

                        CALL logger_error("MPP INIT READ: dimension from file "//&
                        &     TRIM(tl_file%c_name)//" and mpp strcuture "//&
                        &     TRIM(mpp__init_file%c_name)//"differ ")

                     ELSE

                        ! add processor to mpp strcuture
                        CALL mpp__add_proc(mpp__init_file, tl_mpp%t_proc(1))

                     ENDIF
                  ENDIF

                  ! close file
                  CALL iom_close(tl_file)

               ENDDO
               IF( mpp__init_file%i_nproc /= il_nproc )THEN
                  CALL logger_error("MPP INIT READ: some processors can't be added &
                  &               to mpp structure")
               ENDIF

            ELSE
               mpp__init_file=mpp_copy(tl_mpp)
            ENDIF

            ! mpp type
            mpp__init_file%c_type=TRIM(td_file%c_type)

            ! mpp domain type
            CALL mpp_get_dom(mpp__init_file)

            ! create some attributes for domain decomposition (use with dimg file)
            tl_att=att_init( "DOMAIN_number_total", mpp__init_file%i_nproc )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_I_position_first", mpp__init_file%t_proc(:)%i_impp )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_J_position_first", mpp__init_file%t_proc(:)%i_jmpp )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_I_position_last", mpp__init_file%t_proc(:)%i_lci )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_J_position_last", mpp__init_file%t_proc(:)%i_lcj )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_I_halo_size_start", mpp__init_file%t_proc(:)%i_ldi )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_J_halo_size_start", mpp__init_file%t_proc(:)%i_ldj )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_I_halo_size_end", mpp__init_file%t_proc(:)%i_lei )
            CALL mpp_move_att(mpp__init_file, tl_att)

            tl_att=att_init( "DOMAIN_J_halo_size_end", mpp__init_file%t_proc(:)%i_lej )
            CALL mpp_move_att(mpp__init_file, tl_att)
            
            ! clean
            CALL mpp_clean(tl_mpp)
            CALL att_clean(tl_att)

         CASE('dimg')
            ! domain decomposition could be read in one file

            tl_file=file_copy(td_file)
            ! open file
            CALL logger_debug("MPP INIT READ: open file "//TRIM(tl_file%c_name))
            CALL iom_open(tl_file)

            CALL logger_debug("MPP INIT READ: read mpp structure ")
            ! read mpp structure
            mpp__init_file=mpp__init_file_rstdimg(tl_file)

            ! mpp type
            mpp__init_file%c_type=TRIM(td_file%c_type)

            ! mpp domain type
            CALL logger_debug("MPP INIT READ: mpp_get_dom ")
            CALL mpp_get_dom(mpp__init_file)

            ! get processor size
            CALL logger_debug("MPP INIT READ: get processor size ")
            DO ji=1,mpp__init_file%i_nproc

               il_shape(:)=mpp_get_proc_size( mpp__init_file, ji )

               tl_dim=dim_init('X',il_shape(1))
               CALL file_add_dim(mpp__init_file%t_proc(ji), tl_dim)

               tl_dim=dim_init('Y',il_shape(2))
               CALL file_add_dim(mpp__init_file%t_proc(ji), tl_dim)            

               ! clean
               CALL dim_clean(tl_dim)

            ENDDO

            ! close file
            CALL iom_close(tl_file)

         CASE DEFAULT
            CALL logger_error("MPP INIT READ: invalid type for file "//&
            &              TRIM(tl_file%c_name))
      END SELECT

      ! east west overlap
      IF( PRESENT(id_ew) ) mpp__init_file%i_ew=id_ew
      ! NEMO periodicity
      IF( PRESENT(id_perio) )THEN
         mpp__init_file%i_perio= id_perio
         SELECT CASE(id_perio)
         CASE(3,4)
            mpp__init_file%i_pivot=1
         CASE(5,6)
            mpp__init_file%i_pivot=0
         CASE DEFAULT
            mpp__init_file%i_pivot=1
         END SELECT
      ENDIF

      IF( PRESENT(id_pivot) ) mpp__init_file%i_pivot= id_pivot

      ! clean 
      CALL file_clean(tl_file)

   END FUNCTION mpp__init_file
   !-------------------------------------------------------------------
   !> @brief This function initalise a mpp structure, 
   !> reading some netcdf files.
   !
   !> @details 
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !> @date July, 2015 - add only use dimension in MPP structure
   !>
   !> @param[in] td_file   file strcuture
   !> @return mpp structure
   !-------------------------------------------------------------------
   TYPE(TMPP) FUNCTION mpp__init_file_cdf( td_file )
      IMPLICIT NONE

      ! Argument
      TYPE(TFILE), INTENT(IN) :: td_file

      ! local variable
      INTEGER(i4) :: il_attid  ! attribute id
      
      LOGICAL     :: ll_exist
      LOGICAL     :: ll_open

      TYPE(TATT)  :: tl_att

      TYPE(TDIM)  :: tl_dim
      
      TYPE(TFILE) :: tl_proc
      !----------------------------------------------------------------

      CALL logger_trace("MPP INIT READ: netcdf file "//TRIM(td_file%c_name))

      INQUIRE( FILE=TRIM(td_file%c_name), EXIST=ll_exist, OPENED=ll_open )
      ! ll_open do not work for netcdf file, return always FALSE
      IF( ll_exist )THEN

         IF( td_file%i_id == 0 )THEN
            CALL logger_info(" id "//TRIM(fct_str(td_file%i_id))) 
            CALL logger_error("MPP INIT READ: netcdf file "//&
               &  TRIM(td_file%c_name)//" not opened")
         ELSE

            ! get mpp name
            mpp__init_file_cdf%c_name=TRIM( file_rename(td_file%c_name) )

            ! add type
            mpp__init_file_cdf%c_type="cdf"

            ! global domain size
            il_attid = 0
            IF( ASSOCIATED(td_file%t_att) )THEN
               il_attid=att_get_id( td_file%t_att, "DOMAIN_size_global" )
            ENDIF
            IF( il_attid /= 0 )THEN
               tl_dim=dim_init('X',INT(td_file%t_att(il_attid)%d_value(1)))
               CALL mpp_add_dim(mpp__init_file_cdf,tl_dim)

               tl_dim=dim_init('Y',INT(td_file%t_att(il_attid)%d_value(2)))
               CALL mpp_add_dim(mpp__init_file_cdf,tl_dim)
            ELSE ! assume only one file (not mpp)
               tl_dim=dim_init( td_file%t_dim(1)%c_name, td_file%t_dim(1)%i_len)
               CALL mpp_add_dim(mpp__init_file_cdf,tl_dim)

               tl_dim=dim_init( td_file%t_dim(2)%c_name, td_file%t_dim(2)%i_len)
               CALL mpp_add_dim(mpp__init_file_cdf,tl_dim)
            ENDIF

            IF( td_file%t_dim(3)%l_use )THEN
               tl_dim=dim_init( td_file%t_dim(3)%c_name, td_file%t_dim(3)%i_len)
               CALL mpp_add_dim(mpp__init_file_cdf,tl_dim)
            ENDIF

            IF( td_file%t_dim(4)%l_use )THEN
               tl_dim=dim_init( td_file%t_dim(4)%c_name, td_file%t_dim(4)%i_len)
               CALL mpp_add_dim(mpp__init_file_cdf,tl_dim)
            ENDIF

            ! initialise file/processor
            tl_proc=file_copy(td_file)

            ! processor id
            il_attid = 0
            IF( ASSOCIATED(td_file%t_att) )THEN
               il_attid=att_get_id( td_file%t_att, "DOMAIN_number" )
            ENDIF
            IF( il_attid /= 0 )THEN
               tl_proc%i_pid = INT(td_file%t_att(il_attid)%d_value(1))
            ELSE
               tl_proc%i_pid = 1
            ENDIF

            ! processor dimension
            tl_proc%t_dim(:)=dim_copy(td_file%t_dim(:))

            ! DOMAIN_position_first
            il_attid = 0
            IF( ASSOCIATED(td_file%t_att) )THEN
               il_attid=att_get_id( td_file%t_att, "DOMAIN_position_first" )
            ENDIF
            IF( il_attid /= 0 )THEN
               tl_proc%i_impp = INT(td_file%t_att(il_attid)%d_value(1))
               tl_proc%i_jmpp = INT(td_file%t_att(il_attid)%d_value(2))
            ELSE
               tl_proc%i_impp = 1
               tl_proc%i_jmpp = 1
            ENDIF

            ! DOMAIN_position_last
            il_attid = 0
            IF( ASSOCIATED(td_file%t_att) )THEN
               il_attid=att_get_id( td_file%t_att, "DOMAIN_position_last" )
            ENDIF
            IF( il_attid /= 0 )THEN
               tl_proc%i_lci = INT(td_file%t_att(il_attid)%d_value(1)) + tl_proc%i_impp
               tl_proc%i_lcj = INT(td_file%t_att(il_attid)%d_value(2)) + tl_proc%i_jmpp
            ELSE
               tl_proc%i_lci = mpp__init_file_cdf%t_dim(1)%i_len
               tl_proc%i_lcj = mpp__init_file_cdf%t_dim(2)%i_len
            ENDIF

            ! DOMAIN_halo_size_start
            il_attid = 0
            IF( ASSOCIATED(td_file%t_att) )THEN
               il_attid=att_get_id( td_file%t_att, "DOMAIN_halo_size_start" )
            ENDIF
            IF( il_attid /= 0 )THEN
               tl_proc%i_ldi = INT(td_file%t_att(il_attid)%d_value(1))
               tl_proc%i_ldj = INT(td_file%t_att(il_attid)%d_value(2))
            ELSE
               tl_proc%i_ldi = 1
               tl_proc%i_ldj = 1
            ENDIF

            ! DOMAIN_halo_size_end
            il_attid = 0
            IF( ASSOCIATED(td_file%t_att) )THEN
               il_attid=att_get_id( td_file%t_att, "DOMAIN_halo_size_end" )
            ENDIF
            IF( il_attid /= 0 )THEN
               tl_proc%i_lei = INT(td_file%t_att(il_attid)%d_value(1))
               tl_proc%i_lej = INT(td_file%t_att(il_attid)%d_value(2))
            ELSE
               tl_proc%i_lei = mpp__init_file_cdf%t_dim(1)%i_len
               tl_proc%i_lej = mpp__init_file_cdf%t_dim(2)%i_len
            ENDIF

            ! add attributes
            tl_att=att_init( "DOMAIN_size_global", &
            &                mpp__init_file_cdf%t_dim(:)%i_len)
            CALL file_move_att(tl_proc, tl_att)

            tl_att=att_init( "DOMAIN_number", tl_proc%i_pid )
            CALL file_move_att(tl_proc, tl_att)

            tl_att=att_init( "DOMAIN_position_first", &
            &                (/tl_proc%i_impp, tl_proc%i_jmpp /) )
            CALL file_move_att(tl_proc, tl_att)

            tl_att=att_init( "DOMAIN_position_last", &
            &                (/tl_proc%i_lci, tl_proc%i_lcj /) )
            CALL file_move_att(tl_proc, tl_att)

            tl_att=att_init( "DOMAIN_halo_size_start", &
            &                (/tl_proc%i_ldi, tl_proc%i_ldj /) )
            CALL file_move_att(tl_proc, tl_att)

            tl_att=att_init( "DOMAIN_halo_size_end", &
            &                (/tl_proc%i_lei, tl_proc%i_lej /) )
            CALL file_move_att(tl_proc, tl_att)

            ! add processor to mpp structure
            CALL mpp__add_proc(mpp__init_file_cdf, tl_proc)

            ! clean 
            CALL file_clean(tl_proc)
            CALL att_clean(tl_att)
         ENDIF

      ELSE

         CALL logger_error("MPP INIT READ: netcdf file "//TRIM(td_file%c_name)//&
         &  " do not exist")

      ENDIF      
   END FUNCTION mpp__init_file_cdf
   !-------------------------------------------------------------------
   !> @brief This function initalise a mpp structure, 
   !> reading one dimg restart file.
   !
   !> @details 
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !
   !> @param[in] td_file   file strcuture
   !> @return mpp structure
   !-------------------------------------------------------------------
   TYPE(TMPP) FUNCTION mpp__init_file_rstdimg( td_file )
      IMPLICIT NONE

      ! Argument
      TYPE(TFILE), INTENT(IN) :: td_file

      ! local variable
      INTEGER(i4)       :: il_status
      INTEGER(i4)       :: il_recl                          ! record length 
      INTEGER(i4)       :: il_nx, il_ny, il_nz              ! x,y,z dimension 
      INTEGER(i4)       :: il_n0d, il_n1d, il_n2d, il_n3d   ! number of 0/1/2/3D variables 
      INTEGER(i4)       :: il_iglo, il_jglo                 ! domain global size
      INTEGER(i4)       :: il_rhd                           ! record of the header infos
      INTEGER(i4)       :: il_pni, il_pnj, il_pnij          ! domain decomposition
      INTEGER(i4)       :: il_area                          ! domain index

      LOGICAL           ::  ll_exist
      LOGICAL           ::  ll_open

      CHARACTER(LEN=lc) :: cl_file

      TYPE(TDIM)        :: tl_dim ! dimension structure
      TYPE(TATT)        :: tl_att
      TYPE(TFILE)       :: tl_proc

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      INQUIRE( FILE=TRIM(td_file%c_name), EXIST=ll_exist, OPENED=ll_open)
      IF( ll_exist )THEN

         IF( .NOT. ll_open )THEN
            CALL logger_error("MPP INIT READ: dimg file "//TRIM(td_file%c_name)//&
            &  " not opened")
         ELSE

            ! read first record 
            READ( td_file%i_id, IOSTAT=il_status, REC=1 )& 
            &     il_recl,                         &
            &     il_nx, il_ny, il_nz,             &
            &     il_n0d, il_n1d, il_n2d, il_n3d,  &
            &     il_rhd,                          &
            &     il_pni, il_pnj, il_pnij,         &
            &     il_area
            CALL fct_err(il_status)
            IF( il_status /= 0 )THEN
               CALL logger_error("MPP INIT READ: read first line header of "//&
               &              TRIM(td_file%c_name))
            ENDIF

            ! get mpp name
            mpp__init_file_rstdimg%c_name=TRIM( file_rename(td_file%c_name) )

            ! add type
            mpp__init_file_rstdimg%c_type="dimg"

            ! number of processors to be read
            mpp__init_file_rstdimg%i_nproc  = il_pnij
            mpp__init_file_rstdimg%i_niproc = il_pni
            mpp__init_file_rstdimg%i_njproc = il_pnj

            IF( ASSOCIATED(mpp__init_file_rstdimg%t_proc) )THEN
               CALL file_clean(mpp__init_file_rstdimg%t_proc(:))
               DEALLOCATE(mpp__init_file_rstdimg%t_proc)
            ENDIF
            ALLOCATE( mpp__init_file_rstdimg%t_proc(il_pnij) , stat=il_status )

            tl_proc=file_copy(td_file)
            ! remove dimension from file
            CALL dim_clean(tl_proc%t_dim(:))
            ! initialise file/processors
            DO ji=1,mpp__init_file_rstdimg%i_nproc
               mpp__init_file_rstdimg%t_proc(ji)=file_copy(tl_proc)
            ENDDO

            IF( il_status /= 0 )THEN
               CALL logger_error("MPP INIT READ: not enough space to read domain &
               &              decomposition in file "//TRIM(td_file%c_name))
            ENDIF

            ! read first record 
            READ( td_file%i_id, IOSTAT=il_status, REC=1 )& 
            &     il_recl,                         &
            &     il_nx, il_ny, il_nz,             &
            &     il_n0d, il_n1d, il_n2d, il_n3d,  &
            &     il_rhd,                          &
            &     il_pni, il_pnj, il_pnij,         &
            &     il_area,                         &
            &     il_iglo, il_jglo,                &
            &     mpp__init_file_rstdimg%t_proc(:)%i_lci,    &
            &     mpp__init_file_rstdimg%t_proc(:)%i_lcj,    &
            &     mpp__init_file_rstdimg%t_proc(:)%i_ldi,    &
            &     mpp__init_file_rstdimg%t_proc(:)%i_ldj,    &
            &     mpp__init_file_rstdimg%t_proc(:)%i_lei,    &
            &     mpp__init_file_rstdimg%t_proc(:)%i_lej,    &
            &     mpp__init_file_rstdimg%t_proc(:)%i_impp,   &
            &     mpp__init_file_rstdimg%t_proc(:)%i_jmpp
            CALL fct_err(il_status)
            IF( il_status /= 0 )THEN
               CALL logger_error("MPP INIT READ: read first line of "//&
               &              TRIM(td_file%c_name))
            ENDIF

            ! global domain size
            tl_dim=dim_init('X',il_iglo)
            CALL mpp_add_dim(mpp__init_file_rstdimg,tl_dim)
            tl_dim=dim_init('Y',il_jglo)
            CALL mpp_add_dim(mpp__init_file_rstdimg,tl_dim)

            tl_dim=dim_init('Z',il_nz)
            CALL mpp_add_dim(mpp__init_file_rstdimg,tl_dim)

            DO ji=1,mpp__init_file_rstdimg%i_nproc
               ! get file name
               cl_file =  file_rename(td_file%c_name,ji)
               mpp__init_file_rstdimg%t_proc(ji)%c_name = TRIM(cl_file)
               ! update processor id
               mpp__init_file_rstdimg%t_proc(ji)%i_pid=ji

               ! add attributes
               tl_att=att_init( "DOMAIN_number", ji )
               CALL file_move_att(mpp__init_file_rstdimg%t_proc(ji), tl_att) 

               tl_att=att_init( "DOMAIN_position_first", &
               &                (/mpp__init_file_rstdimg%t_proc(ji)%i_impp, &
               &                  mpp__init_file_rstdimg%t_proc(ji)%i_jmpp /) )
               CALL file_move_att(mpp__init_file_rstdimg%t_proc(ji), tl_att)

               tl_att=att_init( "DOMAIN_position_last", &
               &                (/mpp__init_file_rstdimg%t_proc(ji)%i_lci, &
               &                  mpp__init_file_rstdimg%t_proc(ji)%i_lcj /) )
               CALL file_move_att(mpp__init_file_rstdimg%t_proc(ji), tl_att)

               tl_att=att_init( "DOMAIN_halo_size_start", &
               &                (/mpp__init_file_rstdimg%t_proc(ji)%i_ldi, &
               &                  mpp__init_file_rstdimg%t_proc(ji)%i_ldj /) )
               CALL file_move_att(mpp__init_file_rstdimg%t_proc(ji), tl_att)               

               tl_att=att_init( "DOMAIN_halo_size_end", &
               &                (/mpp__init_file_rstdimg%t_proc(ji)%i_lei, &
               &                  mpp__init_file_rstdimg%t_proc(ji)%i_lej /) )
               CALL file_move_att(mpp__init_file_rstdimg%t_proc(ji), tl_att)
            ENDDO
 
            ! add type
            mpp__init_file_rstdimg%t_proc(:)%c_type="dimg"

            ! add attributes
            tl_att=att_init( "DOMAIN_size_global", &
            &                mpp__init_file_rstdimg%t_dim(:)%i_len)
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_number_total", &
            &                 mpp__init_file_rstdimg%i_nproc )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_I_number_total", &
            &                 mpp__init_file_rstdimg%i_niproc )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_J_number_total", &
            &                 mpp__init_file_rstdimg%i_njproc )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_I_position_first", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_impp )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_J_position_first", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_jmpp )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_I_position_last", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_lci )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_J_position_last", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_lcj )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_I_halo_size_start", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_ldi )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_J_halo_size_start", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_ldj )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_I_halo_size_end", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_lei )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            tl_att=att_init( "DOMAIN_J_halo_size_end", &
            &                 mpp__init_file_rstdimg%t_proc(:)%i_lej )
            CALL mpp_move_att(mpp__init_file_rstdimg, tl_att)

            ! clean
            CALL dim_clean(tl_dim)
            CALL att_clean(tl_att)
         ENDIF

      ELSE

         CALL logger_error("MPP INIT READ: dimg file "//TRIM(td_file%c_name)//&
         &  " do not exist")

      ENDIF

   END FUNCTION mpp__init_file_rstdimg
   !-------------------------------------------------------------------
   !> @brief This function check if variable and mpp structure use same
   !> dimension.
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !
   !> @param[in] td_mpp    mpp structure
   !> @param[in] td_proc   processor structure
   !> @return dimension of processor and mpp structure agree (or not)
   !-------------------------------------------------------------------
   LOGICAL FUNCTION mpp__check_proc_dim(td_mpp, td_proc)
      IMPLICIT NONE
      ! Argument      
      TYPE(TMPP),  INTENT(IN) :: td_mpp
      TYPE(TFILE), INTENT(IN) :: td_proc

      ! local variable
      INTEGER(i4) :: il_isize !< i-direction maximum sub domain size 
      INTEGER(i4) :: il_jsize !< j-direction maximum sub domain size

      !----------------------------------------------------------------
      mpp__check_proc_dim=.TRUE.
      ! check used dimension 
      IF( td_mpp%i_niproc /= 0 .AND. td_mpp%i_njproc /= 0 )THEN
         ! check with maximum size of sub domain
         il_isize = ( td_mpp%t_dim(1)%i_len - 2*td_mpp%i_preci + &
         &           (td_mpp%i_niproc-1) ) / td_mpp%i_niproc + 2*td_mpp%i_preci
         il_jsize = ( td_mpp%t_dim(2)%i_len - 2*td_mpp%i_precj + &
         &           (td_mpp%i_njproc-1) ) / td_mpp%i_njproc + 2*td_mpp%i_precj

         IF( il_isize < td_proc%i_lci .OR.                     &
         &   il_jsize < td_proc%i_lcj )THEN

            mpp__check_proc_dim=.FALSE.

            CALL logger_error( "MPP CHECK DIM: processor and mpp dimension differ" )

         ENDIF

      ELSE
         ! check with global domain size
         IF( td_mpp%t_dim(1)%i_len < td_proc%i_lci .OR.                     &
         &   td_mpp%t_dim(2)%i_len < td_proc%i_lcj )THEN

            mpp__check_proc_dim=.FALSE.

            CALL logger_error( "MPP CHECK DIM: processor and mpp dimension differ" )

         ENDIF
      ENDIF

   END FUNCTION mpp__check_proc_dim
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine add variable in all files of mpp structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in]    td_var variable strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp_add_var( td_mpp, td_var )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TVAR), INTENT(IN)    :: td_var

      ! local variable
      INTEGER(i4) :: il_varid
      TYPE(TVAR)  :: tl_var

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      ! check if mpp exist
      IF( .NOT. ASSOCIATED(td_mpp%t_proc) )THEN

         CALL logger_error( "MPP ADD VAR: processor decomposition not "//&
         &  "define for mpp "//TRIM(td_mpp%c_name))

      ELSE
         ! check if variable exist
         IF( TRIM(td_var%c_name) == '' .AND. &
         &   TRIM(td_var%c_stdname) == '' )THEN
            CALL logger_error("MPP ADD VAR: variable not define ")
         ELSE
            ! check if variable already in mpp structure
            il_varid=0
            IF( ASSOCIATED(td_mpp%t_proc(1)%t_var) )THEN
               il_varid=var_get_index( td_mpp%t_proc(1)%t_var(:), &
               &                       td_var%c_name, td_var%c_stdname )
            ENDIF

            IF( il_varid /= 0 )THEN

               DO ji=1,td_mpp%t_proc(1)%i_nvar
                  CALL logger_debug( " MPP ADD VAR: in mpp structure : &
                  &  variable "//TRIM(td_mpp%t_proc(1)%t_var(ji)%c_name)//&
                  &  ", standard name "//&
                  &  TRIM(td_mpp%t_proc(1)%t_var(ji)%c_stdname) )
               ENDDO
               CALL logger_error( " MPP ADD VAR: variable "//TRIM(td_var%c_name)//&
               &  ", standard name "//TRIM(td_var%c_stdname)//&
               &  ", already in mpp "//TRIM(td_mpp%c_name) )

            ELSE
            
               CALL logger_info( &
               &  " MPP ADD VAR: add variable "//TRIM(td_var%c_name)//&
               &  ", standard name "//TRIM(td_var%c_stdname)//&
               &  ", in mpp "//TRIM(td_mpp%c_name) )
               ! check used dimension 
               IF( mpp__check_dim(td_mpp, td_var) )THEN

                  ! update dimension if need be
                  DO ji=1,ip_maxdim
                     IF( td_var%t_dim(ji)%l_use .AND. &
                     &   .NOT. td_mpp%t_dim(ji)%l_use )THEN
                        CALL mpp_add_dim(td_mpp,td_var%t_dim(ji))
                     ENDIF
                  ENDDO

                  ! add variable in each processor
                  DO ji=1,td_mpp%i_nproc

                     ! split variable on domain decomposition
                     tl_var=mpp__split_var(td_mpp, td_var, ji)

                     CALL file_add_var(td_mpp%t_proc(ji), tl_var)

                     ! clean
                     CALL var_clean(tl_var)
                  ENDDO

               ENDIF
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE mpp_add_var
   !-------------------------------------------------------------------
   !> @brief This function extract, from variable structure, part that will 
   !> be written in processor id_procid.<br/>
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !
   !> @param[in] td_mpp    mpp structure
   !> @param[in] td_var    variable structure
   !> @param[in] id_procid processor id
   !> @return variable structure
   !-------------------------------------------------------------------
   TYPE(TVAR) FUNCTION mpp__split_var(td_mpp, td_var, id_procid)
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),  INTENT(IN) :: td_mpp
      TYPE(TVAR),  INTENT(IN) :: td_var
      INTEGER(i4), INTENT(IN) :: id_procid

      ! local variable
      TYPE(TDIM)  :: tl_dim

      INTEGER(i4), DIMENSION(4) :: il_ind
      INTEGER(i4), DIMENSION(2) :: il_size
      INTEGER(i4) :: il_i1
      INTEGER(i4) :: il_i2
      INTEGER(i4) :: il_j1
      INTEGER(i4) :: il_j2
      !----------------------------------------------------------------

      ! copy mpp
      mpp__split_var=var_copy(td_var)

      IF( ASSOCIATED(td_var%d_value) )THEN
         ! remove value over global domain from pointer
         CALL var_del_value( mpp__split_var )

         ! get processor dimension
         il_size(:)=mpp_get_proc_size( td_mpp, id_procid )

         ! define new dimension in variable structure 
         IF( td_var%t_dim(1)%l_use )THEN
            tl_dim=dim_init( TRIM(td_var%t_dim(1)%c_name), il_size(1) )
            CALL var_move_dim( mpp__split_var, tl_dim )
         ENDIF
         IF( td_var%t_dim(2)%l_use )THEN
            tl_dim=dim_init( TRIM(td_var%t_dim(2)%c_name), il_size(2) )
            CALL var_move_dim( mpp__split_var, tl_dim )      
         ENDIF

         ! get processor indices
         il_ind(:)=mpp_get_proc_index( td_mpp, id_procid )
         il_i1 = il_ind(1)
         il_i2 = il_ind(2)
         il_j1 = il_ind(3)
         il_j2 = il_ind(4)

         IF( .NOT. td_var%t_dim(1)%l_use )THEN
            il_i1=1 
            il_i2=1 
         ENDIF

         IF( .NOT. td_var%t_dim(2)%l_use )THEN
            il_j1=1 
            il_j2=1 
         ENDIF      

         ! add variable value on processor
         CALL var_add_value( mpp__split_var, &
         &                   td_var%d_value(il_i1:il_i2, il_j1:il_j2, :, :) )
      ENDIF

   END FUNCTION mpp__split_var
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine delete all variable in mpp strcuture.
   !>
   !> @author J.Paul
   !> @date October, 2014 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_var_mpp( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp

      ! local variable
      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      CALL logger_info( &
      &  "MPP CLEAN VAR: reset all variable "//&
      &  "in mpp strcuture "//TRIM(td_mpp%c_name) )

      IF( ASSOCIATED(td_mpp%t_proc) )THEN
         DO ji=td_mpp%t_proc(1)%i_nvar,1,-1
            CALL mpp_del_var(td_mpp, td_mpp%t_proc(1)%t_var(ji))
         ENDDO
      ENDIF

   END SUBROUTINE mpp__del_var_mpp
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine delete variable in mpp structure, given variable
   !> structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in]    td_var variable strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_var_str( td_mpp, td_var )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TVAR), INTENT(IN)    :: td_var

      ! local variable
      INTEGER(i4)       :: il_varid
      CHARACTER(LEN=lc) :: cl_name

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      ! check if mpp exist
      IF( .NOT. ASSOCIATED(td_mpp%t_proc) )THEN

         CALL logger_error( "MPP DEL VAR: domain decomposition not define "//&
         &               " in mpp strcuture "//TRIM(td_mpp%c_name))

      ELSE

         ! check if variable already in mpp structure
         il_varid = 0
         IF( ASSOCIATED(td_mpp%t_proc(1)%t_var) )THEN
            il_varid=var_get_index( td_mpp%t_proc(1)%t_var(:), &
            &                       td_var%c_name, td_var%c_stdname )
         ENDIF
         IF( il_varid == 0 )THEN
            CALL logger_error( &
            &  "MPP DEL VAR: no variable "//TRIM(td_var%c_name)//&
            &  ", in mpp structure "//TRIM(td_mpp%c_name) )

            DO ji=1,td_mpp%t_proc(1)%i_nvar
               CALL logger_debug( "MPP DEL VAR: in mpp structure : &
               &  variable : "//TRIM(td_mpp%t_proc(1)%t_var(ji)%c_name)//&
               &  ", standard name "//&
               &  TRIM(td_mpp%t_proc(1)%t_var(ji)%c_stdname) )
            ENDDO

         ELSE

            cl_name=TRIM(td_var%c_name)
            DO ji=1,td_mpp%i_nproc
               CALL file_del_var(td_mpp%t_proc(ji), TRIM(cl_name)) 
            ENDDO

         ENDIF

      ENDIF
   END SUBROUTINE mpp__del_var_str
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine delete variable in mpp structure, given variable name.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !> @date February, 2015 
   !> - define local variable structure to avoid mistake with pointer
   !
   !> @param[inout] td_mpp    mpp strcuture
   !> @param[in]    cd_name   variable name
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_var_name( td_mpp, cd_name )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP)      , INTENT(INOUT) :: td_mpp
      CHARACTER(LEN=*), INTENT(IN   ) :: cd_name

      ! local variable
      INTEGER(i4)       :: il_varid
      TYPE(TVAR)        :: tl_var
      !----------------------------------------------------------------
      ! check if mpp exist
      IF( .NOT. ASSOCIATED(td_mpp%t_proc) )THEN

         CALL logger_error( "MPP DEL VAR: domain decomposition not define "//&
         &               " in mpp strcuture "//TRIM(td_mpp%c_name))

      ELSE

         IF( td_mpp%t_proc(1)%i_nvar == 0 )THEN
            CALL logger_debug( "MPP DEL VAR NAME: no variable associated to mpp &
            &                 structure "//TRIM(td_mpp%c_name) )
         ELSE

            ! get the variable id, in file variable structure
            il_varid=0
            IF( ASSOCIATED(td_mpp%t_proc(1)%t_var) )THEN
               il_varid=var_get_index( td_mpp%t_proc(1)%t_var(:), &
               &                       cd_name )
            ENDIF

            IF( il_varid == 0 )THEN

               CALL logger_warn( &
               &  "MPP DEL VAR : there is no variable with name "//&
               &  "or standard name "//TRIM(ADJUSTL(cd_name))//&
               &  " in mpp structure "//TRIM(td_mpp%c_name))

            ELSE

               tl_var=var_copy(td_mpp%t_proc(1)%t_var(il_varid))
               CALL mpp_del_var(td_mpp, tl_var)

            ENDIF
         ENDIF

      ENDIF
   END SUBROUTINE mpp__del_var_name
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine overwrite variable in mpp structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in]    td_var variable structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp_move_var( td_mpp, td_var )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TVAR), INTENT(IN)    :: td_var

      !local variable
      TYPE(TVAR) :: tl_var
      !----------------------------------------------------------------
      ! copy variable
      tl_var=var_copy(td_var)

      ! remove processor
      CALL mpp_del_var(td_mpp, tl_var)

      ! add processor
      CALL mpp_add_var(td_mpp, tl_var)

      ! clean 
      CALL var_clean(tl_var)

   END SUBROUTINE mpp_move_var
   !> @endcode
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine add processor to mpp structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp    mpp strcuture
   !> @param[in]    td_proc   processor strcuture
   !
   !> @todo 
   !> - check proc type
   !-------------------------------------------------------------------
   SUBROUTINE mpp__add_proc( td_mpp, td_proc )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP) , INTENT(INOUT) :: td_mpp
      TYPE(TFILE), INTENT(IN)    :: td_proc

      ! local variable
      INTEGER(i4)                                  :: il_status
      INTEGER(i4)                                  :: il_procid
      INTEGER(i4)      , DIMENSION(1)              :: il_ind

      TYPE(TFILE)      , DIMENSION(:), ALLOCATABLE :: tl_proc

      CHARACTER(LEN=lc)                            :: cl_name
      !----------------------------------------------------------------

      ! check file name
      cl_name=TRIM( file_rename(td_proc%c_name) )
      IF( TRIM(cl_name) /=  TRIM(td_mpp%c_name) )THEN
         CALL logger_warn("MPP ADD PROC: processor name do not match mpp name")
      ENDIF

      il_procid=0
      IF( ASSOCIATED(td_mpp%t_proc) )THEN
         ! check if processor already in mpp structure
         il_ind(:)=MINLOC( td_mpp%t_proc(:)%i_pid, &
                     mask=(td_mpp%t_proc(:)%i_pid==td_proc%i_pid) )
         il_procid=il_ind(1)
      ENDIF

      IF( il_procid /= 0 )THEN

            CALL logger_error( &
            &  "MPP ADD PROC: processor "//TRIM(fct_str(td_proc%i_pid))//&
            &  ", already in mpp structure " )

      ELSE
 
         CALL logger_trace("MPP ADD PROC: add processor "//&
         &               TRIM(fct_str(td_mpp%i_nproc+1))//" in mpp structure")

         IF( td_mpp%i_nproc > 0 )THEN
            ! 
            il_ind(:)=MAXLOC( td_mpp%t_proc(:)%i_pid, &
                        mask=(td_mpp%t_proc(:)%i_pid < td_proc%i_pid) )
            il_procid=il_ind(1)

            ! already other processor in mpp structure
            ALLOCATE( tl_proc(td_mpp%i_nproc), stat=il_status )
            IF(il_status /= 0 )THEN

               CALL logger_error( "MPP ADD PROC: not enough space to put processor &
               &               in mpp structure")

            ELSE
               ! save temporary mpp structure
               tl_proc(:)=file_copy(td_mpp%t_proc(:))

               CALL file_clean( td_mpp%t_proc(:) )
               DEALLOCATE(td_mpp%t_proc)
               ALLOCATE( td_mpp%t_proc(td_mpp%i_nproc+1), stat=il_status)
               IF(il_status /= 0 )THEN

                  CALL logger_error( "MPP ADD PROC: not enough space to put "//&
                  &  "processor in mpp structure ")

               ENDIF

               ! copy processor in mpp before
               ! processor with lower id than new processor
               td_mpp%t_proc( 1:il_procid ) = file_copy(tl_proc( 1:il_procid ))

               ! processor with greater id than new processor
               td_mpp%t_proc( il_procid+1 : td_mpp%i_nproc+1 ) = &
               &              file_copy(tl_proc( il_procid : td_mpp%i_nproc ))

               ! clean
               CALL file_clean(tl_proc(:))
               DEALLOCATE(tl_proc)
            ENDIF

         ELSE
            ! no processor in mpp structure
            IF( ASSOCIATED(td_mpp%t_proc) )THEN
               CALL file_clean(td_mpp%t_proc(:))
               DEALLOCATE(td_mpp%t_proc)
            ENDIF
            ALLOCATE( td_mpp%t_proc(td_mpp%i_nproc+1), stat=il_status )
            IF(il_status /= 0 )THEN

               CALL logger_error( "MPP ADD PROC: not enough space to put "//&
               &  "processor in mpp structure " )

            ENDIF
         ENDIF

         ! check dimension
         IF( ANY(td_mpp%t_dim(1:2)%i_len < td_proc%t_dim(1:2)%i_len) )THEN
            CALL logger_error( "MPP ADD PROC: mpp structure and new processor "//&
            &  " dimension differ. ")
            CALL logger_debug("MPP ADD PROC: mpp dimension ("//&
            &  TRIM(fct_str(td_mpp%t_dim(1)%i_len))//","//&
            &  TRIM(fct_str(td_mpp%t_dim(2)%i_len))//")" )
            CALL logger_debug("MPP ADD PROC: processor dimension ("//&
            &  TRIM(fct_str(td_proc%t_dim(1)%i_len))//","//&
            &  TRIM(fct_str(td_proc%t_dim(2)%i_len))//")" )
         ELSE
            td_mpp%i_nproc=td_mpp%i_nproc+1

            ! add new processor
            td_mpp%t_proc(td_mpp%i_nproc)=file_copy(td_proc)
         ENDIF

      ENDIF
   END SUBROUTINE mpp__add_proc
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine delete processor in mpp structure, given processor id.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp    mpp strcuture
   !> @param[in]    id_procid processor id
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_proc_id( td_mpp, id_procid )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),   INTENT(INOUT) :: td_mpp
      INTEGER(i4),  INTENT(IN)    :: id_procid

      ! local variable
      INTEGER(i4) :: il_status
      INTEGER(i4) :: il_procid
      INTEGER(i4), DIMENSION(1) :: il_ind
      TYPE(TFILE), DIMENSION(:), ALLOCATABLE :: tl_proc

      ! loop indices
      !----------------------------------------------------------------

      il_ind(:)=MINLOC(td_mpp%t_proc(:)%i_pid,td_mpp%t_proc(:)%i_pid==id_procid)
      il_procid=il_ind(1)
      IF( il_procid == 0 )THEN
         CALL logger_error("MPP DEL PROC: no processor "//&
         &                 TRIM(fct_str(id_procid))//&
         &                 " associated to mpp structure")
      ELSE
         CALL logger_trace("DEL PROC: remove processor "//&
         &                 TRIM(fct_str(id_procid)))

         IF( td_mpp%i_nproc > 1 )THEN
            ALLOCATE( tl_proc(td_mpp%i_nproc-1), stat=il_status )
            IF(il_status /= 0 )THEN
               CALL logger_error( "MPP DEL PROC: not enough space to put &
               &  processor in temporary mpp structure")

            ELSE

               ! save temporary processor's mpp structure
               IF( il_procid > 1 )THEN
                  tl_proc(1:il_procid-1)=file_copy(td_mpp%t_proc(1:il_procid-1))
               ENDIF

               IF( il_procid < td_mpp%i_nproc )THEN
                  tl_proc(il_procid:)=file_copy(td_mpp%t_proc(il_procid+1:))
               ENDIF

               ! new number of processor in mpp
               td_mpp%i_nproc=td_mpp%i_nproc-1

               CALL file_clean( td_mpp%t_proc(:) )
               DEALLOCATE(td_mpp%t_proc)
               ALLOCATE( td_mpp%t_proc(td_mpp%i_nproc), stat=il_status )
               IF(il_status /= 0 )THEN

                  CALL logger_error( "MPP DEL PROC: not enough space &
                  &  to put processors in mpp structure " )

               ELSE

                  ! copy processor in mpp before
                  td_mpp%t_proc(:)=file_copy(tl_proc(:))

                  ! update processor id
                  td_mpp%t_proc( il_procid : td_mpp%i_nproc )%i_pid = &
                  &     td_mpp%t_proc( il_procid : td_mpp%i_nproc )%i_pid - 1

               ENDIF
            ENDIF
            ! clean
            CALL file_clean( tl_proc(:) )
            DEALLOCATE(tl_proc)
         ELSE
            CALL file_clean( td_mpp%t_proc(:) )
            DEALLOCATE(td_mpp%t_proc)

            ! new number of processor in mpp
            td_mpp%i_nproc=td_mpp%i_nproc-1
         ENDIF
      ENDIF
   END SUBROUTINE mpp__del_proc_id
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine delete processor in mpp structure, given processor
   !>    structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp : mpp strcuture
   !> @param[in]    td_proc : file/processor structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_proc_str( td_mpp, td_proc )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),   INTENT(INOUT) :: td_mpp
      TYPE(TFILE),  INTENT(IN)    :: td_proc
      !----------------------------------------------------------------

      IF( td_proc%i_pid >= 0 )THEN
         CALL mpp__del_proc( td_mpp, td_proc%i_pid )
      ELSE
         CALL logger_error("MPP DEL PROC: processor not defined")
      ENDIF

   END SUBROUTINE mpp__del_proc_str
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine overwrite processor in mpp structure.
   !>
   !> @detail
   !
   !> @author J.Paul
   !> @date Nov, 2013 - Initial version
   !
   !> @param[inout] td_mpp    mpp strcuture
   !> @param[in]    id_procid processor id
   !-------------------------------------------------------------------
   SUBROUTINE mpp__move_proc( td_mpp, td_proc )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),  INTENT(INOUT) :: td_mpp
      TYPE(TFILE), INTENT(IN)    :: td_proc
      !----------------------------------------------------------------

      ! remove processor
      CALL mpp__del_proc(td_mpp, td_proc)

      ! add processor
      CALL mpp__add_proc(td_mpp, td_proc)

   END SUBROUTINE mpp__move_proc
   !-------------------------------------------------------------------
   !> @brief This subroutine add a dimension structure in a mpp 
   !> structure.
   !> Do not overwrite, if dimension already in mpp structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !> @date July, 2015 
   !> - rewrite the same as way var_add_dim
   !>
   !> @param[inout] td_mpp mpp structure
   !> @param[in] td_dim    dimension structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp_add_dim(td_mpp, td_dim)
      IMPLICIT NONE
      ! Argument      
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TDIM), INTENT(IN)    :: td_dim

      ! local variable
      INTEGER(i4) :: il_ind

      ! loop indices
      !----------------------------------------------------------------
      IF( td_mpp%i_ndim <= ip_maxdim )THEN

         ! check if dimension already used in mpp structure
         il_ind=SCAN(TRIM(cp_dimorder),TRIM(td_dim%c_sname))
         IF( il_ind == 0 )THEN
            CALL logger_warn( &
            &  " MPP ADD DIM: dimension "//TRIM(td_dim%c_name)//&
            &  ", short name "//TRIM(td_dim%c_sname)//&
            &  ", will not be added in mpp "//TRIM(td_mpp%c_name) )
         ELSEIF( td_mpp%t_dim(il_ind)%l_use )THEN
            CALL logger_error( &
            &  " MPP ADD DIM: dimension "//TRIM(td_dim%c_name)//&
            &  ", short name "//TRIM(td_dim%c_sname)//&
            &  ", already used in mpp "//TRIM(td_mpp%c_name) )
         ELSE

            ! back to disorder dimension array 
            CALL dim_disorder(td_mpp%t_dim(:))

            ! add new dimension
            td_mpp%t_dim(td_mpp%i_ndim+1)=dim_copy(td_dim)

            ! update number of attribute
            td_mpp%i_ndim=COUNT(td_mpp%t_dim(:)%l_use)

         ENDIF
         ! reorder dimension to ('x','y','z','t')
         CALL dim_reorder(td_mpp%t_dim(:))

      ELSE
         CALL logger_error( &
         &  "MPP ADD DIM: too much dimension in mpp "//&
         &  TRIM(td_mpp%c_name)//" ("//TRIM(fct_str(td_mpp%i_ndim))//")")
      ENDIF

   END SUBROUTINE mpp_add_dim
   !-------------------------------------------------------------------
   !> @brief This subroutine delete a dimension structure in a mpp 
   !> structure.<br/>
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !> @date July, 2015 
   !> - rewrite the same as way var_del_dim
   !>
   !> @param[inout] td_mpp mpp structure
   !> @param[in] td_dim    dimension structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp_del_dim(td_mpp, td_dim)
      IMPLICIT NONE
      ! Argument      
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TDIM), INTENT(IN)    :: td_dim

      ! local variable
      INTEGER(i4) :: il_ind
      TYPE(TDIM)  :: tl_dim

      ! loop indices
      !----------------------------------------------------------------


      IF( td_mpp%i_ndim <= ip_maxdim )THEN

         CALL logger_trace( &
         &  " MPP DEL DIM: delete dimension "//TRIM(td_dim%c_name)//&
         &  ", short name "//TRIM(td_dim%c_sname)//&
         &  ", in mpp "//TRIM(td_mpp%c_name) )
         
         ! check if dimension already in variable structure
         il_ind=SCAN(TRIM(cp_dimorder),TRIM(td_dim%c_sname))

         ! replace dimension by empty one
         td_mpp%t_dim(il_ind)=dim_copy(tl_dim)

         ! update number of dimension
         td_mpp%i_ndim=COUNT(td_mpp%t_dim(:)%l_use)

         ! reorder dimension to ('x','y','z','t')
         CALL dim_reorder(td_mpp%t_dim)

      ELSE
         CALL logger_error( &
         &  " MPP DEL DIM: too much dimension in mpp "//&
         &  TRIM(td_mpp%c_name)//" ("//TRIM(fct_str(td_mpp%i_ndim))//")")
      ENDIF

   END SUBROUTINE mpp_del_dim
   !-------------------------------------------------------------------
   !> @brief This subroutine move a dimension structure 
   !> in mpp structure.
   !> @warning dimension order may have changed
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !>
   !> @param[inout] td_mpp mpp structure
   !> @param[in] td_dim    dimension structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp_move_dim(td_mpp, td_dim)
      IMPLICIT NONE
      ! Argument      
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TDIM), INTENT(IN)    :: td_dim

      ! local variable
      INTEGER(i4) :: il_ind
      INTEGER(i4) :: il_dimid
      !----------------------------------------------------------------
      IF( td_mpp%i_ndim <= ip_maxdim )THEN

         ! check if dimension already in mpp structure
         il_ind=dim_get_index(td_mpp%t_dim(:), td_dim%c_name, td_dim%c_sname)
         IF( il_ind /= 0 )THEN

            il_dimid=td_mpp%t_dim(il_ind)%i_id
            ! replace dimension
            td_mpp%t_dim(il_ind)=dim_copy(td_dim)
            td_mpp%t_dim(il_ind)%i_id=il_dimid
            td_mpp%t_dim(il_ind)%l_use=.TRUE.

         ELSE
            CALL mpp_add_dim(td_mpp, td_dim)
         ENDIF

      ELSE
         CALL logger_error( &
         &  "MPP MOVE DIM: too much dimension in mpp "//&
         &  TRIM(td_mpp%c_name)//" ("//TRIM(fct_str(td_mpp%i_ndim))//")")
      ENDIF
   END SUBROUTINE mpp_move_dim
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine add global attribute to mpp structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in]    td_att attribute strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp_add_att( td_mpp, td_att )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TATT), INTENT(IN)    :: td_att

      ! local variable
      INTEGER(i4) :: il_attid

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      ! check if mpp exist
      IF( .NOT. ASSOCIATED(td_mpp%t_proc) )THEN

         CALL logger_error( "MPP ADD ATT: domain decomposition not define "//&
         &               "for mpp "//TRIM(td_mpp%c_name))

      ELSE
         ! check if variable exist
         IF( TRIM(td_att%c_name) == '' )THEN
            CALL logger_error("MPP ADD ATT: attribute not define ")
         ELSE
            ! check if attribute already in mpp structure
            il_attid=0
            IF( ASSOCIATED(td_mpp%t_proc(1)%t_att) )THEN
               il_attid=att_get_index( td_mpp%t_proc(1)%t_att(:), &
               &                    td_att%c_name )
            ENDIF
            IF( il_attid /= 0 )THEN

               CALL logger_error( " MPP ADD ATT: attribute "//&
               &                 TRIM(td_att%c_name)//&
               &                 ", already in mpp "//TRIM(td_mpp%c_name) )

               DO ji=1,td_mpp%t_proc(1)%i_natt
                  CALL logger_debug( " MPP ADD ATT: in mpp structure : &
                  &  attribute "//TRIM(td_mpp%t_proc(1)%t_att(ji)%c_name) )
               ENDDO

            ELSE
            
               CALL logger_info( &
               &  " MPP ADD ATT: add attribute "//TRIM(td_att%c_name)//&
               &  ", in mpp "//TRIM(td_mpp%c_name) )

               ! add attribute in each processor
               DO ji=1,td_mpp%i_nproc

                  CALL file_add_att(td_mpp%t_proc(ji), td_att)

               ENDDO

            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE mpp_add_att
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine delete attribute in mpp structure, given attribute
   !> structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in]    td_att attribute strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_att_str( td_mpp, td_att )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TATT), INTENT(IN)    :: td_att

      ! local variable
      INTEGER(i4)       :: il_attid
      CHARACTER(LEN=lc) :: cl_name

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      ! check if mpp exist
      IF( .NOT. ASSOCIATED(td_mpp%t_proc) )THEN

         CALL logger_warn( "MPP DEL VAR: domain decomposition not define "//&
         &               " in mpp strcuture "//TRIM(td_mpp%c_name))

      ELSE

         ! check if attribute already in mpp structure
         il_attid=0
         IF( ASSOCIATED(td_mpp%t_proc(1)%t_att) )THEN
            il_attid=att_get_index( td_mpp%t_proc(1)%t_att(:), &
            &                    td_att%c_name )
         ENDIF
         IF( il_attid == 0 )THEN
            CALL logger_warn( &
            &  "MPP DEL VAR: no attribute "//TRIM(td_att%c_name)//&
            &  ", in mpp structure "//TRIM(td_mpp%c_name) )

            IF( ASSOCIATED(td_mpp%t_proc(1)%t_att) )THEN
               DO ji=1,td_mpp%t_proc(1)%i_natt
                  CALL logger_debug( "MPP DEL ATT: in mpp structure : &
                  &  attribute : "//TRIM(td_mpp%t_proc(1)%t_att(ji)%c_name) )
               ENDDO
            ENDIF

         ELSE

            cl_name=TRIM(td_att%c_name)
            CALL logger_debug( "MPP DEL ATT: delete in mpp structure : &
            &  attribute : "//TRIM(cl_name) )
            DO ji=1,td_mpp%i_nproc
               CALL file_del_att(td_mpp%t_proc(ji), TRIM(cl_name)) 
            ENDDO

         ENDIF

      ENDIF
   END SUBROUTINE mpp__del_att_str
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine delete attribute in mpp structure, given attribute name.
   !>
   !> @detail
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !> @date February, 2015 
   !> - define local attribute structure to avoid mistake with pointer
   !
   !> @param[inout] td_mpp    mpp strcuture
   !> @param[in]    cd_name   attribute name
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_att_name( td_mpp, cd_name )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP)       , INTENT(INOUT) :: td_mpp
      CHARACTER(LEN=*) , INTENT(IN   ) :: cd_name

      ! local variable
      INTEGER(i4) :: il_attid
      TYPE(TATT)  :: tl_att
      !----------------------------------------------------------------
      ! check if mpp exist
      IF( .NOT. ASSOCIATED(td_mpp%t_proc) )THEN

         CALL logger_warn( "MPP DEL ATT: domain decomposition not define "//&
         &               " in mpp strcuture "//TRIM(td_mpp%c_name))

      ELSE

         IF( td_mpp%t_proc(1)%i_natt == 0 )THEN
            CALL logger_debug( "MPP DEL ATT NAME: no attribute associated to mpp &
            &                 structure "//TRIM(td_mpp%c_name) )
         ELSE

            ! get the attribute id, in file variable structure
            il_attid=0
            IF( ASSOCIATED(td_mpp%t_proc(1)%t_att) )THEN
               il_attid=att_get_id( td_mpp%t_proc(1)%t_att(:), &
               &                    cd_name )
            ENDIF

            IF( il_attid == 0 )THEN

               CALL logger_debug( &
               &  "MPP DEL ATT : there is no attribute with "//&
               &  "name "//TRIM(cd_name)//" in mpp structure "//&
               &  TRIM(td_mpp%c_name))

            ELSE

               tl_att=att_copy(td_mpp%t_proc(1)%t_att(il_attid))
               CALL mpp_del_att(td_mpp, tl_att) 

            ENDIF
         ENDIF

      ENDIF
   END SUBROUTINE mpp__del_att_name
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine overwrite attribute in mpp structure.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in]    td_att attribute structure
   !-------------------------------------------------------------------
   SUBROUTINE mpp_move_att( td_mpp, td_att )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp
      TYPE(TATT), INTENT(IN)    :: td_att

      !local variable
      TYPE(TATT)  :: tl_att
      !----------------------------------------------------------------
      ! copy variable
      tl_att=att_copy(td_att)

      ! remove processor
      CALL mpp_del_att(td_mpp, tl_att)

      ! add processor
      CALL mpp_add_att(td_mpp, tl_att)

      ! clean
      CALL att_clean(tl_att)

   END SUBROUTINE mpp_move_att
   !-------------------------------------------------------------------
   !> @brief
   !>    This subroutine compute domain decomposition for niproc and njproc 
   !> processors following I and J.
   !>
   !> @detail
   !> To do so, it need to know :
   !> - global domain dimension
   !> - overlap region length
   !> - number of processors following I and J
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp mpp strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp__compute( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP), INTENT(INOUT) :: td_mpp

      ! local variable
      INTEGER(i4)                              :: il_isize !< i-direction maximum sub domain size 
      INTEGER(i4)                              :: il_jsize !< j-direction maximum sub domain size
      INTEGER(i4)                              :: il_resti !<  
      INTEGER(i4)                              :: il_restj !<  
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_nlci
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_nlcj
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_impp
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: il_jmpp

      CHARACTER(LEN=lc)                        :: cl_file
      TYPE(TFILE)                              :: tl_proc
      TYPE(TATT)                               :: tl_att

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      INTEGER(i4) :: jk
      !----------------------------------------------------------------

      ! intialise
      td_mpp%i_nproc=0

      CALL logger_trace( "MPP COMPUTE: compute domain decomposition with "//&
      &               TRIM(fct_str(td_mpp%i_niproc))//" x "//&
      &               TRIM(fct_str(td_mpp%i_njproc))//" processors")
      ! maximum size of sub domain
      il_isize = ((td_mpp%t_dim(1)%i_len - 2*td_mpp%i_preci + (td_mpp%i_niproc-1))/ &
      &           td_mpp%i_niproc) + 2*td_mpp%i_preci
      il_jsize = ((td_mpp%t_dim(2)%i_len - 2*td_mpp%i_precj + (td_mpp%i_njproc-1))/ &
      &           td_mpp%i_njproc) + 2*td_mpp%i_precj

      il_resti = MOD(td_mpp%t_dim(1)%i_len - 2*td_mpp%i_preci, td_mpp%i_niproc)
      il_restj = MOD(td_mpp%t_dim(2)%i_len - 2*td_mpp%i_precj, td_mpp%i_njproc)
      IF( il_resti == 0 ) il_resti = td_mpp%i_niproc
      IF( il_restj == 0 ) il_restj = td_mpp%i_njproc

      ! compute dimension of each sub domain
      ALLOCATE( il_nlci(td_mpp%i_niproc,td_mpp%i_njproc) )
      ALLOCATE( il_nlcj(td_mpp%i_niproc,td_mpp%i_njproc) )

      il_nlci( 1 : il_resti                , : ) = il_isize
      il_nlci( il_resti+1 : td_mpp%i_niproc, : ) = il_isize-1

      il_nlcj( : , 1 : il_restj                ) = il_jsize
      il_nlcj( : , il_restj+1 : td_mpp%i_njproc) = il_jsize-1

      ! compute first index of each sub domain
      ALLOCATE( il_impp(td_mpp%i_niproc,td_mpp%i_njproc) )
      ALLOCATE( il_jmpp(td_mpp%i_niproc,td_mpp%i_njproc) )

      il_impp(:,:)=1
      il_jmpp(:,:)=1

      DO jj=1,td_mpp%i_njproc
         DO ji=2,td_mpp%i_niproc
            il_impp(ji,jj)=il_impp(ji-1,jj)+il_nlci(ji-1,jj)-2*td_mpp%i_preci
         ENDDO
      ENDDO

      DO jj=2,td_mpp%i_njproc
         DO ji=1,td_mpp%i_niproc
            il_jmpp(ji,jj)=il_jmpp(ji,jj-1)+il_nlcj(ji,jj-1)-2*td_mpp%i_precj
         ENDDO
      ENDDO 

      DO jj=1,td_mpp%i_njproc
         DO ji=1,td_mpp%i_niproc

            jk=ji+(jj-1)*td_mpp%i_niproc

            ! get processor file name
            cl_file=file_rename(td_mpp%c_name,jk)
            ! initialise file structure
            tl_proc=file_init(cl_file,td_mpp%c_type)

            ! procesor id
            tl_proc%i_pid=jk

            tl_att=att_init("DOMAIN_number",tl_proc%i_pid)
            CALL file_add_att(tl_proc, tl_att)

            ! processor indices
            tl_proc%i_iind=ji
            tl_proc%i_jind=jj

            ! fill processor dimension and first indices
            tl_proc%i_impp = il_impp(ji,jj)
            tl_proc%i_jmpp = il_jmpp(ji,jj)

            tl_att=att_init( "DOMAIN_poistion_first", &
            &                (/tl_proc%i_impp, tl_proc%i_jmpp/) )
            CALL file_add_att(tl_proc, tl_att)

            tl_proc%i_lci  = il_nlci(ji,jj)
            tl_proc%i_lcj  = il_nlcj(ji,jj)

            tl_att=att_init( "DOMAIN_poistion_last", &
            &                (/tl_proc%i_lci, tl_proc%i_lcj/) )
            CALL file_add_att(tl_proc, tl_att)

            ! compute first and last indoor indices
            
            ! west boundary
            IF( ji == 1 )THEN
               tl_proc%i_ldi = 1 
               tl_proc%l_ctr = .TRUE.
            ELSE
               tl_proc%i_ldi = 1 + td_mpp%i_preci
            ENDIF

            ! south boundary
            IF( jj == 1 )THEN
               tl_proc%i_ldj = 1 
               tl_proc%l_ctr = .TRUE.
            ELSE
               tl_proc%i_ldj = 1 + td_mpp%i_precj
            ENDIF

            ! east boundary
            IF( ji == td_mpp%i_niproc )THEN
               tl_proc%i_lei = il_nlci(ji,jj)
               tl_proc%l_ctr = .TRUE.
            ELSE
               tl_proc%i_lei = il_nlci(ji,jj) - td_mpp%i_preci
            ENDIF

            ! north boundary
            IF( jj == td_mpp%i_njproc )THEN
               tl_proc%i_lej = il_nlcj(ji,jj)
               tl_proc%l_ctr = .TRUE.
            ELSE
               tl_proc%i_lej = il_nlcj(ji,jj) - td_mpp%i_precj
            ENDIF

            tl_att=att_init( "DOMAIN_halo_size_start", &
            &                (/tl_proc%i_ldi, tl_proc%i_ldj/) )
            CALL file_add_att(tl_proc, tl_att)
            tl_att=att_init( "DOMAIN_halo_size_end", &
            &                (/tl_proc%i_ldi, tl_proc%i_ldj/) )
            CALL file_add_att(tl_proc, tl_att)

            ! add processor to mpp structure
            CALL mpp__add_proc(td_mpp, tl_proc)

            ! clean
            CALL att_clean(tl_att)
            CALL file_clean(tl_proc)

         ENDDO
      ENDDO

      DEALLOCATE( il_impp, il_jmpp )
      DEALLOCATE( il_nlci, il_nlcj )

   END SUBROUTINE mpp__compute
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine remove land processor from domain decomposition.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in] id_mask   sub domain mask (sea=1, land=0)
   !-------------------------------------------------------------------
   SUBROUTINE mpp__del_land( td_mpp, id_mask )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),                  INTENT(INOUT) :: td_mpp
      INTEGER(i4), DIMENSION(:,:), INTENT(IN)    :: id_mask

      ! loop indices
      INTEGER(i4) :: jk
      !----------------------------------------------------------------

      IF( ASSOCIATED(td_mpp%t_proc) )THEN
         jk=1
         DO WHILE( jk <= td_mpp%i_nproc )
            IF( mpp__land_proc(td_mpp, jk, id_mask(:,:)) )THEN
               CALL mpp__del_proc(td_mpp, jk)
            ELSE
               jk=jk+1
            ENDIF
         ENDDO
      ELSE
         CALL logger_error("MPP DEL LAND: domain decomposition not define.")
      ENDIF

   END SUBROUTINE mpp__del_land
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine optimize the number of sub domain to be used, given mask.
   !> @details
   !>  Actually it get the domain decomposition with the most land 
   !>  processor removed.
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in] id_mask   sub domain mask (sea=1, land=0) 
   !-------------------------------------------------------------------
   SUBROUTINE mpp__optimiz( td_mpp, id_mask )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),                  INTENT(INOUT) :: td_mpp
      INTEGER(i4), DIMENSION(:,:), INTENT(IN)    :: id_mask

      ! local variable
      TYPE(TMPP)  :: tl_mpp
      INTEGER(i4) :: il_maxproc

      TYPE(TFILE), DIMENSION(:), ALLOCATABLE :: tl_proc
      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jj
      !----------------------------------------------------------------

      CALL logger_trace("MPP OPTIMIZ: look for best domain decomposition")
      tl_mpp=mpp_copy(td_mpp)

      ! save maximum number of processor to be used
      il_maxproc=td_mpp%i_nproc
      ! 
      td_mpp%i_nproc=0
      DO ji=1,il_maxproc
         DO jj=1,il_maxproc

            ! clean mpp processor
            IF( ASSOCIATED(tl_mpp%t_proc) )THEN
               CALL file_clean(tl_mpp%t_proc(:))
               DEALLOCATE(tl_mpp%t_proc)
            ENDIF

            ! compute domain decomposition
            tl_mpp%i_niproc=ji
            tl_mpp%i_njproc=jj
            
            CALL mpp__compute( tl_mpp )
            
            ! remove land sub domain
            CALL mpp__del_land( tl_mpp, id_mask )

            CALL logger_info("MPP OPTIMIZ: number of processor "//&
            &   TRIM(fct_str(ji))//"x"//TRIM(fct_str(jj))//"="//&
            &   TRIM(fct_str(tl_mpp%i_nproc)) )
            IF( tl_mpp%i_nproc > td_mpp%i_nproc .AND. &
            &   tl_mpp%i_nproc <= il_maxproc )THEN
               ! save optimiz decomposition 

               CALL logger_info("MPP OPTIMIZ:save this decomposition "//&
               &   TRIM(fct_str(ji))//"x"//TRIM(fct_str(jj))//"="//&
               &   TRIM(fct_str(tl_mpp%i_nproc)) )

               ! clean mpp
               CALL mpp_clean(td_mpp)

               ! save processor array
               ALLOCATE( tl_proc(tl_mpp%i_nproc) )
               tl_proc(:)=file_copy(tl_mpp%t_proc(:))

               ! remove pointer on processor array
               CALL file_clean(tl_mpp%t_proc(:))
               DEALLOCATE(tl_mpp%t_proc)
 
               ! save data except processor array
               td_mpp=mpp_copy(tl_mpp)

               ! save processor array
               ALLOCATE( td_mpp%t_proc(td_mpp%i_nproc) )
               td_mpp%t_proc(:)=file_copy(tl_proc(:))

               ! clean
               CALL file_clean( tl_proc(:) )
               DEALLOCATE(tl_proc)

            ENDIF
            
         ENDDO
      ENDDO

      ! clean
      CALL mpp_clean(tl_mpp)

   END SUBROUTINE mpp__optimiz
   !-------------------------------------------------------------------
   !> @brief
   !>    This function check if processor is a land processor.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[in] td_mpp    mpp strcuture
   !> @param[in] id_proc   processor id
   !> @param[in] id_mask   sub domain mask (sea=1, land=0)
   !-------------------------------------------------------------------
   LOGICAL FUNCTION mpp__land_proc( td_mpp , id_proc, id_mask )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),                  INTENT(IN) :: td_mpp
      INTEGER(i4),                 INTENT(IN) :: id_proc
      INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: id_mask

      ! local variable
      INTEGER(i4), DIMENSION(2) :: il_shape
      !----------------------------------------------------------------

      CALL logger_trace("MPP LAND PROC: check processor "//TRIM(fct_str(id_proc))//&
      &  " of mpp "//TRIM(td_mpp%c_name) )
      mpp__land_proc=.FALSE.
      IF( ASSOCIATED(td_mpp%t_proc) )THEN

         il_shape(:)=SHAPE(id_mask)
         IF( il_shape(1) /= td_mpp%t_dim(1)%i_len .OR. &
         &   il_shape(2) /= td_mpp%t_dim(2)%i_len )THEN
             CALL logger_debug("MPP LAND PROC: mask size ("//&
             &                  TRIM(fct_str(il_shape(1)))//","//&
             &                  TRIM(fct_str(il_shape(2)))//")")
             CALL logger_debug("MPP LAND PROC: domain size ("//&
             &                  TRIM(fct_str(td_mpp%t_dim(1)%i_len))//","//&
             &                  TRIM(fct_str(td_mpp%t_dim(2)%i_len))//")")
             CALL logger_error("MPP LAND PROC: mask and domain size differ")
         ELSE
            IF( ALL(id_mask( td_mpp%t_proc(id_proc)%i_impp +            &
            &                       td_mpp%t_proc(id_proc)%i_ldi - 1 : &
            &                td_mpp%t_proc(id_proc)%i_impp +            &
            &                       td_mpp%t_proc(id_proc)%i_lei - 1,  &
            &                td_mpp%t_proc(id_proc)%i_jmpp +            &
            &                       td_mpp%t_proc(id_proc)%i_ldj - 1 : &
            &                td_mpp%t_proc(id_proc)%i_jmpp +            &
            &                       td_mpp%t_proc(id_proc)%i_lej - 1)  &
            &      /= 1 ) )THEN
               ! land domain
               CALL logger_info("MPP LAND PROC: processor "//TRIM(fct_str(id_proc))//&
               &             " is land processor")
               mpp__land_proc=.TRUE.
            ENDIF
         ENDIF

      ELSE
         CALL logger_error("MPP LAND PROC: domain decomposition not define.")
      ENDIF

   END FUNCTION mpp__land_proc
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine clean mpp strcuture.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp__clean_unit( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),  INTENT(INOUT) :: td_mpp

      ! local variable
      TYPE(TMPP) :: tl_mpp ! empty mpp structure

      ! loop indices
      !----------------------------------------------------------------

      CALL logger_info( &
      &  "MPP CLEAN: reset mpp "//TRIM(td_mpp%c_name) )

      ! del dimension
      IF( td_mpp%i_ndim /= 0 )THEN
         CALL dim_clean( td_mpp%t_dim(:) )
      ENDIF

      IF( ASSOCIATED(td_mpp%t_proc) )THEN
         ! clean array of file processor
         CALL file_clean( td_mpp%t_proc(:) )
         DEALLOCATE(td_mpp%t_proc)
      ENDIF

      ! replace by empty structure
      td_mpp=mpp_copy(tl_mpp)

   END SUBROUTINE mpp__clean_unit
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine clean mpp strcuture.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp__clean_arr( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),  DIMENSION(:), INTENT(INOUT) :: td_mpp

      ! local variable
      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------

      DO ji=SIZE(td_mpp(:)),1,-1
         CALL mpp_clean(td_mpp(ji))
      ENDDO

   END SUBROUTINE mpp__clean_arr
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine get sub domains which cover "zoom domain".
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !> @param[in] id_imin   i-direction lower indice
   !> @param[in] id_imax   i-direction upper indice
   !> @param[in] id_jmin   j-direction lower indice
   !> @param[in] id_jmax   j-direction upper indice
   !-------------------------------------------------------------------
   SUBROUTINE mpp__get_use_unit( td_mpp, id_imin, id_imax, &
   &                                     id_jmin, id_jmax )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP) ,  INTENT(INOUT) :: td_mpp
      INTEGER(i4),  INTENT(IN), OPTIONAL :: id_imin
      INTEGER(i4),  INTENT(IN), OPTIONAL :: id_imax
      INTEGER(i4),  INTENT(IN), OPTIONAL :: id_jmin
      INTEGER(i4),  INTENT(IN), OPTIONAL :: id_jmax

      ! local variable
      LOGICAL     :: ll_iuse
      LOGICAL     :: ll_juse

      INTEGER(i4) :: il_imin
      INTEGER(i4) :: il_imax
      INTEGER(i4) :: il_jmin
      INTEGER(i4) :: il_jmax

      ! loop indices
      INTEGER(i4) :: jk
      !----------------------------------------------------------------
      IF( ASSOCIATED(td_mpp%t_proc) )THEN
   
         il_imin=1
         il_imax=td_mpp%t_dim(1)%i_len
         IF( PRESENT(id_imin) ) il_imin=id_imin
         IF( PRESENT(id_imax) ) il_imax=id_imax
         il_jmin=1
         il_jmax=td_mpp%t_dim(2)%i_len
         IF( PRESENT(id_jmin) ) il_jmin=id_jmin
         IF( PRESENT(id_jmax) ) il_jmax=id_jmax

         ! check domain
         IF( il_imin < 1 .OR. il_imin > td_mpp%t_dim(1)%i_len .OR. &
         &   il_imax < 1 .OR. il_imax > td_mpp%t_dim(1)%i_len .OR. &
         &   il_jmin < 1 .OR. il_jmin > td_mpp%t_dim(2)%i_len .OR. &
         &   il_jmax < 1 .OR. il_jmax > td_mpp%t_dim(2)%i_len )THEN
            CALL logger_debug("MPP GET USE: mpp gloabl size "//&
            &        TRIM(fct_str(td_mpp%t_dim(1)%i_len))//","//&
            &        TRIM(fct_str(td_mpp%t_dim(2)%i_len)))
            CALL logger_debug("MPP GET USE: i-indices "//&
            &        TRIM(fct_str(il_imin))//","//TRIM(fct_str(il_imax)))
            CALL logger_debug("MPP GET USE: j-indices "//&
            &        TRIM(fct_str(il_jmin))//","//TRIM(fct_str(il_jmax)))
            CALL logger_error("MPP GET USE: invalid indices ")
         ELSE
            td_mpp%t_proc(:)%l_use=.FALSE.
            DO jk=1,td_mpp%i_nproc

               ! check i-direction
               ll_iuse=.FALSE.
               IF( il_imin < il_imax )THEN

                  ! not overlap east west boundary
                  IF( td_mpp%t_proc(jk)%i_impp + td_mpp%t_proc(jk)%i_lci > &
                  &   il_imin .AND.                                  &
                  &   td_mpp%t_proc(jk)%i_impp < il_imax )THEN
                      ll_iuse=.TRUE.
                  ENDIF

               ELSEIF( il_imin == il_imax )THEN

                  ! east west cyclic
                  ll_iuse=.TRUE.

               ELSE ! il_imin > id_imax

                  ! overlap east west boundary
                  IF( ( td_mpp%t_proc(jk)%i_impp + td_mpp%t_proc(jk)%i_lci >  &
                  &     il_imin )                                             &
                  &   .OR.                                                    &
                  &   ( td_mpp%t_proc(jk)%i_impp < il_imax) )THEN
                     ll_iuse=.TRUE.
                  ENDIF

               ENDIF

               ! check j-direction
               ll_juse=.FALSE.
               IF( il_jmin < il_jmax )THEN

                  ! not overlap north fold
                  IF( td_mpp%t_proc(jk)%i_jmpp + td_mpp%t_proc(jk)%i_lcj > &
                  &   il_jmin .AND.                                  &
                  &   td_mpp%t_proc(jk)%i_jmpp < il_jmax )THEN
                     ll_juse=.TRUE.
                  ENDIF

               ELSE ! id_jmin >= id_jmax

                  IF( td_mpp%t_proc(jk)%i_jmpp + td_mpp%t_proc(jk)%i_lcj > &
                  &  il_jmin )THEN
                     ll_juse=.TRUE.
                  ENDIF

               ENDIF

               IF( ll_iuse .AND. ll_juse ) td_mpp%t_proc(jk)%l_use=.TRUE.

            ENDDO
         ENDIF

      ELSE
         CALL logger_error("MPP GET USE: mpp decomposition not define.")
      ENDIF

   END SUBROUTINE mpp__get_use_unit
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine get sub domains which form global domain border.
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp_get_contour( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),  INTENT(INOUT) :: td_mpp

      ! loop indices
      INTEGER(i4) :: jk
      !----------------------------------------------------------------

      IF( ASSOCIATED(td_mpp%t_proc) )THEN

         td_mpp%t_proc(:)%l_use = .FALSE.
         DO jk=1,td_mpp%i_nproc
            IF( td_mpp%t_proc(jk)%i_ldi == 1 .OR. &
            &   td_mpp%t_proc(jk)%i_ldj == 1 .OR. &
            &   td_mpp%t_proc(jk)%i_lei == td_mpp%t_proc(jk)%i_lci .OR. &
            &   td_mpp%t_proc(jk)%i_lej == td_mpp%t_proc(jk)%i_lcj )THEN

               td_mpp%t_proc(jk)%l_use = .TRUE.
 
            ENDIF
         ENDDO
   
      ELSE
         CALL logger_error("MPP GET CONTOUR: domain decomposition not define.")
      ENDIF

   END SUBROUTINE mpp_get_contour
   !-------------------------------------------------------------------
   !> @brief
   !> This function return processor indices, without overlap boundary,
   !> given processor id. 
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[in] td_mpp    mpp strcuture
   !> @param[in] id_procid processor id
   !> @return array of index (/ i1, i2, j1, j2 /)
   !-------------------------------------------------------------------
   FUNCTION mpp_get_proc_index( td_mpp, id_procid )
      IMPLICIT NONE

      ! Argument
      TYPE(TMPP) , INTENT(IN) :: td_mpp
      INTEGER(i4), INTENT(IN) :: id_procid

      ! function
      INTEGER(i4), DIMENSION(4) :: mpp_get_proc_index

      ! local variable
      INTEGER(i4) :: il_i1, il_i2
      INTEGER(i4) :: il_j1, il_j2
      !----------------------------------------------------------------

      IF( ASSOCIATED(td_mpp%t_proc) )THEN

         IF( TRIM(td_mpp%c_dom) == '' )THEN
            CALL logger_fatal("MPP GET PROC INDEX: decomposition type unknown. "//&
            &                 "you should ahve run mpp_get_dom before.")
         ENDIF

         SELECT CASE(TRIM(td_mpp%c_dom))
            CASE('full')
               il_i1 = 1 ; il_i2 = td_mpp%t_dim(1)%i_len
               il_j1 = 1 ; il_j2 = td_mpp%t_dim(2)%i_len
            CASE('overlap')
                il_i1 = td_mpp%t_proc(id_procid)%i_impp
                il_j1 = td_mpp%t_proc(id_procid)%i_jmpp

                il_i2 = il_i1 + td_mpp%t_proc(id_procid)%i_lci - 1 
                il_j2 = il_j1 + td_mpp%t_proc(id_procid)%i_lcj - 1 
            CASE('nooverlap')
               il_i1 = td_mpp%t_proc(id_procid)%i_impp + &
               &        td_mpp%t_proc(id_procid)%i_ldi - 1
               il_j1 = td_mpp%t_proc(id_procid)%i_jmpp + &
               &        td_mpp%t_proc(id_procid)%i_ldj - 1

               il_i2 = td_mpp%t_proc(id_procid)%i_impp + &
               &        td_mpp%t_proc(id_procid)%i_lei - 1
               il_j2 = td_mpp%t_proc(id_procid)%i_jmpp + &
               &        td_mpp%t_proc(id_procid)%i_lej - 1
            CASE DEFAULT
               CALL logger_error("MPP GET PROC INDEX: invalid decomposition type.")
         END SELECT

         mpp_get_proc_index(:)=(/il_i1, il_i2, il_j1, il_j2/)

      ELSE
         CALL logger_error("MPP GET PROC INDEX: domain decomposition not define.")
      ENDIF

   END FUNCTION mpp_get_proc_index
   !-------------------------------------------------------------------
   !> @brief
   !> This function return processor domain size, depending of domain 
   !> decompisition type, given sub domain id. 
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !
   !> @param[in] td_mpp    mpp strcuture
   !> @param[in] id_procid sub domain id
   !> @return array of index (/ isize, jsize /)
   !-------------------------------------------------------------------
   FUNCTION mpp_get_proc_size( td_mpp, id_procid )
      IMPLICIT NONE

      ! Argument
      TYPE(TMPP),  INTENT(IN) :: td_mpp
      INTEGER(i4), INTENT(IN) :: id_procid

      ! function
      INTEGER(i4), DIMENSION(2) :: mpp_get_proc_size

      ! local variable
      INTEGER(i4) :: il_isize
      INTEGER(i4) :: il_jsize
      !----------------------------------------------------------------

      IF( ASSOCIATED(td_mpp%t_proc) )THEN

         IF( TRIM(td_mpp%c_dom) == '' )THEN
            CALL logger_fatal("MPP GET PROC SIZE: decomposition type unknown. "//&
            &                 "you should ahve run mpp_get_dom before.")
         ENDIF

         SELECT CASE(TRIM(td_mpp%c_dom))
            CASE('full')
               
               il_isize = td_mpp%t_dim(1)%i_len
               il_jsize = td_mpp%t_dim(2)%i_len

            CASE('overlap')

                il_isize = td_mpp%t_proc(id_procid)%i_lci
                il_jsize = td_mpp%t_proc(id_procid)%i_lcj

            CASE('nooverlap')
               il_isize = td_mpp%t_proc(id_procid)%i_lei - &
               &          td_mpp%t_proc(id_procid)%i_ldi + 1
               il_jsize = td_mpp%t_proc(id_procid)%i_lej - &
               &          td_mpp%t_proc(id_procid)%i_ldj + 1
            CASE DEFAULT
               CALL logger_error("MPP GET PROC SIZE: invalid decomposition type : "//&
               &  TRIM(td_mpp%c_dom) )
         END SELECT

         mpp_get_proc_size(:)=(/il_isize, il_jsize/)

      ELSE
         CALL logger_error("MPP GET PROC SIZE: domain decomposition not define.")
      ENDIF

   END FUNCTION mpp_get_proc_size
   !-------------------------------------------------------------------
   !> @brief 
   !>  This subroutine determine domain decomposition type.
   !>  (full, overlap, noverlap)
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial version
   !>
   !> @param[inout] td_mpp mpp strcuture
   !-------------------------------------------------------------------
   SUBROUTINE mpp_get_dom( td_mpp )
      IMPLICIT NONE
      ! Argument
      TYPE(TMPP),  INTENT(INOUT) :: td_mpp

      ! local variable
      INTEGER(i4) :: il_isize
      INTEGER(i4) :: il_jsize
      !----------------------------------------------------------------

      IF( ASSOCIATED(td_mpp%t_proc) )THEN

         IF( td_mpp%i_niproc == 0 .AND. td_mpp%i_niproc == 0 )THEN
            CALL logger_info("MPP GET DOM: use indoor indices to get domain "//&
            &             "decomposition type.")
            IF((td_mpp%t_proc(1)%t_dim(1)%i_len ==                         &
            &   td_mpp%t_proc(1)%i_lei - td_mpp%t_proc(1)%i_ldi + 1) .AND. &
            &  (td_mpp%t_proc(1)%t_dim(2)%i_len ==                         &
            &   td_mpp%t_proc(1)%i_lej - td_mpp%t_proc(1)%i_ldj + 1) )THEN

               td_mpp%c_dom='nooverlap'

            ELSEIF((td_mpp%t_proc(1)%t_dim(1)%i_len ==                     &
            &       td_mpp%t_proc(1)%i_lci                     )     .AND. &
            &      (td_mpp%t_proc(1)%t_dim(2)%i_len ==                     &
            &       td_mpp%t_proc(1)%i_lcj                     )     )THEN

               td_mpp%c_dom='overlap'

            ELSEIF((td_mpp%t_proc(1)%t_dim(1)%i_len ==                     &
            &       td_mpp%t_dim(1)%i_len             )              .AND. &
            &      (td_mpp%t_proc(1)%t_dim(2)%i_len ==                     &
            &       td_mpp%t_dim(2)%i_len )                          )THEN

               td_mpp%c_dom='full'

            ELSE

               CALL logger_error("MPP GET DOM: should have been an impossible case")

               il_isize=td_mpp%t_proc(1)%t_dim(1)%i_len
               il_jsize=td_mpp%t_proc(1)%t_dim(2)%i_len
               CALL logger_debug("MPP GET DOM: proc size "//&
               &  TRIM(fct_str(il_isize))//" x "//TRIM(fct_str(il_jsize)) )

               il_isize=td_mpp%t_proc(1)%i_lei - td_mpp%t_proc(1)%i_ldi + 1
               il_jsize=td_mpp%t_proc(1)%i_lej - td_mpp%t_proc(1)%i_ldj + 1
               CALL logger_debug("MPP GET DOM: no overlap size "//&
               &  TRIM(fct_str(il_isize))//" x "//TRIM(fct_str(il_jsize)) )

               il_isize=td_mpp%t_proc(1)%i_lci
               il_jsize=td_mpp%t_proc(1)%i_lcj
               CALL logger_debug("MPP GET DOM: overlap size "//&
               &  TRIM(fct_str(il_isize))//" x "//TRIM(fct_str(il_jsize)) )

               il_isize=td_mpp%t_dim(1)%i_len
               il_jsize=td_mpp%t_dim(2)%i_len
               CALL logger_debug("MPP GET DOM: full size "//&
               &  TRIM(fct_str(il_isize))//" x "//TRIM(fct_str(il_jsize)) )

            ENDIF

         ELSE

            CALL logger_info("MPP GET DOM: use number of processors following "//&
            &             "I and J to get domain decomposition type.")
            IF( td_mpp%i_niproc*td_mpp%i_njproc==td_mpp%i_nproc )THEN
               IF( td_mpp%i_nproc == 1 )THEN
                  td_mpp%c_dom='full'
               ENDIF
               td_mpp%c_dom='nooverlap'
            ELSE
               td_mpp%c_dom='overlap'
            ENDIF

         ENDIF

      ELSE
         CALL logger_error("MPP GET DOM: domain decomposition not define.")
      ENDIF

   END SUBROUTINE mpp_get_dom
   !-------------------------------------------------------------------
   !> @brief This function check if variable  and mpp structure use same
   !> dimension.
   !>
   !> @details
   !>
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !>
   !> @param[in] td_mpp mpp structure
   !> @param[in] td_var variable structure
   !> @return dimension of variable and mpp structure agree (or not)
   !-------------------------------------------------------------------
   LOGICAL FUNCTION mpp__check_var_dim(td_mpp, td_var)
      IMPLICIT NONE
      ! Argument      
      TYPE(TMPP), INTENT(IN) :: td_mpp
      TYPE(TVAR), INTENT(IN) :: td_var

      ! local variable

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      mpp__check_var_dim=.TRUE.
      ! check used dimension 
      IF( ANY( td_var%t_dim(:)%l_use .AND. &
      &        td_var%t_dim(:)%i_len /= td_mpp%t_dim(:)%i_len) )THEN

         mpp__check_var_dim=.FALSE.

         CALL logger_debug( &
         &  " mpp dimension: "//TRIM(fct_str(td_mpp%i_ndim))//&
         &  " variable dimension: "//TRIM(fct_str(td_var%i_ndim)) )
         DO ji = 1, ip_maxdim
            CALL logger_debug( &
            &  "MPP CHECK DIM: for dimension "//&
            &  TRIM(td_mpp%t_dim(ji)%c_name)//&
            &  ", mpp length: "//&
            &  TRIM(fct_str(td_mpp%t_dim(ji)%i_len))//&
            &  ", variable length: "//&
            &  TRIM(fct_str(td_var%t_dim(ji)%i_len))//&
            &  ", variable used "//TRIM(fct_str(td_var%t_dim(ji)%l_use)))
         ENDDO

         CALL logger_error( &
         &  "MPP CHECK DIM: variable and mpp dimension differ"//&
         &  " for variable "//TRIM(td_var%c_name)//&
         &  " and mpp "//TRIM(td_mpp%c_name))

      ENDIF

   END FUNCTION mpp__check_var_dim
   !-------------------------------------------------------------------
   !> @brief This function return the mpp id, in a array of mpp
   !> structure,  given mpp base name. 
   !
   !> @author J.Paul
   !> @date November, 2013 - Initial Version
   !
   !> @param[in] td_file   array of file structure
   !> @param[in] cd_name   file name
   !> @return file id in array of file structure (0 if not found)
   !-------------------------------------------------------------------
   INTEGER(i4) FUNCTION mpp_get_index(td_mpp, cd_name)
      IMPLICIT NONE
      ! Argument      
      TYPE(TMPP)      , DIMENSION(:), INTENT(IN) :: td_mpp
      CHARACTER(LEN=*),               INTENT(IN) :: cd_name

      ! local variable
      CHARACTER(LEN=lc) :: cl_name
      INTEGER(i4)       :: il_size

      ! loop indices
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      mpp_get_index=0
      il_size=SIZE(td_mpp(:))

      cl_name=TRIM( file_rename(cd_name) )

      ! check if mpp is in array of mpp structure
      DO ji=1,il_size
         ! look for file name
         IF( TRIM(fct_lower(td_mpp(ji)%c_name)) == TRIM(fct_lower(cd_name)) )THEN
 
            mpp_get_index=ji
            EXIT

         ENDIF
      ENDDO

   END FUNCTION mpp_get_index
   !-------------------------------------------------------------------
   !> @brief This function recombine variable splitted mpp structure. 
   !
   !> @author J.Paul
   !> @date Ocotber, 2014 - Initial Version
   !
   !> @param[in] td_mpp   mpp file structure
   !> @param[in] cd_name  variable name
   !> @return variable strucutre
   !-------------------------------------------------------------------
   TYPE(TVAR) FUNCTION mpp_recombine_var(td_mpp, cd_name) 
   IMPLICIT NONE
      ! Argument      
      TYPE(TMPP)      , INTENT(IN) :: td_mpp
      CHARACTER(LEN=*), INTENT(IN) :: cd_name

      ! local variable
      INTEGER(i4)                       :: il_varid
      INTEGER(i4)                       :: il_status
      INTEGER(i4)                       :: il_i1p
      INTEGER(i4)                       :: il_i2p
      INTEGER(i4)                       :: il_j1p
      INTEGER(i4)                       :: il_j2p
      INTEGER(i4), DIMENSION(4)         :: il_ind

      INTEGER(i4), DIMENSION(ip_maxdim) :: il_strt
      INTEGER(i4), DIMENSION(ip_maxdim) :: il_cnt

      TYPE(TVAR)                        :: tl_tmp
      TYPE(TVAR)                        :: tl_var

      ! loop indices
      INTEGER(i4) :: ji
      INTEGER(i4) :: jk
      !----------------------------------------------------------------

      il_varid=var_get_index( td_mpp%t_proc(1)%t_var(:), cd_name)
      IF( il_varid /= 0 )THEN
      
         tl_var=var_copy(td_mpp%t_proc(1)%t_var(il_varid))
         ! Allocate space to hold variable value in structure 
         IF( ASSOCIATED(tl_var%d_value) )THEN
            DEALLOCATE(tl_var%d_value)   
         ENDIF
         ! 
         DO ji=1,ip_maxdim
            IF( tl_var%t_dim(ji)%l_use )THEN
               tl_var%t_dim(ji)%i_len=td_mpp%t_dim(ji)%i_len
            ENDIF
         ENDDO

         ALLOCATE(tl_var%d_value( tl_var%t_dim(1)%i_len, &
         &                        tl_var%t_dim(2)%i_len, &
         &                        tl_var%t_dim(3)%i_len, &
         &                        tl_var%t_dim(4)%i_len),&
         &        stat=il_status)
         IF(il_status /= 0 )THEN

           CALL logger_error( &
            &  " MPP RECOMBINE VAR: not enough space to put variable "//&
            &  TRIM(tl_var%c_name)//" in variable structure")

         ENDIF

         ! FillValue by default
         tl_var%d_value(:,:,:,:)=tl_var%d_fill

         ! read processor 
         DO jk=1,td_mpp%i_nproc
            IF( td_mpp%t_proc(jk)%l_use )THEN
               ! get processor indices
               il_ind(:)=mpp_get_proc_index( td_mpp, jk )
               il_i1p = il_ind(1)
               il_i2p = il_ind(2)
               il_j1p = il_ind(3)
               il_j2p = il_ind(4)
 
               il_strt(:)=(/ 1,1,1,1 /)

               il_cnt(:)=(/ il_i2p-il_i1p+1,         &
               &            il_j2p-il_j1p+1,         &
               &            tl_var%t_dim(3)%i_len, &
               &            tl_var%t_dim(4)%i_len /)

               tl_tmp=iom_read_var( td_mpp%t_proc(jk), tl_var%c_name,&
               &                    il_strt(:), il_cnt(:) )
               
               ! replace value in output variable structure
               tl_var%d_value( il_i1p : il_i2p,  &
               &               il_j1p : il_j2p,  &
               &               :,:) = tl_tmp%d_value(:,:,:,:)

               ! clean
               CALL var_clean(tl_tmp)

            ENDIF
         ENDDO

         mpp_recombine_var=var_copy(tl_var)

         ! clean
         CALL var_clean(tl_var)

      ELSE

         CALL logger_error( &
         &  " MPP RECOMBINE VAR: there is no variable with "//&
         &  "name or standard name"//TRIM(cd_name)//&
         &  " in mpp file "//TRIM(td_mpp%c_name))
      ENDIF
   END FUNCTION mpp_recombine_var
END MODULE mpp

