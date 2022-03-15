#include "messy_main_ppd_bi.inc" ! mz_rj_20140410

  MODULE messy_oracle_si
  !
  ! DESCRIPTION
  ! -----------
  ! oracle INTERFACE LAYER FOR ECHAM5/MESSY
  !
  ! AUTHOR
  ! ------
  ! Alexandra Tsimpidi,  Forschungzentrum Juelich, IEK-8:Troposphere, Germany
  ! Vlassis Karydis,    Forschungzentrum Juelich, IEK-8:Troposphere, Germany
  ! questions/suggestions: a.tsimpidi@fz-juelich.de
  !
  ! LAST MODIFICATIONS - 
  !*****************************************************************************

  USE messy_oracle
  USE messy_main_tracer,        ONLY: t_ident
  USE messy_main_constants_mem, ONLY: DP,STRLEN_MEDIUM, STRLEN_ULONG 
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  ! POINTER TO STREAM ELEMENTS
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  ! USE SPECIAL NCREGRID EVENT TRIGGER
  USE messy_main_blather_bi,    ONLY: info_bi, warning_bi, error_bi
  IMPLICIT NONE
  PRIVATE
 
  ! SUBROUTINES
  PUBLIC :: oracle_initialize              ! initialization
  PUBLIC :: oracle_new_tracer              ! define tracers
  PUBLIC :: oracle_init_memory             ! allocate memory
  PUBLIC :: oracle_init_tracer             ! initialize tracers
  PUBLIC :: oracle_init_coupling           ! coupling to echam5
  PUBLIC :: oracle_vdiff                   ! distribute online emissions
  PUBLIC :: oracle_physc                   ! calls oracle core layer
  PUBLIC :: oracle_radiation               ! 
  PUBLIC :: oracle_free_memory             ! deallocate memory 

  INTRINSIC ABS, ASSOCIATED, ALLOCATED, MAX, TRIM

  INTEGER, PUBLIC,  SAVE :: npre !,nmode,tmode(3)
  INTEGER, PUBLIC,  SAVE :: trac_id_gp(1) = 0  !location of gas phase speies in tracer array
!idt

   INTEGER :: idt_LfPOG01  = 0 ! created here as non reactive (no reaction in MECCA!)
   INTEGER :: idt_LbbPOG01 = 0 ! created here as non reactive (no reaction in MECCA!)  
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_POA
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_POG
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_SOA
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_SOG
!carbon emission mass fluxes insoluble part in kg/kg
!  REAL(dp) ,DIMENSION(:,:), POINTER :: OC_sum_insol => NULL()
!carbon emission mass fluxes soluble part in kg/kg
!  REAL(dp) ,DIMENSION(:,:), POINTER :: OC_sum_sol => NULL()

! carbon emissions
  CHARACTER (LEN=STRLEN_MEDIUM) :: Cemis_channel  = ''
! organic carbon
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_soa_sol   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_ff_sol    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_bb_sol    = ''
   !--- CPL namelist fields:
 
 !-------------------------------------------------------------------------------------------
 !   for emissions
 
  TYPE emspec
    CHARACTER(LEN=STRLEN_MEDIUM) :: name      ! tracer name
    INTEGER                      :: trac_idx  ! tracer index
    REAL(dp)                     :: molarmass ! molar mass
    REAL(dp)                     :: frac      ! fraction of the total emission flux used
    INTEGER                      :: mode      ! mode of the target species
    LOGICAL                      :: L_numb    ! species is a number concentration and not a mass
  END TYPE emspec

   TYPE emflux
    ! name for each flux (fix, helps for identification)
    CHARACTER(LEN=STRLEN_MEDIUM)        :: name
    ! 3D-array for the mass flux
    REAL(dp), DIMENSION(:,:,:), POINTER :: flux          => NULL()
    ! 2D-array for the mass flux
    REAL(dp), DIMENSION(:,:),   POINTER :: flux_2D       => NULL()
    ! 3D-array for the corresponding number flux (if it exists)
    REAL(dp), DIMENSION(:,:,:), POINTER :: nflux         => NULL()
    ! 2D-array for the corresponding number flux (if it exists)
    REAL(dp), DIMENSION(:,:),   POINTER :: nflux_2D      => NULL()
    ! 3D-array for the V(ertical)IND(ex) for NxD emissions
    REAL(dp), DIMENSION(:,:,:), POINTER :: vind          => NULL()
    ! density = native density of the emission flux
    REAL(dp)                            :: density
    ! total_frac = total scaling factor for the emission flux
    ! controlled via the parameters.inc file
    REAL(dp)                            :: total_frac
    ! num_spec_emis = number of species which get a tendency from this flux
    INTEGER                             :: num_spec_emis
    ! dim = used dimension of the emission flux array (2D,3D)
    INTEGER                             :: dim
    ! dim_orig = native dimension of the emission flux array (2D,3D)
    INTEGER                             :: dim_orig
    ! NxD = logical whether a 3D flux originates from a NxD flux
    LOGICAL                             :: NxD
    ! mode = native mode of the flux
    INTEGER                             :: mode
    ! scal_fac = scaling factor if for a flux another flux is used and scaled 
    REAL(dp)                            :: scal_fac
    ! diameter = value for the aerosol diameter associated with this flux
    !            used in case of determining the number from the mass flux
    REAL(dp)                            :: diameter
    ! fac_num_emis = conversion factor (depending on the emission flux) to convert
    !                mass mean to count median (also including density if required)
    REAL(dp)                            :: fac_num_emis
    ! unit = unit of the emission flux -> determines conversion of emission
    CHARACTER(LEN=STRLEN_MEDIUM)        :: unit
    ! flux_name = name of the corresponding channel element
    CHARACTER(LEN=STRLEN_OBJECT)        :: flux_name
    ! nflux_name = name of the corresponding number flux channel element
    CHARACTER(LEN=STRLEN_OBJECT)        :: nflux_name
    ! channel_name = name of the corresponding channel
    CHARACTER(LEN=STRLEN_CHANNEL)       :: channel_name
    TYPE(emspec), DIMENSION(:), POINTER :: specs         => NULL()
   END TYPE emflux

   TYPE(emflux), DIMENSION(:),   POINTER, SAVE :: emis_flux_array => NULL()
   TYPE(emflux), SAVE                          :: emis_flux_list(100)
   INTEGER, SAVE                               :: num_fluxes

   REAL(dp), POINTER, DIMENSION(:,:,:,:) :: dryradius   => NULL() !gmxelink
 
   ! EMIS_CASK = Character array which contains for each of maximum 500 emission 
   !              fields the necessary information for emission assignment:
   !          1 = name of the emission object (for identification)
   !          2 = total scaling factor for incoming flux
   !          3 = channel / channel name of emission flux
   !          4 = channel object name of mass emission flux
   !          5 = corresponding name emission flux (if exists)
   !          6 = list of tracers which should get a value from this emission
   !              ";" separated list of tracers (fullname, CASE SENSITIVE)
   !          7 = list of scaling factors for each tracer
   !              ";" separated list of REAL values
   CHARACTER(LEN=STRLEN_ULONG), DIMENSION(50,7), SAVE :: EMIS_CASK


  NAMELIST /CPL/ Cemis_channel,  EMIS_CASK

CONTAINS
!==============================================================================

  SUBROUTINE oracle_initialize
!
   USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi
   USE messy_main_tools,        ONLY: find_next_free_unit
   USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, p_pe
 
   IMPLICIT NONE
   ! LOCAL
   CHARACTER(LEN=*), PARAMETER     :: substr = 'oracle_initialize' ! name of subroutine
   INTEGER                         :: iou    ! I/O unit
   INTEGER                         :: status ! error status
   INTEGER                         :: j,jm,jc
  
   IF (p_parallel_io) &
   CALL start_message_bi(modstr,'INITIALIZATION',substr)
 
   !--- Read namelist and control variables:
  
   ! INITIALIZE MAIN-CTRL
   IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      ! *** CALL CORE ROUTINE:
      CALL oracle_read_nml_ctrl(status, iou)
      IF (status /= 0)  CALL error_bi('error oracle_read_nml_ctrl', substr)
   END IF

   !--- Read CPL namelist
   IF (p_parallel_io) THEN
      EMIS_CASK(:,:) =''
      iou = find_next_free_unit(100,200)
      CALL oracle_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi("error in coupling namelist", substr)
   END IF

   !--- Broadcast over processors:
   CALL p_bcast (NPOA,       p_io)
   CALL p_bcast (NSOA,       p_io)
   CALL p_bcast (NSOAP,        p_io)
   CALL p_bcast (aermod,       p_io)
   CALL p_bcast (nmode,        p_io)
   CALL p_bcast (tmode,        p_io)
   CALL p_bcast (Cemis_channel,p_io)

   DO jm=1,50
    DO jc=1,7
     CALL p_bcast(EMIS_CASK(jm,jc),p_io)
    END DO
   END DO

 !write(200+p_pe,*) 'emis_OC_insol=',emis_OC_insol

!   call flush (200+p_pe)

   DO j=1,NSOAP
    CALL p_bcast (mwsoap(j), p_io)
    CALL p_bcast (csat(j), p_io)
    CALL p_bcast (cstemp(j), p_io)
    CALL p_bcast (deltah(j), p_io)
    CALL p_bcast (flagsoap(j), p_io)
!   call flush (200+p_pe)
   END DO

   !--- Initialize core:

!   CALL oracle_initialize_core
 
   IF (p_parallel_io) &
   CALL end_message_bi(modstr,'INITIALIZATION',substr)
  
  END SUBROUTINE oracle_initialize
 
!***************************************************************************

  SUBROUTINE oracle_read_nml_cpl(status, iou)
 
! oracle MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
! read namelist for 'coupling' to channel containing organic emissions
! Authors: Alexandra Tsimpidi, FZJ, IEK-8, 2021
!          Vlassis Karydis,    FZJ, IEK-8, 2021
 
! MESSy
   USE messy_main_tools,          ONLY: read_nml_open, read_nml_check &
                                       , read_nml_close
   USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
   USE messy_main_mpi_bi,        ONLY: p_parallel_io

  IMPLICIT NONE


! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    INTRINSIC TRIM

   IF(p_parallel_io) &
     CALL start_message_bi(modstr,'Reading coupling namelist',substr)
 
    status = 1
 
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)

    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

     IF(p_parallel_io) &
     CALL end_message_bi(modstr,'Reading coupling namelist',substr)

  END SUBROUTINE oracle_read_nml_cpl
 
!*****************************************************************************



!*****************************************************************************

  SUBROUTINE oracle_new_tracer

   USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
   USE messy_main_tracer_bi,     ONLY: tracer_halt 
   USE messy_main_tracer,        ONLY: new_tracer, get_tracer, set_tracer,   &
                                       AIR, ON, OFF, MODAL, AEROSOL,         &
                                       AMOUNTFRACTION, NUMBERDENSITY,        &
                                       I_ADVECT, I_CONVECT,                  &
                                       I_VDIFF, I_WETDEP,                    &
                                       I_DRYDEP, I_SEDI,                     &
                                       I_SCAV, I_MIX,                        &
                                       I_AEROSOL_METHOD, I_AEROSOL_MODE,     &
                                       I_AEROSOL_SOL, S_AEROSOL_MODEL,       &
                                       R_MOLARMASS, R_AEROSOL_DENSITY,       &
                                       R_henry,     R_dryreac_sf
   USE messy_main_tracer_bi,     ONLY: tracer_halt
   USE messy_main_tracer_mem_bi, ONLY: ti_gp, GPTRSTR
   USE messy_main_mpi_bi,        ONLY: p_parallel_io
   USE MESSY_MAIN_TOOLS,         ONLY: strcrack

   IMPLICIT NONE
! LOCAL
   INTEGER :: i, j, status
 
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_new_tracer'
   CHARACTER(LEN=2)            :: str_soa 
   CHARACTER(LEN=2)            :: str_mode(3)

     IF(p_parallel_io) &
   CALL start_message_bi(modstr, 'TRACER DEFINITION', substr)

   ALLOCATE (idt_POA(NPOA,nmode))
   ALLOCATE (idt_SOA(NSOA,nmode))
!mz_ap_20150311+
! Add tracers non reactive  
      CALL new_tracer(status, GPTRSTR, "LPOG01",  &
           modstr, quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AIR, idx = idt_LfPOG01)
      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, R_molarmass       , 250.)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, R_henry           , 1.0d5)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, R_dryreac_sf      , 0._dp)
      CALL tracer_halt(substr, status)

!mz_ap_20150311-

  str_mode(1)= 'ks'
  str_mode(2)= 'as'
  str_mode(3)= 'cs'
   DO j=1,nmode
  npre = 0
! Add aerosols for POA species 
    DO i = 1,NPOA
    idt_POA(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "POA"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_POA(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_wetdep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), R_molarmass       , mwsoap(i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_POA(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
     npre=NPOA
! Add aerosols for SOA species 
    DO i = 1,NSOA
      idt_SOA(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "SOA"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_aSOAv(i,j))
      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_wetdep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), R_molarmass       , mwsoap(npre+i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_SOA(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)

    END DO

     IF(p_parallel_io) &
    CALL end_message_bi(modstr, 'TRACER DEFINITION', substr)

  END SUBROUTINE oracle_new_tracer
!***************************************************************************

  SUBROUTINE oracle_init_memory

!   ! oracle MODULE ROUTINE (ECHAM-5 INTERFACE)
!   !
!   ! define oracle specific channel(s) and allocate memory for
!   ! global fields
!   !
! Authors: Alexandra Tsimpidi, FZJ, IEK-8, 2021
!          Vlassis Karydis,    FZJ, IEK-8, 2021

    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp     
    USE messy_main_blather_bi,    ONLY: error_bi, info_bi
! 
!   ! ECHAM5/MESSy
   USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
   USE messy_main_mpi_bi,     ONLY: p_parallel_io
!
   IMPLICIT NONE

   INTEGER :: jt,i

!   ! LOCAL
   CHARACTER(LEN=2)            :: str_soa 
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_init_memory'
! 
     IF(p_parallel_io) &
   CALL start_message_bi(modstr, 'Get gasses from MECCA', substr)
!
!CHECK FOR CGs from MECCA

    DO i = 1,NPOA
    idt_POG (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LPOG"//str_soa) THEN
       idt_POG(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "POG"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF  (idt_fPOG(i) == 0) THEN
    CALL error_bi('POG not present in chemical mechanism!', substr)
   ENDIF
   END DO 

    DO i = 1,NSOA
    idt_SOG (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LSOG"//str_soa) THEN
       idt_SOG(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "SOG"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF (idt_SOG(i) == 0) THEN
    CALL error_bi('SOG not present in chemical mechanism!', substr)
   ENDIF
   END DO

  END SUBROUTINE oracle_init_memory
!***************************************************************************

!***************************************************************************
  SUBROUTINE oracle_init_coupling

     USE messy_main_channel,       ONLY: get_channel_object, get_channel_info
     USE messy_main_channel_bi,    ONLY: channel_halt
     USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
     USE messy_main_mpi_bi,       ONLY: p_parallel_io

   IMPLICIT NONE
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_init_coupling'
   INTEGER  :: ierr
   INTEGER  :: status ! error status         !gmxelink


     IF(p_parallel_io) &
     CALL start_message_bi(modstr,'INIT COUPLING',substr)

!gmxelink
        write(6,*) TRIM(aermod//'_gp')
         CALL get_channel_object(status, TRIM(aermod//'_gp'),&
              'dryradius', p4=dryradius)
         IF (status /= 0) & 
              CALL error_bi('dry radius channel object not found !',substr)
!gmxelink

! This part deals with the emission fluxes
   CALL oracle_emis_init_si

 
    IF(p_parallel_io) &
    CALL end_message_bi(modstr,'INIT COUPLING',substr)

  END SUBROUTINE oracle_init_coupling

!***************************************************************************

   SUBROUTINE oracle_emis_init_si


    USE MESSY_MAIN_TRACER,              ONLY: r_molarmass, i_aerosol_mode, &
                                              numberdensity
    USE MESSY_MAIN_TRACER_MEM_BI,       ONLY: ti_gp, ntrac => ntrac_gp
    USE messy_main_tools,               ONLY: strcrack, str2num

    USE messy_main_channel,             ONLY: get_channel_object,  &
                                              get_channel_object_info
    USE messy_main_channel_bi,          ONLY: channel_halt, GP_3D_MID, &
                                              GP_3D_1LEV, GP_2D_HORIZONTAL
    IMPLICIT NONE
    INTEGER :: jm, jc, jt, dummy, counter, status, id_repr
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring(:) => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring2(:) => NULL()

    LOGICAL                                   :: found
    CHARACTER(LEN=*), PARAMETER               :: substr='oracle_emis_si'
    REAL(dp)                                  :: val

    DO jt=1,100
      emis_flux_list(jt)%name         = ""
      emis_flux_list(jt)%flux_name    = ""
      emis_flux_list(jt)%nflux_name   = ""
      emis_flux_list(jt)%channel_name = ""
      emis_flux_list(jt)%density      = 0._dp
      emis_flux_list(jt)%diameter     = 0._dp
      emis_flux_list(jt)%mode         = 0
      emis_flux_list(jt)%scal_fac     = 1._dp
      emis_flux_list(jt)%total_frac   = 1._dp
      emis_flux_list(jt)%fac_num_emis = 1._dp
    END DO
! organic carbon
    emis_flux_list(1)%name         = "oc_mass_ff_ks"
    emis_flux_list(1)%density      = 2.0_dp
    emis_flux_list(1)%diameter     = 2._dp * 0.03e-6_dp
    emis_flux_list(1)%mode         = 2
    emis_flux_list(1)%unit         = "kg/(m^2 s)"

    emis_flux_list(2)%name         = "oc_mass_ff_as"
    emis_flux_list(2)%density      = 2.0_dp
    emis_flux_list(2)%diameter     = 2._dp * 0.3e-6_dp 
    emis_flux_list(2)%mode         = 3
    emis_flux_list(2)%unit         = "kg/(m^2 s)"

    emis_flux_list(3)%name         = "oc_mass_bb_ks"
    emis_flux_list(3)%density      = 2.0_dp
    emis_flux_list(3)%diameter     = 2._dp * 0.075e-6_dp
    emis_flux_list(3)%mode         = 2
!mz_ap_20140921+
!    emis_flux_list(3)%unit         = "kg/(m^2 s)"
    emis_flux_list(3)%unit         = "kg/(m^3 s)"
!mz_ap_20140921-

    emis_flux_list(4)%name         = "oc_mass_bb_as"
    emis_flux_list(4)%density      = 2.0_dp
    emis_flux_list(4)%diameter     = 2._dp * 0.75e-6_dp
    emis_flux_list(4)%mode         = 3
!mz_ap_20140921+
!    emis_flux_list(4)%unit         = "kg/(m^2 s)"
    emis_flux_list(4)%unit         = "kg/(m^3 s)"
!mz_ap_20140921-

    emis_flux_list(5)%name         = "oc_mass_ks"
    emis_flux_list(5)%density      = 2.0_dp
    emis_flux_list(5)%diameter     = 2._dp * 0.258E-7_dp
    emis_flux_list(5)%mode         = 2
    emis_flux_list(5)%unit         = "kg/(m^2 s)"

    emis_flux_list(6)%name         = "oc_mass_as"
    emis_flux_list(6)%density      = 2.0_dp
    emis_flux_list(6)%diameter     = 2._dp * 0.258E-6_dp
    emis_flux_list(6)%mode         = 3
    emis_flux_list(6)%unit         = "kg/(m^2 s)"

    emis_flux_list(7)%name         = "oc_ss_mass_as"
    emis_flux_list(7)%density      = 2.0_dp
    emis_flux_list(7)%diameter     = 2._dp * 0.258E-6_dp
    emis_flux_list(7)%mode         = 3
    emis_flux_list(7)%unit         = "kg/(m^2 s)"
 
    num_fluxes = 0
     DO jm=1,50
      found = .false.
      IF ( ADJUSTL(TRIM(EMIS_CASK(jm,1))) == "") CYCLE
      DO jt=1,100
        IF (TRIM(emis_flux_list(jt)%name) == TRIM(EMIS_CASK(jm,1)) ) THEN
          num_fluxes = num_fluxes + 1
          found = .TRUE.
        ENDIF
      END DO
      IF (.NOT. FOUND) CALL warning_bi("EMIS_CASK named "//&
        TRIM(EMIS_CASK(jm,1))//&
        " not found in the list of fluxes and is therefore ignored!", substr)
    END DO

    counter = 0
    ALLOCATE(emis_flux_array(num_fluxes))
    DO jm=1,50
      IF ( TRIM(EMIS_CASK(jm,1)) == "") CYCLE
      DO jt=1,100
        IF (TRIM(emis_flux_list(jt)%name) == TRIM(EMIS_CASK(jm,1)) ) THEN
          counter = counter + 1
          emis_flux_array(counter)%name         = emis_flux_list(jt)%name
          emis_flux_array(counter)%channel_name = TRIM(EMIS_CASK(jm,3))
          emis_flux_array(counter)%flux_name    = TRIM(EMIS_CASK(jm,4))
          emis_flux_array(counter)%nflux_name   = TRIM(EMIS_CASK(jm,5))
          emis_flux_array(counter)%density      = emis_flux_list(jt)%density
          emis_flux_array(counter)%diameter     = emis_flux_list(jt)%diameter
          emis_flux_array(counter)%mode         = emis_flux_list(jt)%mode
          emis_flux_array(counter)%unit         = emis_flux_list(jt)%unit
          emis_flux_array(counter)%fac_num_emis = emis_flux_list(jt)%fac_num_emis
           IF ( TRIM(EMIS_CASK(jm,2)) == "") THEN 
          emis_flux_array(counter)%total_frac   = emis_flux_list(jt)%total_frac
           ELSE
            call str2num(TRIM(EMIS_CASK(jm,2)), val)
             emis_flux_array(counter)%total_frac = val
           END IF

          call strcrack(EMIS_CASK(jm,6), ';', outstring,  &
            emis_flux_array(counter)%num_spec_emis)
          call strcrack(EMIS_CASK(jm,7), ';', outstring2, dummy)


          ALLOCATE(emis_flux_array(counter)%specs(emis_flux_array(counter)%num_spec_emis))

          DO jc = 1, emis_flux_array(counter)%num_spec_emis
            val = 0._dp
            emis_flux_array(counter)%specs(jc)%name = TRIM(outstring(jc))
            call str2num(TRIM(outstring2(jc)), val)
             emis_flux_array(counter)%specs(jc)%frac = val
          ENDDO

        END IF
      END DO
    END DO
    DO jm = 1,num_fluxes
     DO jc = 1,emis_flux_array(jm)%num_spec_emis
         emis_flux_array(jm)%specs(jc)%trac_idx = 0
         emis_flux_array(jm)%specs(jc)%molarmass = 1._dp
         emis_flux_array(jm)%specs(jc)%mode = 0
         emis_flux_array(jm)%specs(jc)%l_numb = .FALSE.
        DO jt = 1,ntrac
          IF (TRIM(emis_flux_array(jm)%specs(jc)%name) == &
            ti_gp(jt)%tp%ident%fullname ) THEN
            emis_flux_array(jm)%specs(jc)%trac_idx = jt

            emis_flux_array(jm)%specs(jc)%molarmass = &
             ti_gp(jt)%tp%meta%cask_r(R_molarmass)
            emis_flux_array(jm)%specs(jc)%mode = &
              ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE)
             IF (ti_gp(jt)%tp%ident%quantity == numberdensity) &
              emis_flux_array(jm)%specs(jc)%l_numb = .TRUE.
          END IF
        END DO
      END DO
    END DO

    DO jm = 1,num_fluxes
      CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
        TRIM(emis_flux_array(jm)%flux_name), p3=emis_flux_array(jm)%flux)
      IF (status /= 0) &
       CALL error_bi(&
       'requested object name for element '//TRIM(emis_flux_array(jm)%name)//&
       ' not found in channel '//TRIM(emis_flux_array(jm)%channel_name), substr)

       CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
         TRIM(emis_flux_array(jm)%nflux_name), p3=emis_flux_array(jm)%nflux)

      IF (status /= 0) THEN
        CALL info_bi( &
        'requested object name for number flux for element '// &
        TRIM(emis_flux_array(jm)%name)//&
        ' not found; calculating number from mass!', substr)
      ENDIF
 
      CALL get_channel_object_info(status,                     &
        TRIM(emis_flux_array(jm)%channel_name),                &
        TRIM(emis_flux_array(jm)%flux_name), reprid = id_repr )
      emis_flux_array(jm)%dim = 0
      emis_flux_array(jm)%NxD = .FALSE.
      IF (id_repr == GP_2D_HORIZONTAL) emis_flux_array(jm)%dim      = 2
      IF (id_repr == GP_2D_HORIZONTAL) emis_flux_array(jm)%dim_orig = 2
      IF (id_repr == GP_3D_MID)        emis_flux_array(jm)%dim      = 3
      IF (id_repr == GP_3D_MID)        emis_flux_array(jm)%dim_orig = 3
      IF (id_repr == GP_3D_1LEV)       emis_flux_array(jm)%dim_orig = 3
      IF (id_repr == GP_3D_1LEV)       emis_flux_array(jm)%dim      = 2
      IF (emis_flux_array(jm)%dim_orig == 2) THEN
        emis_flux_array(jm)%flux_2D  => emis_flux_array(jm)%flux(:,:,1)
        IF (ASSOCIATED(emis_flux_array(jm)%nflux)) &
          emis_flux_array(jm)%nflux_2D => emis_flux_array(jm)%nflux(:,:,1)
      ELSE IF (emis_flux_array(jm)%dim_orig == 3) THEN
        IF (id_repr == GP_3D_1LEV) THEN
          emis_flux_array(jm)%dim = 2
          emis_flux_array(jm)%flux_2D  => emis_flux_array(jm)%flux(:,1,:)
         IF (ASSOCIATED(emis_flux_array(jm)%nflux)) &
            emis_flux_array(jm)%nflux_2D => emis_flux_array(jm)%nflux(:,1,:)
        ENDIF
      END IF

      ! NxD emissions 
      IF ( emis_flux_array(jm)%dim == 0 ) THEN
        CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
          TRIM(emis_flux_array(jm)%flux_name)//'_vind', &
          p3=emis_flux_array(jm)%vind)
         IF (status /= 0) &
           CALL error_bi(&
           'requested object index array for element '//TRIM(emis_flux_array(jm)%name)//&
           ' not found in channel '//TRIM(emis_flux_array(jm)%channel_name), substr)
        emis_flux_array(jm)%dim = 3
        emis_flux_array(jm)%NxD = .TRUE.
      END IF

    END DO

  END SUBROUTINE oracle_emis_init_si


!***************************************************************************

  SUBROUTINE oracle_init_tracer
  
   IMPLICIT NONE

  END SUBROUTINE oracle_init_tracer

!***************************************************************************
   SUBROUTINE oracle_global_start
 
     IMPLICIT NONE
 
   END SUBROUTINE oracle_global_start
!***************************************************************************
 
  SUBROUTINE oracle_vdiff
! oracle MODULE ROUTINE (ECHAM-5 INTERFACE)
!
! distribute online organic emissions
!
! Authors: Alexandra Tsimpidi, FZJ, IEK-8, 2021
!          Vlassis Karydis,    FZJ, IEK-8, 2021
 
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, ti_gp
    USE messy_main_constants_mem, ONLY: g,M_air,STRLEN_MEDIUM, avo => N_A
    USE messy_main_data_bi,       ONLY: nlev, philat_2d &
                                      , nproma, jrow, kproma     &
! mz_rj_20140918
#if defined(ECHAM5)
                                      , pxtems &
#endif
                                      , tm1_3d, tte_3d   &
                                      , qte_3d, qm1_3d, press_3d &
                                      , grmass, grvol, pressi_3d, geopoti_3d
#ifndef MESSYTIMER
    USE messy_main_data_bi,       ONLY: time_step_len &
                                      , nstep=>current_time_step, lstart
#else
    USE messy_main_timer,         ONLY: time_step_len &
                                      , nstep=>current_time_step, lstart
#endif
    USE messy_main_mpi_bi,        ONLY: p_io, p_bcast, p_parallel_io
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi

   IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='oracle_vdiff'
    INTEGER                     :: jc,jt,jk,idx,jl, mlev, ji, idt
    REAL(dp), POINTER           :: zxtems(:,:)
    REAL(dp)                    :: zdp    (nproma,nlev)
    REAL(dp)                    :: conv_unit(nproma,nlev)
    REAL(dp)                    :: zdz    (nproma,nlev)
    REAL(dp)                    :: flux(nproma,nlev)
! mz_rj_20140918+
#if defined (ECHAM5)
    zxtems => pxtems(:,1,:,jrow)
#endif
! mz_rj_20140918-
     DO jk=1,nlev
       zdp(1:kproma,jk) = pressi_3d(1:kproma,_RI_YZp1_) - &
                          pressi_3d(1:kproma,_RI_YZ_)
     END DO
     zdz(:,:) = 0.0_dp
     zdz(1:kproma,2:nlev) = ( geopoti_3d(1:kproma,_RI_YZcm1_)     &
                          -   geopoti_3d(1:kproma,_RI_YZc2_) )/g
 ! Loop over emission fluxes

    DO jc = 1,num_fluxes
      flux(:,:)  = 0._dp

      conv_unit(1:kproma,1:nlev) = 1._dp
 
      DO jt = 1, emis_flux_array(jc)%num_spec_emis
        flux(:,:)  = 0._dp
        conv_unit(:,:) = 1._dp
        idx = emis_flux_array(jc)%specs(jt)%trac_idx
        IF (idx == 0) CYCLE
        SELECT CASE (TRIM(emis_flux_array(jc)%unit))
        CASE("kg/(m^2 s)")
          conv_unit(1:kproma,1:nlev) = &
            M_air / emis_flux_array(jc)%specs(jt)%molarmass
!mz_ap_20140921+
         CASE("kg/(m^3 s)")
           conv_unit(1:kproma,2:nlev) = &
             M_air / emis_flux_array(jc)%specs(jt)%molarmass*zdz(1:kproma,2:nlev) 
!mz_ap_20140921-
        CASE("molecules/(m^2 s)")
          conv_unit(1:kproma,2:nlev) = M_air / avo / 1.e3_dp
        CASE("molecules/(m^3 s)")
          conv_unit(1:kproma,2:nlev) = M_air / avo / 1.e3_dp * &
                                       zdz(1:kproma,2:nlev)
        END SELECT

        IF (emis_flux_array(jc)%dim == 3) THEN
          IF (emis_flux_array(jc)%NxD) THEN
            mlev = SIZE(emis_flux_array(jc)%VIND,2)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%flux(1:kproma,_RI_YZcm_)
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%flux_2D(1:kproma,jrow)
        ENDIF
 
          idt = idx
          IF (emis_flux_array(jc)%NxD) THEN
            DO ji=1,mlev
              DO jl=1,kproma
                jk= NINT(emis_flux_array(jc)%VIND(jl,ji,jrow))
                pxtte(jl,_RI_ZD_) = pxtte(jl,_RI_ZD_)        &
                  + conv_unit(jl,jk)                       &
                  * flux(jl,ji)                            &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ELSE
! mz_rj_20140410+
            DO jk=1,nlev
! mz_rj_20140410-
              DO jl=1,kproma
                pxtte(jl,_RI_ZD_) = pxtte(jl,_RI_ZD_)        &
                  + conv_unit(jl,jk)                       &
                  * flux(jl,jk)                            &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ENDIF
      END DO
    END DO

 END SUBROUTINE oracle_vdiff
!***************************************************************************

!***************************************************************************
  SUBROUTINE oracle_physc 
   USE messy_main_data_bi,        ONLY:           tm1       &    ! dry air temperature (at time step -1)  [K]
                                        ,         tte       &    ! dry air temperature (tendency)         [K/s]
                                        ,         press_3d  &    ! air pressure                           [Pa]
                                        ,         pressi_3d &    ! air pressure (interface)               [Pa]
                                        ,grvol, grmass      &    ! grid volume [m3] and mass [kg]
                                        ,klev => nlev       &
                                        ,nproma, kproma     &
                                        ,nrow => ngpblks    &
                                        ,jrow             
   USE messy_main_constants_mem, ONLY:   M_air,dp
#ifndef MESSYTIMER
   USE messy_main_data_bi,  ONLY: delta_time,time_step_len &
                                ,nstep=>current_time_step
#else
   USE messy_main_timer,    ONLY: delta_time,time_step_len &
                                ,nstep=>current_time_step
#endif
   USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi
   USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
   USE messy_main_mpi_bi,         ONLY: p_io, p_bcast, &
                                        p_parallel_io, p_pe

   IMPLICIT NONE
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_physc'
   real(dp) caer(NSOAP), cgas(NSOAP), csatT(NSOAP),qaer(NSOAP,nmode),qgas(NSOAP)
   real(dp) tempk,xconvert(kproma,klev),idt_tes(4)    !gmxelink
   real(dp), DIMENSION(nmode) :: rsec
   REAL(dp), DIMENSION(kproma,klev) :: ztemp,zpress
   REAL(dp), DIMENSION(kproma,klev,nmode) :: zdryrad              !gmxelink
   integer i,j,jk,jl 
   integer idt,nlev
   LOGICAL,  SAVE              :: entered = .FALSE.

      IF(p_parallel_io.AND. .NOT. entered) &
     CALL start_message_bi(modstr,'CALL oracle core layer',substr)

  nlev = klev ! mz_rj_20140909
  ztemp       (:,:) = 0.0_dp
  zpress      (:,:) = 0.0_dp

   !--- Assign conversion factor to convert mol mol-1 to umol m-3:
   xconvert(1:kproma,1:klev)=(1.E+9_dp/M_air)* &                        !  [mol mol-1] => [umol g-1]
          (grmass(1:kproma,_RI_YZc_)/grvol(1:kproma,_RI_YZc_)) !  Air density [g m-3]: [umol g-1] => [umol m-3]
          

   !--- Ambient properties: ----------------------------------
   !--- Temperature [K]:
 
   ztemp  (1:kproma,1:klev)       = tm1  (1:kproma,_RI_YZc_)   +       &
                                    tte  (1:kproma,1:klev) * time_step_len
   !--- Pressure [Pa]:
 
   zpress  (1:kproma,1:klev)      = press_3d(1:kproma,_RI_YZc_)     ! [Pa]

   !--- Dry radius [m]:

   zdryrad (1:kproma,1:klev,1:nmode) = dryradius(1:kproma,1:klev,tmode(1):tmode(nmode),jrow) ! [m] gmxelink

   do jk = 1,klev
    do jl = 1,kproma

     tempk = ztemp  (jl,jk)
     caer=0.0_dp
     qaer=0.0_dp
     cgas=0.0_dp
     qgas=0.0_dp
     csatT=0.0_dp
     do j=1,nmode
      rsec(j) = zdryrad(jl,jk,j)     !gmxelink
      npre=0
      do i = 1,NPOA      
      idt = idt_POA(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(jl,_RI_ZD_)    &
                 + pxtte(jl,_RI_ZD_) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                              & ! conversion factor umol m-3 
                 * mwsoap(i)                                      !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      
      idt = idt_POG(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(jl,_RI_ZD_)     &
                 + pxtte(jl,_RI_ZD_) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                              & ! conversion factor umol m-3
                 * mwsoap(i)                                      !convert umol m-3 to ug m-3
      end do

      npre=npre+NPOA 
      do i = 1,NSOA      
      idt = idt_SOA(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(jl,_RI_ZD_)      &
                  + pxtte(jl,_RI_ZD_) * time_step_len))   & ! tracer mass mol mol-1
                  * xconvert(jl,jk)                               & ! conversion factor umol m-3
                  * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      

      idt = idt_SOG(i)
      if(j.eq.1) cgas(npre+i) = MAX(0._dp, (pxtm1(jl,_RI_ZD_)               &
                            + pxtte(jl,_RI_ZD_) * time_step_len )) & ! tracer mass mol mol-1
                            * xconvert(jl,jk)                             & ! conversion factor umol m-3
                            * mwsoap(npre+i)                                !convert umol m-3 to ug m-3

      end do
     end do

    do i=1,NSOAP
     qgas(i)=cgas(i)
    end do
!!!!!!!!!!!!!!!!!!!!!CALL gas/aerosol partitioning subroutine!!!!!!!!!!!!!!!!!!!!!!!!
     CALL oracle_soap(caer,cgas,tempk,p_pe,jl,jk,csatT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!CALL mode distribution subroutine!!!!!!!!!!!!!!!!!!!!!!!!
     if(nmode.ne.1) then
       CALL oracle_mode(caer,cgas,qaer,qgas,csatT,rsec,p_pe,jl,jk) !gmxelink
     else
       do i=1,NSOAP
        qaer(i,1)=caer(i)
       end do
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j=1,nmode
      npre=0
      do i = 1,NPOA      
      idt = idt_POA(i,j)
       pxtte(jl,_RI_ZD_) = pxtte(jl,_RI_ZD_)              &
                           + (qaer(i,j)/xconvert(jl,jk)/mwsoap(i)           &
                           - MAX(0._dp, (pxtm1(jl,_RI_ZD_)      &
                           + pxtte(jl,_RI_ZD_)*time_step_len)))   &
                           / time_step_len 

       idt = idt_POG(i)                    
       if(j.eq.1)        pxtte(jl,_RI_ZD_) = pxtte(jl,_RI_ZD_)                & 
                           + (cgas(i)/xconvert(jl,jk)/mwsoap(i)           &
                           - MAX(0._dp, (pxtm1(jl,_RI_ZD_)       &
                           + pxtte(jl,_RI_ZD_)*time_step_len)))    &
                           / time_step_len 
      end do

      npre=npre+NPOA
      do i = 1,NSOA      
       idt = idt_SOA(i,j)
       pxtte(jl,_RI_ZD_) = pxtte(jl,_RI_ZD_)                &
                           + (qaer(npre+i,j)/xconvert(jl,jk)/mwsoap(npre+i)   &
                           - MAX(0._dp, (pxtm1(jl,_RI_ZD_)        &
                           + pxtte(jl,_RI_ZD_)*time_step_len)))     &
                           / time_step_len 

       idt = idt_SOG(i)                    
       if(j.eq.1)        pxtte(jl,_RI_ZD_) = pxtte(jl,_RI_ZD_)                  &
                           + (cgas(npre+i)/xconvert(jl,jk)/mwsoap(npre+i)   &
                           - MAX(0._dp, (pxtm1(jl,_RI_ZD_)         &
                           + pxtte(jl,_RI_ZD_)*time_step_len)))      &
                           / time_step_len
      end do
     end do
    end do
   end do

      IF(p_parallel_io.AND. .NOT. entered) &
     CALL end_message_bi(modstr,'CALL oracle core layer',substr)
   entered = .TRUE.

  END SUBROUTINE oracle_physc 
!***************************************************************************
   SUBROUTINE oracle_radiation

    IMPLICIT NONE

   END SUBROUTINE oracle_radiation
!***************************************************************************
   SUBROUTINE oracle_free_memory

    IMPLICIT NONE

   IF (ALLOCATED (idt_POA)) DEALLOCATE (idt_POA)
   IF (ALLOCATED (idt_POG))  DEALLOCATE (idt_POG)
   IF (ALLOCATED (idt_SOA)) DEALLOCATE (idt_SOA)
   IF (ALLOCATED (idt_SOG))  DEALLOCATE (idt_SOG)

   END SUBROUTINE oracle_free_memory
!***************************************************************************
!***************************************************************************
  END MODULE messy_oracle_si
!***************************************************************************

