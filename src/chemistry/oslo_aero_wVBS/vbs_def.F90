module vbs_def
    use shr_kind_mod,          only: r8 => shr_kind_r8

    implicit none
    save
    public

    integer, public, parameter ::  NPOA  = 3
    integer, public, parameter ::  NSOAv  = 6
    integer, public, parameter ::  NSOAP   = 9
    character(len=5), public, dimension(NSOAv), parameter    :: vbs_soag_names = (/'SOG01','SOG02','SOG03', 'SOG04','SOG05','SOG06'/)

    character(len=5), public, dimension(NPOA), parameter     :: vbs_pogpoa_names = (/'POA01','POG02','POG03'/)
    !character(len=5), public, dimension(NSOAP), parameter    :: og_all_names = (/'POG01','POG02','POG03', 'SOG01','SOG02','SOG03', 'SOG04','SOG05','SOG06'/)
    !character(len=5), public, dimension(NSOAP), parameter    :: og_all_names = (/'POG01','POG02','POG03', 'SOG01','SOG02','SOG03', 'SOG04','SOG05','SOG06' /)

    character(len=5), public, dimension(NSOAP), parameter    :: oa_all_names = (/'POA01','POA02','POA03',  'SOA01','SOA02','SOA03', 'SOA04','SOA05','SOA06'/)
    character(len=5), public, dimension(NSOAP), parameter    :: og_all_names = (/'POG01','POG02','POG03',  'SOG01','SOG02','SOG03', 'SOG04','SOG05','SOG06'/)
    character(len=5), public, dimension(NSOAP+1), parameter  :: oa_all_names_p00 = (/'POA01','POA02','POA03', &
            & 'SOA01','SOA02','SOA03', 'SOA04','SOA05','SOA06','SOA00'/)


    real(r8), public, dimension(NPOA)   ::  mw_pogpoa = (/250.d0,250.d0,250.d0 /)                               ! mwsoap(n11:n12)   : molecular weights of POG/POA species (g/mol)
    real(r8), public, dimension(NSOAv)  ::  mw_soa = (/250.d0 , 250.d0,   150.d0, 150.d0, 180.d0, 180.d0 /)     ! mwsoap(n11:n12)   : molecular weights of SOG/SOA species (g/mol)
    real(r8), public, dimension(NPOA)   ::  csat_pogpoa = (/1.e-2, 1.e1, 1.e4/)                                 ! saturation concentrations of POG/POA species (ug/m3)
    real(r8), public, dimension(NSOAv)  ::  csat_soa = (/  1.e-2,   1.e1,     1.e0,   1.e3,   1.e0,   1.e3 /)   ! saturation concentrations of SOG/SOA species (ug/m3)
    real(r8), public, dimension(NPOA)   ::  cstemp_pogpoa = (/300.0,300.0,300.0/)                               ! temperatures corresponding to csat
    real(r8), public, dimension(NSOAv)  ::  cstemp_soa = (/300.0 , 300.0,   300.0, 300.0, 300.0, 300.0/)        ! temperatures corresponding to csat
    real(r8), public, dimension(NPOA)   ::  deltah_pogpoa = (/112.3,94.3, 76.3/)                                ! deltah(1:NPOA)   : enthalpy of vaporization of POG/POA species (kJ/mol)
    real(r8), public, dimension(NSOAv)  ::  deltah_soa = (/112.3 ,  94.3,    30.3,  30.3,  30.3,  30.3 /)       ! deltah(n11:n12)   : enthalpy of vaporization of SOG/SOA species (kJ/mol)
    integer,  public, dimension(NPOA)   ::  floagsoap_pogpoa = (/1, 1,1 /)                                      ! flagsoap(n11:n12) : set to 1 if POG/POA species forms solutions; 0 if not
    integer,  public, dimension(NSOAv)  ::  floagsoap_soa = (/1, 1, 1, 1, 1,1 /)                                ! flagsoap(n11:n12) : set to 1 if SOG/SOA species forms solutions; 0 if not


    real(r8), public, dimension(NSOAP)  ::  mw_all_vbs = (/250.d0,250.d0,250.d0, 250.d0 , 250.d0,   150.d0, 150.d0, 180.d0, 180.d0 /)
    real(r8), public, dimension(NSOAP)   ::  csat_all_vbs = (/1.e-2, 1.e1, 1.e4,  1.e-2,   1.e1,     1.e0,   1.e3,   1.e0,   1.e3 /)
    real(r8), public, dimension(NSOAP)   ::  cstemp_all_vbs = (/300.0,300.0,300.0, 300.0 , 300.0,   300.0, 300.0, 300.0, 300.0/)     ! temperatures corresponding to csat

    real(r8), public, dimension(NSOAP)   ::  deltah_all_vbs = (/112.3, 94.3, 76.3,  112.3 ,  94.3,    30.3,  30.3,  30.3,  30.3 /)    ! deltah(1:NPOA)   : enthalpy of vaporization of POG/POA species (kJ/mol)
    integer,  public, dimension(NSOAP)   ::  floagsoap_all_vbs = (/1, 1,1, 1, 1, 1, 1, 1,1 /)                                        ! flagsoap(n11:n12) : set to 1 if POG/POA species forms solutions; 0 if not



end module vbs_def