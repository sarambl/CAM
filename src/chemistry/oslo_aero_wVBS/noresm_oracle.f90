 MODULE noresm_oracle
!_______________________________________________________________________________
! DESCRIPTION
! -----------
! oracle INTERFACE LAYER FOR ECHAM5/MESSY
!
! AUTHOR
! -----
! Sara Blichner, Stockholm University
!     based on work by
! Alexandra Tsimpidi,  Forschungzentrum Juelich, IEK-8:Troposphere, Germany

! -----------

     use vbs_def, only: csat_all_vbs, mw_all_vbs, cstemp_soa,cstemp_all_vbs, deltah_all_vbs, &
             floagsoap_pogpoa, NPOA, NSOAP, NSOAv
     use shr_kind_mod, only: r8 => shr_kind_r8
     use const, only : smallNumber

     IMPLICIT NONE


     ! SUBMODEL
     CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr  = 'oracle'     ! name of module
     CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver  = '1.0'    ! module version





     ! SUBROUTINES

     !PUBLIC :: oracle_read_nml_ctrl
     PUBLIC :: oracle_soap
     PUBLIC :: oracle_spfcn
     PUBLIC :: calc_loss
     !PUBLIC :: oracle_mode



     !
     !-----------------------------------------------------------------------
     !     Variables that are initialized in soapdat.f
     !
     !     mw_all_vbs  -- molecular weights of CG/SOA species (g/mol)
     !     csat_all_vbs    -- saturation concentrations of CG/SOA species (ug/m3)
     !     cstemp_all_vbs  -- temperatures corresponding to saturation concentrations
     !                of CG/SOA species (K)
     !    deltah_all_vbs  -- enthalpy of vaporization of CG/SOA species (kJ/mol)
     !     flagsoap-- set to 1 if CG/SOA species forms solutions; 0 if not
     !-----------------------------------------------------------------------
     !
     !DOUBLE PRECISION, PUBLIC    :: mw_all_vbs(40)
     !DOUBLE PRECISION, PUBLIC    :: csat(40)
     !DOUBLE PRECISION, PUBLIC    :: cstemp(40)
     !OUBLE PRECISION, PUBLIC    :: deltah(40)

     INTEGER, PUBLIC :: flagsoap(40)
!
!
!-----------------------------------------------------------------------


 CONTAINS



     subroutine oracle_soap(caer,cgas,tempk,kproma,klev,csatT)

         !USE messy_main_constants_mem,   ONLY: dp


         implicit none


         ! DEFINITION OF VARIABLES:
         !
         !  INPUTS
         !
         !     ntot    - total number of CG/SOA species pairs
         !     caer    - aerosol-phase concentrations of SOA species (ug/m3)
         !     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
         !     tempk   - temperature (K)
         !     convfac - conversion factor: umol/m3 = ppm * convfac
         !     iout    - standard output file unit
         !     igrdchm - index for grid containing the grid cell
         !     ichm    - i index of grid cell
         !     jchm    - j index of grid cell
         !     kchm    - k index of grid cell
         !     mwpre   - molecular weight of pre-existing organic aerosol (g/mol)
         !
         !   OUTPUTS
         !
         !     caer    - aerosol-phase concentrations of SOA species (ug/m3)
         !     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
         !     csatT   - saturation concentrations of CG/SOA species at current T
         !               (ug/m3)
         !
         !   VARIABLES USED WITHIN oracle_SOAP
         !
         !     i       - counter
         !     icont   - counter
         !     sumsc     - counter
         !     nsol    - total number of solution-forming SOA species
         !     cstemp_all_vbs  - temperatures corresponding to saturation concentrations
         !               of CG/SOA species (K)
         !     csat_all_vbs    - saturation concentrations of CG/SOA species (ug/m3)
         !    deltah_all_vbs  - enthalpy of vaporization of CG/SOA species (J/mol)
         !     flagsoap- set to 1 if CG/SOA species forms solutions; 0 if not
         !     mw_all_vbs  - molecular weights of CG/SOA species (g/mol)
         !     scaer   - aerosol-phase concentrations of solution-forming
         !               SOA species (ug/m3)
         !     scgas   - gas-phase concentrations of solution-forming
         !               SOA species (ug/m3)
         !     scsat   - saturation concentrations of solution-forming
         !               SOA species (ug/m3)
         !     sctot   - total concentrations of solution-forming SOA species
         !               (ug/m3)
         !     smw     - molecular weights of solution-forming SOA species
         !               (g/mol)
         !     znum    - counter for number of iterations
         !     conmin  - use simple solution for species below this level (ug/m3)
         !     cpremin - no pre-existing organics if cpre < cpremin (ug/m3)
         !     xtol    - error tolerance for bi-section method
         !
         !***********************************************************************
         !
         ! VARIABLE DECLARATION
         !
         real(r8), parameter       ::  conmin  = 1.d-6
         real(r8),parameter         ::  cpremin = 1.01d-9
         real(r8),parameter         ::  xtol    = 5.0d-5
         !conmin, cpremin, xtol
         !real (r8)      ::  conmin, cpremin, xtol
         !real (r8)      ::  conmin, cpremin, xtol
         !parameter ( conmin  = 1.d-6 )
         !parameter ( cpremin = 1.01d-9 )
         !parameter ( xtol    = 5.0d-5 )
!
         !integer     :: ntot
         real(r8),  intent(inout)   :: caer(NSOAP), cgas(NSOAP)
         real(r8),  intent(in)      :: tempk
         integer,   intent(in)      :: kproma, klev
         real(r8),  intent(out)     :: csatT(NSOAP)
         real(r8)                   :: ctot(NSOAP)
         real(r8)                   :: csat(NSOAP)
         !real (r8)   :: ctot(NSOAP)  ! caer(NSOAP), cgas(NSOAP),
         real (r8)   :: smw(NSOAP), scsat(NSOAP)
         real (r8)   :: sctot(NSOAP), scaer(NSOAP), scgas(NSOAP)
         integer     :: idx(NSOAP)
!
         real (r8)   :: mwpre, cpre, sumsc,convfac
         integer     :: iout
         integer     :: i, icont, nsol, znum

         real (r8)   :: cpx, bb, cc, xend, fend, xmid, fmid, dx
         integer     :: ntot
         integer     :: j,k
         !
         !***********************************************************************
         !
         ! Entry point
         !
         iout=6
         ntot=NSOAP
         do i=1,ntot
           ctot(i) = caer(i) + cgas(i)
         enddo
         !
         !
         mwpre = 220.d0
         cpre = 0.d0
         cpx = cpre/mwpre
         !
         ! CHANGE SATURATION CONCENTRATIONS ACCORDING TO CURRENT TEMPERATURE
         !
         do i=1,ntot
             csatT(i)=csat(i)*(cstemp_all_vbs(i)/tempk)*dexp((deltah_all_vbs(i)/8.314d0)&
                     &                *(1.d0/cstemp_all_vbs(i)-1.d0/tempk))
         enddo
         !
         ! CALCULATE AEROSOL-PHASE CONCENTRATION (CAER) AND GAS-PHASE
         ! CONCENTRATION (CGAS) FOR NON-SOLUTION-FORMING COMPOUNDS
         ! COMPOUNDS THAT HAVE A CONCENTRATION OF LESS THAN conmin ARE IGN0RED
         ! MAP COMPOUNDS THAT FORM SOLUTIONS ONTO ARRAYS
         !
         icont=0
         do i=1,ntot
             if (flagsoap(i).eq.0) then
                 cgas(i) = dmin1(ctot(i), csatT(i))
                 caer(i) = ctot(i) - cgas(i)
             elseif (ctot(i).lt.conmin) then
                 cgas(i) = ctot(i)
                 caer(i) = 0.d0
             else
                 icont=icont+1
                 idx(icont) = i
                 smw(icont)=mw_all_vbs(i)
                 scsat(icont)=csatT(i)
                 sctot(icont)=ctot(i)
                 scaer(icont)=caer(i)
             endif
         enddo
         nsol=icont
         !
         ! Check for a trivial solution
         !
         if (nsol.eq.0) goto 1000
         if (nsol.eq.1) then
             if (cpre.lt.cpremin) then
                 scgas(1) = dmin1(sctot(1), scsat(1))
                 scaer(1) = sctot(1) - scgas(1)
             else ! This case has an analytical solution
                 bb = scsat(1)-sctot(1)+cpx*smw(1)
                 cc = -sctot(1)*cpx*smw(1)
                 scaer(1) = dmin1( sctot(1), .5*(-bb+DSQRT(bb*bb-4.d0*cc)) )
                 scgas(1) = sctot(1) - scaer(1)
             endif
             goto 900
         endif
         sumsc=0.d0
         do i=1,nsol
             sumsc = sumsc + sctot(i)/scsat(i)
         enddo
         if (cpre.lt.cpremin .and. sumsc.le.1.d0) then
             do i=1,nsol
                 scgas(i)=sctot(i)
                 scaer(i)=0.d0
             enddo
             goto 900
         endif
         !
         ! Find the solution using a bi-section method (approach from max)
         !
         xend = 0.d0
         do i = 1, nsol
           xend = xend + sctot(i)/smw(i)
         enddo
         xend = xend + cpx
         call oracle_spfcn (nsol,sctot,scsat,scaer,smw,cpx,xend,fend)
         if (dabs(fend).le.xtol*xend) goto 99
         if (fend.gt.0.d0) then
             write (iout,'(//,a)') ' ERROR in oracle_SOAP:'
             write (iout,'(/,a)') ' ERROR: positive end point'
             goto 50
         endif
         dx = xend - cpx
         do znum = 1, 200
             dx = 0.5d0 * dx
             xmid = xend - dx
             call oracle_spfcn (nsol,sctot,scsat,scaer,smw,cpx,xmid,fmid)
             if (dabs(fmid).le.xtol*xmid .or. dx.le.xtol*xmid) goto 100
             if (fmid.lt.0.d0) xend = xmid
         enddo
         write (iout,'(//,a)') ' ERROR in oracle_SOAP:'
         write (iout,'(/,a)') ' ERROR: max number of iterations reached'
50       write (iout,'(a,i3,i4)') &
              &                 ' cell(kproma,klev) = ', kproma,klev
         write (iout,'(a5,2a15)') ' spec','total [ug/m3]','c* [ug/m3]'
         write (iout,'(i5,1p2e15.6)') (idx(i),sctot(i),scsat(i),i=1,nsol)
         write (iout,'(a5,e15.6)') ' cpre',cpre
         STOP
         !
         ! Converged
         !
99       xmid = xend
100      continue
         do i=1,nsol
             scaer(i) = dmin1( sctot(i), scaer(i) )
             scgas(i) = sctot(i) - scaer(i)
         enddo

         !
         ! REMAP COMPOUNDS THAT FORM SOLUTIONS BACK ONTO ORIGINAL ARRAYS
         !
900      continue
         do i=1,nsol
            caer(idx(i))=scaer(i)
            cgas(idx(i))=scgas(i)
         enddo
         !
         ! Convert to ppm if inputs in ppm
         !
1000     continue

         !



         return





     end subroutine oracle_soap





     subroutine oracle_spfcn (n,ct,cs,ca,mw,cpx,tom,fval)

         !USE messy_main_constants_mem,   ONLY: dp

         implicit none
         !
         ! oracle_SPFCN calculates the objective function for the bi-section solver in oracle_SOAP
         !     Total Organics in Mole (TOM) = sumsc_i(C_i(aer)/MW_i) + C_pre/MW_pre
         !     C_i(aer) = C_i(tot) - x_i * Cstar_i
         !              = C_i(tot) - (C_i(aer)/MW_i/TOM) * Cstar_i
         !  => C_i(aer) = C_i(tot) * TOM / (TOM + Cstar_i/MW_i)
         !  => sumsc_i(C_i(tot) * TOM / (TOM*MW_i + Cstar_i)) + C_pre/MW_pre - TOM = 0
         !
         ! Called by oracle_SOAP
         !
         integer     n,i
         real (r8)       ct(n),cs(n),ca(n),mw(n),cpx,tom,fval
         !
         fval = 0.d0
         do i = 1, n
             ca(i) = ct(i) * tom / ( tom + cs(i) / mw(i) )
             fval  = fval + ca(i) / mw(i)
         enddo
         fval = fval + cpx - tom
         !
         return

     end subroutine oracle_spfcn





     subroutine calc_loss(soa_a1_before, q, ncol, totalLoss)

         use constituents,  only: pcnst, cnst_name, cnst_get_ind
         use ppgrid,           only : pcols, pver

         use physconst,     only: avogad, rair
         use chem_mods,     only: gas_pcnst

         use vbs_def,       only: NPOA, NSOAv, NSOAP, oa_all_names_p00
         use aerosoldef,    only:l_soa_a1, chemistryIndex
         implicit none
         real(r8), intent(inout)    :: q(ncol,pver,gas_pcnst)             ! TMR [kg/kg] including moisture
         real(r8), intent(inout)    :: totalLoss(ncol,pver,gas_pcnst)     ! TMR [kg/kg] including moisture
         real(r8), intent(in)       :: soa_a1_before(ncol,pver)    ! TMR [kg/kg] including moisture
         integer,  intent(in)       :: ncol                       ! number of columns

         real(r8)                   :: soa_a1_after
         real(r8)                   :: lossR
         real(r8)                   :: ndxs_all_oap1(NSOAP+1)
         integer                    :: ndx_dummy, m, k, i, tracerIndex
        WRITE(*,*) 'Hey smb', pcols, pver, gas_pcnst
         do m=1,NSOAP+1

             call cnst_get_ind ( oa_all_names_p00(m), ndx_dummy,abort=.true.)
             ndxs_all_oap1(m) = chemistryIndex(ndx_dummy)
         end do



         do k=1,pver
             do i=1,ncol
                 !WRITE(*,*) q(i,k,chemistryIndex(l_soa_a1)), soa_a1_before(i,k)
                 soa_a1_after = q(i,k,chemistryIndex(l_soa_a1))
                 lossR = (soa_a1_after-soa_a1_before(i,k))/(soa_a1_before(i,k)+smallNumber)
                 lossR = max(lossR,0._r8) !make sure it's positive.

                 ! if
                 do m = 1,NSOAP+1
                     tracerIndex = ndxs_all_oap1(m)
                     ! reduce each tracer by the same factor
                     totalLoss(i,k,tracerIndex)= q(i,k,tracerIndex)*lossR ! kg_air/m3*ug_aer/kg_aer*kg_aer/kg_air
                     q(i,k,tracerIndex)= q(i,k,tracerIndex)*(1-lossR) ! kg_air/m3*ug_aer/kg_aer*kg_aer/kg_air
                 end do

             end do
         end do

         return

     end subroutine calc_loss

 END MODULE noresm_oracle


