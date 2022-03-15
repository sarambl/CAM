# TODO

- Add SOA gas phase tracers
  - from BVOCs
  - from anthropogenics?
    - Use emissions from e.g. cam/bld/namelist_files/use_cases/2010_trop_strat_vbs_cam6.xml
    - cam/bld/namelist_files/namelist_defaults_cam.xml
  - Call oracle_soap from messy_oracle.F90 from condtend? Or before? 
    - INPUT:
      - ntot: total number of CG/SOA species pairs
      - caer (ntot): aerosol-phase concentrations of SOA species (ug/m3)
      - cgas(ntot): gas-phase concentrations of CG species (ppm or ug/m3)
      - tempk: temperature (K)
      - p_pe: not used
      - kproma: not really used
      - klev: not really used
      - csatT(ntot): saturation concentrations of CG/SOA species at current T  (ug/m3)
    - OUTPUT: 
      - caer (ntot): aerosol-phase concentrations of SOA species (ug/m3)
      - cgas (ntot): gas-phase concentrations of CG species (ppm or ug/m3)
      - csatT(ntot): saturation concentrations of CG/SOA species at current T  (ug/m3)


Variables? 
- NfPOA   - fossil fuel combustion POA 
- NbbPOA  - biomass burning POA
- NfSOAsv - fossil fuel SOA oxidation prod of SVOCs
- NbbSOAsv- biomass burdning SOA from oxidation of SVOCs
- NfSOAiv - fossil fuel SOA from oxidation of IVOCs and VOCs? 
- NbbSOAiv- biomass burdning SOA from oxidation of IVOCs and VOCs?
- NaSOAv  -  SOA from oxidation of anthropogenic VOCs?
- NbSOAv  -  SOA from oxidation of anthropogenic VOCs?
- NaOSOAv - 
- NbOSOAv





Emissions: 

Needed:  
- Benzene: VOC13-benzene-em-speciated-VOC-anthro
- Toluene: VOC14-toluene-em-speciated-VOC-anthro 
- Trimethyl_Benzenes: VOC16-trimethylb-em-speciated-VOC-anthro 
- Xylene: VOC15-xylene-em-speciated-VOC-anthro
- Other_Aromatics: VOC17-other-arom-em-speciated-VOC-anthro 
- Pentanes: VOC05-pentanes-em-speciated-VOC-anthro 
- Hexanes_and_higher_alkanes: VOC06-hexanes-pl-em-speciated-VOC-anthro 
- Propene: VOC08-propene-em-speciated-VOC-anthro  
- Other_alkenes_and_alkynes: VOC12-other-alke-em-speciated-VOC-anthro





LSOG01 from S/IVOCs with C*=0.01 μg/m3
LSOG02 from S/IVOCs with C*=10 μg/m3
LSOG03 from anthropogenic VOCs with C*=1 μg/m3
LSOG04 from anthropogenic VOCs with C*=1000 μg/m3  
LSOG05 from biogenic VOCs with C*=1 μg/m3
LSOG06 from biogenic VOCs with C*=1000 μg/m3




- How to deal with SOA_A1 vs 
  - OM_AI (OM coexisting with BC and coated, Aitken mode)
  - OM_AC (OM and SOA particles coagulated with other aerosols)
  - OM_NI? (OM emitted internally mixed with BC, Aitken mode)
Maybe we cannot deal with primary aerosols? 
  - OBS: nice detail SOA_A1 does not change with coagulation. 
  - Maybe we can get away with a fixed OM_AI/OM_AC ratio and leave OM_NI out of it? Will have to ask people at met how the currec


primary OM lifecylce in NorESM:
- OM_AI coagulate --> OM_AC
- OM_NI coated then transferred --> OM_AI (mode 04)
- OM_NI is emitted. 

Alternativs: 
1) don't fuck with this, it's annoying. Ignore primary organics and just go with it. 
2) Will anyways need to move some emissions to the anthropogenic. Therefore we will need to do something about it. Maybe we can treat these as a split more? And then remove the life-cycling...
3) OK: for all the real SOAs, we do what we've previously done --> all into SOA_A1
   1) for OM





Okey, how to deal: 

1) chemistry deals gives conc of SOG01-SOG06 and POG02, POG03
2) run condensation thingy and get total aerosol conc in particle phase. 
3) new dist between SOA01 vs SOG01 etc. 
4) SOA_A1 = SUM(SOA01:SOA06)
5) OM_AI+OM_AC = POA02 + POA03 ?? 
6) 
7) run rest of the code
8) resulting SOA_A1 and OM_AI,OM_AC:
   1) lossR = SOA_A1(start timestep)/SOA_A1(end timestep)
   2) let SOA01 = lossR*SOA01
   3) same for POM
9) 



Questions left: 
- How to deal with OD? 