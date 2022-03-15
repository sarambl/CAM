# Implementation:


- chemistry handles production of SOG
- SOG then goes into the condtend and oracle handles how much SOG goes to SOA
- SOA_A1 = sum(SOG) + cond_SOA_LV
- coagulation and dep see only SOA_A1
- Afterwards we get:
  - SOA_A1 pre loss process = SOA_A1_1
  - SOA_A1 after loss process = SOA_A1_2
  - lossR = (SOA_A1_2-SOA_A1_1)/SOA_A1_1
  - SOGx_2 = lossR*SOGx_1
  - Need to track also SOA_LV cond